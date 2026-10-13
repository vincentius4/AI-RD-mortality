* deepcars_tables.do - Stata 15 code for every table of the paper: Table 1-2 and eTable 1-14 -> results/JAMA_tables.xlsx (one sheet per table).
* Part 1 builds the landmark cohort (results/work/) and defines the estimation and formatting programs; parts 2-17 are one table each, in order.
* Each part was a separate do file in the previous release; the code is unchanged. Requires rdrobust (eTable 2) and rddensity 3.0 (eTable 4).
* Put the data files in this folder and run: do deepcars_tables.do

* PART 1. 01_data_and_programs
* Data and programs used by every table and figure file
* Output: results/JAMA_tables.xlsx (one sheet per manuscript table, cells as printed) and results/figures/ (figures)
version 15
clear all
set more off
set seed 7
cap mkdir results
cap mkdir results/work
cap mkdir results/figures
cap log close
log using "results/run_log.log", replace text
global OUTX "results/JAMA_tables.xlsx"
global W    "results/work"
global COVS age_z male cci htn sofa surgical_dept summer fall winter
global CUT  89.5

* Landmark cohort: AI score in the first 24 h and in hospital at 24 h
use "paper_data_v9.dta", clear
keep if rv_dec < .
keep if h_stay > 24
gen double s    = rv_dec - $CUT
gen byte   D    = (s >= 0)
gen byte   post = treat
gen int    disp = round(rv_dec)
gen byte   arr_nondnr = (arr_any == 1) | (death == 1 & dnr_case == 0)
gen byte   disc3 = (total_los <= 3 & death == 0)
save "$W/base.dta", replace

* ---------- estimation programs ----------
* seadj: small-sample factor applied to cluster-robust standard errors (as in the Python code that produced the paper)
program define seadj, rclass
    local adj = 1
    if "`e(cmd)'" == "regress" {
        if "`e(clustvar)'" != "" local adj = sqrt((e(N_clust)-1)/e(N_clust) * (e(N)-e(rank))/(e(N)-1))
        else local adj = sqrt((e(N)-e(rank))/e(N))
    }
    else {
        if "`e(clustvar)'" != "" local adj = sqrt((e(N_clust)-1)/e(N_clust))
        else local adj = sqrt((e(N)-1)/e(N))
    }
    return scalar adj = `adj'
end

* mkw: kernel weight within window h (triangular unless stated)
program define mkw
    args h kern
    cap drop kw
    if "`kern'" == "uni"      gen double kw = 1
    else if "`kern'" == "epa" gen double kw = max(0, 0.75*(1 - (s/`h')^2))
    else                      gen double kw = max(0, 1 - abs(s)/`h')
end

* didisc: difference-in-discontinuities (pooled local logistic model; coefficient on post x above threshold)
program define didisc, rclass
    syntax , h(real) [outcome(string) covs(string)]
    if "`outcome'" == "" local outcome death
    if "`covs'" == ""    local covs $COVS
    if "`covs'" == "none" local covs
    preserve
    qui keep if abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui gen double pD = post*D
    qui gen double ps = post*s
    qui gen double pDs = post*D*s
    qui gen double hB_post = hB*post
    qui gen double hC_post = hC*post
    qui gen double hB_D = hB*D
    qui gen double hC_D = hC*D
    qui count if `outcome' == 1
    local ev = r(N)
    qui count
    local n = r(N)
    cap qui logit `outcome' D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D `covs' [pw=kw], vce(cluster clust) iter(100)
    if _rc == 0 & `ev' >= 15 {
        seadj
        return scalar b = _b[pD]
        return scalar se = _se[pD]*r(adj)
        return scalar p = 2*normal(-abs(_b[pD]/(_se[pD]*r(adj))))
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar p = .
    }
    return scalar ev = `ev'
    return scalar n = `n'
    restore
end

* srd: single-period regression discontinuity (post = 1 after, 0 before implementation)
program define srd, rclass
    syntax , h(real) post(integer) [outcome(string) covs(string)]
    if "`outcome'" == "" local outcome death
    if "`covs'" == ""    local covs $COVS
    if "`covs'" == "none" local covs
    preserve
    qui keep if post == `post' & abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui count if `outcome' == 1
    local ev = r(N)
    qui count
    local n = r(N)
    cap qui logit `outcome' D s Ds hB hC `covs' [pw=kw], vce(cluster clust) iter(100)
    if _rc == 0 & `ev' >= 8 {
        seadj
        return scalar b = _b[D]
        return scalar se = _se[D]*r(adj)
        return scalar p = 2*normal(-abs(_b[D]/(_se[D]*r(adj))))
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar p = .
    }
    return scalar ev = `ev'
    return scalar n = `n'
    restore
end

* trio: difference-in-discontinuities and the two single-period discontinuities for one outcome (odds ratios, 95% CI with 1.96)
program define trio, rclass
    syntax , h(real) outcome(string) [covs(string)]
    didisc, h(`h') outcome(`outcome') covs(`covs')
    local dd = exp(r(b))
    local ddlo = exp(r(b) - 1.96*r(se))
    local ddhi = exp(r(b) + 1.96*r(se))
    local ddp = r(p)
    local ev = r(ev)
    local n = r(n)
    srd, h(`h') post(1) outcome(`outcome') covs(`covs')
    local po = exp(r(b))
    local polo = exp(r(b) - 1.96*r(se))
    local pohi = exp(r(b) + 1.96*r(se))
    local pop = r(p)
    local evpo = r(ev)
    local npo = r(n)
    srd, h(`h') post(0) outcome(`outcome') covs(`covs')
    local pr = exp(r(b))
    local prlo = exp(r(b) - 1.96*r(se))
    local prhi = exp(r(b) + 1.96*r(se))
    local prp = r(p)
    local evpr = r(ev)
    local npr = r(n)
    foreach x in dd ddlo ddhi ddp ev n po polo pohi pop evpo npo pr prlo prhi prp evpr npr {
        return scalar `x' = ``x''
    }
end

* rdpp: model-based risk difference at the threshold (percentage points, delta method; est = dd, post or pre)
program define rdpp, rclass
    syntax , h(real) est(string) [outcome(string)]
    if "`outcome'" == "" local outcome death
    preserve
    qui keep if abs(s) < `h'
    if "`est'" == "post" qui keep if post == 1
    if "`est'" == "pre"  qui keep if post == 0
    mkw `h' tri
    qui gen double Ds = D*s
    if "`est'" == "dd" {
        qui gen double pD = post*D
        qui gen double ps = post*s
        qui gen double pDs = post*D*s
        qui gen double hB_post = hB*post
        qui gen double hC_post = hC*post
        qui gen double hB_D = hB*D
        qui gen double hC_D = hC*D
        qui logit `outcome' D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust) iter(100)
        seadj
        local adj = r(adj)
        * post-implementation admissions at the threshold (D = 1, s = 0); remove only the period-by-threshold term
        qui replace hB_D = hB
        qui replace hC_D = hC
        qui margins if post == 1, at(D = 1 s = 0 Ds = 0 post = 1 pD = 1 ps = 0 pDs = 0) ///
            expression(invlogit(predict(xb)) - invlogit(predict(xb) - _b[pD]))
    }
    else {
        qui logit `outcome' D s Ds hB hC $COVS [pw=kw], vce(cluster clust) iter(100)
        seadj
        local adj = r(adj)
        qui margins, at(D = 1 s = 0 Ds = 0) expression(invlogit(predict(xb)) - invlogit(predict(xb) - _b[D]))
    }
    tempname B V
    matrix `B' = r(b)
    matrix `V' = r(V)
    local rd = `B'[1,1]
    local se = sqrt(`V'[1,1])*`adj'
    return scalar rd = 100*`rd'
    return scalar lo = 100*(`rd' - 1.959964*`se')
    return scalar hi = 100*(`rd' + 1.959964*`se')
    return scalar p = 2*normal(-abs(`rd'/`se'))
    restore
end

* ---------- formatting programs (JAMA style as printed) ----------
* fnum: x with dp decimals, minus sign as in the paper
program define fnum, rclass
    args x dp
    if `x' >= . {
        return local s "NA"
        exit
    }
    local s = strtrim(string(`x', "%20.`dp'f"))
    if substr("`s'", 1, 1) == "-" local s = "−" + substr("`s'", 2, .)
    return local s "`s'"
end
* fp: P value (<.001; 3 decimals below .01; otherwise 2 decimals)
program define fp, rclass
    args p
    if `p' >= . {
        return local s "NA"
        exit
    }
    if `p' < 0.001 local s "<.001"
    else {
        if `p' < 0.01 {
            local s = string(`p', "%5.3f")
            if "`s'" == "0.010" local s "0.01"
        }
        else local s = string(`p', "%4.2f")
        local s = substr("`s'", 2, .)
    }
    return local s "`s'"
end
* fn: count, space as thousands separator from 10 000
program define fn, rclass
    args x
    if `x' >= 10000 local s = subinstr(strtrim(string(`x', "%15.0fc")), ",", " ", .)
    else local s = string(`x', "%15.0f")
    return local s "`s'"
end
* orci: OR (lo to hi) or OR (lo-hi); sep = to or dash
program define orci, rclass
    args or lo hi sep
    if `or' >= . {
        return local s "NA"
        exit
    }
    local a = strtrim(string(`or', "%20.2f"))
    local b = strtrim(string(`lo', "%20.2f"))
    local c = strtrim(string(`hi', "%20.2f"))
    if "`sep'" == "to" return local s "`a' (`b' to `c')"
    else return local s "`a' (`b'-`c')"
end
* npc: n (%) of v == 1 among observations meeting cond
program define npc, rclass
    args v cond
    qui count if `cond'
    local n = r(N)
    qui count if `cond' & `v' == 1
    local k = r(N)
    fn `k'
    return local s = "`r(s)' (" + string(100*`k'/`n', "%4.1f") + ")"
end
* miqr: median (IQR) with linear interpolation between order statistics; dp decimals
program define miqr, rclass
    args v cond dp
    preserve
    qui keep if `cond' & `v' < .
    sort `v'
    local n = _N
    foreach q in 25 50 75 {
        local pos = (`n' - 1)*`q'/100 + 1
        local lo = floor(`pos')
        local hi = ceil(`pos')
        local q`q' = `v'[`lo'] + (`pos' - `lo')*(`v'[`hi'] - `v'[`lo'])
    }
    restore
    return local s = string(`q50', "%9.`dp'f") + " (" + string(`q25', "%9.`dp'f") + " to " + string(`q75', "%9.`dp'f") + ")"
end
* msd: mean (SD), 2 decimals
program define msd, rclass
    args v cond
    qui su `v' if `cond'
    return local s = string(r(mean), "%9.2f") + " (" + string(r(sd), "%9.2f") + ")"
end

* ---------- table writer: one row of up to 12 cells per call; tsave writes the sheet ----------
program define topen
    global TBL tbl_handle
    cap postclose $TBL
    postfile $TBL str244(c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12) using "$W/_table.dta", replace
end
program define prow
    local nw : word count `0'
    local args
    forvalues k = 1/12 {
        local c
        if `k' <= `nw' local c : word `k' of `0'
        local args `"`args' (`"`c'"')"'
    }
    post $TBL `args'
end
program define tsave
    args sheet
    postclose $TBL
    preserve
    use "$W/_table.dta", clear
    export excel using "$OUTX", sheet("`sheet'") sheetreplace
    restore
end
* prow_mark: one row of eTables 9 and 12 for indicator y in data file `data' (hospitals ABC, AC or BC; window 5):
* n (%) below and above the threshold in each period, then the three estimates with 95% CIs and P values
* (na: estimates shown as NA, for treatments rarely recorded before admission, as in the paper)
program define prow_mark
    args y hs lab data na
    use "`data'", clear
    if "`hs'" == "AC" qui keep if hB == 0
    if "`hs'" == "BC" qui keep if hB == 1 | hC == 1
    local c
    foreach p in 0 1 {
        foreach sd in 0 1 {
            npc `y' "post == `p' & D == `sd' & abs(s) < 5"
            local c `"`c' "`r(s)'""'
        }
    }
    trio, h(5) outcome(`y')
    foreach x in dd ddlo ddhi ddp po polo pohi pop pr prlo prhi prp {
        local `x' = r(`x')
    }
    foreach k in dd po pr {
        if "`na'" == "na" {
            local c `"`c' "NA" "NA""'
            continue
        }
        orci ``k'' ``k'lo' ``k'hi' to
        local c `"`c' "`r(s)'""'
        fp ``k'p'
        local c `"`c' "`r(s)'""'
    }
    prow "`lab'" `c'
end

* PART 2. 02_Table1
* Table 1. Characteristics of admissions within 5 points of the alert threshold
use "$W/base.dta", clear
merge 1:1 adm_id using "process_flags_v1.dta", keep(master match) keepusing(b_any) nogen
merge 1:1 adm_id using "process_t90_v2.dta", keep(master match) keepusing(vaso_all_Bd12_any o2_escalation_Bd12_any blood_culture_Bd12_any abx_iv_Bd12_any) nogen
keep if abs(s) < 5
gen byte female = 1 - male
gen byte hA = (hB == 0 & hC == 0)
local G1 "post == 0 & D == 0"
local G2 "post == 0 & D == 1"
local G3 "post == 1 & D == 0"
local G4 "post == 1 & D == 1"
topen
prow "Characteristic" "Before implementation, No. (%)a" "After implementation, No. (%)a"
local c
forvalues g = 1/4 {
    qui count if `G`g''
    local lab = cond(inlist(`g', 1, 3), "Score 85-89", "Score 90-94")
    local c `"`c' "`lab' (n = `r(N)')""'
}
prow "" `c'
prow "Characteristics at admission"
local R1 "Age, median (IQR), y|m0|age"
local R2 "Female sex|b|female"
local R3 "Charlson comorbidity index, mean (SD)|sd|cci"
local R4 "Hypertension|b|htn"
local R5 "SOFA score, mean (SD)|sd|sofa"
local R6 "Surgical specialty|b|surgical_dept"
local R7 "Hospital|h|"
local R8 "A|b|hA"
local R9 "B|b|hB"
local R10 "C|b|hC"
local R11 "Highest AI score in the first 24 h, median (IQR)|m1|rv_dec"
local R12 "Any measured treatment in the 24 h before admissionb|b|b_any"
local R13 "Treatment on days 1-2 after admissionc|h|"
local R14 "Vasopressor or inotrope|b|vaso_all_Bd12_any"
local R15 "Oxygen escalation|b|o2_escalation_Bd12_any"
local R16 "Blood culture|b|blood_culture_Bd12_any"
local R17 "Intravenous antibiotic|b|abx_iv_Bd12_any"
forvalues r = 1/17 {
    tokenize "`R`r''", parse("|")
    local lab "`1'"
    local kind "`3'"
    local var "`5'"
    if "`kind'" == "h" {
        prow "`lab'"
        continue
    }
    local c
    forvalues g = 1/4 {
        if "`kind'" == "b"       npc `var' "`G`g''"
        else if "`kind'" == "sd" msd `var' "`G`g''"
        else                     miqr `var' "`G`g''" `=cond("`kind'" == "m0", 0, 1)'
        local c `"`c' "`r(s)'""'
    }
    prow "`lab'" `c'
}
tsave "Table 1"

* PART 3. 03_Table2
* Table 2. In-hospital death and death or ward cardiac arrest at the alert threshold
use "$W/base.dta", clear
topen
prow "Outcome and estimator" "Window, points" "Events/No." "Unadjusted OR (95% CI)" "Adjusted OR (95% CI)a" "P valueb" "Risk difference, percentage points (95% CI)c" "P valueb"
foreach y in death comp_gw {
    if "`y'" == "death" prow "In-hospital death (primary outcome)"
    else                prow "Death or ward cardiac arrest"
    local ests = cond("`y'" == "death", "dd post pre", "dd post")
    foreach e of local ests {
        foreach h in 5 8 {
            * unadjusted (hospital terms only) and adjusted odds ratios
            foreach cv in none adj {
                local co = cond("`cv'" == "none", "none", "")
                if "`e'" == "dd" didisc, h(`h') outcome(`y') covs(`co')
                else             srd, h(`h') post(`=cond("`e'" == "post", 1, 0)') outcome(`y') covs(`co')
                local p_`cv' = r(p)
                local evn "`r(ev)'/`r(n)'"
                orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' to
                local or_`cv' "`r(s)'"
            }
            fp `p_adj'
            local pa "`r(s)'"
            * risk difference at the threshold, percentage points (0.0 when |x| < 0.05)
            rdpp, h(`h') est(`e') outcome(`y')
            local rp = r(p)
            foreach x in rd lo hi {
                local r_`x' = r(`x')
            }
            foreach x in rd lo hi {
                local v = `r_`x''
                if abs(`v') < 0.05 local v = 0
                fnum `v' 1
                local s_`x' "`r(s)'"
            }
            local rdc "`s_rd' (`s_lo' to `s_hi')"
            fp `rp'
            local lab
            if `h' == 5 {
                if "`e'" == "dd"   local lab "Difference-in-discontinuities (primary)"
                if "`e'" == "post" local lab "Postimplementation discontinuity"
                if "`e'" == "pre"  local lab "Preimplementation discontinuityd"
            }
            prow "`lab'" "±`h'" "`evn'" "`or_none'" "`or_adj'" "`pa'" "`rdc'" "`r(s)'"
        }
    }
}
tsave "Table 2"

* PART 4. 04_eTable1
* eTable 1. Admissions, deaths and cardiac arrests by hospital and period; derivation of the analysis sample
use "paper_data_v9.dta", clear
gen byte lm = (rv_dec < . & h_stay > 24)
topen
prow "" "Before implementation" "After implementation"
prow "Hospital" "Implementation date" "Population" "Admissions" "Deaths" "Cardiac arrests" "Admissions" "Deaths" "Cardiac arrests"
forvalues k = 1/3 {
    local hl : word `k' of A B C
    local dt : word `k' of "April 22, 2023" "October 5, 2023" "March 1, 2024"
    foreach pop in all landmark {
        local c
        foreach p in 0 1 {
            local cnd "hosp_num == `k' & treat == `p'"
            if "`pop'" == "landmark" local cnd "`cnd' & lm == 1"
            qui count if `cnd'
            fn `r(N)'
            local a "`r(s)'"
            qui count if `cnd' & death == 1
            fn `r(N)'
            local b "`r(s)'"
            qui count if `cnd' & arr_gw == 1
            local w = r(N)
            qui count if `cnd' & arr_any == 1
            local c `"`c' "`a'" "`b'" "`r(N)' (`w')""'
        }
        if "`pop'" == "all" prow "Hospital `hl'" "`dt'" "All admissions" `c'
        else                prow "" "" "Landmark cohort" `c'
    }
}
prow "Step" "All" "Before implementation" "After implementation"
gen double s = rv_dec - $CUT
local S1 "Admissions in the analysis data (adults on general wards, January 2022 to June 2025)|1"
local S2 "No AI score in the first 24 h (excluded)|rv_dec >= ."
local S3 "AI score in the first 24 h|rv_dec < ."
local S4 "Died or discharged before 24 h (excluded)|rv_dec < . & h_stay <= 24"
local S5 "Landmark cohort (in hospital at 24 h)|lm == 1"
local S6 "In-hospital deaths|lm == 1 & death == 1"
local S7 "Within ±5 points of the threshold (primary window)|lm == 1 & abs(s) < 5"
local S8 "In-hospital deaths|lm == 1 & abs(s) < 5 & death == 1"
local S9 "Within ±8 points of the threshold (secondary window)|lm == 1 & abs(s) < 8"
local S10 "In-hospital deaths|lm == 1 & abs(s) < 8 & death == 1"
forvalues r = 1/10 {
    tokenize "`S`r''", parse("|")
    local c
    foreach p in all 0 1 {
        if "`p'" == "all" qui count if `3'
        else              qui count if `3' & treat == `p'
        fn `r(N)'
        local c `"`c' "`r(s)'""'
    }
    prow "`1'" `c'
}
tsave "eTable 1"

* PART 5. 05_eTable2
* eTable 2. Probability-scale estimates, robust bias-corrected inference, curvature and bias bandwidth (in-hospital death)
* Requires rdrobust (ssc install rdrobust)
use "$W/base.dta", clear
topen
* A and B: conventional and robust bias-corrected estimates, windows 5, 8 and 10
foreach h in 5 8 10 {
    foreach q in 1 0 {
        cap qui rdrobust death rv_dec if post == `q', c($CUT) h(`h') b(`h') p(1) kernel(triangular) covs($COVS) vce(cluster clust)
        local cl`q' = 100*e(tau_cl)
        local clse`q' = 100*e(se_tau_cl)
        local bc`q' = 100*e(tau_bc)
        local rbse`q' = 100*e(se_tau_rb)
    }
    local A`h'
    foreach q in 1 0 {
        fnum `cl`q'' 2
        local e "`r(s)'"
        fnum `clse`q'' 2
        local A`h' `"`A`h'' "`e' (`r(s)')""'
        fp `=2*normal(-abs(`cl`q''/`clse`q''))'
        local A`h' `"`A`h'' "`r(s)'""'
        fnum `bc`q'' 2
        local e "`r(s)'"
        fnum `rbse`q'' 2
        local A`h' `"`A`h'' "`e' (`r(s)')""'
        fp `=2*normal(-abs(`bc`q''/`rbse`q''))'
        local A`h' `"`A`h'' "`r(s)'""'
    }
    local B`h'
    local d = `cl1' - `cl0'
    local dse = sqrt(`clse1'^2 + `clse0'^2)
    fnum `d' 2
    local e "`r(s)'"
    fnum `dse' 2
    local B`h' `"`B`h'' "`e' (`r(s)')""'
    fp `=2*normal(-abs(`d'/`dse'))'
    local B`h' `"`B`h'' "`r(s)'""'
    local d = `bc1' - `bc0'
    local dse = sqrt(`rbse1'^2 + `rbse0'^2)
    fnum `d' 2
    local e "`r(s)'"
    fnum `dse' 2
    local B`h' `"`B`h'' "`e' (`r(s)')""'
    fp `=2*normal(-abs(`d'/`dse'))'
    local B`h' `"`B`h'' "`r(s)'""'
}
prow "A. In-hospital death: discontinuity at the threshold in each period on the probability scale (percentage points), under conventional and robust bias-corrected inference"
prow "" "After implementation, conventional" "After implementation, robust bias-corrected" "Before implementation, conventional" "Before implementation, robust bias-corrected"
prow "Window" "Estimate (SE)" "P value" "Estimate (SE)" "P value" "Estimate (SE)" "P value" "Estimate (SE)" "P value"
foreach h in 5 8 10 {
    prow "±`h'" `A`h''
}
prow "B. In-hospital death: difference-in-discontinuities (after minus before) on the probability scale (percentage points), under conventional and robust bias-corrected inference"
prow "" "Conventional" "Robust bias-corrected"
prow "Window" "Estimate (SE)" "P value" "Estimate (SE)" "P value"
foreach h in 5 8 10 {
    prow "±`h'" `B`h''
}

* C: local quadratic fit on each side of the threshold; relative curvature = |quadratic term| x h^2 / fitted level at the threshold
prow "C. In-hospital death: curvature of mortality in the score, by scale (local quadratic fits on each side of the threshold)"
prow "" "Quadratic term P value" ""
prow "Period" "Window" "Side" "Deaths" "Probability scale" "Logit scale" "Relative curvature"
foreach p in 1 0 {
    local first = 1
    foreach h in 5 8 10 {
        foreach side in 0 1 {
            preserve
            qui keep if post == `p' & abs(s) < `h' & D == `side'
            mkw `h' tri
            qui gen double s2 = s^2
            qui count if death == 1
            local dth = r(N)
            qui regress death s s2 [pw=kw], vce(cluster clust)
            seadj
            local pq1 = 2*normal(-abs(_b[s2]/(_se[s2]*r(adj))))
            local rel = abs(_b[s2])*`h'^2/abs(_b[_cons])
            qui logit death s s2 [pw=kw], vce(cluster clust) iter(100)
            seadj
            local pq2 = 2*normal(-abs(_b[s2]/(_se[s2]*r(adj))))
            restore
            local per
            if `first' local per = cond(`p' == 1, "After implementation", "Before implementation")
            local first = 0
            fp `pq1'
            local a "`r(s)'"
            fp `pq2'
            local b "`r(s)'"
            local sl = cond(`side' == 1, "At or above", "Below")
            fnum `rel' 2
            prow "`per'" "±`h'" "`sl'" "`dth'" "`a'" "`b'" "`r(s)'"
        }
    }
}

* D: bias bandwidth (after implementation)
prow "D. In-hospital death: sensitivity of robust bias-corrected inference to the bias bandwidth (after implementation)"
prow "" "Conventional" "Robust bias-corrected"
prow "Window" "Bias bandwidth" "Estimate, pp" "P value" "Estimate, pp" "Robust SE, pp" "Robust P value"
qui rdbwselect death rv_dec if post == 1, c($CUT) p(1) kernel(triangular) bwselect(mserd) covs($COVS) vce(cluster clust)
local bsel = e(b_mserd)
foreach h in 5 8 {
    foreach m in 1 1.5 2 0 {
        if `m' > 0 {
            local b = `h'*`m'
            local bl = cond(`m' == 1, "b = h", cond(`m' == 1.5, "b = 1.5h", "b = 2h"))
        }
        else {
            local b = `bsel'
            local bl = "b = data-driven (" + string(`bsel', "%4.2f") + ")"
        }
        qui rdrobust death rv_dec if post == 1, c($CUT) h(`h') b(`b') p(1) kernel(triangular) covs($COVS) vce(cluster clust)
        local t1 = 100*e(tau_cl)
        local p1 = e(pv_cl)
        local t2 = 100*e(tau_bc)
        local s2 = 100*e(se_tau_rb)
        local p2 = e(pv_rb)
        fnum `t1' 2
        local c1 "`r(s)'"
        fp `p1'
        local c2 "`r(s)'"
        fnum `t2' 2
        local c3 "`r(s)'"
        fnum `s2' 2
        local c4 "`r(s)'"
        fp `p2'
        prow "±`h'" "`bl'" "`c1'" "`c2'" "`c3'" "`c4'" "`r(s)'"
    }
}
qui rdrobust death rv_dec if post == 1, c($CUT) p(1) kernel(triangular) bwselect(mserd) covs($COVS) vce(cluster clust)
local hl = "±" + string(e(h_l), "%3.1f") + "a"
local bl = "both data-driven (b = " + string(e(b_l), "%4.2f") + ")"
local t1 = 100*e(tau_cl)
local p1 = e(pv_cl)
local t2 = 100*e(tau_bc)
local s2 = 100*e(se_tau_rb)
local p2 = e(pv_rb)
fnum `t1' 2
local c1 "`r(s)'"
fp `p1'
local c2 "`r(s)'"
fnum `t2' 2
local c3 "`r(s)'"
fnum `s2' 2
local c4 "`r(s)'"
fp `p2'
prow "`hl'" "`bl'" "`c1'" "`c2'" "`c3'" "`c4'" "`r(s)'"
tsave "eTable 2"

* PART 6. 06_eTable3
* eTable 3. Events during the admission within 5 points of the alert threshold
use "$W/base.dta", clear
keep if abs(s) < 5
local G1 "post == 0 & D == 0"
local G2 "post == 0 & D == 1"
local G3 "post == 1 & D == 0"
local G4 "post == 1 & D == 1"
topen
prow "Event" "Before implementation, No. (%)" "After implementation, No. (%)"
local c
forvalues g = 1/4 {
    qui count if `G`g''
    local lab = cond(inlist(`g', 1, 3), "Score 85-89", "Score 90-94")
    local c `"`c' "`lab' (n = `r(N)')""'
}
prow "" `c'
local R1 "Sepsisa|sepsis"
local R2 "DNR order|dnr_case"
local R3 "Unplanned transfer to intensive care after 24 h|icu_uit_24h"
local R4 "Ward cardiac arrest|arr_gw"
local R5 "In-hospital death|death"
forvalues r = 1/5 {
    tokenize "`R`r''", parse("|")
    local c
    forvalues g = 1/4 {
        npc `3' "`G`g''"
        local c `"`c' "`r(s)'""'
    }
    prow "`1'" `c'
}
local c
forvalues g = 1/4 {
    miqr total_los "`G`g''" 1
    local c `"`c' "`r(s)'""'
}
prow "Length of stay, median (IQR), d" `c'
tsave "eTable 3"

* PART 7. 07_eTable4
* eTable 4. Density, covariate continuity, covariate discontinuities by period, and entry into the landmark cohort
* Requires rddensity (Cattaneo, Jansson and Ma), version 3.0 from github.com/rdpackages/rddensity
use "$W/base.dta", clear
topen
* A: local polynomial density test (rddensity, default options) in each period; the change between periods is the difference of
* the two density discontinuities (density just above minus just below the threshold) with summed variances
foreach p in 0 1 {
    rddensity s if post == `p', c(0)
    local fl`p' = e(f_ql)
    local fr`p' = e(f_qr)
    local se`p' = e(se_q)
    local P`p' = e(pv_q)
    local N`p' = e(N_h_l) + e(N_h_r)
    local J`p' = ln(e(f_qr)/e(f_ql))
}
local Jd = `J1' - `J0'
local Nd = `N0' + `N1'
local Pd = 2*normal(-abs(((`fr1' - `fl1') - (`fr0' - `fl0'))/sqrt(`se1'^2 + `se0'^2)))
prow "A. Density of the assignment variable in each period and the change between periods"
prow "Quantity" "Log-density jump" "Density ratio (above/below)" "P value" "Admissions"
foreach k in 0 1 d {
    local lab = cond("`k'" == "0", "Before implementation", cond("`k'" == "1", "After implementation", "Difference (after minus before)"))
    local n = `N`k''
    fnum `J`k'' 3
    local a "`r(s)'"
    fnum `=exp(`J`k'')' 2
    local b "`r(s)'"
    fp `P`k''
    prow "`lab'" "`a'" "`b'" "`r(s)'" "`n'"
}

* B: covariate continuity after implementation (binary covariates in percentage points)
prow "B. Continuity of covariates at the threshold, after implementation"
prow "" "±5" "±8"
prow "Covariate" "Discontinuity" "P value" "Discontinuity" "P value"
local L_age_z "Age (standardized)"
local L_male "Male sex"
local L_cci "Charlson comorbidity index"
local L_htn "Hypertension"
local L_sofa "SOFA score at admission"
local L_surgical_dept "Surgical specialty"
local L_summer "Admitted in summer"
local L_fall "Admitted in autumn"
local L_winter "Admitted in winter"
local L_rv_news "National Early Warning Score, first 24 h"
local BIN male htn surgical_dept summer fall winter sepsis hB hC
foreach v in $COVS rv_news {
    local c
    foreach h in 5 8 {
        preserve
        qui keep if post == 1 & abs(s) < `h'
        mkw `h' tri
        qui gen double Ds = D*s
        qui regress `v' D s Ds hB hC [pw=kw], vce(cluster clust)
        seadj
        local m = cond(strpos(" `BIN' ", " `v' "), 100, 1)
        local est = `m'*_b[D]
        local pv = 2*normal(-abs(_b[D]/(_se[D]*r(adj))))
        restore
        fnum `est' 2
        local c `"`c' "`r(s)'""'
        fp `pv'
        local c `"`c' "`r(s)'""'
    }
    prow "`L_`v''" `c'
}

* C: covariate discontinuities before and after implementation and their difference (window 5)
prow "C. Covariate discontinuities at the threshold before and after implementation, and their difference (window ±5)"
prow "" "Before implementation" "After implementation" "Between periods"
prow "Covariate" "Discontinuity" "P value" "Discontinuity" "P value" "Difference" "P value"
local L_age "Age, years"
local L_sepsis "Sepsis"
local L_hB "Hospital B"
local L_hC "Hospital C"
foreach v in age male cci htn sofa surgical_dept sepsis summer fall winter hB hC rv_news {
    preserve
    qui keep if abs(s) < 5
    mkw 5 tri
    qui gen double Ds = D*s
    qui gen double pD = post*D
    qui gen double ps = post*s
    qui gen double pDs = post*D*s
    qui gen double hB_post = hB*post
    qui gen double hC_post = hC*post
    local hosp hB hC hB_post hC_post
    if "`v'" == "hB" | "`v'" == "hC" local hosp
    qui regress `v' D s Ds post pD ps pDs `hosp' [pw=kw], vce(cluster clust)
    seadj
    local adj = r(adj)
    local m = cond(strpos(" `BIN' ", " `v' "), 100, 1)
    local e0 = `m'*_b[D]
    local p0 = 2*normal(-abs(_b[D]/(_se[D]*`adj')))
    local ed = `m'*_b[pD]
    local pd = 2*normal(-abs(_b[pD]/(_se[pD]*`adj')))
    qui lincom D + pD
    local e1 = `m'*r(estimate)
    local p1 = 2*normal(-abs(r(estimate)/(r(se)*`adj')))
    restore
    local c
    foreach k in 0 1 d {
        fnum `e`k'' 2
        local c `"`c' "`r(s)'""'
        fp `p`k''
        local c `"`c' "`r(s)'""'
    }
    prow "`L_`v''" `c'
}

* D: entry into the landmark cohort (all admissions with a score in the first 24 h)
prow "D. Entry into the landmark cohort at the threshold"
prow "" "Difference-in-discontinuities" "After implementation" "Before implementation"
prow "Check" "Window" "Events" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
use "paper_data_v9.dta", clear
keep if rv_dec < .
gen double s    = rv_dec - $CUT
gen byte   D    = (s >= 0)
gen byte   post = treat
gen byte   excl24  = (h_stay <= 24)
gen byte   died24  = (h_stay <= 24 & death == 1)
gen byte   disch24 = (h_stay <= 24 & death == 0)
foreach v in excl24 died24 disch24 {
    local lab = cond("`v'" == "excl24", "Not in cohort at 24 h (died or discharged)", cond("`v'" == "died24", "Died before 24 h", "Discharged alive before 24 h"))
    foreach h in 5 8 {
        trio, h(`h') outcome(`v')
        foreach x in dd ddlo ddhi ddp ev po polo pohi pop pr prlo prhi prp {
            local `x' = r(`x')
        }
        orci `dd' `ddlo' `ddhi' -
        local a "`r(s)'"
        fp `ddp'
        local b "`r(s)'"
        orci `po' `polo' `pohi' -
        local c "`r(s)'"
        fp `pop'
        local d "`r(s)'"
        orci `pr' `prlo' `prhi' -
        local e "`r(s)'"
        fp `prp'
        local rl
        if `h' == 5 local rl "`lab'"
        prow "`rl'" "±`h'" "`ev'" "`a'" "`b'" "`c'" "`d'" "`e'" "`r(s)'"
    }
}
tsave "eTable 4"

* PART 8. 08_eTable5
* eTable 5. Placebo thresholds, negative-control assignment variable (NEWS) and full-range placebo scan
use "$W/base.dta", clear
topen
* A: placebo thresholds on the AI score scale, admissions on one side of the true threshold only (window 5)
prow "A. Placebo thresholds on the AI score scale (difference-in-discontinuities, ±5)"
prow "Placebo threshold (displayed score)" "Deaths" "OR (95% CI)" "P value"
foreach c in 80 83 85 87 92 94 96 90 {
    preserve
    if `c' < 90 qui keep if rv_dec < $CUT
    if `c' > 90 qui keep if rv_dec >= $CUT
    qui replace s = rv_dec - (`c' - 0.5)
    qui replace D = (s >= 0)
    didisc, h(5)
    local ev = r(ev)
    local p = r(p)
    orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
    local o "`r(s)'"
    fp `p'
    prow "`c'" "`ev'" "`o'" "`r(s)'"
    restore
}
* B: National Early Warning Score as the assignment variable
prow "B. Negative-control assignment variable: National Early Warning Score (difference-in-discontinuities)"
prow "NEWS threshold" "Window" "Deaths" "OR (95% CI)" "P value"
foreach c in 5 6 7 {
    foreach h in 3 4 {
        preserve
        qui replace s = rv_news - `c'
        qui replace D = (s >= 0)
        didisc, h(`h')
        local ev = r(ev)
        local p = r(p)
        orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
        local o "`r(s)'"
        fp `p'
        prow "`c'" "`h'" "`ev'" "`o'" "`r(s)'"
        restore
    }
}
* C: every integer threshold from 45 to 99; estimable when both periods have at least 120 admissions, 6 deaths and 15%-85% on each side
tempname SC
postfile `SC' h displayed dd_or dd_p using "$W/scan.dta", replace
foreach hh in 3 5 10 {
    forvalues c = 45/99 {
        preserve
        if `c' < 90 qui keep if rv_dec < $CUT
        if `c' > 90 qui keep if rv_dec >= $CUT
        qui replace s = rv_dec - (`c' - 0.5)
        qui replace D = (s >= 0)
        local ok = 1
        foreach p in 0 1 {
            qui count if post == `p' & abs(s) < `hh'
            local nn = r(N)
            qui count if post == `p' & abs(s) < `hh' & death == 1
            local ee = r(N)
            qui count if post == `p' & abs(s) < `hh' & D == 1
            local sh = cond(`nn' > 0, r(N)/`nn', 0)
            if `nn' < 120 | `ee' < 6 | `sh' < 0.15 | `sh' > 0.85 local ok = 0
        }
        if `ok' {
            didisc, h(`hh')
            if r(b) < . post `SC' (`hh') (`c') (exp(r(b))) (r(p))
        }
        restore
    }
}
postclose `SC'
prow "C. Full-range placebo scan: every integer threshold from 45 to 99"
prow "" "Placebo thresholds" ""
prow "Window" "True threshold, OR" "P value" "Estimable, n" "P < .05, n" "More extreme than true threshold, n" "Randomization P value" "Range of ORs"
preserve
use "$W/scan.dta", clear
gen double absl = abs(ln(dd_or))
foreach hh in 3 5 10 {
    qui su absl if displayed == 90 & h == `hh'
    local tr = r(mean)
    qui su dd_or if displayed == 90 & h == `hh'
    local o90 = r(mean)
    qui su dd_p if displayed == 90 & h == `hh'
    local p90 = r(mean)
    qui count if displayed != 90 & h == `hh'
    local npl = r(N)
    qui count if displayed != 90 & h == `hh' & dd_p < 0.05
    local nsig = r(N)
    qui count if displayed != 90 & h == `hh' & absl >= `tr' - 1e-12
    local nex = r(N)
    qui su dd_or if displayed != 90 & h == `hh'
    local rng = string(r(min), "%4.2f") + "-" + string(r(max), "%4.2f")
    fnum `o90' 2
    local a "`r(s)'"
    fp `p90'
    local b "`r(s)'"
    fp `=(1 + `nex')/(`npl' + 1)'
    prow "±`hh'" "`a'" "`b'" "`npl'" "`nsig'" "`nex'" "`r(s)'" "`rng'"
}
restore
tsave "eTable 5"

* PART 9. 09_eTable6
* eTable 6. Pandemic period and calendar-year stability of the preimplementation discontinuity
use "$W/base.dta", clear
topen
prow "A. Sensitivity to the pandemic period in the preimplementation data"
prow "" "Single RD before implementation" "Difference-in-discontinuities"
prow "Sample" "Window" "Admissions" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
foreach k in 0 6 12 {
    local lab = cond(`k' == 0, "Full preimplementation period", cond(`k' == 6, "Excluding January to June 2022", "Excluding all of 2022"))
    foreach h in 5 8 {
        preserve
        if `k' > 0 qui drop if post == 0 & t <= `k'
        qui count if post == 0
        fn `r(N)'
        local n "`r(s)'"
        srd, h(`h') post(0)
        local p = r(p)
        orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
        local a "`r(s)'"
        fp `p'
        local b "`r(s)'"
        didisc, h(`h')
        local p = r(p)
        orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
        local c "`r(s)'"
        fp `p'
        local rl
        if `h' == 5 local rl "`lab'"
        prow "`rl'" "±`h'" "`n'" "`a'" "`b'" "`c'" "`r(s)'"
        restore
    }
}
prow "B. Stability of the preimplementation discontinuity across calendar years"
prow "" "Preimplementation discontinuity"
prow "Preimplementation period" "Window" "Deaths" "Admissions" "OR (95% CI)" "P value" "P value for equality"
foreach h in 5 8 {
    preserve
    qui keep if post == 0 & abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui gen byte y22 = (t <= 12)
    foreach y in 1 0 {
        qui count if y22 == `y' & death == 1
        local dth`y' = r(N)
        qui count if y22 == `y'
        local nn`y' = r(N)
        qui logit death D s Ds hB hC $COVS [pw=kw] if y22 == `y', vce(cluster clust) iter(100)
        seadj
        local b`y' = _b[D]
        local se`y' = _se[D]*r(adj)
    }
    restore
    local peq = 2*normal(-abs((`b1' - `b0')/sqrt(`se1'^2 + `se0'^2)))
    foreach y in 1 0 {
        orci `=exp(`b`y'')' `=exp(`b`y'' - 1.96*`se`y'')' `=exp(`b`y'' + 1.96*`se`y'')' -
        local a "`r(s)'"
        fp `=2*normal(-abs(`b`y''/`se`y''))'
        local b "`r(s)'"
        if `y' == 1 {
            fp `peq'
            local e "`r(s)'"
        }
        else local e "NA"
        local yl = cond(`y' == 1, "2022", "2023 onwards")
        prow "`yl'" "±`h'" "`dth`y''" "`nn`y''" "`a'" "`b'" "`e'"
    }
}
tsave "eTable 6"

* PART 10. 10_eTable7
* eTable 7. Every window from 1 to 10 points (internal value) with bias-aware inference; displayed integer score as the assignment variable
* Mata functions for the bias-aware (Armstrong-Kolesar) interval; defined here because rdrobust clears Mata
cap mata: mata drop wbias()
cap mata: mata drop hci()
mata:
real scalar wbias(string scalar xvars, string scalar wvar, string scalar svar, real scalar k)
{
    real matrix X, H, Hi
    real colvector w, s, b, p, W, A, u, f, sR, aR, sL, aL
    real scalar i, tot
    X = st_data(., tokens(xvars))
    X = X, J(rows(X), 1, 1)
    w = st_data(., wvar)
    s = st_data(., svar)
    b = st_matrix("e(b)")'
    p = invlogit(X*b)
    W = w :* p :* (1 :- p)
    H = quadcross(X, W, X) + 1e-6*I(cols(X))
    Hi = luinv(H)
    A = ((Hi*(X :* W)')[k, .])'
    tot = 0
    sR = select(s, s :>= 0)
    aR = select(A, s :>= 0)
    if (rows(sR) > 0) {
        u = rangen(0, max(sR), 401)
        f = J(401, 1, 0)
        for (i = 1; i <= 401; i++) f[i] = abs(sum(aR :* ((sR :- u[i]) :* ((sR :- u[i]) :> 0))))
        tot = tot + sum((u[|2 \ 401|] - u[|1 \ 400|]) :* (f[|2 \ 401|] + f[|1 \ 400|]))/2
    }
    sL = select(s, s :< 0)
    aL = select(A, s :< 0)
    if (rows(sL) > 0) {
        u = rangen(min(sL), 0, 401)
        f = J(401, 1, 0)
        for (i = 1; i <= 401; i++) f[i] = abs(sum(aL :* ((u[i] :- sL) :* ((u[i] :- sL) :> 0))))
        tot = tot + sum((u[|2 \ 401|] - u[|1 \ 400|]) :* (f[|2 \ 401|] + f[|1 \ 400|]))/2
    }
    return(tot)
}
void hci(real scalar est, real scalar se, real scalar B)
{
    real scalar t, z, lo, hi, mid, i
    t = B/se
    z = abs(est)/se
    lo = 0
    hi = 500
    for (i = 1; i <= 200; i++) {
        mid = (lo + hi)/2
        if (normal(mid - t) - normal(-mid - t) > 0.95) hi = mid
        else lo = mid
    }
    st_numscalar("r_hp", 1 - normal(z - t) + normal(-z - t))
}
end
* curvM: curvature bound M from the preimplementation data (local quadratic on the logit of displayed-score means)
cap program drop curvM
program define curvM, rclass
    preserve
    qui keep if post == 0 & s > -10 & s < 10
    qui gen double sr = round(s)
    collapse (mean) m = death, by(D sr)
    qui gen double g = min(max(m, 1e-4), 1 - 1e-4)
    qui gen double lg = ln(g/(1 - g))
    qui gen double sr2 = sr^2
    local M = 0
    foreach side in 0 1 {
        qui regress lg sr sr2 if D == `side'
        local M = max(`M', 2*abs(_b[sr2]))
    }
    restore
    return scalar M = `M'
end
* hp: bias-aware P value for the postimplementation discontinuity at window h
cap program drop hp
program define hp, rclass
    syntax , h(real) m(real)
    preserve
    qui keep if post == 1 & abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui logit death D s Ds hB hC $COVS [pw=kw], vce(cluster clust) iter(100)
    seadj
    local b = _b[D]
    local se = _se[D]*r(adj)
    mata: st_numscalar("r_Bu", wbias("D s Ds hB hC $COVS", "kw", "s", 1))
    mata: hci(`b', `se', `=`m'*scalar(r_Bu)')
    return scalar p = scalar(r_hp)
    restore
end

use "$W/base.dta", clear
curvM
local M = r(M)
topen
prow "A. Internal value as the assignment variable (threshold 89.5), every window from 1 to 10 points"
prow "" "Difference-in-discontinuities" "Single RD after implementation" "Single RD before implementation"
prow "Window" "Deaths (admissions)" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "Bias-aware P value" "OR (95% CI)" "P value"
forvalues h = 1/10 {
    trio, h(`h') outcome(death)
    foreach x in dd ddlo ddhi ddp ev n po polo pohi pop pr prlo prhi prp {
        local `x' = r(`x')
    }
    hp, h(`h') m(`M')
    fp `r(p)'
    local hpv "`r(s)'"
    orci `dd' `ddlo' `ddhi' -
    local a "`r(s)'"
    fp `ddp'
    local b "`r(s)'"
    orci `po' `polo' `pohi' -
    local c "`r(s)'"
    fp `pop'
    local d "`r(s)'"
    orci `pr' `prlo' `prhi' -
    local e "`r(s)'"
    fp `prp'
    prow "±`h'" "`ev' (`n')" "`a'" "`b'" "`c'" "`d'" "`hpv'" "`e'" "`r(s)'"
}
prow "B. Displayed integer score as the assignment variable (threshold 90)"
prow "" "Difference-in-discontinuities" "Single RD after implementation" "Single RD before implementation"
prow "Window" "Deaths (admissions)" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR" "P value"
replace s = rv - 90
replace D = (s >= 0)
foreach h in 3 5 8 10 {
    trio, h(`h') outcome(death)
    foreach x in dd ddlo ddhi ddp ev n po polo pohi pop pr prp {
        local `x' = r(`x')
    }
    orci `dd' `ddlo' `ddhi' -
    local a "`r(s)'"
    fp `ddp'
    local b "`r(s)'"
    orci `po' `polo' `pohi' -
    local c "`r(s)'"
    fp `pop'
    local d "`r(s)'"
    fnum `pr' 2
    local e "`r(s)'"
    fp `prp'
    prow "±`h'" "`ev' (`n')" "`a'" "`b'" "`c'" "`d'" "`e'" "`r(s)'"
}
tsave "eTable 7"

* PART 11. 11_eTable8
* eTable 8. Secondary outcomes and sensitivity analyses (window 5)
use "$W/base.dta", clear
gen byte o6 = (death == 1) | (arr_any == 1) | (icu_uit_24h == 1)
foreach hh in B C {
    foreach ss in summer fall winter {
        gen double h`hh'_`ss' = h`hh'*`ss'
    }
}
save "$W/base_e8.dta", replace
topen
prow "" "Difference-in-discontinuities" "Single RD after implementation"
prow "Outcome or analysis" "Events" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
prow "Secondary outcomes" "" "" "" "" ""
local R1 "Death or ward cardiac arrest|comp_gw"
local R2 "Death or any cardiac arrest|comp_any"
local R3 "Any cardiac arrest or death without a DNR order|arr_nondnr"
local R4 "Death, any cardiac arrest or unplanned intensive care transfer|o6"
local R5 "Ward cardiac arrest|arr_gw"
local R6 "Unplanned transfer to intensive care after 24 h|icu_uit_24h"
forvalues r = 1/6 {
    tokenize "`R`r''", parse("|")
    local lab "`1'"
    local y "`3'"
    didisc, h(5) outcome(`y')
    local ev = r(ev)
    local p = r(p)
    orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
    local a "`r(s)'"
    fp `p'
    local b "`r(s)'"
    srd, h(5) post(1) outcome(`y')
    local p = r(p)
    orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' -
    local c "`r(s)'"
    fp `p'
    prow "`lab'" "`ev'" "`a'" "`b'" "`c'" "`r(s)'"
}
prow "Sensitivity analyses" "" "" "" "" ""
* Events column: admissions removed below / above the threshold (NA = none removed)
foreach k in 1 2 3 4 {
    use "$W/base_e8.dta", clear
    local cv "$COVS"
    local rem "NA"
    if `k' == 1 {
        local lab "Excluding any DNR order during stay"
        qui count if abs(s) < 5 & D == 0 & dnr_case == 1
        local r0 = r(N)
        qui count if abs(s) < 5 & D == 1 & dnr_case == 1
        local rem "−`r0' / −`r(N)'"
        qui drop if dnr_case == 1
    }
    if `k' == 2 {
        local lab "Season terms removed"
        local cv "age_z male cci htn sofa surgical_dept"
    }
    if `k' == 3 {
        local lab "Hospital-by-season terms added"
        local cv "$COVS hB_summer hB_fall hB_winter hC_summer hC_fall hC_winter"
    }
    if `k' == 4 {
        local lab "Summer admissions excluded"
        qui count if abs(s) < 5 & D == 0 & summer == 1
        local r0 = r(N)
        qui count if abs(s) < 5 & D == 1 & summer == 1
        local rem "−`r0' / −`r(N)'"
        qui drop if summer == 1
        local cv "age_z male cci htn sofa surgical_dept fall winter"
    }
    didisc, h(5) covs(`cv')
    local p = r(p)
    fnum `=exp(r(b))' 2
    local a "`r(s)'"
    fp `p'
    local b "`r(s)'"
    srd, h(5) post(1) covs(`cv')
    local p = r(p)
    fnum `=exp(r(b))' 2
    local c "`r(s)'"
    fp `p'
    prow "`lab'" "`rem'" "`a'" "`b'" "`c'" "`r(s)'"
}
tsave "eTable 8"

* PART 12. 12_eTable9
* eTable 9. Treatments recorded in the 24 h before admission at the threshold (window 5), and in-hospital death adjusted for them
use "$W/base.dta", clear
merge 1:1 adm_id using "process_flags_v1.dta", keep(master match) nogen
save "$W/base_flags.dta", replace
topen
prow "Outcome (hospitals)" "Before implementation, No. (%)" "After implementation, No. (%)" "DiD-RD" "Post-RD" "Pre-RD"
prow "" "Below" "Above" "Below" "Above" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
prow_mark b_any  ABC "Any measured treatment (A-C)"   "$W/base_flags.dta"
prow_mark b_abx  AC  "Intravenous antibiotic (A, C)"  "$W/base_flags.dta"
prow_mark b_bc   ABC "Blood culture (A-C)"            "$W/base_flags.dta" na
prow_mark b_lac  BC  "Lactate measurement (B, C)"     "$W/base_flags.dta" na
prow_mark b_o2   ABC "Oxygen escalation (A-C)"        "$W/base_flags.dta"
prow_mark b_vaso ABC "Vasopressor or inotrope (A-C)"  "$W/base_flags.dta"
prow "Window" "Adjustment (treatment in the 24 h before admission)" "OR (95% CI)" "P value"
use "$W/base_flags.dta", clear
foreach h in 5 8 {
    foreach a in none b_abx b_bc b_lac b_o2 b_vaso all {
        if "`a'" == "none" local c "$COVS"
        else if "`a'" == "all" local c "$COVS b_abx b_bc b_lac b_o2 b_vaso"
        else local c "$COVS `a'"
        if "`a'" == "none"   local lab "None (primary analysis)"
        if "`a'" == "b_abx"  local lab "+ intravenous antibiotic"
        if "`a'" == "b_bc"   local lab "+ blood culture"
        if "`a'" == "b_lac"  local lab "+ lactate measurement"
        if "`a'" == "b_o2"   local lab "+ oxygen escalation"
        if "`a'" == "b_vaso" local lab "+ vasopressor or inotrope"
        if "`a'" == "all"    local lab "+ all five"
        didisc, h(`h') covs(`c')
        local p = r(p)
        orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' to
        local o "`r(s)'"
        fp `p'
        prow "Window ±`h'" "`lab'" "`o'" "`r(s)'"
    }
}
tsave "eTable 9"

* PART 13. 13_eTable10
* eTable 10. Alternative assignment variables (in-hospital death)
use "$W/base.dta", clear
merge 1:1 adm_id using "alt_assign_v1.dta", keep(match) nogen
save "$W/base_alt.dta", replace
topen
prow "Assignment variable and window" "Admissions, No." "Deaths, No." "DiD-RD" "Post-RD" "Pre-RD"
prow "" "" "" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
foreach k in 1 2 3 4 {
    use "$W/base_alt.dta", clear
    local h = cond(`k' == 3, 8, 5)
    if `k' == 1 local lab "Highest value in the first 24 h (primary), ±5"
    if inlist(`k', 2, 3) {
        local lab "First value at or above the threshold within 24 h, ±`h'"
        qui replace s = rv_first24 - $CUT
    }
    if `k' == 4 {
        local lab "Highest value over the whole stay, ±5"
        qui replace s = rv_whole - $CUT
    }
    qui replace D = (s >= 0)
    trio, h(`h') outcome(death)
    foreach x in dd ddlo ddhi ddp ev n po polo pohi pop pr prlo prhi prp {
        local `x' = r(`x')
    }
    fn `n'
    local c `""`r(s)'" "`ev'""'
    foreach e in dd po pr {
        orci ``e'' ``e'lo' ``e'hi' to
        local c `"`c' "`r(s)'""'
        fp ``e'p'
        local c `"`c' "`r(s)'""'
    }
    prow "`lab'" `c'
}
tsave "eTable 10"

* PART 14. 14_eTable11
* eTable 11. Charting frequency between 24 and 72 h and discharge alive within 3 days
use "$W/base.dta", clear
gen double hrs = min(max(h_stay - 24, 0), 48)
gen double rate = n_24_72/hrs if hrs > 0
topen
prow "Outcome" "Window" "Period" "Below threshold" "At or above" "Discontinuity" "P value"
foreach y in rate disc3 {
    local m = cond("`y'" == "disc3", 100, 1)
    local first = 1
    foreach h in 5 8 {
        foreach p in 0 1 2 {
            preserve
            qui keep if abs(s) < `h' & `y' < .
            if `p' < 2 qui keep if post == `p'
            mkw `h' tri
            qui gen double Ds = D*s
            if `p' < 2 {
                qui su `y' if D == 0
                local mb = `m'*r(mean)
                qui su `y' if D == 1
                local ma = `m'*r(mean)
                qui regress `y' D s Ds hB hC $COVS [pw=kw], vce(cluster clust)
                seadj
                local est = `m'*_b[D]
                local pv = 2*normal(-abs(_b[D]/(_se[D]*r(adj))))
            }
            else {
                qui gen double pD = post*D
                qui gen double ps = post*s
                qui gen double pDs = post*D*s
                qui gen double hB_post = hB*post
                qui gen double hC_post = hC*post
                qui gen double hB_D = hB*D
                qui gen double hC_D = hC*D
                qui regress `y' D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust)
                seadj
                local est = `m'*_b[pD]
                local pv = 2*normal(-abs(_b[pD]/(_se[pD]*r(adj))))
                local mb = .
                local ma = .
            }
            restore
            local lab
            if `first' local lab = cond("`y'" == "rate", "Vital-sign charting, entries per hour", "Discharged alive within 3 days, %")
            local first = 0
            local per = cond(`p' == 0, "Before implementation", cond(`p' == 1, "After implementation", "Difference-in-discontinuities"))
            fnum `mb' 2
            local a "`r(s)'"
            fnum `ma' 2
            local b "`r(s)'"
            fnum `est' 2
            local c "`r(s)'"
            fp `pv'
            prow "`lab'" "±`h'" "`per'" "`a'" "`b'" "`c'" "`r(s)'"
        }
    }
}
tsave "eTable 11"

* PART 15. 15_eTable12
* eTable 12. Treatment within 72 h of admission at the threshold (part 1) and vasopressor or inotrope at hospitals B and C by day and by time after the alert (part 2); window 5
* Part 1 uses $W/base_flags.dta (made in 12_eTable9.do); part 2 the markers in process_t90_v2.dta
use "$W/base.dta", clear
merge 1:1 adm_id using "process_t90_v2.dta", keep(match) nogen
save "$W/base_pm.dta", replace
topen
prow "Outcome (hospitals)" "Before implementation, No. (%)" "After implementation, No. (%)" "DiD-RD" "Post-RD" "Pre-RD"
prow "" "Below" "Above" "Below" "Above" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
prow_mark a_any  ABC "Any measured treatment (A-C)"   "$W/base_flags.dta"
prow_mark a_abx  AC  "Intravenous antibiotic (A, C)"  "$W/base_flags.dta"
prow_mark a_bc   ABC "Blood culture (A-C)"            "$W/base_flags.dta"
prow_mark a_lac  BC  "Lactate measurement (B, C)"     "$W/base_flags.dta"
prow_mark a_o2   ABC "Oxygen escalation (A-C)"        "$W/base_flags.dta"
prow_mark a_vaso ABC "Vasopressor or inotrope (A-C)"  "$W/base_flags.dta"
prow "Outcome (hospitals)" "Before implementation, No. (%)" "After implementation, No. (%)" "DiD-RD" "Post-RD" "Pre-RD"
prow "" "Below" "Above" "Below" "Above" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR (95% CI)" "P value"
prow_mark vaso_all_Bd0_any  BC "Day 0 (admission day)"   "$W/base_pm.dta"
prow_mark vaso_all_Bd2_any  BC "Day 2"                   "$W/base_pm.dta"
prow_mark vaso_all_Bd3_any  BC "Day 3"                   "$W/base_pm.dta"
prow_mark vaso_all_A24h_any BC "0-24 h after the alert"  "$W/base_pm.dta"
prow_mark vaso_all_A48h_any BC "0-48 h after the alert"  "$W/base_pm.dta"
prow_mark vaso_all_A72h_any BC "0-72 h after the alert"  "$W/base_pm.dta"
tsave "eTable 12"

* PART 16. 16_eTable13
* eTable 13. Process-of-care markers at the alert threshold: specifications (part 1), window 8, linear probability models and
* placebo thresholds (part 2), hospital-level estimates (part 3). Also saves the placebo-threshold estimates used by eFigure 6.
* lpm: the primary design as a linear probability model (percentage points); est = dd or post
cap program drop lpm
program define lpm, rclass
    syntax , h(real) outcome(string) est(string)
    preserve
    qui keep if abs(s) < `h'
    if "`est'" == "post" qui keep if post == 1
    mkw `h' tri
    qui gen double Ds = D*s
    if "`est'" == "dd" {
        qui gen double pD = post*D
        qui gen double ps = post*s
        qui gen double pDs = post*D*s
        qui gen double hB_post = hB*post
        qui gen double hC_post = hC*post
        qui gen double hB_D = hB*D
        qui gen double hC_D = hC*D
        qui regress `outcome' D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust)
        local b = _b[pD]
        local se = _se[pD]
    }
    else {
        qui regress `outcome' D s Ds hB hC $COVS [pw=kw], vce(cluster clust)
        local b = _b[D]
        local se = _se[D]
    }
    seadj
    local se = `se'*r(adj)
    return scalar pp = 100*`b'
    return scalar lo = 100*(`b' - 1.96*`se')
    return scalar hi = 100*(`b' + 1.96*`se')
    return scalar p = 2*normal(-abs(`b'/`se'))
    restore
end
local MARK "vital_any_A1h_any vital_any_A2h_any vital_any_A3h_any vaso_all_Bd1_any vaso_all_Bd12_any o2_escalation_A1h_new o2_escalation_A24h_new o2_escalation_A72h_new blood_culture_Bd01_any"
local L_vital_any_A1h_any "VS ≤1 h"
local L_vital_any_A2h_any "VS ≤2 h (Figure 3)"
local L_vital_any_A3h_any "VS ≤3 h"
local L_vaso_all_Bd1_any "VI, day 1"
local L_vaso_all_Bd12_any "VI, days 1-2 (Figure 3)"
local L_o2_escalation_A1h_new "O2 ≤1 h, new"
local L_o2_escalation_A24h_new "O2 ≤24 h, new"
local L_o2_escalation_A72h_new "O2 ≤72 h, new"
local L_blood_culture_Bd01_any "BC, days 0-1"

* placebo thresholds (displayed scores 80, 78, ..., 62): markers re-timed at each threshold; specifications with at least 15 events and
* 15 non-events within 5 points (all hospitals); the two Figure 3 markers also at hospitals B and C (eFigure 6)
tempname PP PB
postfile `PP' cut str32 spec dd ddlo ddhi ddp po polo pohi pop using "$W/placebo_markers.dta", replace
postfile `PB' cut str32 spec dd ddlo ddhi ddp po polo pohi pop using "$W/placebo_markers_BC.dta", replace
foreach c in 80 78 76 74 72 70 68 66 64 62 {
    use "$W/base.dta", clear
    qui replace s = rv_dec - (`c' - 0.5)
    qui replace D = (s >= 0)
    merge 1:1 adm_id using "process_t`c'_v1.dta", keep(match) nogen
    qui keep if abs(s) < 5
    save "$W/base_plac.dta", replace
    foreach y of local MARK {
        qui count if `y' == 1
        local ev = r(N)
        qui count
        if `ev' >= 15 & r(N) - `ev' >= 15 {
            didisc, h(5) outcome(`y')
            local dd = exp(r(b))
            local ddlo = exp(r(b) - 1.96*r(se))
            local ddhi = exp(r(b) + 1.96*r(se))
            local ddp = r(p)
            srd, h(5) post(1) outcome(`y')
            post `PP' (`c') ("`y'") (`dd') (`ddlo') (`ddhi') (`ddp') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
        }
    }
    qui keep if hB == 1 | hC == 1
    foreach y in vital_any_A2h_any vaso_all_Bd12_any {
        didisc, h(5) outcome(`y')
        local dd = exp(r(b))
        local ddlo = exp(r(b) - 1.96*r(se))
        local ddhi = exp(r(b) + 1.96*r(se))
        local ddp = r(p)
        srd, h(5) post(1) outcome(`y')
        post `PB' (`c') ("`y'") (`dd') (`ddlo') (`ddhi') (`ddp') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    }
}
postclose `PP'
postclose `PB'
use "$W/placebo_markers.dta", clear
gen byte hit = (dd > 1 & ddp < .05)
gen byte one = 1
collapse (sum) n_eval = one n_hit = hit, by(spec)
foreach i of numlist 1/`=_N' {
    local PL_`=spec[`i']' = string(n_hit[`i']) + "/" + string(n_eval[`i'])
}

use "process_t90_v2.dta", clear
merge 1:1 adm_id using "$W/base.dta", keep(match) nogen
save "$W/base_pm.dta", replace
topen
* part 1: window 5, hospitals A-C and B-C
prow "Marker (hospitals)" "Events/No." "DiD-RD" "Post-RD" "Pre-RD"
prow "" "" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR" "P value"
foreach y of local MARK {
    foreach hs in ABC BC {
        use "$W/base_pm.dta", clear
        if "`hs'" == "BC" qui keep if hB == 1 | hC == 1
        trio, h(5) outcome(`y')
        foreach x in dd ddlo ddhi ddp ev n po polo pohi pop pr prp {
            local `x' = r(`x')
        }
        orci `dd' `ddlo' `ddhi' to
        local a "`r(s)'"
        fp `ddp'
        local b "`r(s)'"
        orci `po' `polo' `pohi' to
        local c "`r(s)'"
        fp `pop'
        local d "`r(s)'"
        fnum `pr' 2
        local e "`r(s)'"
        fp `prp'
        local hl = cond("`hs'" == "ABC", "(A-C)", "(B, C)")
        prow "`L_`y'' `hl'" "`ev'/`n'" "`a'" "`b'" "`c'" "`d'" "`e'" "`r(s)'"
    }
}
* part 2: window 8, linear probability models (window 5) and placebo thresholds; hospitals A-C
prow "Marker" "DiD-RD, ±8" "Post-RD, ±8" "DiD-RD, linear, percentage points" "Post-RD, linear, percentage points" "Placebo thresholds"
prow "" "OR (95% CI)" "P value" "OR" "P value" "Difference (95% CI)" "P value" "Difference" "P value" "Significant, No./10"
foreach y of local MARK {
    use "$W/base_pm.dta", clear
    trio, h(8) outcome(`y')
    foreach x in dd ddlo ddhi ddp po pop {
        local `x' = r(`x')
    }
    lpm, h(5) outcome(`y') est(dd)
    foreach x in pp lo hi p {
        local l_`x' = r(`x')
    }
    lpm, h(5) outcome(`y') est(post)
    foreach x in pp p {
        local m_`x' = r(`x')
    }
    orci `dd' `ddlo' `ddhi' to
    local a "`r(s)'"
    fp `ddp'
    local b "`r(s)'"
    fnum `po' 2
    local c "`r(s)'"
    fp `pop'
    local d "`r(s)'"
    fnum `l_pp' 1
    local e1 "`r(s)'"
    fnum `l_lo' 1
    local e2 "`r(s)'"
    fnum `l_hi' 1
    local e "`e1' (`e2' to `r(s)')"
    fp `l_p'
    local f "`r(s)'"
    fnum `m_pp' 1
    local g "`r(s)'"
    fp `m_p'
    prow "`L_`y''" "`a'" "`b'" "`c'" "`d'" "`e'" "`f'" "`g'" "`r(s)'" "`PL_`y''"
}
* part 3: hospital-level difference-in-discontinuities (window 5)
prow "Marker" "DiD-RD, OR (95% CI); P value"
prow "" "Hospital A" "Hospital B" "Hospital C"
foreach y in vital_any_A1h_any vital_any_A2h_any vaso_all_Bd1_any vaso_all_Bd12_any o2_escalation_A1h_new {
    local c
    forvalues k = 1/3 {
        use "$W/base_pm.dta", clear
        if `k' == 1 qui keep if hB == 0 & hC == 0
        if `k' == 2 qui keep if hB == 1
        if `k' == 3 qui keep if hC == 1
        didisc, h(5) outcome(`y')
        local p = r(p)
        orci `=exp(r(b))' `=exp(r(b) - 1.96*r(se))' `=exp(r(b) + 1.96*r(se))' to
        local o "`r(s)'"
        fp `p'
        local c `"`c' "`o'; P = `r(s)'""'
    }
    prow "`L_`y''" `c'
}
tsave "eTable 13"

* PART 17. 17_eTable14
* eTable 14. Reference-time rules R0-R9 for patients below the threshold (window 5, hospitals A-C)
topen
prow "Marker" "Rule for unalerted patients" "Median reference time, h" "DiD-RD" "Post-RD" "Pre-RD"
prow "" "" "" "OR (95% CI)" "P value" "OR (95% CI)" "P value" "OR" "P value"
local RL0 "R0 (primary): 24-h maximum"
local RL1 "R1: first value ≥84.5"
local RL2 "R2: first value ≥79.5"
local RL3 "R3: within 1 point of own maximum"
local RL4 "R4: within 2.5 points of own maximum"
local RL5 "R5: within 5 points of own maximum"
local RL6 "R6: 24-h maximum, both sides"
local RL7 "R7: within 1 point, both sides"
local RL8 "R8: matched above-threshold patient"
local RL9 "R9: hospital-period median"
foreach y in vital_any_A2h_any o2_escalation_A1h_new {
    forvalues r = 0/9 {
        use "rules_t90_v1.dta", clear
        qui keep if rule == `r'
        drop rule
        merge 1:1 adm_id using "$W/base.dta", keep(match) nogen
        qui summarize t0hr if s < 0 & abs(s) < 5, detail
        local med = string(r(p50), "%4.1f")
        trio, h(5) outcome(`y')
        foreach x in dd ddlo ddhi ddp po polo pohi pop pr prp {
            local `x' = r(`x')
        }
        orci `dd' `ddlo' `ddhi' to
        local a "`r(s)'"
        fp `ddp'
        local b "`r(s)'"
        orci `po' `polo' `pohi' to
        local c "`r(s)'"
        fp `pop'
        local d "`r(s)'"
        fnum `pr' 2
        local e "`r(s)'"
        fp `prp'
        local ml
        if `r' == 0 local ml = cond("`y'" == "vital_any_A2h_any", "VS ≤2 h (Figure 3)", "O2 ≤1 h, new")
        prow "`ml'" "`RL`r''" "`med'" "`a'" "`b'" "`c'" "`d'" "`e'" "`r(s)'"
    }
}
tsave "eTable 14"

