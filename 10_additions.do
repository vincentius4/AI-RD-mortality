* JAMA additions: absolute risk differences with and without covariates (Table 2), process of care (eTable 10),
* treatments before admission and mortality adjusted for them (eTable 8)
* Requires 01_data_and_programs.do (programs and replicate_v12/base.dta) and process_flags_v1.dta
version 15
set more off
global OUTJ "replicate_v12\jama_additions.xlsx"
cap log close
log using "replicate_v12\jama_additions.log", replace text

use "replicate_v12/base.dta", clear
merge 1:1 adm_id using "process_flags_v1.dta", keep(master match) nogen
gen byte H = cond(hB == 1, 2, cond(hC == 1, 3, 1))
label define HL 1 "A" 2 "B" 3 "C"
label values H HL
save "replicate_v12/base_jama.dta", replace

* Model-based risk difference at the threshold (percentage points), delta method
* est = dd (difference-in-discontinuities) or post / pre (single-period discontinuity)
capture program drop rdpp
program define rdpp, rclass
    syntax , h(real) est(string) [outcome(string) NOCOVs]
    if "`outcome'" == "" local outcome death
    local cv "$COVS"
    if "`nocovs'" != "" local cv ""
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
        qui logit `outcome' D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D `cv' [pw=kw], vce(cluster clust) iter(100)
        seadj
        local adj = r(adj)
        local or = exp(_b[pD])
        local orp = 2*normal(-abs(_b[pD]/(_se[pD]*`adj')))
        * counterfactual at the threshold for post-implementation admissions: D = 1, s = 0; remove only the period-by-threshold term
        qui replace hB_D = hB
        qui replace hC_D = hC
        qui margins if post == 1, at(D = 1 s = 0 Ds = 0 post = 1 pD = 1 ps = 0 pDs = 0) ///
            expression(invlogit(predict(xb)) - invlogit(predict(xb) - _b[pD]))
    }
    else {
        qui logit `outcome' D s Ds hB hC `cv' [pw=kw], vce(cluster clust) iter(100)
        seadj
        local adj = r(adj)
        local or = exp(_b[D])
        local orp = 2*normal(-abs(_b[D]/(_se[D]*`adj')))
        qui margins, at(D = 1 s = 0 Ds = 0) expression(invlogit(predict(xb)) - invlogit(predict(xb) - _b[D]))
    }
    tempname B V
    matrix `B' = r(b)
    matrix `V' = r(V)
    local rd = `B'[1,1]
    local se = sqrt(`V'[1,1])*`adj'
    return scalar or = `or'
    return scalar orp = `orp'
    return scalar rd = 100*`rd'
    return scalar lo = 100*(`rd' - 1.959964*`se')
    return scalar hi = 100*(`rd' + 1.959964*`se')
    return scalar p = 2*normal(-abs(`rd'/`se'))
    restore
end

tempname PF
tempfile R1
postfile `PF' str32 outcome h str8 est or orp rd lo hi p using `R1', replace
foreach y in death comp_gw {
    foreach h in 3 5 8 10 {
        foreach e in dd post pre {
            rdpp, h(`h') est(`e') outcome(`y')
            post `PF' ("`y'") (`h') ("`e'") (r(or)) (r(orp)) (r(rd)) (r(lo)) (r(hi)) (r(p))
        }
    }
}
postclose `PF'
use `R1', clear
export excel using "$OUTJ", sheet("RD_pp") sheetreplace firstrow(variables)
list, clean noobs

* Table 2, rows without covariates (hospital terms kept)
use "replicate_v12/base_jama.dta", clear
tempname PU
tempfile RU
postfile `PU' str32 outcome h str8 est or orp rd lo hi p using `RU', replace
foreach y in death comp_gw {
    foreach h in 5 8 {
        foreach e in dd post pre {
            rdpp, h(`h') est(`e') outcome(`y') nocovs
            post `PU' ("`y'") (`h') ("`e'") (r(or)) (r(orp)) (r(rd)) (r(lo)) (r(hi)) (r(p))
        }
    }
}
postclose `PU'
use `RU', clear
export excel using "$OUTJ", sheet("RD_pp_nocovs") sheetreplace firstrow(variables)
list, clean noobs

* Process of care and treatments before admission: DiD-RD and single-period RDs (odds ratios), window 5
tempname PG
tempfile R2
postfile `PG' str24 outcome str4 hosp h dd ddp post postp pre prep using `R2', replace
local specs "a_any ABC 5 a_abx AC 5 a_bc ABC 5 a_lac BC 5 a_o2 ABC 5 a_vaso ABC 5 newpress BC 5 newpress BC 8 newpress ABC 5 newpress ABC 8 b_any ABC 5 b_abx AC 5 b_o2 ABC 5 b_vaso ABC 5"
local n : word count `specs'
forvalues i = 1(3)`n' {
    local y : word `i' of `specs'
    local hs : word `=`i'+1' of `specs'
    local h : word `=`i'+2' of `specs'
    use "replicate_v12/base_jama.dta", clear
    if "`hs'" == "AC" keep if inlist(H, 1, 3)
    if "`hs'" == "BC" keep if inlist(H, 2, 3)
    didisc, h(`h') outcome(`y')
    local dd = exp(r(b))
    local ddp = r(p)
    srd, h(`h') post(1) outcome(`y')
    local po = exp(r(b))
    local pop = r(p)
    srd, h(`h') post(0) outcome(`y')
    local pr = exp(r(b))
    local prp = r(p)
    post `PG' ("`y'") ("`hs'") (`h') (`dd') (`ddp') (`po') (`pop') (`pr') (`prp')
}
postclose `PG'
use `R2', clear
export excel using "$OUTJ", sheet("process") sheetreplace firstrow(variables)
list, clean noobs

* In-hospital death adjusted for treatments recorded before admission (primary estimator)
tempname PH
tempfile R3
postfile `PH' h str16 adj or p using `R3', replace
foreach h in 5 8 {
    foreach a in none b_abx b_bc b_lac b_o2 b_vaso all {
        use "replicate_v12/base_jama.dta", clear
        if "`a'" == "none" local c "$COVS"
        else if "`a'" == "all" local c "$COVS b_abx b_bc b_lac b_o2 b_vaso"
        else local c "$COVS `a'"
        didisc, h(`h') covs(`c')
        post `PH' (`h') ("`a'") (exp(r(b))) (r(p))
    }
}
postclose `PH'
use `R3', clear
export excel using "$OUTJ", sheet("adj_mort") sheetreplace firstrow(variables)
list, clean noobs
log close
