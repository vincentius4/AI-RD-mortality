* Process-of-care markers at the alert threshold: Figure 3, eTable 10 (alternative assignment variables), eTable 12 panel B,
* eTable 13, eTable 14, eFigure 6 (JAMA manuscript; eFigure 6 panels C and D: hospitals B and C)
* Requires 01_data_and_programs.do (programs mkw, seadj, didisc, srd; replicate_v12/base.dta) and the marker data files built by
* data_preparation/build_process_markers.py (one row per admission, indicators 0/1):
*   process_t90_v1.dta   admissions within 8 points of the alert threshold; t0hr = reference time in hours after admission
*   process_t80_v1.dta ... process_t62_v1.dta   the same markers at the placebo thresholds (displayed scores 80, 78, ..., 62)
*   rules_t90_v1.dta     markers recomputed under the reference-time rules R0-R9 (admissions within 5 points)
*   alt_assign_v1.dta    alternative assignment variables (first value at or above the threshold within 24 h; whole stay)
* Marker variable names: measure_timing_definition (e.g. vital_any_A2h_any = vital-sign measurement within 2 h of the reference time;
*   vaso_all_Bd12_any = vasopressor or inotrope on calendar days 1-2; o2_escalation_A1h_new = new oxygen escalation within 1 h).
* Confidence intervals use 1.96 as in the Python code that produced the paper's numbers.
version 15
set more off
global OUTP "replicate_v12\process_markers.xlsx"
cap log close
log using "replicate_v12\process_markers.log", replace text

* ---------- programs ----------
* trio: difference-in-discontinuities and single-period discontinuities (odds ratios) for one marker, window h
capture program drop trio
program define trio, rclass
    syntax , h(real) outcome(string)
    didisc, h(`h') outcome(`outcome')
    local dd = exp(r(b))
    local ddlo = exp(r(b) - 1.96*r(se))
    local ddhi = exp(r(b) + 1.96*r(se))
    local ddp = r(p)
    local ev = r(ev)
    local n = r(n)
    srd, h(`h') post(1) outcome(`outcome')
    local po = exp(r(b))
    local polo = exp(r(b) - 1.96*r(se))
    local pohi = exp(r(b) + 1.96*r(se))
    local pop = r(p)
    srd, h(`h') post(0) outcome(`outcome')
    local pr = exp(r(b))
    local prlo = exp(r(b) - 1.96*r(se))
    local prhi = exp(r(b) + 1.96*r(se))
    local prp = r(p)
    foreach x in dd ddlo ddhi ddp ev n po polo pohi pop pr prlo prhi prp {
        return scalar `x' = ``x''
    }
end

* lpm: the same design as the primary model as a linear probability model (percentage points), est = dd or post
capture program drop lpm
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

* locfit: unadjusted local logistic fit on one side of the threshold (triangular kernel, window 5, hospital-month clusters);
* returns the fitted proportion at the threshold and posts the fitted curve (80 points) to the open postfile `curve'
capture program drop locfit
program define locfit, rclass
    syntax , outcome(string) post(integer) side(integer) curve(string) panel(string)
    preserve
    qui keep if post == `post' & abs(s) < 5 & D == `side'
    mkw 5 tri
    qui logit `outcome' s [pw=kw], vce(cluster clust) iter(100)
    seadj
    local adj = r(adj)
    tempname V
    matrix `V' = e(V)*(`adj')^2
    local b1 = _b[s]
    local b0 = _b[_cons]
    local vss = `V'[1,1]
    local vsc = `V'[1,2]
    local vcc = `V'[2,2]
    * at the threshold (s = 0 above; s -> 0 from below, evaluated at -0.001 as in the Python curve)
    local x0 = cond(`side' == 1, 0, -0.001)
    local eta = `b0' + `b1'*`x0'
    return scalar p0 = invlogit(`eta')
    forvalues j = 0/79 {
        if `side' == 0 local x = -5 + (5 - 0.001)*`j'/79
        else local x = 0 + (5 - 0.001)*`j'/79
        local eta = `b0' + `b1'*`x'
        local se = sqrt(`x'^2*`vss' + 2*`x'*`vsc' + `vcc')
        post `curve' ("`panel'") (`post') (`side') (`x' + $CUT + 0.5) (invlogit(`eta')) (invlogit(`eta' - 1.96*`se')) (invlogit(`eta' + 1.96*`se'))
    }
    restore
end

* ---------- data ----------
use "replicate_v12/base.dta", clear
merge 1:1 adm_id using "process_t90_v1.dta", keep(match) nogen
gen byte H = cond(hB == 1, 2, cond(hC == 1, 3, 1))
cap label drop HL
label define HL 1 "A" 2 "B" 3 "C"
label values H HL
count
save "replicate_v12/base_pm.dta", replace

* ---------- Figure 3: vital-sign measurement within 2 h, vasopressor or inotrope on days 1-2 (hospitals A-C) ----------
tempname CV
tempfile CURVE
postfile `CV' str40 panel post side x p lo hi using `CURVE', replace
tempname F3
tempfile RF3
postfile `F3' str40 panel post n below above jump using `RF3', replace
foreach y in vital_any_A2h_any vaso_all_Bd12_any {
    foreach t in 0 1 {
        use "replicate_v12/base_pm.dta", clear
        qui count if post == `t' & abs(s) < 5
        local n = r(N)
        locfit, outcome(`y') post(`t') side(0) curve(`CV') panel(`y')
        local pb = 100*r(p0)
        locfit, outcome(`y') post(`t') side(1) curve(`CV') panel(`y')
        local pa = 100*r(p0)
        post `F3' ("`y'") (`t') (`n') (`pb') (`pa') (`pa' - `pb')
    }
}
postclose `F3'
postclose `CV'
* observed proportions in 0.5-point bins
use "replicate_v12/base_pm.dta", clear
keep if abs(s) < 5
gen double bin = floor((rv_dec + 0.5)/0.5)*0.5 + 0.25
collapse (mean) vital_any_A2h_any vaso_all_Bd12_any (count) n = s, by(post bin)
save "replicate_v12/F3_bins.dta", replace
export excel using "$OUTP", sheet("F3_bins") sheetreplace firstrow(variables)
use `RF3', clear
list, clean noobs
export excel using "$OUTP", sheet("F3_fit") sheetreplace firstrow(variables)
use `CURVE', clear
save "replicate_v12/F3_curves.dta", replace
export excel using "$OUTP", sheet("F3_curves") sheetreplace firstrow(variables)
* figure (2 x 2): A/B vital-sign measurement, C/D vasopressor or inotrope; before (gray) and after (green) implementation
local k = 0
foreach y in vital_any_A2h_any vaso_all_Bd12_any {
    local ymax = cond("`y'" == "vital_any_A2h_any", 100, 25)
    local ylab = cond("`y'" == "vital_any_A2h_any", "Vital-sign measurement within 2 h, %", "Vasopressor or inotrope on days 1-2, %")
    foreach t in 0 1 {
        local k = `k' + 1
        local L : word `k' of A B C D
        local per = cond(`t' == 0, "Before implementation", "After implementation")
        local col = cond(`t' == 0, "gs9", "dkgreen")
        local fil = cond(`t' == 0, "gs13", "eltgreen")
        use "replicate_v12/F3_curves.dta", clear
        keep if panel == "`y'" & post == `t'
        foreach v in p lo hi {
            qui replace `v' = 100*`v'
        }
        append using "replicate_v12/F3_bins.dta"
        qui gen double pb = 100*`y' if bin < . & post == `t'
        twoway (rarea lo hi x if side == 0, color(`fil') lwidth(none)) (rarea lo hi x if side == 1, color(`fil') lwidth(none)) ///
            (line p x if side == 0, lcolor(`col') lwidth(medthick)) (line p x if side == 1, lcolor(`col') lwidth(medthick)) ///
            (scatter pb bin if post == `t', mcolor(`col') msize(small)), ///
            xline(90, lpattern(dash) lcolor(gs7)) ylabel(0(`=`ymax'/5')`ymax', angle(0)) xlabel(85(1)95) ///
            ytitle("`ylab'") xtitle("Highest AI score in first 24 h") title("`L'  `per'", position(11) size(medsmall)) legend(off) ///
            name(F3`L', replace) nodraw
    }
}
graph combine F3A F3B F3C F3D, cols(2) ysize(6) xsize(7)
graph export "replicate_v12/F3_process_stata.png", replace width(2100)

* ---------- eTable 10: alternative assignment variables (in-hospital death, primary model) ----------
tempname PA
tempfile RA
postfile `PA' str60 assignment h n deaths dd ddlo ddhi ddp po polo pohi pop pr prlo prhi prp using `RA', replace
use "replicate_v12/base.dta", clear
merge 1:1 adm_id using "alt_assign_v1.dta", keep(match) nogen
count
save "replicate_v12/base_alt.dta", replace
trio, h(5) outcome(death)
post `PA' ("Highest value in the first 24 h (primary)") (5) (r(n)) (r(ev)) (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
foreach h in 5 8 {
    use "replicate_v12/base_alt.dta", clear
    replace s = rv_first24 - $CUT
    replace D = (s >= 0)
    trio, h(`h') outcome(death)
    post `PA' ("First value at or above the threshold within 24 h") (`h') (r(n)) (r(ev)) (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
}
use "replicate_v12/base_alt.dta", clear
replace s = rv_whole - $CUT
replace D = (s >= 0)
trio, h(5) outcome(death)
post `PA' ("Highest value over the whole stay") (5) (r(n)) (r(ev)) (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
postclose `PA'
use `RA', clear
list, clean noobs
export excel using "$OUTP", sheet("eT10_alt_assign") sheetreplace firstrow(variables)

* ---------- eTable 12 panel B: vasopressor or inotrope at hospitals B and C by calendar day and by window after the alert ----------
tempname PB
tempfile RB
postfile `PB' str32 outcome n ev pre_below pre_above post_below post_above dd ddlo ddhi ddp po polo pohi pop pr prlo prhi prp using `RB', replace
foreach y in vaso_all_Bd0_any vaso_all_Bd1_any vaso_all_Bd2_any vaso_all_Bd3_any vaso_all_A24h_any vaso_all_A48h_any vaso_all_A72h_any {
    use "replicate_v12/base_pm.dta", clear
    keep if inlist(H, 2, 3)
    foreach t in 0 1 {
        foreach sd in 0 1 {
            qui count if post == `t' & D == `sd' & abs(s) < 5 & `y' == 1
            local e`t'`sd' = r(N)
        }
    }
    trio, h(5) outcome(`y')
    post `PB' ("`y'") (r(n)) (r(ev)) (`e00') (`e01') (`e10') (`e11') (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
}
postclose `PB'
use `RB', clear
list, clean noobs
export excel using "$OUTP", sheet("eT12B_vaso_BC") sheetreplace firstrow(variables)

* ---------- eTable 13 A: nine markers, hospitals A-C and B-C, windows 5 and 8 ----------
local MARK "vital_any_A1h_any vital_any_A2h_any vital_any_A3h_any vaso_all_Bd1_any vaso_all_Bd12_any o2_escalation_A1h_new o2_escalation_A24h_new o2_escalation_A72h_new blood_culture_Bd01_any"
tempname P13
tempfile R13
postfile `P13' str32 marker str4 hosp h n ev dd ddlo ddhi ddp po polo pohi pop pr prlo prhi prp using `R13', replace
foreach y of local MARK {
    foreach hs in ABC BC {
        foreach h in 5 8 {
            use "replicate_v12/base_pm.dta", clear
            if "`hs'" == "BC" keep if inlist(H, 2, 3)
            trio, h(`h') outcome(`y')
            post `P13' ("`y'") ("`hs'") (`h') (r(n)) (r(ev)) (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
        }
    }
}
postclose `P13'
use `R13', clear
list, clean noobs
export excel using "$OUTP", sheet("eT13A_markers") sheetreplace firstrow(variables)

* ---------- eTable 13 B: linear probability models (window 5, hospitals A-C) ----------
tempname PL
tempfile RL
postfile `PL' str32 marker dd_pp dd_lo dd_hi dd_p post_pp post_lo post_hi post_p using `RL', replace
foreach y of local MARK {
    use "replicate_v12/base_pm.dta", clear
    lpm, h(5) outcome(`y') est(dd)
    local a1 = r(pp)
    local a2 = r(lo)
    local a3 = r(hi)
    local a4 = r(p)
    lpm, h(5) outcome(`y') est(post)
    post `PL' ("`y'") (`a1') (`a2') (`a3') (`a4') (r(pp)) (r(lo)) (r(hi)) (r(p))
}
postclose `PL'
use `RL', clear
list, clean noobs
export excel using "$OUTP", sheet("eT13B_lpm") sheetreplace firstrow(variables)

* ---------- eTable 13 C: hospital-level difference-in-discontinuities (window 5) ----------
tempname PH
tempfile RH
postfile `PH' str32 marker str4 hosp n ev dd ddlo ddhi ddp po pop using `RH', replace
foreach y of local MARK {
    forvalues hh = 1/3 {
        use "replicate_v12/base_pm.dta", clear
        keep if H == `hh'
        local hn : word `hh' of A B C
        trio, h(5) outcome(`y')
        post `PH' ("`y'") ("`hn'") (r(n)) (r(ev)) (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(pop))
    }
}
postclose `PH'
use `RH', clear
list, clean noobs
export excel using "$OUTP", sheet("eT13C_hospitals") sheetreplace firstrow(variables)

* ---------- eTable 14: reference-time rules R0-R9 for patients below the threshold (window 5, hospitals A-C) ----------
tempname PR
tempfile RR
postfile `PR' str32 marker rule n ev t0_median_below dd ddlo ddhi ddp po polo pohi pop pr prlo prhi prp using `RR', replace
forvalues r = 0/9 {
    foreach y in vital_any_A1h_any vital_any_A2h_any o2_escalation_A1h_any o2_escalation_A1h_new vaso_all_Bd12_any vaso_all_Bd1_any {
        use "rules_t90_v1.dta", clear
        keep if rule == `r'
        drop rule
        merge 1:1 adm_id using "replicate_v12/base.dta", keep(match) nogen
        qui summarize t0hr if s < 0 & abs(s) < 5, detail
        local med = r(p50)
        trio, h(5) outcome(`y')
        post `PR' ("`y'") (`r') (r(n)) (r(ev)) (`med') (r(dd)) (r(ddlo)) (r(ddhi)) (r(ddp)) (r(po)) (r(polo)) (r(pohi)) (r(pop)) (r(pr)) (r(prlo)) (r(prhi)) (r(prp))
    }
}
postclose `PR'
use `RR', clear
list, clean noobs
export excel using "$OUTP", sheet("eT14_rules") sheetreplace firstrow(variables)

* ---------- placebo thresholds: the nine markers at displayed scores 80 ... 62 (eTable 13 B placebo column; eFigure 6) ----------
* specifications with at least 15 events and 15 non-events within 5 points of the placebo threshold;
* the reference time is applied to that threshold (first value at or above it within 24 h; 24-h maximum below it)
tempname PP
tempfile RP
postfile `PP' cut str32 spec n ev dd ddlo ddhi ddp po polo pohi pop using `RP', replace
foreach c in 80 78 76 74 72 70 68 66 64 62 {
    use "replicate_v12/base.dta", clear
    replace s = rv_dec - (`c' - 0.5)
    replace D = (s >= 0)
    merge 1:1 adm_id using "process_t`c'_v1.dta", keep(match) nogen
    keep if abs(s) < 5
    save "replicate_v12/base_plac.dta", replace
    foreach y of local MARK {
        use "replicate_v12/base_plac.dta", clear
        qui count if `y' == 1
        local ev = r(N)
        qui count
        local n = r(N)
        if `ev' >= 15 & `n' - `ev' >= 15 {
            didisc, h(5) outcome(`y')
            local dd = exp(r(b))
            local ddlo = exp(r(b) - 1.96*r(se))
            local ddhi = exp(r(b) + 1.96*r(se))
            local ddp = r(p)
            srd, h(5) post(1) outcome(`y')
            local po = exp(r(b))
            local polo = exp(r(b) - 1.96*r(se))
            local pohi = exp(r(b) + 1.96*r(se))
            local pop = r(p)
            post `PP' (`c') ("`y'") (`n') (`ev') (`dd') (`ddlo') (`ddhi') (`ddp') (`po') (`polo') (`pohi') (`pop')
        }
    }
    di "placebo threshold `c' done"
}
postclose `PP'

* placebo thresholds at hospitals B and C for the two Figure 3 markers (eFigure 6, panels C and D)
tempname PB
tempfile RB
postfile `PB' cut str32 spec str4 hosp n ev dd ddlo ddhi ddp po polo pohi pop using `RB', replace
foreach c in 80 78 76 74 72 70 68 66 64 62 {
    use "replicate_v12/base.dta", clear
    replace s = rv_dec - (`c' - 0.5)
    replace D = (s >= 0)
    merge 1:1 adm_id using "process_t`c'_v1.dta", keep(match) nogen
    keep if abs(s) < 5 & (hB == 1 | hC == 1)
    save "replicate_v12/base_plac_bc.dta", replace
    foreach y in vital_any_A2h_any vaso_all_Bd12_any {
        use "replicate_v12/base_plac_bc.dta", clear
        qui count if `y' == 1
        local ev = r(N)
        qui count
        local n = r(N)
        didisc, h(5) outcome(`y')
        local dd = exp(r(b))
        local ddlo = exp(r(b) - 1.96*r(se))
        local ddhi = exp(r(b) + 1.96*r(se))
        local ddp = r(p)
        srd, h(5) post(1) outcome(`y')
        local po = exp(r(b))
        local polo = exp(r(b) - 1.96*r(se))
        local pohi = exp(r(b) + 1.96*r(se))
        local pop = r(p)
        post `PB' (`c') ("`y'") ("BC") (`n') (`ev') (`dd') (`ddlo') (`ddhi') (`ddp') (`po') (`polo') (`pohi') (`pop')
    }
    di "placebo threshold `c' (hospitals B and C) done"
}
postclose `PB'
use `RB', clear
save "replicate_v12/placebo_markers_BC.dta", replace
export excel using "$OUTP", sheet("placebo_markers_BC") sheetreplace firstrow(variables)
use `RP', clear
save "replicate_v12/placebo_markers.dta", replace
export excel using "$OUTP", sheet("placebo_markers") sheetreplace firstrow(variables)

* eTable 13 B placebo column: thresholds (of 10) at which the difference-in-discontinuities OR exceeded 1 with P < .05
gen byte hit = (dd > 1 & ddp < .05) if ddp < .
replace hit = 0 if hit == .
gen byte one = 1
collapse (sum) placebo_eval = one placebo_hits = hit, by(spec)
rename spec marker
merge 1:1 marker using `RL', keep(match) nogen
list, clean noobs
export excel using "$OUTP", sheet("eT13B_placebo") sheetreplace firstrow(variables)

* ---------- eFigure 6: Figure 3 markers at the alert threshold and the 10 placebo thresholds, hospitals A-C (panels A and B) and hospitals B and C (panels C and D) ----------
use "replicate_v12/placebo_markers.dta", clear
keep if inlist(spec, "vital_any_A2h_any", "vaso_all_Bd12_any")
keep cut spec n ev dd ddlo ddhi ddp po polo pohi pop
gen str4 hosp = "ABC"
append using "replicate_v12/placebo_markers_BC.dta"
tempfile PLF
save `PLF'
use `R13', clear
keep if h == 5 & inlist(marker, "vital_any_A2h_any", "vaso_all_Bd12_any")
gen cut = 90
rename marker spec
keep cut spec hosp n ev dd ddlo ddhi ddp po polo pohi pop
append using `PLF'
sort hosp spec cut
list, clean noobs
save "replicate_v12/eF6_placebo.dta", replace
export excel using "$OUTP", sheet("eF6_placebo") sheetreplace firstrow(variables)
gen x = (cut - 60)/2
replace x = 11 if cut == 90
gen x1 = x - 0.15
gen x2 = x + 0.15
local k = 0
foreach hs in ABC BC {
    foreach y in vital_any_A2h_any vaso_all_Bd12_any {
        local k = `k' + 1
        local L : word `k' of A B C D
        local ttl = cond("`y'" == "vital_any_A2h_any", "Vital-sign measurement within 2 h", "Vasopressor or inotrope on days 1-2")
        local hl = cond("`hs'" == "ABC", "hospitals A-C", "hospitals B and C")
        twoway (rcap ddlo ddhi x1 if spec == "`y'" & hosp == "`hs'", lcolor(navy)) (scatter dd x1 if spec == "`y'" & hosp == "`hs'", mcolor(navy) msymbol(O)) ///
            (rcap polo pohi x2 if spec == "`y'" & hosp == "`hs'", lcolor(dkgreen)) (scatter po x2 if spec == "`y'" & hosp == "`hs'", mcolor(dkgreen) msymbol(O)), ///
            yscale(log range(0.1 40)) ylabel(0.1 0.25 0.5 1 2 5 10 20, angle(0)) yline(1, lcolor(gs8)) ///
            xlabel(1 "62" 2 "64" 3 "66" 4 "68" 5 "70" 6 "72" 7 "74" 8 "76" 9 "78" 10 "80" 11 "90") ///
            xtitle("Threshold (displayed score)", size(small)) ytitle("OR (95% CI), log scale", size(small)) ///
            title("`L'  `ttl'" "`hl'", position(11) size(small)) legend(order(2 "Difference-in-discontinuities" 4 "Postimplementation discontinuity") rows(2) size(vsmall)) ///
            name(EF6`L', replace) nodraw
    }
}
graph combine EF6A EF6B EF6C EF6D, cols(2) ysize(7) xsize(8)
graph export "replicate_v12/eF6_placebo_stata.png", replace width(2400)
log close
