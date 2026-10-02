* Tables (manuscript layout)
version 15
clear all
set more off
global R "replicate_v12"
global X "$R/tables.xlsx"
global CUT 89.5
cap log close
log using "$R/tables_v2.log", replace text
cap erase "$X"

program define fmt_or
    args new or lo hi
    if "`lo'" == "" {
        qui gen str40 `new' = string(`or', "%9.2f") if `or' < .
    }
    else {
        qui gen str40 `new' = string(`or', "%9.2f") + " (" + string(`lo', "%9.2f") + "–" + string(`hi', "%9.2f") + ")" if `or' < .
    }
    qui replace `new' = "not estimable" if `or' >= .
    qui replace `new' = subinstr(subinstr(`new', ".", "·", .), "-", "−", .)
end
program define fmt_p
    args new p
    qui gen str12 `new' = cond(`p' < 0.0001, "<0.0001", cond(`p' >= 0.0995, string(`p', "%4.2f"), cond(`p' >= 0.00995, string(`p', "%5.3f"), string(`p', "%6.4f")))) if `p' < .
    qui replace `new' = "··" if `p' >= .
    qui replace `new' = subinstr(`new', ".", "·", .)
end
program define fmt_x
    args new x dec
    qui gen str20 `new' = string(`x', "%12.`dec'f") if `x' < .
    qui replace `new' = "··" if `x' >= .
    qui replace `new' = subinstr(subinstr(`new', ".", "·", .), "-", "−", .)
end
program define fmt_es
    args new x se
    qui gen str30 `new' = string(`x', "%12.2f") + " (" + string(`se', "%12.2f") + ")" if `x' < .
    qui replace `new' = subinstr(subinstr(`new', ".", "·", .), "-", "−", .)
end
program define fmt_n
    args new x
    qui gen str14 `new' = string(`x', "%12.0f") if `x' < .
    qui replace `new' = substr(`new', 1, length(`new') - 3) + " " + substr(`new', -3, 3) if `x' >= 10000 & `x' < .
    qui replace `new' = "··" if `x' >= .
end
program define tosheet
    args sheet
    export excel using "$X", sheet("`sheet'") sheetreplace firstrow(varlabels)
end
program define mkw
    args h
    cap drop kw
    gen double kw = max(0, 1 - abs(s)/`h')
end

di as txt _n "===== Table 1 ====="
use "$R/res_T1.dta", clear
drop if row == "DNR order at admission"
replace row = "AI score, median (IQR)" if row == "DeepCARS value, median (IQR)"
replace row = "Surgical specialty" if row == "Surgical service"
forvalues j = 1/8 {
    replace c`j' = ustrregexra(c`j', "(?<![0-9·])([0-9]{2})([0-9]{3})(?![0-9·])", "$1 $2")
}
local heads `" "Before: score 85–89" "Before: score 90–94" "Before: score <90" "Before: score ≥90" "After: score 85–89" "After: score 90–94" "After: score <90" "After: score ≥90" "'
label var row "Characteristic"
forvalues j = 1/8 {
    local hd : word `j' of `heads'
    local nn = c`j'[1]
    label var c`j' "`hd' (n=`nn')"
}
list, noobs clean
tosheet "Table 1"

di as txt _n "===== Table 2 ====="
program define t2panel
    args file panel
    use "$R/`file'.dta", clear
    keep if inlist(h, 3, 5, 8, 10)
    tempfile a b
    preserve
    gen str60 estimator = "Difference-in-discontinuities (primary)"
    fmt_or orci dd_or dd_lo dd_hi
    fmt_p pv dd_p
    fmt_n events ev
    fmt_n admissions n
    gen byte ord = 1
    keep h estimator orci pv events admissions ord
    save `a'
    restore
    preserve
    gen str60 estimator = "Single RD, after implementation"
    fmt_or orci rd_or rd_lo rd_hi
    fmt_p p1 rd_p
    fmt_p p2 h_p
    gen str30 pv = p1 + " (" + p2 + ")†"
    fmt_n events evp
    fmt_n admissions n_post
    gen byte ord = 2
    keep h estimator orci pv events admissions ord
    save `b'
    restore
    gen str60 estimator = "Single RD, before implementation"
    fmt_or orci pre_or pre_lo pre_hi
    fmt_p pv pre_p
    fmt_n events ev_pre
    fmt_n admissions n_pre
    gen byte ord = 3
    keep h estimator orci pv events admissions ord
    append using `a' `b'
    gen str60 panel = "`panel'"
    gen str6 window = "±" + string(h)
    sort ord h
end
t2panel res_windows "Panel A. In-hospital death (primary outcome)"
tempfile pa
save `pa'
t2panel res_windows_comp_gw "Panel B. Death or ward cardiac arrest"
gen byte pn = 2
append using `pa'
replace pn = 1 if pn >= .
sort pn ord h
keep panel estimator window orci pv events admissions
order panel estimator window orci pv events admissions
label var panel "Panel"
label var estimator "Estimator"
label var window "Window"
label var orci "Odds ratio (95% CI)"
label var pv "p value (bias-aware p)†"
label var events "Events"
label var admissions "Admissions"
list, noobs clean
tosheet "Table 2"

use "$R/res_apptable1.dta", clear
fmt_n adm admissions
fmt_n dth deaths
gen str20 arr = string(arrests_any) + " (" + string(arrests_ward) + ")"
keep population hosp_num post adm dth arr
reshape wide adm dth arr, i(population hosp_num) j(post)
gen str12 hospital = "Hospital " + char(64 + hosp_num)
gen str20 pop = cond(trim(population) == "all", "All admissions", "Landmark cohort")
gen byte po = cond(trim(population) == "all", 1, 2)
sort hosp_num po
keep hospital pop adm0 dth0 arr0 adm1 dth1 arr1
order hospital pop adm0 dth0 arr0 adm1 dth1 arr1
label var hospital "Hospital"
label var pop "Population"
label var adm0 "Before: admissions"
label var dth0 "Before: deaths"
label var arr0 "Before: cardiac arrests, all (ward)"
label var adm1 "After: admissions"
label var dth1 "After: deaths"
label var arr1 "After: cardiac arrests, all (ward)"
list, noobs clean
tosheet "App T1"

use "$R/res_rdrobust.dta", clear
gen str6 window = "±" + string(h)
fmt_es a1 post_conv post_conv_se
fmt_p a2 post_conv_p
fmt_es a3 post_bc post_rb_se
fmt_p a4 post_rb_p
fmt_es a5 pre_conv pre_conv_se
fmt_p a6 pre_conv_p
fmt_es a7 pre_bc pre_rb_se
fmt_p a8 pre_rb_p
fmt_es b1 dd_conv dd_conv_se
fmt_p b2 dd_conv_p
fmt_es b3 dd_rb dd_rb_se
fmt_p b4 dd_rb_p
preserve
keep window a1-a8
label var window "Window"
label var a1 "After, conventional: estimate (SE)"
label var a2 "p value"
label var a3 "After, robust bias-corrected: estimate (SE)"
label var a4 "p value"
label var a5 "Before, conventional: estimate (SE)"
label var a6 "p value"
label var a7 "Before, robust bias-corrected: estimate (SE)"
label var a8 "p value"
list, noobs clean
tosheet "App T2a"
restore
keep window b1-b4
label var window "Window"
label var b1 "Conventional: estimate (SE)"
label var b2 "p value"
label var b3 "Robust bias-corrected: estimate (SE)"
label var b4 "p value"
list, noobs clean
tosheet "App T2b"

use "paper_data_v9.dta", clear
keep if rv_dec < . & h_stay > 24 & h_stay < .
gen double s = rv_dec - $CUT
gen byte D = (s >= 0)
gen byte post = treat
tempname RC
tempfile relc
postfile `RC' double(post h side rel) using `relc', replace
foreach p in 1 0 {
    foreach h in 5 8 10 {
        foreach side in 0 1 {
            preserve
            qui keep if post == `p' & abs(s) < `h' & D == `side'
            mkw `h'
            qui gen double s2 = s^2
            qui regress death s s2 [pw=kw]
            post `RC' (`p') (`h') (`side') (abs(_b[s2])*`h'^2/max(abs(_b[_cons]), 1e-6))
            restore
        }
    }
}
postclose `RC'
use "$R/res_curvature.dta", clear
merge 1:1 post h side using `relc', nogen
gen byte po = cond(post == 1, 1, 2)
sort po h side
gen str24 period = cond(post == 1, "After implementation", "Before implementation")
gen str6 window = "±" + string(h)
gen str12 sidel = cond(side == 1, "At or above", "Below")
fmt_n dth deaths
fmt_p pp p_quad_prob
fmt_p pl p_quad_logit
fmt_x rc rel 2
keep period window sidel dth pp pl rc
label var period "Period"
label var window "Window"
label var sidel "Side"
label var dth "Deaths"
label var pp "Quadratic term p value: probability scale"
label var pl "Quadratic term p value: logit scale"
label var rc "Relative curvature"
list, noobs clean
tosheet "App T2c"

use "$R/res_biasbw.dta", clear
gen str12 window = "±" + string(h) if kind < 5
replace window = "±" + subinstr(string(h, "%3.1f"), ".", "·", .) + "†" if kind == 5
gen str40 bw = cond(kind == 1, "b = h", cond(kind == 2, "b = 1·5h", cond(kind == 3, "b = 2h", "")))
replace bw = "b = data-driven (" + subinstr(string(b_used, "%4.2f"), ".", "·", .) + ")" if kind == 4
replace bw = "both data-driven (b=" + subinstr(string(b_used, "%4.2f"), ".", "·", .) + ")" if kind == 5
fmt_x d1 tau_cl 2
fmt_p d2 p_cl
fmt_x d3 tau_bc 2
fmt_x d4 se_rb 2
fmt_p d5 p_rb
keep window bw d1-d5
label var window "Window"
label var bw "Bias bandwidth"
label var d1 "Conventional: estimate, pp"
label var d2 "Conventional: p value"
label var d3 "Robust bias-corrected: estimate, pp"
label var d4 "Robust SE, pp"
label var d5 "Robust p value"
list, noobs clean
tosheet "App T2d"

use "$R/res_windows.dta", clear
gen str6 window = "±" + string(h)
gen str20 dn = string(ev) + " (" + string(n) + ")"
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
fmt_or c3 rd_or rd_lo rd_hi
fmt_p c4 rd_p
fmt_p c5 h_p
fmt_or c6 pre_or pre_lo pre_hi
fmt_p c7 pre_p
keep window dn c1-c7
label var window "Window"
label var dn "Deaths (admissions)"
label var c1 "Difference-in-discontinuities: OR (95% CI)"
label var c2 "p value"
label var c3 "Single RD after implementation: OR (95% CI)"
label var c4 "p value"
label var c5 "Bias-aware p value"
label var c6 "Single RD before implementation: OR (95% CI)"
label var c7 "p value"
list, noobs clean
tosheet "App T3a"
use "$R/res_integer.dta", clear
gen str6 window = "±" + string(h)
gen str20 dn = string(ev) + " (" + string(n) + ")"
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
fmt_or c3 rd_or rd_lo rd_hi
fmt_p c4 rd_p
fmt_or c5 pre_or
fmt_p c6 pre_p
keep window dn c1-c6
label var window "Window"
label var dn "Deaths (admissions)"
label var c1 "Difference-in-discontinuities: OR (95% CI)"
label var c2 "p value"
label var c3 "Single RD after implementation: OR (95% CI)"
label var c4 "p value"
label var c5 "Single RD before implementation: OR"
label var c6 "p value"
list, noobs clean
tosheet "App T3b"

use "$R/res_placebo_onesided.dta", clear
gen byte last = (displayed == 90)
sort last displayed
gen str6 thr = string(displayed)
fmt_n dth ev
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
keep thr dth c1 c2
label var thr "Placebo threshold (displayed score); 90 = true threshold"
label var dth "Deaths"
label var c1 "OR (95% CI)"
label var c2 "p value"
list, noobs clean
tosheet "App T4a"
use "$R/res_news.dta", clear
sort news_cut h
gen str6 thr = string(news_cut)
gen str6 window = string(h)
fmt_n dth ev
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
keep thr window dth c1 c2
label var thr "NEWS threshold"
label var window "Window"
label var dth "Deaths"
label var c1 "OR (95% CI)"
label var c2 "p value"
list, noobs clean
tosheet "App T4b"
use "$R/res_scan_summary.dta", clear
sort h
gen str6 window = "±" + string(h)
fmt_or c1 or_90
fmt_p c2 p_90
fmt_n c3 n_placebo
fmt_n c4 n_placebo_sig
fmt_n c5 n_more_extreme
fmt_p c6 rand_p
gen str20 c7 = subinstr(string(or_min, "%9.2f") + "–" + string(or_max, "%9.2f"), ".", "·", .)
keep window c1-c7
label var window "Window"
label var c1 "True threshold, OR"
label var c2 "p value"
label var c3 "Placebo thresholds: estimable, n"
label var c4 "p<0·05, n"
label var c5 "More extreme than true threshold, n"
label var c6 "Randomisation p value"
label var c7 "Range of ORs"
list, noobs clean
tosheet "App T4c"

use "$R/res_windows.dta", clear
keep if h == 5
local nall = n[1]
local npost = n_post[1]
local npre = n_pre[1]
use "$R/res_density.dta", clear
gen byte ord = cond(period == "Before", 1, cond(period == "After", 2, 3))
sort ord
gen str40 q = cond(ord == 1, "Before implementation", cond(ord == 2, "After implementation", "Difference (after minus before)"))
fmt_x c1 jump 3
gen double ratio = exp(jump)
fmt_x c2 ratio 2
fmt_p c3 jump_p
replace c3 = c3 + " (Stata)"
gen double nn = cond(ord == 1, `npre', cond(ord == 2, `npost', `nall'))
fmt_n c4 nn
keep q c1 c2 c3 c4
order q c1 c2 c3 c4
label var q "Quantity"
label var c1 "Log-density jump"
label var c2 "Density ratio (above/below)"
label var c3 "Bootstrap p value (random numbers; Stata)"
label var c4 "Admissions"
list, noobs clean
tosheet "App T5a"

use "$R/res_cov_after.dta", clear
gen byte bin = inlist(covariate, "male", "htn", "surgical_dept", "summer", "fall", "winter")
gen double dsc = cond(bin, 100*disc, disc)
fmt_x c1 dsc 2
fmt_p c2 p
keep covariate h c1 c2
reshape wide c1 c2, i(covariate) j(h)
gen byte ord = .
local k = 0
foreach v in age_z male cci htn sofa surgical_dept summer fall winter rv_news {
    local ++k
    replace ord = `k' if covariate == "`v'"
}
sort ord
gen str50 lab = ""
replace lab = "Age (standardised)" if covariate == "age_z"
replace lab = "Male sex" if covariate == "male"
replace lab = "Charlson comorbidity index" if covariate == "cci"
replace lab = "Hypertension" if covariate == "htn"
replace lab = "SOFA score at admission" if covariate == "sofa"
replace lab = "Surgical specialty" if covariate == "surgical_dept"
replace lab = "Admitted in summer" if covariate == "summer"
replace lab = "Admitted in autumn" if covariate == "fall"
replace lab = "Admitted in winter" if covariate == "winter"
replace lab = "National Early Warning Score, first 24 h" if covariate == "rv_news"
keep lab c15 c25 c18 c28
order lab c15 c25 c18 c28
label var lab "Covariate"
label var c15 "±5: discontinuity"
label var c25 "±5: p value"
label var c18 "±8: discontinuity"
label var c28 "±8: p value"
list, noobs clean
tosheet "App T5b"

use "$R/res_cov_didisc.dta", clear
keep if h == 5 & covariate != "dnr_at_admit"
gen byte bin = inlist(covariate, "male", "htn", "surgical_dept", "sepsis", "summer", "fall", "winter", "hB", "hC")
foreach v in disc_before disc_after diff {
    replace `v' = 100*`v' if bin
}
fmt_x c1 disc_before 2
fmt_p c2 p_before
fmt_x c3 disc_after 2
fmt_p c4 p_after
fmt_x c5 diff 2
fmt_p c6 p_diff
gen byte ord = .
local k = 0
foreach v in age male cci htn sofa surgical_dept sepsis summer fall winter hB hC rv_news {
    local ++k
    replace ord = `k' if covariate == "`v'"
}
sort ord
gen str50 lab = ""
replace lab = "Age, years" if covariate == "age"
replace lab = "Male sex" if covariate == "male"
replace lab = "Charlson comorbidity index" if covariate == "cci"
replace lab = "Hypertension" if covariate == "htn"
replace lab = "SOFA score at admission" if covariate == "sofa"
replace lab = "Surgical specialty" if covariate == "surgical_dept"
replace lab = "Sepsis" if covariate == "sepsis"
replace lab = "Admitted in summer" if covariate == "summer"
replace lab = "Admitted in autumn" if covariate == "fall"
replace lab = "Admitted in winter" if covariate == "winter"
replace lab = "Hospital B" if covariate == "hB"
replace lab = "Hospital C" if covariate == "hC"
replace lab = "National Early Warning Score, first 24 h" if covariate == "rv_news"
keep lab c1 c2 c3 c4 c5 c6
order lab c1 c2 c3 c4 c5 c6
label var lab "Covariate"
label var c1 "Before implementation: discontinuity"
label var c2 "p value"
label var c3 "After implementation: discontinuity"
label var c4 "p value"
label var c5 "Between periods: difference"
label var c6 "p value"
list, noobs clean
tosheet "App T5c"

use "$R/res_entry.dta", clear
gen byte ord = cond(check == "excl24", 1, cond(check == "died24", 2, 3))
sort ord h
gen str50 lab = cond(ord == 1, "Not in cohort at 24 h (died or discharged)", cond(ord == 2, "Died before 24 h", "Discharged alive before 24 h"))
gen str6 window = "±" + string(h)
fmt_n c0 ev
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
fmt_or c3 post_or post_lo post_hi
fmt_p c4 post_p
fmt_or c5 pre_or pre_lo pre_hi
fmt_p c6 pre_p
keep lab window c0-c6
label var lab "Check"
label var window "Window"
label var c0 "Events"
label var c1 "Difference-in-discontinuities: OR (95% CI)"
label var c2 "p value"
label var c3 "After implementation: OR (95% CI)"
label var c4 "p value"
label var c5 "Before implementation: OR (95% CI)"
label var c6 "p value"
list, noobs clean
tosheet "App T5d"

use "$R/res_secondary.dta", clear
append using "$R/res_secondary_add.dta"
keep if h == 5
preserve
gen byte ord = .
local k = 0
foreach v in comp_gw comp_any arr_nondnr "Death, any arrest, or unplanned ICU" arr_gw icu_uit_24h {
    local ++k
    replace ord = `k' if analysis == "`v'"
}
keep if ord < .
sort ord
gen str70 lab = ""
replace lab = "Death or ward cardiac arrest" if ord == 1
replace lab = "Death or any cardiac arrest" if ord == 2
replace lab = "Any cardiac arrest or death without a DNR order" if ord == 3
replace lab = "Death, any cardiac arrest, or unplanned intensive care transfer" if ord == 4
replace lab = "Ward cardiac arrest" if ord == 5
replace lab = "Unplanned transfer to intensive care after 24 h" if ord == 6
fmt_n c0 ev
fmt_or c1 dd_or dd_lo dd_hi
fmt_p c2 dd_p
fmt_or c3 rd_or rd_lo rd_hi
fmt_p c4 rd_p
keep lab c0-c4
label var lab "Outcome"
label var c0 "Events"
label var c1 "Difference-in-discontinuities: OR (95% CI)"
label var c2 "p value"
label var c3 "Single RD after implementation: OR (95% CI)"
label var c4 "p value"
list, noobs clean
tosheet "App T6A"
restore
preserve
use "paper_data_v9.dta", clear
keep if rv_dec < . & h_stay > 24 & h_stay < .
gen double s = rv_dec - $CUT
qui count if abs(s) < 5 & s < 0 & dnr_case == 1
local dL = r(N)
qui count if abs(s) < 5 & s >= 0 & dnr_case == 1
local dR = r(N)
qui count if abs(s) < 5 & s < 0 & summer == 1
local sL = r(N)
qui count if abs(s) < 5 & s >= 0 & summer == 1
local sR = r(N)
restore
gen byte ord = .
replace ord = 1 if analysis == "Excluding dnr_case"
replace ord = 2 if analysis == "Season terms removed"
replace ord = 3 if analysis == "Hospital-by-season terms added"
replace ord = 4 if analysis == "Summer admissions excluded"
keep if ord < .
sort ord
gen str50 lab = cond(ord == 1, "Excluding any DNR order during stay", analysis)
gen str20 c0 = cond(ord == 1, "−`dL' / −`dR'", cond(ord == 4, "−`sL' / −`sR'", "··"))
fmt_or c1 dd_or
fmt_p c2 dd_p
fmt_or c3 rd_or
fmt_p c4 rd_p
keep lab c0-c4
label var lab "Analysis"
label var c0 "Admissions removed (below / above)"
label var c1 "Difference-in-discontinuities: OR"
label var c2 "p value"
label var c3 "Single RD after implementation: OR"
label var c4 "p value"
list, noobs clean
tosheet "App T6B"

use "$R/res_charting.dta", clear
gen byte o = 1
append using "$R/res_early_discharge.dta"
replace o = 2 if o >= .
gen byte po = cond(period == "Before", 1, cond(period == "After", 2, 3))
sort o h po
gen str50 outcome = cond(o == 1, "Vital-sign charting, entries per hour", "Discharged alive within 3 days, %")
gen str6 window = "±" + string(h)
gen str40 per = cond(po == 1, "Before implementation", cond(po == 2, "After implementation", "Difference-in-discontinuities"))
fmt_x c1 below 2
fmt_x c2 above 2
fmt_x c3 disc 2
fmt_p c4 p
keep outcome window per c1-c4
label var outcome "Outcome"
label var window "Window"
label var per "Period"
label var c1 "Below threshold"
label var c2 "At or above"
label var c3 "Discontinuity"
label var c4 "p value"
list, noobs clean
tosheet "App T7"

use "$R/res_pandemic.dta", clear
gen byte ord = cond(sample == "Full pre-implementation period", 1, cond(sample == "Excluding Jan-Jun 2022", 2, 3))
sort ord h
gen str40 lab = cond(ord == 1, "Full pre-implementation period", cond(ord == 2, "Excluding Jan–June 2022", "Excluding all of 2022"))
gen str6 window = "±" + string(h)
fmt_n c0 n_pre
fmt_or c1 pre_or
fmt_p c2 pre_p
fmt_or c3 dd_or dd_lo dd_hi
fmt_p c4 dd_p
keep lab window c0-c4
label var lab "Sample"
label var window "Window"
label var c0 "Admissions"
label var c1 "Single RD before implementation: OR"
label var c2 "p value"
label var c3 "Difference-in-discontinuities: OR (95% CI)"
label var c4 "p value"
list, noobs clean
tosheet "App T8a"
use "$R/res_pre_by_year.dta", clear
gen byte ord = cond(years == "2022", 1, 2)
sort h ord
gen str20 lab = cond(ord == 1, "2022", "2023 onwards")
gen str6 window = "±" + string(h)
fmt_n c0 deaths
fmt_n c1 n
fmt_or c2 pre_or pre_lo pre_hi
fmt_p c3 pre_p
fmt_p c4 p_equal
keep lab window c0-c4
label var lab "Pre-implementation period"
label var window "Window"
label var c0 "Deaths"
label var c1 "Admissions"
label var c2 "OR (95% CI)"
label var c3 "p value"
label var c4 "p value for equality"
list, noobs clean
tosheet "App T8b"

use "$R/res_text.dta", clear
label var item "Item"
label var sub "Detail"
label var value "Value"
tosheet "Text and note numbers"
use "$R/res_honestM.dta", clear
export excel using "$X", sheet("Bias-aware M sensitivity") sheetreplace firstrow(variables)
use "$R/res_loglinear.dta", clear
export excel using "$X", sheet("Log linearity") sheetreplace firstrow(variables)
use "$R/res_localrand.dta", clear
export excel using "$X", sheet("LPM and local randomization") sheetreplace firstrow(variables)
use "$R/res_power.dta", clear
export excel using "$X", sheet("Detectable OR") sheetreplace firstrow(variables)
use "$R/res_hosp_deaths.dta", clear
export excel using "$X", sheet("Deaths by hospital") sheetreplace firstrow(variables)
use "$R/res_death_timing.dta", clear
export excel using "$X", sheet("Deaths by timing") sheetreplace firstrow(variables)

log close
