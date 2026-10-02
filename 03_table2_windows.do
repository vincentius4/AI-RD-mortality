* Table 2 (panel A); Appendix Table 3; Figures 2-3
curvM, outcome(death)
global M_death = r(M)
tempname M
postfile `M' double(h n ev evp dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_lo pre_hi pre_p dd_p_stata rd_p_stata n_post ev_pre n_pre h_lo h_hi h_p M) using "replicate_v12\res_windows.dta", replace
forvalues h = 1/10 {
    didisc, h(`h')
    local n = r(n)
    local ev = r(ev)
    local evp = r(evp)
    local a = exp(r(b))
    local al = exp(r(b) - 1.96*r(se))
    local ah = exp(r(b) + 1.96*r(se))
    local ap = r(p)
    local aps = r(p_stata)
    srd, h(`h') post(1)
    local r1 = exp(r(b))
    local r1l = exp(r(b) - 1.96*r(se))
    local r1h = exp(r(b) + 1.96*r(se))
    local r1p = r(p)
    local r1ps = r(p_stata)
    local n1 = r(n)
    srdh, h(`h') post(1) m($M_death)
    local hl = exp(r(lo))
    local hh = exp(r(hi))
    local hp = r(p)
    srd, h(`h') post(0)
    post `M' (`h') (`n') (`ev') (`evp') (`a') (`al') (`ah') (`ap') (`r1') (`r1l') (`r1h') (`r1p') ///
        (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p)) (`aps') (`r1ps') (`n1') (r(ev)) (r(n)) (`hl') (`hh') (`hp') ($M_death)
}
postclose `M'
xout "replicate_v12\res_windows.dta" "windows"

tempname HM
postfile `HM' double(h M_mult or lo hi honest_p) using "replicate_v12\res_honestM.dta", replace
foreach h in 5 8 {
    foreach mu in 0.5 1 2 {
        srdh, h(`h') post(1) m(`=`mu'*$M_death')
        post `HM' (`h') (`mu') (exp(r(b))) (exp(r(lo))) (exp(r(hi))) (r(p))
    }
}
postclose `HM'
xout "replicate_v12\res_honestM.dta" "honestM"
preserve
use "replicate_v12\res_windows.dta", clear
keep if h >= 2
gen double h1 = h - 0.2
gen double h3 = h + 0.2
cap noi twoway (rcap dd_lo dd_hi h1, lcolor(black)) (scatter dd_or h1, mcolor(black)) ///
    (rcap rd_lo rd_hi h, lcolor(gs8)) (scatter rd_or h, mcolor(gs8) msymbol(D)) ///
    (rcap pre_lo pre_hi h3, lcolor(gs11)) (scatter pre_or h3, mcolor(gs11) msymbol(T)), ///
    yline(1, lpattern(dash)) yscale(log) ylabel(0.05 0.1 0.25 0.5 1 2 5 10) xlabel(2(1)10) ///
    xtitle("Window (points)") ytitle("Odds ratio for in-hospital death") ///
    legend(order(2 "Difference-in-discontinuities" 4 "Post-implementation" 6 "Pre-implementation (placebo)") rows(1) size(vsmall)) ///
    graphregion(color(white))
cap noi graph export "replicate_v12\F3_windows.png", replace width(1800)
restore

preserve
gen double bin = floor(rv_dec/0.5)*0.5 + 0.25
collapse (mean) death (count) n = death, by(post bin)
keep if n >= 20 & bin >= 40 & bin <= 100
gen double pct = 100*death
save "replicate_v12\res_F2_bins.dta", replace
export excel using "$OUTX", sheet("F2_bins") sheetreplace firstrow(variables)
restore

