* Appendix Figures 2-6
tempname SP
postfile `SP' str40 spec double(h dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "replicate_v12\res_spec.dta", replace
foreach h in 5 8 {
    local i = 0
    foreach spec in "Reference" "No covariates" "Local quadratic" "Uniform kernel" "Epanechnikov kernel" "Donut 0.5 point" "Donut 1 point" "Cluster by displayed score" "No clustering (HC)" "No hospital-specific terms" {
        local ++i
        local opt
        local optr
        if `i' == 2 {
            local opt covs(none)
            local optr covs(none)
        }
        if `i' == 3 {
            local opt poly(2)
            local optr poly(2)
        }
        if `i' == 4 {
            local opt kern(uni)
            local optr kern(uni)
        }
        if `i' == 5 {
            local opt kern(epa)
            local optr kern(epa)
        }
        if `i' == 6 {
            local opt donut(0.5)
            local optr donut(0.5)
        }
        if `i' == 7 {
            local opt donut(1)
            local optr donut(1)
        }
        if `i' == 8 {
            local opt vce(cluster spt)
            local optr vce(cluster spt)
        }
        if `i' == 9 {
            local opt vce(robust)
            local optr vce(robust)
        }
        if `i' == 10 local opt nosat
        didisc, h(`h') `opt'
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        if `i' < 10 {
            srd, h(`h') post(1) `optr'
            post `SP' ("`spec'") (`h') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
        }
        else post `SP' ("`spec'") (`h') (`a') (`al') (`ah') (`ap') (.) (.) (.) (.)
    }
}
postclose `SP'
xout "replicate_v12\res_spec.dta" "spec"

tempname SG
postfile `SG' str24 subgroup double(level ev dd_or dd_lo dd_hi dd_p rd_or rd_p p_interaction) using "replicate_v12\res_subgroup.dta", replace
foreach g in age75 male cci2 sepsis surgical_dept news5 {
    local b1 = .
    local s1 = .
    local b0 = .
    local s0 = .
    foreach lv in 1 0 {
        preserve
        qui keep if `g' == `lv'
        didisc, h(5)
        local b`lv' = r(b)
        local s`lv' = r(se)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        local ev = r(ev)
        srd, h(5) post(1)
        post `SG' ("`g'") (`lv') (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (r(p)) (.)
        restore
    }
    post `SG' ("`g' interaction") (.) (.) (.) (.) (.) (.) (.) (.) (2*normal(-abs((`b1'-`b0')/sqrt(`s1'^2+`s0'^2))))
}
postclose `SG'
xout "replicate_v12\res_subgroup.dta" "subgroup"

tempname HO
postfile `HO' str20 unit double(h ev dd_or dd_lo dd_hi dd_p rd_or rd_p pre_or pre_p) using "replicate_v12\res_hospital.dta", replace
foreach h in 5 8 {
    forvalues k = 1/3 {
        preserve
        qui keep if hosp_num == `k'
        didisc, h(`h') nosat
        local ev = r(ev)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1)
        local r1 = exp(r(b))
        local r1p = r(p)
        srd, h(`h') post(0)
        post `HO' ("Hospital `k'") (`h') (`ev') (`a') (`al') (`ah') (`ap') (`r1') (`r1p') (exp(r(b))) (r(p))
        restore
    }
    forvalues k = 1/3 {
        preserve
        qui drop if hosp_num == `k'
        didisc, h(`h')
        local ev = r(ev)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1)
        local r1 = exp(r(b))
        local r1p = r(p)
        srd, h(`h') post(0)
        post `HO' ("Omitting hospital `k'") (`h') (`ev') (`a') (`al') (`ah') (`ap') (`r1') (`r1p') (exp(r(b))) (r(p))
        restore
    }
}
postclose `HO'
xout "replicate_v12\res_hospital.dta" "hospital"

cap confirm variable h_death
if _rc {
}
else {
    tempname GR
    postfile `GR' double(L H h n ev dd_or dd_p) using "replicate_v12\res_grid.dta", replace
    foreach L in 24 48 72 120 {
        foreach H in 168 336 720 1000000 {
            foreach h in 5 8 {
                preserve
                qui keep if h_stay > `L'
                qui replace death = (death == 1 & h_death <= `L' + `H')
                didisc, h(`h')
                post `GR' (`L') (`H') (`h') (r(n)) (r(ev)) (exp(r(b))) (r(p))
                restore
            }
        }
    }
    postclose `GR'
    xout "replicate_v12\res_grid.dta" "grid_L_by_H"
}

tempname TP
postfile `TP' str70 scenario double(n_added dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "replicate_v12\res_tipping.dta", replace
local S1 "0"
local S2 "abs(s)<5 & post==1 & D==1 & dnr_by_24h==1 & death==0"
local S3 "abs(s)<5 & post==1 & dnr_by_24h==1 & death==0"
local S4 "abs(s)<5 & dnr_by_24h==1 & death==0"
local S5 "abs(s)<5 & post==1 & D==1 & dnr_case==1 & death==0"
local S6 "abs(s)<5 & dnr_case==1 & death==0"
local S7 "abs(s)<5 & post==1 & D==1 & death==0 & pick==1"
local L1 "Observed data"
local L2 "Order within 24 h, above the threshold, after implementation"
local L3 "Order within 24 h, both sides, after implementation"
local L4 "Order within 24 h, all four groups"
local L5 "Order at any time, above the threshold, after implementation"
local L6 "Order at any time, all four groups"
local L7 "Survivors at random, above the threshold, after implementation"
set seed 7
cap drop u pick
gen double u = runiform() if abs(s)<5 & post==1 & D==1 & death==0
gen byte pick = 0
sort u
replace pick = 1 if u < . & _n <= 12
forvalues i = 1/7 {
    preserve
    local nadd = 0
    if `i' > 1 {
        qui count if `S`i''
        local nadd = r(N)
        qui replace death = 1 if `S`i''
    }
    didisc, h(5)
    local a = exp(r(b))
    local al = exp(r(b) - 1.96*r(se))
    local ah = exp(r(b) + 1.96*r(se))
    local ap = r(p)
    srd, h(5) post(1)
    post `TP' ("`L`i''") (`nadd') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    restore
}
postclose `TP'
xout "replicate_v12\res_tipping.dta" "tipping"
preserve
use "replicate_v12\res_tipping.dta", clear
gen int row = _n
gen double row_dd = row - 0.15
gen double row_rd = row + 0.15
cap noi twoway (rcap dd_lo dd_hi row_dd, horizontal lcolor(black)) (scatter row_dd dd_or, mcolor(black) msymbol(O)) ///
    (rcap rd_lo rd_hi row_rd, horizontal lcolor(gs9)) (scatter row_rd rd_or, mcolor(gs9) msymbol(S)), ///
    xline(1, lpattern(dash)) xscale(log) xlabel(0.1 0.25 0.5 1 2) yscale(reverse) ///
    ylabel(1 "Observed data" 2 "Order within 24 h, above, after" 3 "Order within 24 h, both sides, after" 4 "Order within 24 h, all four groups" 5 "Order at any time, above, after" 6 "Order at any time, all four groups" 7 "Survivors at random (12), above, after", angle(0) labsize(vsmall)) ///
    ytitle("") xtitle("Odds ratio for in-hospital death") ///
    legend(order(2 "Difference-in-discontinuities" 4 "Post-implementation discontinuity") rows(1) size(vsmall)) graphregion(color(white))
cap noi graph export "replicate_v12\FA_tipping.png", replace width(1800)
restore

