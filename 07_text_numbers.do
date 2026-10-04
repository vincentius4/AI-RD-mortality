* Numbers in text and legends
tempname TX
postfile `TX' str40 item str30 sub double(value) using "replicate_v12\res_text.dta", replace
use "paper_data_v9.dta", clear
forvalues k = 1/3 {
    qui su death if hosp_num == `k'
    post `TX' ("crude_mortality_pct") ("hospital `k'") (100*r(mean))
    qui count if hosp_num == `k'
    post `TX' ("admissions_all") ("hospital `k'") (r(N))
}
qui count if death == 1
post `TX' ("events_all") ("deaths") (r(N))
qui count if arr_any == 1
post `TX' ("events_all") ("arrests any") (r(N))
qui count if arr_gw == 1
post `TX' ("events_all") ("arrests ward") (r(N))
qui su dcars_over90 if treat == 1
post `TX' ("listed_pct_after") ("all admissions after") (100*r(mean))
preserve
qui keep if hosp_num == 1 & treat == 0 & h_stay > 24 & h_stay < .
gen double rv_w24 = rv_dec
foreach W in 24 6 12 48 {
    cap drop sw
    gen double sw = rv_w`W' - $CUT
    qui su death if sw < 0 & sw > -5
    post `TX' ("hospA_pre_W`W'") ("below pct") (100*r(mean))
    post `TX' ("hospA_pre_W`W'") ("below deaths") (r(sum))
    qui su death if sw >= 0 & sw < 5
    post `TX' ("hospA_pre_W`W'") ("above pct") (100*r(mean))
    post `TX' ("hospA_pre_W`W'") ("above deaths") (r(sum))
}
restore
use "replicate_v12/base.dta", clear
qui count if death == 1
local dth = r(N)
post `TX' ("landmark") ("deaths") (`dth')
qui count if arr_any == 1
post `TX' ("landmark") ("arrests any") (r(N))
qui count if arr_gw == 1
post `TX' ("landmark") ("arrests ward") (r(N))
qui count if death == 1 & arr_any == 0
post `TX' ("deaths_without_arrest") ("n") (r(N))
post `TX' ("deaths_without_arrest") ("pct") (100*r(N)/`dth')
preserve
qui keep if post == 1 & abs(s) < 5
qui count
post `TX' ("distinct_internal_values") ("admissions") (r(N))
contract rv_dec
qui count
post `TX' ("distinct_internal_values") ("distinct") (r(N))
restore
foreach h in 5 8 {
    qui count if post == 1 & s < 0 & s > -`h'
    local nn = r(N)
    qui count if post == 1 & s < 0 & s > -`h' & dcars_over90 == 1
    post `TX' ("leak_h`h'") ("below n") (`nn')
    post `TX' ("leak_h`h'") ("crossed") (r(N))
    post `TX' ("leak_h`h'") ("pct") (100*r(N)/`nn')
}
qui count if abs(s) < 5 & arr_gw == 1
post `TX' ("ward_arrest_h5") ("all") (r(N))
qui count if abs(s) < 5 & arr_gw == 1 & post == 1
post `TX' ("ward_arrest_h5") ("after") (r(N))
qui count if abs(s) < 5 & arr_gw == 1 & post == 0
post `TX' ("ward_arrest_h5") ("before") (r(N))
forvalues k = 1/3 {
    qui su dnr_time_rec if hosp_num == `k' & dnr_case == 1
    post `TX' ("dnr_time_recorded_pct") ("hospital `k'") (100*r(mean))
}
preserve
qui keep if post == 1 & abs(s) < 5
qui su summer if D == 0
post `TX' ("summer_pct") ("below") (100*r(mean))
qui su summer if D == 1
post `TX' ("summer_pct") ("above") (100*r(mean))
qui count if D == 0 & summer == 0
local a = r(N)
qui count if D == 0 & summer == 1
local b = r(N)
qui count if D == 1 & summer == 0
local c = r(N)
qui count if D == 1 & summer == 1
local d = r(N)
local N = `a' + `b' + `c' + `d'
local chi = `N'*(abs(`a'*`d' - `b'*`c') - `N'/2)^2/((`a' + `b')*(`c' + `d')*(`a' + `c')*(`b' + `d'))
post `TX' ("summer_pct") ("p chi2 corrected") (chi2tail(1, `chi'))
foreach dd in 88 89 {
    qui su summer if disp == `dd'
    post `TX' ("summer_pct") ("displayed `dd'") (100*r(mean))
}
restore
foreach h in 5 8 {
    foreach p in 1 0 {
        preserve
        qui keep if post == `p' & abs(s) < `h'
        mkw `h' tri
        qui gen double Ds = D*s
        qui logit death D s Ds [pw=kw], vce(cluster clust) iter(100)
        seadj
        local b = _b[D]
        local se = _se[D]*r(adj)
        local per = cond(`p' == 1, "after", "before")
        post `TX' ("fig2_unadjusted_h`h'") ("`per' OR") (exp(`b'))
        post `TX' ("fig2_unadjusted_h`h'") ("`per' lo") (exp(`b' - 1.96*`se'))
        post `TX' ("fig2_unadjusted_h`h'") ("`per' hi") (exp(`b' + 1.96*`se'))
        post `TX' ("fig2_unadjusted_h`h'") ("`per' p") (2*normal(-abs(`b'/`se')))
        restore
    }
}
post `TX' ("curvature_bound_M") ("death") ($M_death)
postclose `TX'
xout "replicate_v12\res_text.dta" "text_numbers"
preserve
use "replicate_v12\res_text.dta", clear
list, noobs clean
restore
use "replicate_v12/base.dta", clear
tempname HD
postfile `HD' str1 hospital str6 period str5 side double(n deaths los_median_deaths pct_within3d pct_dnr) using "replicate_v12\res_hosp_deaths.dta", replace
forvalues k = 1/3 {
    foreach p in 0 1 {
        foreach sd in 0 1 {
            local hl = char(64 + `k')
            local per = cond(`p' == 1, "after", "before")
            local sl = cond(`sd' == 1, "above", "below")
            qui count if abs(s) < 5 & hosp_num == `k' & post == `p' & D == `sd'
            local n = r(N)
            qui su total_los if abs(s) < 5 & hosp_num == `k' & post == `p' & D == `sd' & death == 1, detail
            local dth = r(N)
            local med = r(p50)
            qui count if abs(s) < 5 & hosp_num == `k' & post == `p' & D == `sd' & death == 1 & total_los <= 3
            local w3 = 100*r(N)/`dth'
            qui su dnr_case if abs(s) < 5 & hosp_num == `k' & post == `p' & D == `sd' & death == 1
            post `HD' ("`hl'") ("`per'") ("`sl'") (`n') (`dth') (`med') (`w3') (100*r(mean))
        }
    }
}
postclose `HD'
xout "replicate_v12\res_hosp_deaths.dta" "hosp_deaths"
tempname HG
postfile `HG' str1 hospital str6 period str5 group double(n deaths) using "replicate_v12\res_hospital_disp.dta", replace
forvalues k = 1/3 {
    foreach p in 0 1 {
        local hl = char(64 + `k')
        local per = cond(`p' == 1, "after", "before")
        forvalues dd = 85/94 {
            qui count if hosp_num == `k' & post == `p' & disp == `dd' & abs(s) < 5
            local n = r(N)
            qui count if hosp_num == `k' & post == `p' & disp == `dd' & abs(s) < 5 & death == 1
            post `HG' ("`hl'") ("`per'") ("`dd'") (`n') (r(N))
        }
        qui count if hosp_num == `k' & post == `p' & inrange(disp, 87, 89) & abs(s) < 5
        local n = r(N)
        qui count if hosp_num == `k' & post == `p' & inrange(disp, 87, 89) & abs(s) < 5 & death == 1
        post `HG' ("`hl'") ("`per'") ("87-89") (`n') (r(N))
        qui count if hosp_num == `k' & post == `p' & inrange(disp, 90, 91) & abs(s) < 5
        local n = r(N)
        qui count if hosp_num == `k' & post == `p' & inrange(disp, 90, 91) & abs(s) < 5 & death == 1
        post `HG' ("`hl'") ("`per'") ("90-91") (`n') (r(N))
    }
}
postclose `HG'
xout "replicate_v12\res_hospital_disp.dta" "hospital_disp"

preserve
qui keep if abs(s) < 5
gen byte tcat = cond(h_death <= 72, 1, cond(h_death <= 168, 2, cond(h_death <= 720, 3, 4))) if death == 1 & h_death > 0 & h_death < .
tempname DT
postfile `DT' str8 period str20 timing double(deaths above below at_risk_above at_risk_below rate_above rate_below dnr_any dnr_24h icu sepsis age sofa) using "replicate_v12\res_death_timing.dta", replace
foreach p in 0 1 {
    local per = cond(`p' == 1, "After", "Before")
    qui count if post == `p' & D == 1
    local na = r(N)
    qui count if post == `p' & D == 0
    local nb = r(N)
    forvalues t = 1/4 {
        local tl : word `t' of "Death on days 2–3" "Death on days 4–7" "Death on days 8–30" "Death after day 30"
        qui count if post == `p' & tcat == `t'
        local dth = r(N)
        qui count if post == `p' & tcat == `t' & D == 1
        local a = r(N)
        local b = `dth' - `a'
        foreach v in dnr_case dnr_by_24h icu_admission sepsis age sofa {
            qui su `v' if post == `p' & tcat == `t'
            local m_`v' = r(mean)
        }
        post `DT' ("`per'") ("`tl'") (`dth') (`a') (`b') (`na') (`nb') (100*`a'/`na') (100*`b'/`nb') (100*`m_dnr_case') (100*`m_dnr_by_24h') (100*`m_icu_admission') (100*`m_sepsis') (`m_age') (`m_sofa')
    }
}
postclose `DT'
restore
xout "replicate_v12\res_death_timing.dta" "death_timing"

tempname PW
postfile `PW' str50 outcome double(h events events_above events_below mde_or) using "replicate_v12\res_power.dta", replace
foreach v in death comp_gw arr_nondnr arr_any arr_gw {
    local lab = cond("`v'" == "death", "In-hospital death", cond("`v'" == "comp_gw", "Death or ward cardiac arrest", cond("`v'" == "arr_nondnr", "Cardiac arrest or death without a DNR order", cond("`v'" == "arr_any", "Any cardiac arrest", "Ward cardiac arrest"))))
    foreach h in 5 8 {
        qui count if post == 1 & abs(s) < `h' & D == 1 & `v' == 1
        local a = r(N)
        qui count if post == 1 & abs(s) < `h' & D == 1 & `v' == 0
        local b = r(N)
        qui count if post == 1 & abs(s) < `h' & D == 0 & `v' == 1
        local c = r(N)
        qui count if post == 1 & abs(s) < `h' & D == 0 & `v' == 0
        local d = r(N)
        local se = sqrt(1/max(`a', .5) + 1/max(`b', .5) + 1/max(`c', .5) + 1/max(`d', .5))
        post `PW' ("`lab'") (`h') (`a' + `c') (`a') (`c') (exp(-(1.96 + 0.84)*`se'))
    }
}
postclose `PW'
xout "replicate_v12\res_power.dta" "power"

log close
