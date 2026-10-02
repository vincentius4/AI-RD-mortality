* Appendix Tables 2, 4-8
di as txt _n "===== 4. rdrobust ====="
tempname R
postfile `R' double(h post_conv post_conv_se post_conv_p post_bc post_rb_se post_rb_p pre_conv pre_conv_se pre_conv_p pre_bc pre_rb_se pre_rb_p ///
    dd_conv dd_conv_se dd_conv_p dd_rb dd_rb_se dd_rb_p) using "replicate_v12\res_rdrobust.dta", replace
forvalues h = 1/10 {
    foreach q in 1 0 {
        local nm = cond(`q' == 1, "post", "pre")
        foreach k in conv conv_se conv_p bc rb_se rb_p {
            local `nm'_`k' = .
        }
        cap qui rdrobust death rv_dec if post == `q', c($CUT) h(`h') b(`h') p(1) kernel(triangular) covs($COVS) vce(cluster clust)
        if _rc == 0 {
            local `nm'_conv = 100*e(tau_cl)
            local `nm'_conv_se = 100*e(se_tau_cl)
            local `nm'_conv_p = 2*normal(-abs(e(tau_cl)/e(se_tau_cl)))
            local `nm'_bc = 100*e(tau_bc)
            local `nm'_rb_se = 100*e(se_tau_rb)
            local `nm'_rb_p = 2*normal(-abs(e(tau_bc)/e(se_tau_rb)))
        }
    }
    local dd_conv = `post_conv' - `pre_conv'
    local dd_conv_se = sqrt(`post_conv_se'^2 + `pre_conv_se'^2)
    local dd_rb = `post_bc' - `pre_bc'
    local dd_rb_se = sqrt(`post_rb_se'^2 + `pre_rb_se'^2)
    post `R' (`h') (`post_conv') (`post_conv_se') (`post_conv_p') (`post_bc') (`post_rb_se') (`post_rb_p') ///
        (`pre_conv') (`pre_conv_se') (`pre_conv_p') (`pre_bc') (`pre_rb_se') (`pre_rb_p') ///
        (`dd_conv') (`dd_conv_se') (2*normal(-abs(`dd_conv'/`dd_conv_se'))) (`dd_rb') (`dd_rb_se') (2*normal(-abs(`dd_rb'/`dd_rb_se')))
}
postclose `R'
xout "replicate_v12\res_rdrobust.dta" "rdrobust"

tempname B
postfile `B' double(kind h b_used tau_cl se_cl p_cl tau_bc se_rb p_rb) using "replicate_v12\res_biasbw.dta", replace
foreach h in 5 8 {
    local k = 0
    foreach m in 1 1.5 2 {
        local ++k
        local b = `h'*`m'
        cap qui rdrobust death rv_dec if post == 1, c($CUT) h(`h') b(`b') p(1) kernel(triangular) covs($COVS) vce(cluster clust)
        if _rc == 0 post `B' (`k') (`h') (`b') (100*e(tau_cl)) (100*e(se_tau_cl)) (e(pv_cl)) (100*e(tau_bc)) (100*e(se_tau_rb)) (e(pv_rb))
    }
    cap qui rdbwselect death rv_dec if post == 1, c($CUT) p(1) kernel(triangular) bwselect(mserd) covs($COVS) vce(cluster clust)
    if _rc == 0 {
        local bsel = e(b_mserd)
        cap qui rdrobust death rv_dec if post == 1, c($CUT) h(`h') b(`bsel') p(1) kernel(triangular) covs($COVS) vce(cluster clust)
        if _rc == 0 post `B' (4) (`h') (`bsel') (100*e(tau_cl)) (100*e(se_tau_cl)) (e(pv_cl)) (100*e(tau_bc)) (100*e(se_tau_rb)) (e(pv_rb))
    }
}
cap qui rdrobust death rv_dec if post == 1, c($CUT) p(1) kernel(triangular) bwselect(mserd) covs($COVS) vce(cluster clust)
if _rc == 0 {
    post `B' (5) (e(h_l)) (e(b_l)) (100*e(tau_cl)) (100*e(se_tau_cl)) (e(pv_cl)) (100*e(tau_bc)) (100*e(se_tau_rb)) (e(pv_rb))
}
postclose `B'
xout "replicate_v12\res_biasbw.dta" "biasbw"

tempname C
postfile `C' double(post h side n deaths quad_prob p_quad_prob quad_logit p_quad_logit) using "replicate_v12\res_curvature.dta", replace
foreach p in 1 0 {
    foreach h in 5 8 10 {
        foreach side in 0 1 {
            preserve
            qui keep if post == `p' & abs(s) < `h' & D == `side'
            mkw `h' tri
            qui gen double s2 = s^2
            qui count if death == 1
            local dth = r(N)
            qui count
            local nn = r(N)
            local q1 = .
            local pq1 = .
            local q2 = .
            local pq2 = .
            cap qui regress death s s2 [pw=kw], vce(cluster clust)
            if _rc == 0 {
                seadj
                local q1 = _b[s2]
                local pq1 = 2*normal(-abs(_b[s2]/(_se[s2]*r(adj))))
            }
            cap qui logit death s s2 [pw=kw], vce(cluster clust) iter(100)
            if _rc == 0 {
                seadj
                local q2 = _b[s2]
                local pq2 = 2*normal(-abs(_b[s2]/(_se[s2]*r(adj))))
            }
            post `C' (`p') (`h') (`side') (`nn') (`dth') (`q1') (`pq1') (`q2') (`pq2')
            restore
        }
    }
}
postclose `C'
xout "replicate_v12\res_curvature.dta" "curvature"

tempname LL
postfile `LL' str8 period str16 scale double(points r2_linear slope) using "replicate_v12\res_loglinear.dta", replace
foreach p in 1 0 2 {
    preserve
    if `p' < 2 qui keep if post == `p'
    qui keep if disp >= 70 & disp <= 99
    collapse (mean) m = death (count) n = death, by(disp)
    qui keep if n >= 20
    gen double lm = ln(max(m, 1e-4))
    local per = cond(`p' == 1, "After", cond(`p' == 0, "Before", "Both"))
    qui regress lm disp
    post `LL' ("`per'") ("log mortality") (e(N)) (e(r2)) (_b[disp])
    qui regress m disp
    post `LL' ("`per'") ("mortality") (e(N)) (e(r2)) (_b[disp])
    restore
}
postclose `LL'
xout "replicate_v12\res_loglinear.dta" "loglinear"

tempname DN
postfile `DN' str24 period double(rdd_p jump jump_se jump_p) using "replicate_v12\res_density.dta", replace
foreach p in 1 0 {
    local lab = cond(`p' == 1, "After", "Before")
    local rdp = .
    cap qui rddensity s if post == `p', c(0)
    if _rc == 0 local rdp = e(pv_q)
    preserve
    qui keep if post == `p'
    cap noi bootstrap r(jump), reps(200) seed(7) nowarn: densjump, h(5)
    if _rc == 0 {
        matrix b = e(b)
        matrix V = e(V)
        local j = b[1,1]
        local jse = sqrt(V[1,1])
    }
    else {
        densjump, h(5)
        local j = r(jump)
        local jse = .
    }
    restore
    post `DN' ("`lab'") (`rdp') (`j') (`jse') (2*normal(-abs(`j'/`jse')))
    local J`p' = `j'
    local S`p' = `jse'
}
local dj = `J1' - `J0'
local dse = sqrt(`S1'^2 + `S0'^2)
post `DN' ("Change (after-before)") (.) (`dj') (`dse') (2*normal(-abs(`dj'/`dse')))
postclose `DN'
xout "replicate_v12\res_density.dta" "density"

tempname CV
postfile `CV' str16 covariate double(h disc se p) using "replicate_v12\res_cov_after.dta", replace
foreach v in $COVS rv_news {
    foreach h in 5 8 {
        preserve
        qui keep if post == 1 & abs(s) < `h'
        mkw `h' tri
        qui gen double Ds = D*s
        cap qui regress `v' D s Ds hB hC [pw=kw], vce(cluster clust)
        if _rc == 0 {
            seadj
            post `CV' ("`v'") (`h') (_b[D]) (_se[D]*r(adj)) (2*normal(-abs(_b[D]/(_se[D]*r(adj)))))
        }
        restore
    }
}
postclose `CV'
xout "replicate_v12\res_cov_after.dta" "cov_after"

tempname CD
postfile `CD' str16 covariate double(h disc_before p_before disc_after p_after diff p_diff) using "replicate_v12\res_cov_didisc.dta", replace
foreach v in age male cci htn sofa surgical_dept summer fall winter sepsis dnr_at_admit hB hC rv_news {
    foreach h in 5 8 {
        preserve
        qui keep if abs(s) < `h'
        mkw `h' tri
        qui gen double Ds = D*s
        qui gen double pD = post*D
        qui gen double ps = post*s
        qui gen double pDs = post*D*s
        qui gen double hB_post = hB*post
        qui gen double hC_post = hC*post
        local hosp hB hC hB_post hC_post
        if "`v'" == "hB" | "`v'" == "hC" local hosp
        cap qui regress `v' D s Ds post pD ps pDs `hosp' [pw=kw], vce(cluster clust)
        if _rc == 0 {
            seadj
            local adj = r(adj)
            local b0 = _b[D]
            local p0 = 2*normal(-abs(_b[D]/(_se[D]*`adj')))
            local bd = _b[pD]
            local pd = 2*normal(-abs(_b[pD]/(_se[pD]*`adj')))
            qui lincom D + pD
            post `CD' ("`v'") (`h') (`b0') (`p0') (r(estimate)) (2*normal(-abs(r(estimate)/(r(se)*`adj')))) (`bd') (`pd')
        }
        restore
    }
}
postclose `CD'
xout "replicate_v12\res_cov_didisc.dta" "cov_didisc"

tempname S
postfile `S' str40 analysis double(h ev dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "replicate_v12\res_secondary.dta", replace
foreach h in 5 8 {
    foreach v in death comp_gw comp_any arr_nondnr arr_gw icu_uit_24h {
        didisc, h(`h') outcome(`v')
        local ev = r(ev)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1) outcome(`v')
        post `S' ("`v'") (`h') (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    }
    foreach ex in dnr_at_admit dnr_by_24h dnr_case {
        preserve
        qui drop if `ex' == 1
        didisc, h(`h')
        local ev = r(ev)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1)
        post `S' ("Excluding `ex'") (`h') (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
        restore
    }
    didisc, h(`h') covs(age_z male cci htn sofa surgical_dept)
    local ev = r(ev)
    local a = exp(r(b))
    local al = exp(r(b) - 1.96*r(se))
    local ah = exp(r(b) + 1.96*r(se))
    local ap = r(p)
    srd, h(`h') post(1) covs(age_z male cci htn sofa surgical_dept)
    post `S' ("Season terms removed") (`h') (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    preserve
    qui drop if summer == 1
    didisc, h(`h') covs(age_z male cci htn sofa surgical_dept fall winter)
    local ev = r(ev)
    local a = exp(r(b))
    local al = exp(r(b) - 1.96*r(se))
    local ah = exp(r(b) + 1.96*r(se))
    local ap = r(p)
    srd, h(`h') post(1) covs(age_z male cci htn sofa surgical_dept fall winter)
    post `S' ("Summer admissions excluded") (`h') (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    restore
}
postclose `S'
xout "replicate_v12\res_secondary.dta" "secondary"

tempname P
postfile `P' double(displayed cutoff ev dd_or dd_lo dd_hi dd_p) using "replicate_v12\res_placebo.dta", replace
foreach c in 80 83 85 87 90 92 94 96 {
    preserve
    qui replace s = rv_dec - (`c' - 0.5)
    qui replace D = (s >= 0)
    didisc, h(5)
    post `P' (`c') (`c'-0.5) (r(ev)) (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
    restore
}
postclose `P'
xout "replicate_v12\res_placebo.dta" "placebo"

tempname P1
postfile `P1' double(displayed cutoff ev dd_or dd_lo dd_hi dd_p) using "replicate_v12\res_placebo_onesided.dta", replace
foreach c in 80 83 85 87 90 92 94 96 {
    preserve
    if `c' < 90 qui keep if rv_dec < 89.5
    if `c' > 90 qui keep if rv_dec >= 89.5
    qui replace s = rv_dec - (`c' - 0.5)
    qui replace D = (s >= 0)
    didisc, h(5)
    post `P1' (`c') (`c'-0.5) (r(ev)) (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
    restore
}
postclose `P1'
xout "replicate_v12\res_placebo_onesided.dta" "placebo_onesided"

tempname NW
postfile `NW' double(news_cut h ev dd_or dd_lo dd_hi dd_p) using "replicate_v12\res_news.dta", replace
foreach c in 5 6 7 {
    foreach h in 3 4 {
        preserve
        qui replace s = rv_news - `c'
        qui replace D = (s >= 0)
        didisc, h(`h')
        post `NW' (`c') (`h') (r(ev)) (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
        restore
    }
}
postclose `NW'
xout "replicate_v12\res_news.dta" "news"

tempname SC
postfile `SC' double(h displayed ev n dd_or dd_lo dd_hi dd_p) using "replicate_v12\res_scan.dta", replace
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
            if r(b) < . post `SC' (`hh') (`c') (r(ev)) (r(n)) (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
        }
        restore
    }
}
postclose `SC'
xout "replicate_v12\res_scan.dta" "scan"
preserve
use "replicate_v12\res_scan.dta", clear
gen double absl = abs(ln(dd_or))
tempname SS
postfile `SS' double(h or_90 p_90 n_placebo n_placebo_sig n_more_extreme rand_p or_min or_max) using "replicate_v12\res_scan_summary.dta", replace
foreach hh in 3 5 10 {
    qui su absl if displayed == 90 & h == `hh'
    local tr = r(mean)
    qui su dd_or if displayed == 90 & h == `hh'
    local o90 = r(mean)
    qui su dd_p if displayed == 90 & h == `hh'
    local p90 = r(mean)
    qui count if displayed != 90 & h == `hh' & absl < .
    local npl = r(N)
    qui count if displayed != 90 & h == `hh' & dd_p < 0.05
    local nsig = r(N)
    qui count if displayed != 90 & h == `hh' & absl >= `tr' - 1e-12 & absl < .
    local nex = r(N)
    qui su dd_or if displayed != 90 & h == `hh'
    post `SS' (`hh') (`o90') (`p90') (`npl') (`nsig') (`nex') ((1 + `nex')/(`npl' + 1)) (r(min)) (r(max))
}
postclose `SS'
restore
xout "replicate_v12\res_scan_summary.dta" "scan_summary"

tempname PD
postfile `PD' str40 sample double(h n_pre pre_or pre_p dd_or dd_lo dd_hi dd_p) using "replicate_v12\res_pandemic.dta", replace
foreach h in 5 8 {
    foreach k in 0 6 12 {
        preserve
        if `k' > 0 qui drop if post == 0 & t <= `k'
        qui count if post == 0
        local npre = r(N)
        srd, h(`h') post(0)
        local por = exp(r(b))
        local pp = r(p)
        didisc, h(`h')
        local lab = cond(`k' == 0, "Full pre-implementation period", cond(`k' == 6, "Excluding Jan-Jun 2022", "Excluding 2022"))
        post `PD' ("`lab'") (`h') (`npre') (`por') (`pp') (exp(r(b))) (exp(r(b)-1.96*r(se))) (exp(r(b)+1.96*r(se))) (r(p))
        restore
    }
}
postclose `PD'
xout "replicate_v12\res_pandemic.dta" "pandemic"

tempname PY
postfile `PY' str12 years double(h deaths n pre_or pre_lo pre_hi pre_p p_equal) using "replicate_v12\res_pre_by_year.dta", replace
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
        local b`y' = .
        local se`y' = .
        cap qui logit death D s Ds hB hC $COVS [pw=kw] if y22 == `y', vce(cluster clust) iter(100)
        if _rc == 0 {
            seadj
            local b`y' = _b[D]
            local se`y' = _se[D]*r(adj)
        }
    }
    local peq = 2*normal(-abs((`b1' - `b0')/sqrt(`se1'^2 + `se0'^2)))
    foreach y in 1 0 {
        post `PY' (cond(`y'==1,"2022","2023 onward")) (`h') (`dth`y'') (`nn`y'') (exp(`b`y'')) (exp(`b`y''-1.96*`se`y'')) (exp(`b`y''+1.96*`se`y'')) (2*normal(-abs(`b`y''/`se`y''))) (cond(`y'==1, `peq', .))
    }
    restore
}
postclose `PY\'
xout "replicate_v12\res_pre_by_year.dta" "pre_by_year"

tempname ED
postfile `ED' str8 period double(h below above disc p) using "replicate_v12\res_early_discharge.dta", replace
foreach h in 5 8 {
    foreach p in 0 1 {
        preserve
        qui keep if post == `p' & abs(s) < `h'
        mkw `h' tri
        qui gen double Ds = D*s
        qui su disc3 if D == 0
        local mb = 100*r(mean)
        qui su disc3 if D == 1
        local ma = 100*r(mean)
        cap qui regress disc3 D s Ds hB hC $COVS [pw=kw], vce(cluster clust)
        if _rc == 0 {
            seadj
            post `ED' (cond(`p' == 1, "After", "Before")) (`h') (`mb') (`ma') (100*_b[D]) (2*normal(-abs(_b[D]/(_se[D]*r(adj)))))
        }
        restore
    }
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
    cap qui regress disc3 D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust)
    if _rc == 0 {
        seadj
        post `ED' ("DiDisc") (`h') (.) (.) (100*_b[pD]) (2*normal(-abs(_b[pD]/(_se[pD]*r(adj)))))
    }
    restore
}
postclose `ED'
xout "replicate_v12\res_early_discharge.dta" "early_discharge"

