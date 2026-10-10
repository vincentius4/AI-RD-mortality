* Appendix Tables 1, 5d, 6, 7; Table 2 (panel B); Appendix Figures 1F, 2
use "paper_data_v9.dta", clear
qui count
local n_all = r(N)
qui count if rv_dec < .
local n_sc = r(N)
qui count if rv_dec < . & h_stay > 24
local n_lm = r(N)
tempname A1
postfile `A1' str12 population double(hosp_num post admissions deaths arrests_any arrests_ward) using "replicate_v12\res_apptable1.dta", replace
foreach pop in all landmark {
    forvalues k = 1/3 {
        foreach p in 0 1 {
            local cnd "hosp_num == `k' & treat == `p'"
            if "`pop'" == "landmark" local cnd "`cnd' & rv_dec < . & h_stay > 24"
            qui count if `cnd'
            local n = r(N)
            qui count if `cnd' & death == 1
            local dd = r(N)
            qui count if `cnd' & arr_any == 1
            local aa = r(N)
            qui count if `cnd' & arr_gw == 1
            post `A1' ("`pop'") (`k') (`p') (`n') (`dd') (`aa') (r(N))
        }
    }
}
postclose `A1'
preserve
clear
set obs 1
gen double n_all = `n_all'
gen double n_scored = `n_sc'
gen double n_landmark = `n_lm'
save "replicate_v12\res_cohort_n.dta", replace
restore
xout "replicate_v12\res_apptable1.dta" "apptable1"
xout "replicate_v12\res_cohort_n.dta" "cohort_n"

cap mata: mata drop wbias()
cap mata: mata drop hci()
cap mata: mata drop permp()
cap mata: mata drop permp2()
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
    st_numscalar("r_hlo", est - mid*se)
    st_numscalar("r_hhi", est + mid*se)
    st_numscalar("r_hp", 1 - normal(z - t) + normal(-z - t))
}
real scalar permp(real colvector y, real colvector D, real scalar reps)
{
    real scalar stat, cnt, i
    real colvector Dp
    stat = mean(select(y, D)) - mean(select(y, !D))
    cnt = 0
    for (i = 1; i <= reps; i++) {
        Dp = jumble(D)
        cnt = cnt + (abs(mean(select(y, Dp)) - mean(select(y, !Dp))) >= abs(stat) - 1e-12)
    }
    return((cnt + 1)/(reps + 1))
}
real scalar permp2(real colvector y1, real colvector D1, real colvector y0, real colvector D0, real scalar reps)
{
    real scalar stat, cnt, i, d
    real colvector P1, P0
    stat = (mean(select(y1, D1)) - mean(select(y1, !D1))) - (mean(select(y0, D0)) - mean(select(y0, !D0)))
    cnt = 0
    for (i = 1; i <= reps; i++) {
        P1 = jumble(D1)
        P0 = jumble(D0)
        d = (mean(select(y1, P1)) - mean(select(y1, !P1))) - (mean(select(y0, P0)) - mean(select(y0, !P0)))
        cnt = cnt + (abs(d) >= abs(stat) - 1e-12)
    }
    return((cnt + 1)/(reps + 1))
}
end
use "replicate_v12/base.dta", clear
foreach v in comp_gw comp_any {
    curvM, outcome(`v')
    local Mc = r(M)
    tempname MC
    postfile `MC' double(h n ev evp dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_lo pre_hi pre_p n_post ev_pre n_pre h_lo h_hi h_p M) using "replicate_v12\res_windows_`v'.dta", replace
    forvalues h = 1/10 {
        didisc, h(`h') outcome(`v')
        local n = r(n)
        local ev = r(ev)
        local evp = r(evp)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1) outcome(`v')
        local r1 = exp(r(b))
        local r1l = exp(r(b) - 1.96*r(se))
        local r1h = exp(r(b) + 1.96*r(se))
        local r1p = r(p)
        local n1 = r(n)
        srdh, h(`h') post(1) m(`Mc') outcome(`v')
        local hl = exp(r(lo))
        local hh = exp(r(hi))
        local hp = r(p)
        srd, h(`h') post(0) outcome(`v')
        post `MC' (`h') (`n') (`ev') (`evp') (`a') (`al') (`ah') (`ap') (`r1') (`r1l') (`r1h') (`r1p') ///
            (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p)) (`n1') (r(ev)) (r(n)) (`hl') (`hh') (`hp') (`Mc')
    }
    postclose `MC'
    xout "replicate_v12\res_windows_`v'.dta" "windows_`v'"
}

use "paper_data_v9.dta", clear
keep if rv_dec < .
gen double s    = rv_dec - $CUT
gen byte   D    = (s >= 0)
gen byte   post = treat
gen byte   excl24  = (h_stay <= 24)
gen byte   died24  = (h_stay <= 24 & death == 1)
gen byte   disch24 = (h_stay <= 24 & death == 0)
tempname EN
postfile `EN' str12 check double(h ev n dd_or dd_lo dd_hi dd_p post_or post_lo post_hi post_p pre_or pre_lo pre_hi pre_p) using "replicate_v12\res_entry.dta", replace
foreach v in excl24 died24 disch24 {
    foreach h in 5 8 {
        didisc, h(`h') outcome(`v')
        local ev = r(ev)
        local n = r(n)
        local a = exp(r(b))
        local al = exp(r(b) - 1.96*r(se))
        local ah = exp(r(b) + 1.96*r(se))
        local ap = r(p)
        srd, h(`h') post(1) outcome(`v')
        local r1 = exp(r(b))
        local r1l = exp(r(b) - 1.96*r(se))
        local r1h = exp(r(b) + 1.96*r(se))
        local r1p = r(p)
        srd, h(`h') post(0) outcome(`v')
        post `EN' ("`v'") (`h') (`ev') (`n') (`a') (`al') (`ah') (`ap') (`r1') (`r1l') (`r1h') (`r1p') ///
            (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
    }
}
postclose `EN'
xout "replicate_v12\res_entry.dta" "entry"

use "replicate_v12/base.dta", clear
gen byte o6 = (death == 1) | (arr_any == 1) | (icu_uit_24h == 1)
foreach hh in B C {
    foreach ss in summer fall winter {
        gen double h`hh'_`ss' = h`hh'*`ss'
    }
}
tempname S6
postfile `S6' str40 analysis double(h ev dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "replicate_v12\res_secondary_add.dta", replace
didisc, h(5) outcome(o6)
local ev = r(ev)
local a = exp(r(b))
local al = exp(r(b) - 1.96*r(se))
local ah = exp(r(b) + 1.96*r(se))
local ap = r(p)
srd, h(5) post(1) outcome(o6)
post `S6' ("Death, any arrest, or unplanned ICU") (5) (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
local hs "$COVS hB_summer hB_fall hB_winter hC_summer hC_fall hC_winter"
didisc, h(5) covs(`hs')
local ev = r(ev)
local a = exp(r(b))
local al = exp(r(b) - 1.96*r(se))
local ah = exp(r(b) + 1.96*r(se))
local ap = r(p)
srd, h(5) post(1) covs(`hs')
post `S6' ("Hospital-by-season terms added") (5) (`ev') (`a') (`al') (`ah') (`ap') (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p))
postclose `S6'
xout "replicate_v12\res_secondary_add.dta" "secondary_add"

use "replicate_v12/base.dta", clear
replace s = rv - 90
replace D = (s >= 0)
tempname IN
postfile `IN' double(h n ev dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_p) using "replicate_v12\res_integer.dta", replace
foreach h in 3 5 8 10 {
    didisc, h(`h')
    local n = r(n)
    local ev = r(ev)
    local a = exp(r(b))
    local al = exp(r(b) - 1.96*r(se))
    local ah = exp(r(b) + 1.96*r(se))
    local ap = r(p)
    srd, h(`h') post(1)
    local r1 = exp(r(b))
    local r1l = exp(r(b) - 1.96*r(se))
    local r1h = exp(r(b) + 1.96*r(se))
    local r1p = r(p)
    srd, h(`h') post(0)
    post `IN' (`h') (`n') (`ev') (`a') (`al') (`ah') (`ap') (`r1') (`r1l') (`r1h') (`r1p') (exp(r(b))) (r(p))
}
postclose `IN'
xout "replicate_v12\res_integer.dta" "integer"

use "replicate_v12/base.dta", clear
keep if abs(s) < 5
gen str5 group = cond(D == 1, "above", "below")
collapse (count) n = death (sum) deaths = death, by(hosp_num post group)
save "replicate_v12\res_hospital_window.dta", replace
export excel using "$OUTX", sheet("hospital_window") sheetreplace firstrow(variables)

use "replicate_v12/base.dta", clear
qui logit death $COVS hB hC if abs(s) >= 5, iter(100)
predict double risk, pr
replace risk = 100*risk
tempname CI
postfile `CI' double(h disc_before p_before disc_after p_after dd p_dd) using "replicate_v12\res_covindex.dta", replace
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
    qui regress risk D s Ds post pD ps pDs hB hC hB_post hC_post [pw=kw], vce(cluster clust)
    seadj
    local adj = r(adj)
    local b0 = _b[D]
    local p0 = 2*normal(-abs(_b[D]/(_se[D]*`adj')))
    local bd = _b[pD]
    local pd = 2*normal(-abs(_b[pD]/(_se[pD]*`adj')))
    qui lincom D + pD
    post `CI' (`h') (`b0') (`p0') (r(estimate)) (2*normal(-abs(r(estimate)/(r(se)*`adj')))) (`bd') (`pd')
    restore
}
postclose `CI'
xout "replicate_v12\res_covindex.dta" "covindex"

cap mata: mata drop wbias()
cap mata: mata drop hci()
cap mata: mata drop permp()
cap mata: mata drop permp2()
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
    st_numscalar("r_hlo", est - mid*se)
    st_numscalar("r_hhi", est + mid*se)
    st_numscalar("r_hp", 1 - normal(z - t) + normal(-z - t))
}
real scalar permp(real colvector y, real colvector D, real scalar reps)
{
    real scalar stat, cnt, i
    real colvector Dp
    stat = mean(select(y, D)) - mean(select(y, !D))
    cnt = 0
    for (i = 1; i <= reps; i++) {
        Dp = jumble(D)
        cnt = cnt + (abs(mean(select(y, Dp)) - mean(select(y, !Dp))) >= abs(stat) - 1e-12)
    }
    return((cnt + 1)/(reps + 1))
}
real scalar permp2(real colvector y1, real colvector D1, real colvector y0, real colvector D0, real scalar reps)
{
    real scalar stat, cnt, i, d
    real colvector P1, P0
    stat = (mean(select(y1, D1)) - mean(select(y1, !D1))) - (mean(select(y0, D0)) - mean(select(y0, !D0)))
    cnt = 0
    for (i = 1; i <= reps; i++) {
        P1 = jumble(D1)
        P0 = jumble(D0)
        d = (mean(select(y1, P1)) - mean(select(y1, !P1))) - (mean(select(y0, P0)) - mean(select(y0, !P0)))
        cnt = cnt + (abs(d) >= abs(stat) - 1e-12)
    }
    return((cnt + 1)/(reps + 1))
}
end
use "replicate_v12/base.dta", clear
gen double hrs = min(max(h_stay - 24, 0), 48)
gen double rate = n_24_72/hrs if hrs > 0
tempname CH
postfile `CH' str8 period double(h below above disc p) using "replicate_v12\res_charting.dta", replace
foreach h in 5 8 {
    foreach p in 0 1 {
        preserve
        qui keep if post == `p' & abs(s) < `h' & rate < .
        mkw `h' tri
        qui gen double Ds = D*s
        qui su rate if D == 0
        local mb = r(mean)
        qui su rate if D == 1
        local ma = r(mean)
        qui regress rate D s Ds hB hC $COVS [pw=kw], vce(cluster clust)
        seadj
        post `CH' (cond(`p' == 1, "After", "Before")) (`h') (`mb') (`ma') (_b[D]) (2*normal(-abs(_b[D]/(_se[D]*r(adj)))))
        restore
    }
    preserve
    qui keep if abs(s) < `h' & rate < .
    mkw `h' tri
    qui gen double Ds = D*s
    qui gen double pD = post*D
    qui gen double ps = post*s
    qui gen double pDs = post*D*s
    qui gen double hB_post = hB*post
    qui gen double hC_post = hC*post
    qui gen double hB_D = hB*D
    qui gen double hC_D = hC*D
    qui regress rate D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust)
    seadj
    post `CH' ("DiDisc") (`h') (.) (.) (_b[pD]) (2*normal(-abs(_b[pD]/(_se[pD]*r(adj)))))
    restore
}
postclose `CH'
xout "replicate_v12\res_charting.dta" "charting"

use "replicate_v12/base.dta", clear
tempname LR
postfile `LR' str90 analysis double(h est lo hi p n ev) using "replicate_v12\res_localrand.dta", replace
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
    qui gen double hB_D = hB*D
    qui gen double hC_D = hC*D
    qui regress death D s Ds post pD ps pDs hB hC hB_post hC_post hB_D hC_D $COVS [pw=kw], vce(cluster clust)
    seadj
    local se = _se[pD]*r(adj)
    post `LR' ("Difference-in-discontinuities, linear probability (percentage points)") (`h') (100*_b[pD]) (100*(_b[pD] - 1.96*`se')) (100*(_b[pD] + 1.96*`se')) (2*normal(-abs(_b[pD]/`se'))) (e(N)) (.)
    restore
    preserve
    qui keep if post == 1 & abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui regress t_of_max D s Ds hB hC $COVS [pw=kw], vce(cluster clust)
    seadj
    local se = _se[D]*r(adj)
    post `LR' ("Timing of the maximum score (hours after admission), after implementation") (`h') (_b[D]) (_b[D] - 1.96*`se') (_b[D] + 1.96*`se') (2*normal(-abs(_b[D]/`se'))) (e(N)) (.)
    restore
}
set seed 20260907
foreach h in 2 5 {
    foreach p in 1 0 {
        preserve
        qui keep if post == `p' & abs(s) < `h'
        qui su death if D == 1
        local m1 = r(mean)
        qui su death if D == 0
        local m0 = r(mean)
        qui count if death == 1
        local ev = r(N)
        qui count
        local n = r(N)
        mata: st_numscalar("r_pp", permp(st_data(., "death"), st_data(., "D"), 2000))
        local d`p' = 100*(`m1' - `m0')
        post `LR' (cond(`p' == 1, "Local randomisation, after implementation (raw difference, 2000 permutations)", "Local randomisation, before implementation (placebo)")) (`h') (`d`p'') (.) (.) (scalar(r_pp)) (`n') (`ev')
        restore
    }
    preserve
    qui keep if abs(s) < `h'
    qui gen byte pre0 = (post == 0)
    mata: st_numscalar("r_pp", permp2(st_data(., "death", "post"), st_data(., "D", "post"), st_data(., "death", "pre0"), st_data(., "D", "pre0"), 2000))
    post `LR' ("Local randomisation, difference-in-discontinuities") (`h') (`d1' - `d0') (.) (.) (scalar(r_pp)) (.) (.)
    restore
}
postclose `LR'
xout "replicate_v12\res_localrand.dta" "localrand"

