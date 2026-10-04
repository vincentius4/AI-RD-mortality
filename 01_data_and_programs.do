* Data and programs
version 15
clear all
set more off
set seed 7
cap mkdir replicate_v12
cap log close
log using "replicate_v12\replicate_v12.log", replace text
global OUTX "replicate_v12\replicate_v12.xlsx"
global COVS age_z male cci htn sofa surgical_dept summer fall winter
global CUT 89.5

foreach p in rdrobust rddensity {
    cap which `p'
}

use "paper_data_v9.dta", clear
keep if rv_dec < .
keep if h_stay > 24
gen double s    = rv_dec - $CUT
gen byte   D    = (s >= 0)
gen byte   post = treat
gen int    spt  = round(rv_dec)
gen int    disp = round(rv_dec)
gen byte   arr_nondnr = (arr_any == 1) | (death == 1 & dnr_case == 0)
gen byte   disc3 = (total_los <= 3 & death == 0)
gen byte   news5 = (rv_news >= 5) if rv_news < .
gen byte   age75 = (age >= 75)
gen byte   cci2  = (cci >= 2)
count
save "replicate_v12/base.dta", replace

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

program define mkw
    args h kern
    cap drop kw
    if "`kern'" == "uni"      gen double kw = 1
    else if "`kern'" == "epa" gen double kw = max(0, 0.75*(1 - (s/`h')^2))
    else                      gen double kw = max(0, 1 - abs(s)/`h')
end

program define didisc, rclass
    syntax , h(real) [outcome(string) kern(string) poly(integer 1) donut(real 0) covs(string) vce(string) nosat]
    if "`outcome'" == "" local outcome death
    if "`kern'" == ""    local kern tri
    if "`vce'" == ""     local vce cluster clust
    if "`covs'" == ""    local covs $COVS
    if "`covs'" == "none" local covs
    preserve
    qui keep if abs(s) < `h' & abs(s) >= `donut'
    mkw `h' `kern'
    qui gen double Ds = D*s
    qui gen double pD = post*D
    qui gen double ps = post*s
    qui gen double pDs = post*D*s
    local q
    if `poly' == 2 {
        qui gen double s2 = s^2
        qui gen double Ds2 = D*s2
        qui gen double ps2 = post*s2
        qui gen double pDs2 = post*D*s2
        local q s2 Ds2 ps2 pDs2
    }
    local hosp hB hC
    if "`sat'" == "" {
        qui gen double hB_post = hB*post
        qui gen double hC_post = hC*post
        qui gen double hB_D = hB*D
        qui gen double hC_D = hC*D
        local hosp hB hC hB_post hC_post hB_D hC_D
    }
    qui count if `outcome' == 1
    local ev = r(N)
    qui count if `outcome' == 1 & post == 1
    local evp = r(N)
    qui count
    local n = r(N)
    cap qui logit `outcome' D s Ds post pD ps pDs `q' `hosp' `covs' [pw=kw], vce(`vce') iter(100)
    if _rc == 0 & `ev' >= 15 {
        seadj
        local adj = r(adj)
        return scalar b = _b[pD]
        return scalar se = _se[pD]*`adj'
        return scalar p = 2*normal(-abs(_b[pD]/(_se[pD]*`adj')))
        return scalar p_stata = 2*normal(-abs(_b[pD]/_se[pD]))
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar p = .
        return scalar p_stata = .
    }
    return scalar ev = `ev'
    return scalar evp = `evp'
    return scalar n = `n'
    restore
end

program define srd, rclass
    syntax , h(real) post(integer) [outcome(string) kern(string) poly(integer 1) donut(real 0) covs(string) vce(string)]
    if "`outcome'" == "" local outcome death
    if "`kern'" == ""    local kern tri
    if "`vce'" == ""     local vce cluster clust
    if "`covs'" == ""    local covs $COVS
    if "`covs'" == "none" local covs
    preserve
    qui keep if post == `post' & abs(s) < `h' & abs(s) >= `donut'
    mkw `h' `kern'
    qui gen double Ds = D*s
    local q
    if `poly' == 2 {
        qui gen double s2 = s^2
        qui gen double Ds2 = D*s2
        local q s2 Ds2
    }
    qui count if `outcome' == 1
    local ev = r(N)
    qui count
    local n = r(N)
    cap qui logit `outcome' D s Ds `q' hB hC `covs' [pw=kw], vce(`vce') iter(100)
    if _rc == 0 & `ev' >= 8 {
        seadj
        local adj = r(adj)
        return scalar b = _b[D]
        return scalar se = _se[D]*`adj'
        return scalar p = 2*normal(-abs(_b[D]/(_se[D]*`adj')))
        return scalar p_stata = 2*normal(-abs(_b[D]/_se[D]))
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar p = .
        return scalar p_stata = .
    }
    return scalar ev = `ev'
    return scalar n = `n'
    restore
end

program define curvM, rclass
    syntax , outcome(string)
    preserve
    qui keep if post == 0 & s > -10 & s < 10
    qui gen double sr = round(s)
    collapse (mean) m = `outcome', by(D sr)
    qui gen double g = min(max(m, 1e-4), 1 - 1e-4)
    qui gen double lg = ln(g/(1 - g))
    qui gen double sr2 = sr^2
    local M = 0
    foreach side in 0 1 {
        qui count if D == `side'
        if r(N) < 4 local M = .
        else {
            qui regress lg sr sr2 if D == `side'
            if `M' < . local M = max(`M', 2*abs(_b[sr2]))
        }
    }
    restore
    return scalar M = `M'
end

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

program define srdh, rclass
    syntax , h(real) post(integer) m(real) [outcome(string)]
    if "`outcome'" == "" local outcome death
    preserve
    qui keep if post == `post' & abs(s) < `h'
    mkw `h' tri
    qui gen double Ds = D*s
    qui count if `outcome' == 1
    local ev = r(N)
    qui count
    local n = r(N)
    cap qui logit `outcome' D s Ds hB hC $COVS [pw=kw], vce(cluster clust) iter(100)
    if _rc == 0 & `ev' >= 8 {
        seadj
        local b = _b[D]
        local se = _se[D]*r(adj)
        mata: st_numscalar("r_Bu", wbias("D s Ds hB hC $COVS", "kw", "s", 1))
        local B = `m'*scalar(r_Bu)
        mata: hci(`b', `se', `B')
        return scalar b = `b'
        return scalar se = `se'
        return scalar Bunit = scalar(r_Bu)
        return scalar lo = scalar(r_hlo)
        return scalar hi = scalar(r_hhi)
        return scalar p = scalar(r_hp)
    }
    else {
        return scalar b = .
        return scalar se = .
        return scalar Bunit = .
        return scalar lo = .
        return scalar hi = .
        return scalar p = .
    }
    return scalar ev = `ev'
    return scalar n = `n'
    restore
end
program define densjump, rclass
    syntax , h(real) [bw(real 0.5)]
    preserve
    qui keep if abs(s) < `h'
    qui gen double bin = floor(s/`bw')*`bw' + `bw'/2
    qui contract bin, freq(cnt)
    qui gen double lc = ln(cnt)
    qui regress lc bin if bin < 0
    local aL = _b[_cons]
    qui regress lc bin if bin > 0
    local aR = _b[_cons]
    return scalar jump = `aR' - `aL'
    restore
end

program define xout
    args dta sheet
    preserve
    use "`dta'", clear
    export excel using "$OUTX", sheet("`sheet'") sheetreplace firstrow(variables)
    restore
end

