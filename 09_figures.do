* Figures 1-3; Appendix Figures 1-6
version 15
clear all
set more off
set seed 7

global DATA  "paper_data_v9.dta"
global OUT   "figures_stata_v3"
global WHICH "F1 F2 F3 A1 A2 A3 A4 A5 A6"
if `"`0'"' != "" global WHICH `0'
global TAG   ""
global FMT   "png"
global PNGW  3000

global CUT   89.5
global SH    0.5
global COVS  age_z male cci htn sofa surgical_dept summer fall winter

global F2_H      5
global F2_H2     8
global F2_BIN    0.5
global F2_POLY   1
global F2_KERN   tri
global F2_YSCALE prob
global F2_MINN   20
global F2_YTOP   50
global F2_Y1     22
global F2_Y2     55

global C_MAIN   "181 0 66"
global C_AFTER  "64 112 117"
global C_BEFORE "148 178 181"
global F_AFTER  "207 239 234"
global F_BEFORE "230 237 237"

cap mkdir "$OUT"
cap log close
log using "$OUT/figures_stata_v3_run$TAG.log", replace text
set scheme s1color
cap graph set window fontface "Arial"
global XLSX "$OUT/figures_stata_v3$TAG.xlsx"

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
    syntax , h(real) [outcome(string) kern(string) poly(integer 1) donut(real 0) covs(string) vce(string) minside(integer 0) nosat]
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
    qui count if `outcome' == 1 & post == 1 & D == 1
    local e11 = r(N)
    qui count if `outcome' == 1 & post == 1 & D == 0
    local e10 = r(N)
    qui count
    local n = r(N)
    local b = .
    local se = .
    local p = .
    if `ev' >= 15 & `e11' >= `minside' & `e10' >= `minside' {
        cap qui logit `outcome' D s Ds post pD ps pDs `q' `hosp' `covs' [pw=kw], vce(`vce') iter(100)
        if _rc == 0 {
            seadj
            local b = _b[pD]
            local se = _se[pD]*r(adj)
            local p = 2*normal(-abs(`b'/`se'))
        }
    }
    return scalar b = `b'
    return scalar se = `se'
    return scalar p = `p'
    return scalar ev = `ev'
    return scalar evp = `evp'
    return scalar n = `n'
    restore
end

program define srd, rclass
    syntax , h(real) post(integer) [outcome(string) kern(string) poly(integer 1) donut(real 0) covs(string) vce(string) minside(integer 0)]
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
    qui count if `outcome' == 1 & D == 1
    local e1 = r(N)
    qui count if `outcome' == 1 & D == 0
    local e0 = r(N)
    qui count
    local n = r(N)
    local b = .
    local se = .
    local p = .
    if `ev' >= 8 & `e1' >= `minside' & `e0' >= `minside' {
        cap qui logit `outcome' D s Ds `q' hB hC `covs' [pw=kw], vce(`vce') iter(100)
        if _rc == 0 {
            seadj
            local b = _b[D]
            local se = _se[D]*r(adj)
            local p = 2*normal(-abs(`b'/`se'))
        }
    }
    return scalar b = `b'
    return scalar se = `se'
    return scalar p = `p'
    return scalar ev = `ev'
    return scalar n = `n'
    restore
end

program define trio, rclass
    syntax , h(real) [outcome(string) minside(integer 0) nosat *]
    if "`outcome'" == "" local outcome death
    didisc, h(`h') outcome(`outcome') minside(`minside') `sat' `options'
    foreach k in b se p ev evp n {
        local dd_`k' = r(`k')
    }
    srd, h(`h') post(1) outcome(`outcome') minside(`minside') `options'
    foreach k in b se p {
        local rd_`k' = r(`k')
    }
    srd, h(`h') post(0) outcome(`outcome') minside(`minside') `options'
    foreach k in b se p {
        local pre_`k' = r(`k')
    }
    foreach e in dd rd pre {
        return scalar `e'_or = exp(``e'_b')
        return scalar `e'_lo = exp(``e'_b' - 1.96*``e'_se')
        return scalar `e'_hi = exp(``e'_b' + 1.96*``e'_se')
        return scalar `e'_p  = ``e'_p'
        return scalar `e'_b  = ``e'_b'
        return scalar `e'_se = ``e'_se'
    }
    return scalar ev = `dd_ev'
    return scalar evp = `dd_evp'
    return scalar n = `dd_n'
end

program define poolre, rclass
    args b1 s1 b2 s2 b3 s3
    local sw = 0
    local swb = 0
    local sw2 = 0
    forvalues k = 1/3 {
        local w`k' = 1/((`s`k'')^2)
        local sw  = `sw' + `w`k''
        local swb = `swb' + `w`k''*(`b`k'')
        local sw2 = `sw2' + (`w`k'')^2
    }
    local bF = `swb'/`sw'
    local Q = 0
    forvalues k = 1/3 {
        local Q = `Q' + `w`k''*((`b`k'') - (`bF'))^2
    }
    local tau2 = max(0, (`Q' - 2)/(`sw' - `sw2'/`sw'))
    local swr = 0
    local swrb = 0
    forvalues k = 1/3 {
        local wr = 1/((`s`k'')^2 + `tau2')
        local swr = `swr' + `wr'
        local swrb = `swrb' + `wr'*(`b`k'')
    }
    return scalar bF = `bF'
    return scalar seF = 1/sqrt(`sw')
    return scalar bR = `swrb'/`swr'
    return scalar seR = 1/sqrt(`swr')
    return scalar Q = `Q'
    return scalar tau2 = `tau2'
end

program define gexp
    args nm
    foreach f in $FMT {
        if "`f'" == "png" cap noi graph export "$OUT/`nm'$TAG.png", replace width($PNGW)
        else              cap noi graph export "$OUT/`nm'$TAG.`f'", replace
    }
end

program define xsheet
    args sheet
    cap noi export excel using "$XLSX", sheet("`sheet'") sheetreplace firstrow(variables)
end

program define useest
    args pre
    cap drop est lo hi txt
    gen double est = `pre'_or
    gen double lo = `pre'_lo
    gen double hi = `pre'_hi
    gen str60 txt = subinstr(string(est, "%9.2f") + " (" + string(lo, "%9.2f") + "–" + string(hi, "%9.2f") + ")", ".", "·", .) if est < .
    replace txt = "not estimable" if est >= . & lab != ""
end

program define forest
    syntax , name(string) title(string) xmin(real) xmax(real) xlab(string asis) color(string) [textcol noylab xtitle(string) rmargin(real 3) head(string)]
    tempvar loc hic xt
    qui gen double `loc' = max(lo, `xmin') if lo < . & hi < .
    qui gen double `hic' = min(hi, `xmax') if lo < . & hi < .
    qui gen double `xt' = `xmax'
    qui su row
    local n = r(max)
    local yl
    forvalues i = 1/`=_N' {
        local l = lab[`i']
        local r = row[`i']
        if `"`l'"' != "" local yl `yl' `r' `"`l'"'
    }
    local ylopt ylabel(`yl', angle(0) labsize(small) noticks nogrid)
    if "`ylab'" != "" local ylopt ylabel(1(1)`n', nolabels noticks nogrid)
    if "`head'" == "" local head OR (95% CI)
    local tx
    local hd
    if "`textcol'" != "" {
        local tx (scatter row `xt' if txt != "", msymbol(none) mlabel(txt) mlabposition(3) mlabgap(*4) mlabsize(small) mlabcolor(black))
        local hd text(0.1 `xmax' "{bf:`head'}", placement(e) size(small) margin(l=2))
    }
    twoway (rspike `loc' `hic' row, horizontal lcolor("`color'") lwidth(medthin)) ///
           (scatter row est if inrange(est, `xmin', `xmax'), mcolor("`color'") msymbol(O) msize(small)) ///
           `tx', ///
           xline(1, lpattern(dash) lcolor(gs7) lwidth(thin)) ///
           xscale(log range(`xmin' `xmax')) xlabel(`xlab', labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           yscale(reverse range(0 `=`n' + 0.6') noline) `ylopt' ///
           ytitle("") xtitle("`xtitle'", size(small)) ///
           title("{bf:`title'}", position(11) size(medsmall) color(black)) ///
           legend(off) `hd' ///
           graphregion(color(white) margin(r=`rmargin')) plotregion(lstyle(none)) ///
           name(`name', replace) nodraw
end

program define f2fit
    syntax , post(integer) h(real) saving(string)
    tempname PF
    postfile `PF' double(h post side x p lo hi eta eta_lo eta_hi b_cons b_s b_s2) using "`saving'", replace
    forvalues side = 0/1 {
        preserve
        qui keep if post == `post' & D == `side' & abs(s) < `h'
        mkw `h' $F2_KERN
        local rhs s
        if $F2_POLY == 2 {
            qui gen double s2 = s^2
            local rhs s s2
        }
        qui logit death `rhs' [pw=kw], vce(cluster clust) iter(100)
        seadj
        matrix b = e(b)
        matrix V = e(V)*(r(adj)^2)
        local bc = _b[_cons]
        local bs = _b[s]
        local bs2 = .
        if $F2_POLY == 2 local bs2 = _b[s2]
        restore
        forvalues i = 0/79 {
            if `side' == 0 local x = -`h' + `i'*(`h' - 0.001)/79
            else           local x = `i'*(`h' - 0.001)/79
            local x2 = (`x')^2
            if $F2_POLY == 2 matrix g = (`x', `x2', 1)
            else             matrix g = (`x', 1)
            matrix e1 = g*b'
            matrix v1 = g*V*g'
            local eta = e1[1,1]
            local se = sqrt(v1[1,1])
            post `PF' (`h') (`post') (`side') (`x' + $CUT + $SH) (100*invlogit(`eta')) (100*invlogit(`eta' - 1.96*`se')) (100*invlogit(`eta' + 1.96*`se')) ///
                (`eta') (`eta' - 1.96*`se') (`eta' + 1.96*`se') (`bc') (`bs') (`bs2')
        }
    }
    postclose `PF'
end

use "$DATA", clear
keep if rv_dec < . & h_stay > 24 & h_stay < .
gen double s    = rv_dec - $CUT
gen byte   D    = (s >= 0)
gen byte   post = treat
gen int    disp = round(rv_dec)
count
tempfile base win
save `base'
global BASE "`base'"
keep if abs(s) < 10
save `win'
global WIN "`win'"

if strpos(" $WHICH ", " F1 ") {
    clear
    set obs 2
    gen double x = cond(_n == 1, 0, 24)
    gen double ylo = -0.55
    gen double yhi = 2.75
    twoway (rarea ylo yhi x, fcolor("$F_AFTER") lwidth(none)) ///
           (pci -0.55 24 3.02 24, lcolor(black) lpattern(dash) lwidth(medthick)) ///
           (pci 2 0 2 14,  lcolor("$C_BEFORE") lwidth(medthick)) ///
           (pci 1 0 1 240, lcolor("$C_AFTER") lwidth(thick)) ///
           (pci 0 0 0 168, lcolor("$C_AFTER") lwidth(medthick)) ///
           (scatteri 2 14,  msymbol(X) mcolor("$C_BEFORE") msize(vlarge) mlwidth(medthick)) ///
           (scatteri 1 240, msymbol(O) mcolor("$C_AFTER") msize(large)) ///
           (scatteri 0 168, msymbol(S) mfcolor(white) mlcolor("$C_AFTER") mlwidth(medthick) msize(large)) ///
           (pcarrowi 3.08 32 3.08 252, lcolor(gs7) mcolor(gs7) lwidth(thin)), ///
           text(2 -8 "Patient 1", placement(w) size(small)) text(1 -8 "Patient 2", placement(w) size(small)) text(0 -8 "Patient 3", placement(w) size(small)) ///
           text(4.15 12 "Score observation window" "W = 24 h", size(small)) ///
           text(3.62 160 "Outcome observation:  in-hospital death", size(small)) ///
           text(2.55 28 "Landmark  L = 24 h", placement(e) size(small)) ///
           xscale(range(-52 262)) yscale(range(-0.95 4.9) off) ///
           xlabel(0 `""Admission" "0 h""' 24 "24 h" 48 "48 h" 72 "72 h" 120 "day 5" 168 "day 7" 240 "day 10", labsize(small)) ///
           xtitle("Time from admission", size(small)) ytitle("") legend(off) ///
           title("{bf:A  Cohort and timing}", position(11) size(small) color(black)) ///
           graphregion(color(white) margin(t=0 b=0)) plotregion(lstyle(none)) name(f1a, replace) nodraw

    clear
    set obs 15
    gen int d0 = 87 + ceil(_n/3)
    gen double x = d0 - 0.5 if mod(_n, 3) == 1
    replace x = d0 + 0.5 if mod(_n, 3) == 2
    gen double ylo = 0.42 if x < .
    gen double yhi = 0.88 if x < .
    twoway (rarea ylo yhi x if d0 < 90,  cmissing(n) fcolor("$F_BEFORE") lcolor(white) lwidth(medium)) ///
           (rarea ylo yhi x if d0 >= 90, cmissing(n) fcolor("$F_AFTER") lcolor(white) lwidth(medium)), ///
           xline(89.5, lpattern(dash) lcolor(black) lwidth(medthick)) ///
           text(0.65 88 "displayed 88", size(small)) text(0.65 89 "displayed 89", size(small)) text(0.65 90 "displayed 90", size(small)) ///
           text(0.65 91 "displayed 91", size(small)) text(0.65 92 "displayed 92", size(small)) ///
           text(1.12 89.58 "Alert boundary: internal value 89·5 = start of displayed 90", placement(e) size(small)) ///
           text(0.22 89.42 "89·46 → 89, not listed", placement(w) size(vsmall)) text(0.22 89.58 "89·53 → 90, listed", placement(e) size(vsmall)) ///
           xscale(range(87.4 92.6)) yscale(range(0 1.32) off) ///
           xlabel(87.5 "87·5" 88.5 "88·5" 89.5 "89·5" 90.5 "90·5" 91.5 "91·5" 92.5 "92·5", labsize(small)) ///
           xtitle("Internal value of the AI score (continuous)", size(small)) ytitle("") legend(off) ///
           title("{bf:B  Displayed score = internal value rounded to the nearest integer; the alert follows the displayed score}", position(11) size(small) color(black)) ///
           graphregion(color(white) margin(t=0 b=0)) plotregion(lstyle(none)) name(f1b, replace) nodraw

    use "$BASE", clear
    keep if disp >= 78 & disp <= 99
    collapse (mean) m = death (count) n = death (sum) deaths = death, by(disp)
    gen double pct = 100*m
    gen double x = disp + $SH
    save "$OUT/figdata_F1C$TAG.dta", replace
    xsheet "F1C"
    gen double rx8 = cond(_n == 1, 82, 98) in 1/2
    gen double rx5 = cond(_n == 1, 85, 95) in 1/2
    gen double rlo = 1.5 in 1/2
    gen double rhi = 75 in 1/2
    twoway (rarea rlo rhi rx8, fcolor("$F_BEFORE") lwidth(none)) (rarea rlo rhi rx5, fcolor("$F_AFTER") lwidth(none)) ///
           (connected pct x, lcolor("$C_AFTER") mcolor("$C_AFTER") msize(small) lwidth(medium)), ///
           xline(90, lpattern(dash) lcolor(black) lwidth(medthick)) ///
           text(55 90.25 "Alert threshold (90)", placement(e) size(small)) ///
           yscale(log range(1.5 75)) ylabel(2 5 10 20 50, angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           xscale(range(77.4 100.6)) xlabel(78 82 85 90 95 98, labsize(small)) ///
           xtitle("Highest AI score within 24 h of admission (displayed integer)", size(small)) ytitle("In-hospital mortality (%), log scale", size(small)) ///
           legend(order(2 "Primary window  ±5 points (85 to 94)" 1 "Secondary window  ±8 points (82 to 97)") cols(1) ring(0) position(11) region(lstyle(none) fcolor(none)) size(small) symxsize(*0.4)) ///
           title("{bf:C  Mortality by assignment value and the analysis windows}", position(11) size(small) color(black)) ///
           graphregion(color(white) margin(t=0 b=0)) plotregion(lstyle(none)) name(f1c, replace) nodraw

    clear
    set obs 1
    gen x = 0
    twoway (function y = 5 + 0.55*x,                   range(-5 0) lcolor("$C_BEFORE") lwidth(thick)) ///
           (function y = 5 + 0.55*x + 1.6,             range(0 5)  lcolor("$C_BEFORE") lwidth(thick)) ///
           (function y = 5 + 0.55*x - 0.9,             range(-5 0) lcolor("$C_AFTER") lwidth(thick)) ///
           (function y = 5 + 0.55*x - 0.9 + 1.6 - 3.0, range(0 5)  lcolor("$C_AFTER") lwidth(thick)) ///
           (pci 5 0 6.6 0, lcolor("$C_BEFORE") lpattern(dot)) (pci 4.1 0 2.7 0, lcolor("$C_AFTER") lpattern(dot)), ///
           xline(0, lpattern(dash) lcolor(black) lwidth(medthick)) ///
           xscale(range(-5.4 5.4)) yscale(range(0.5 9.7) noline) ylabel(none) xlabel(-5 "−5" 0 "Threshold" 5 "+5", labsize(small)) ///
           xtitle("Assignment value relative to the alert threshold (points)", size(small)) ytitle("Mortality (schematic)", size(small)) ///
           legend(order(2 "Before implementation (no alert): discontinuity = confounding component" 4 "After implementation: discontinuity = confounding component + effect of the alert") ///
                  cols(1) position(6) region(lstyle(none)) size(small)) ///
           title("{bf:D  Identification: the difference-in-discontinuities subtracts the pre-implementation discontinuity}", position(11) size(small) color(black)) ///
           graphregion(color(white) margin(t=0 b=0)) plotregion(lstyle(none)) name(f1d, replace) nodraw

    graph combine f1a f1b f1c f1d, cols(1) xsize(7.4) ysize(10.3) iscale(*0.75) imargin(zero) graphregion(color(white)) name(F1, replace)
    gexp "F1_design"
}

if strpos(" $WHICH ", " F2 ") {
    use "$BASE", clear
    gen double bin = floor((rv_dec - $CUT)/$F2_BIN)*$F2_BIN + $CUT + $F2_BIN/2 + $SH
    gen byte one = 1
    collapse (mean) m = death (sum) n = one (sum) deaths = death, by(post bin)
    gen double pct = 100*m
    gen byte inwin  = (bin - $SH > $CUT - $F2_H)  & (bin - $SH < $CUT + $F2_H)
    gen byte inwin2 = (bin - $SH > $CUT - $F2_H2) & (bin - $SH < $CUT + $F2_H2)
    gen byte inall = (bin - $SH > 39.5) & (bin - $SH < 100.5) & n >= $F2_MINN
    save "$OUT/figdata_F2_bins$TAG.dta", replace
    xsheet "F2_bins"
    use "$BASE", clear
    local k = 0
    foreach hh in $F2_H $F2_H2 {
        foreach t in 0 1 {
            local ++k
            f2fit, post(`t') h(`hh') saving("$OUT/figdata_F2_part`k'$TAG.dta")
        }
    }
    use "$OUT/figdata_F2_part1$TAG.dta", clear
    forvalues j = 2/4 {
        append using "$OUT/figdata_F2_part`j'$TAG.dta"
    }
    save "$OUT/figdata_F2_fit$TAG.dta", replace
    xsheet "F2_fit"
    forvalues j = 1/4 {
        erase "$OUT/figdata_F2_part`j'$TAG.dta"
    }
    foreach hh in $F2_H $F2_H2 {
        foreach t in 0 1 {
            foreach sd in 0 1 {
                qui su b_s if h == `hh' & post == `t' & side == `sd'
            }
        }
    }
    append using "$OUT/figdata_F2_bins$TAG.dta"
    local binlab = subinstr("$F2_BIN", ".", "·", .)
    local xcut = $CUT + $SH
    local gopt angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)
    if "$F2_YSCALE" == "logit" {
        replace pct = logit(pct/100)
        replace p = eta
        replace lo = eta_lo
        replace hi = eta_hi
        local ytop ylabel(`=logit(.002)' "0·2" `=logit(.005)' "0·5" `=logit(.01)' "1" `=logit(.02)' "2" `=logit(.05)' "5" `=logit(.1)' "10" `=logit(.2)' "20" `=logit(.4)' "40", `gopt')
        local topif
        local ytl "In-hospital mortality (%), log-odds scale"
    }
    else {
        local ytop yscale(range(0 $F2_YTOP)) ylabel(0(10)$F2_YTOP, `gopt')
        local topif & pct <= $F2_YTOP
        local ytl "In-hospital mortality (%)"
    }
    local w = 0
    foreach hh in $F2_H $F2_H2 {
        local ++w
        local inw = cond(`w' == 1, "inwin", "inwin2")
        qui su lo if h == `hh'
        local mn = r(min)
        qui su hi if h == `hh'
        local mx = r(max)
        qui su pct if `inw' == 1
        local mn = min(`mn', r(min))
        local mx = max(`mx', r(max))
        if "$F2_YSCALE" == "logit" {
            local yl
            foreach v in 0.5 1 2 3 5 7 10 15 20 30 40 50 {
                if logit(`v'/100) >= `mn' & logit(`v'/100) <= `mx' local yl `yl' `=logit(`v'/100)' "`=subinstr("`v'", ".", "·", .)'"
            }
            local ybot`w' yscale(range(`mn' `mx')) ylabel(`yl', `gopt')
        }
        else {
            local set = cond(`w' == 1, $F2_Y1, $F2_Y2)
            local top = max(`set', ceil(`mx'))
            local ystep = cond(`top' <= 25, 5, 10)
            local ybot`w' yscale(range(0 `top')) ylabel(0(`ystep')`top', `gopt')
        }
    }
    foreach t in 0 1 {
        local col = cond(`t' == 0, "$C_BEFORE", "$C_AFTER")
        local fil = cond(`t' == 0, "$F_BEFORE", "$F_AFTER")
        local per = cond(`t' == 0, "Before implementation", "After implementation")
        local la = cond(`t' == 0, "A", "B")
        local yt = cond(`t' == 0, "`ytl'", "")
        twoway (scatter pct bin if post == `t' & inall == 1 `topif', mcolor("`col'") msymbol(O) msize(vsmall)) ///
               (line p x if h == $F2_H & post == `t' & side == 0, lcolor("`col'") lwidth(medthick)) (line p x if h == $F2_H & post == `t' & side == 1, lcolor("`col'") lwidth(medthick)), ///
               xline(`xcut', lpattern(dash) lcolor(gs7) lwidth(thin)) xscale(range(39 101)) xlabel(40(10)100, labsize(small)) `ytop' ///
               xtitle("AI score (`binlab'-point bins)", size(small)) ytitle("`yt'", size(small)) legend(off) ///
               title("{bf:`la'  `per', all scores}", position(11) size(medsmall) color(black)) ///
               graphregion(color(white)) plotregion(lstyle(none)) name(f2_r0_`t', replace) nodraw
        local w = 0
        foreach hh in $F2_H $F2_H2 {
            local ++w
            local inw = cond(`w' == 1, "inwin", "inwin2")
            local hlab = subinstr("`hh'", ".", "·", .)
            local xl = `xcut' - `hh'
            local xr = `xcut' + `hh'
            local step = cond(`hh' <= 6, 1, 2)
            local let : word `=2*`w' + `t' + 1' of A B C D E F
            twoway (rarea lo hi x if h == `hh' & post == `t' & side == 0, fcolor("`fil'") lwidth(none)) (rarea lo hi x if h == `hh' & post == `t' & side == 1, fcolor("`fil'") lwidth(none)) ///
                   (line p x if h == `hh' & post == `t' & side == 0, lcolor("`col'") lwidth(medthick)) (line p x if h == `hh' & post == `t' & side == 1, lcolor("`col'") lwidth(medthick)) ///
                   (scatter pct bin if post == `t' & `inw' == 1, mcolor("`col'") msymbol(O) msize(small)), ///
                   xline(`xcut', lpattern(dash) lcolor(gs7) lwidth(thin)) xscale(range(`=`xl' - 0.2' `=`xr' + 0.2')) xlabel(`xl'(`step')`xr', labsize(small)) `ybot`w'' ///
                   xtitle("AI score (`binlab'-point bins)", size(small)) ytitle("`yt'", size(small)) legend(off) ///
                   title("{bf:`let'  `per', ±`hlab' window}", position(11) size(medsmall) color(black)) ///
                   graphregion(color(white)) plotregion(lstyle(none)) name(f2_r`w'_`t', replace) nodraw
        }
    }
    graph combine f2_r0_0 f2_r0_1 f2_r1_0 f2_r1_1 f2_r2_0 f2_r2_1, cols(2) xsize(7.2) ysize(7.4) iscale(*0.8) imargin(small) graphregion(color(white)) name(F2, replace)
    gexp "F2_score_mortality"
}

if strpos(" $WHICH ", " F3 ") {
    use "$WIN", clear
    tempname M
    postfile `M' double(h n ev evp dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_lo pre_hi pre_p) using "$OUT/figdata_F3$TAG.dta", replace
    forvalues h = 2/10 {
        trio, h(`h')
        post `M' (`h') (r(n)) (r(ev)) (r(evp)) (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p)) (r(pre_or)) (r(pre_lo)) (r(pre_hi)) (r(pre_p))
    }
    postclose `M'
    use "$OUT/figdata_F3$TAG.dta", clear
    xsheet "F3"
    gen double h1 = h - 0.16
    gen double h3 = h + 0.16
    foreach v in dd rd pre {
        gen double `v'_loc = max(`v'_lo, 0.02)
        gen double `v'_hic = min(`v'_hi, 12)
    }
    twoway (rcap dd_loc dd_hic h1, lcolor("$C_MAIN") lwidth(medium)) (scatter dd_or h1, mcolor("$C_MAIN") msymbol(O) msize(medium)) ///
           (rcap rd_loc rd_hic h, lcolor("$C_AFTER") lwidth(medthin)) (scatter rd_or h, mcolor("$C_AFTER") msymbol(S) msize(medium)) ///
           (rcap pre_loc pre_hic h3, lcolor("$C_BEFORE") lwidth(medthin)) (scatter pre_or h3, mcolor("$C_BEFORE") msymbol(D) msize(medium)), ///
           yline(1, lpattern(dash) lcolor(gs7) lwidth(thin)) yscale(log range(0.02 12)) ///
           ylabel(0.02 "0·02" 0.05 "0·05" 0.1 "0·1" 0.25 "0·25" 0.5 "0·5" 1 "1" 2 "2" 5 "5" 10 "10", angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           xscale(range(1.4 10.6)) xlabel(2(1)10, labsize(small)) ///
           xtitle("Window (points on each side of the threshold)", size(small)) ytitle("Odds ratio for in-hospital death at the threshold", size(small)) ///
           legend(order(2 "Difference-in-discontinuities (primary)" 4 "Single regression discontinuity, after implementation" 6 "Single regression discontinuity, before implementation") ///
                  cols(1) position(6) region(lstyle(none)) size(small)) ///
           graphregion(color(white)) plotregion(lstyle(none)) xsize(6.8) ysize(4.4) name(F3, replace)
    gexp "F3_bandwidth"
}

if strpos(" $WHICH ", " A1 ") {
    use "$BASE", clear
    keep if post == 1 & disp >= 80 & disp <= 100
    collapse (mean) m = dcars_over90 (count) n = dcars_over90, by(disp)
    gen double pct = 100*m
    gen double x = disp + $SH
    save "$OUT/figdata_A1_A$TAG.dta", replace
    xsheet "A1_A_listentry"
    twoway (connected pct x if disp < 90, lcolor("$C_AFTER") mcolor("$C_AFTER") msize(small)) (connected pct x if disp >= 90, lcolor("$C_AFTER") mcolor("$C_AFTER") msize(small)), ///
           xline(90, lpattern(dash) lcolor(gs7) lwidth(thin)) yscale(range(-3 103)) ylabel(0(25)100, angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           xscale(range(79.5 101)) xlabel(80(5)100, labsize(small)) xtitle("") ytitle("Ever on the high-risk list (%)", size(small)) legend(off) ///
           title("{bf:A  List entry, after implementation}", position(11) size(medsmall) color(black)) graphregion(color(white)) plotregion(lstyle(none)) name(a1_a, replace) nodraw
    use "$BASE", clear
    keep if disp >= 80 & disp <= 99
    gen byte one = 1
    collapse (sum) cnt = one, by(post disp)
    gen double x = disp + $SH
    gen double lc = ln(max(cnt, 1))
    gen double xs = disp - $CUT
    gen double fit = .
    foreach p in 0 1 {
        foreach sd in 0 1 {
            local cond = cond(`sd' == 1, "disp >= 90 & disp <= 94", "disp >= 85 & disp <= 89")
            qui regress lc xs if post == `p' & `cond'
            tempvar f
            qui predict double `f' if e(sample)
            qui replace fit = exp(`f') if e(sample)
            drop `f'
        }
    }
    save "$OUT/figdata_A1_B$TAG.dta", replace
    xsheet "A1_B_density"
    qui su cnt
    local ymax = r(max)*1.8
    local ylb
    foreach v in 10 20 50 100 200 500 1000 {
        if `v' < `ymax' local ylb `ylb' `v'
    }
    twoway (scatter cnt x if post == 0, mcolor("$C_BEFORE") msize(small)) (scatter cnt x if post == 1, mcolor("$C_AFTER") msize(small)) ///
           (line fit x if post == 0 & disp < 90, lcolor("$C_BEFORE") lwidth(medthick)) (line fit x if post == 0 & disp >= 90, lcolor("$C_BEFORE") lwidth(medthick)) ///
           (line fit x if post == 1 & disp < 90, lcolor("$C_AFTER") lwidth(medthick)) (line fit x if post == 1 & disp >= 90, lcolor("$C_AFTER") lwidth(medthick)), ///
           xline(90, lpattern(dash) lcolor(gs7) lwidth(thin)) yscale(log range(8 `ymax')) ylabel(`ylb', angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           xscale(range(79.5 101)) xlabel(80(5)100, labsize(small)) xtitle("") ytitle("Admissions per displayed score", size(small)) ///
           legend(order(1 "Before" 2 "After") cols(1) ring(0) position(7) region(lstyle(none) fcolor(none)) size(small)) ///
           title("{bf:B  Density of the assignment variable}", position(11) size(medsmall) color(black)) graphregion(color(white)) plotregion(lstyle(none)) name(a1_b, replace) nodraw
    use "$BASE", clear
    qui logit death $COVS hB hC if abs(s) >= 5, iter(100)
    predict double risk, pr
    replace risk = 100*risk
    keep if abs(s) < 10
    tempfile a1
    save `a1'
    local T_age      "C  Age (years)"
    local T_cci      "D  Charlson comorbidity index"
    local T_t_of_max "E  Time of maximum score (h after admission)"
    local T_risk     "F  Mortality predicted from covariates (%)"
    tempname DS
    postfile `DS' str12 variable str8 period double(disc se p) using "$OUT/figdata_A1_disc$TAG.dta", replace
    foreach v in age cci t_of_max risk {
        use `a1', clear
        foreach p in 0 1 {
            foreach sd in 0 1 {
                qui regress `v' s [aw = max(0, 1 - abs(s)/5)] if post == `p' & D == `sd' & abs(s) < 5
                local a`p'`sd' = _b[_cons]
                local b`p'`sd' = _b[s]
            }
            preserve
            qui keep if post == `p' & abs(s) < 5
            mkw 5 tri
            qui gen double Ds = D*s
            qui regress `v' D s Ds hB hC [pw=kw], vce(cluster clust)
            seadj
            post `DS' ("`v'") (cond(`p' == 0, "Before", "After")) (_b[D]) (_se[D]*r(adj)) (2*normal(-abs(_b[D]/(_se[D]*r(adj)))))
            restore
        }
        collapse (mean) m = `v', by(post disp)
        gen double x = disp + $SH
        gen str12 variable = "`v'"
        foreach p in 0 1 {
            foreach sd in 0 1 {
                gen double a`p'`sd' = `a`p'`sd''
                gen double b`p'`sd' = `b`p'`sd''
            }
        }
        save "$OUT/figdata_A1_`v'$TAG.dta", replace
        local xt = cond("`v'" == "t_of_max" | "`v'" == "risk", "Displayed AI score", "")
        twoway (scatter m x if post == 0, mcolor("$C_BEFORE") msize(small)) (scatter m x if post == 1, mcolor("$C_AFTER") msize(small)) ///
               (function y = `a00' + `b00'*(x - 90), range(85 90) lcolor("$C_BEFORE") lwidth(medthick)) (function y = `a01' + `b01'*(x - 90), range(90 95) lcolor("$C_BEFORE") lwidth(medthick)) ///
               (function y = `a10' + `b10'*(x - 90), range(85 90) lcolor("$C_AFTER") lwidth(medthick))  (function y = `a11' + `b11'*(x - 90), range(90 95) lcolor("$C_AFTER") lwidth(medthick)), ///
               xline(90, lpattern(dash) lcolor(gs7) lwidth(thin)) ylabel(, angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
               xscale(range(79.5 101)) xlabel(80(5)100, labsize(small)) xtitle("`xt'", size(small)) ytitle("") legend(off) ///
               title("{bf:`T_`v''}", position(11) size(medsmall) color(black)) graphregion(color(white)) plotregion(lstyle(none)) name(a1_`v', replace) nodraw
    }
    postclose `DS'
    use "$OUT/figdata_A1_age$TAG.dta", clear
    foreach v in cci t_of_max risk {
        append using "$OUT/figdata_A1_`v'$TAG.dta"
    }
    save "$OUT/figdata_A1_CF$TAG.dta", replace
    xsheet "A1_CF_means"
    foreach v in age cci t_of_max risk {
        erase "$OUT/figdata_A1_`v'$TAG.dta"
    }
    use "$OUT/figdata_A1_disc$TAG.dta", clear
    xsheet "A1_CF_disc"
    list, noobs clean
    graph combine a1_a a1_b a1_age a1_cci a1_t_of_max a1_risk, cols(2) xsize(7.2) ysize(8.6) iscale(*0.85) imargin(small) graphregion(color(white)) name(A1, replace)
    gexp "FA1_validity"
}

if strpos(" $WHICH ", " A2 ") {
    use "$WIN", clear
    local s1
    local s2 covs(none)
    local s3 poly(2)
    local s4 kern(uni)
    local s5 kern(epa)
    local s6 donut(0.5)
    local s7 donut(1)
    local s8 vce(cluster disp)
    local s9 vce(robust)
    local n1 "Reference specification"
    local n2 "No covariates"
    local n3 "Local quadratic"
    local n4 "Uniform kernel"
    local n5 "Epanechnikov kernel"
    local n6 "Donut, 0·5 point"
    local n7 "Donut, 1 point"
    local n8 "Clustered by displayed score"
    local n9 "No clustering (HC)"
    tempname SP
    postfile `SP' int row str60 lab double(dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "$OUT/figdata_A2$TAG.dta", replace
    forvalues i = 1/9 {
        trio, h(5) `s`i''
        post `SP' (`i') ("`n`i''") (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p))
    }
    didisc, h(5) nosat
    post `SP' (10) ("No hospital-specific terms") (exp(r(b))) (exp(r(b) - 1.96*r(se))) (exp(r(b) + 1.96*r(se))) (r(p)) (.) (.) (.) (.)
    postclose `SP'
    use "$OUT/figdata_A2$TAG.dta", clear
    xsheet "A2_spec"
    list row lab dd_or dd_p rd_or rd_p, noobs clean
    preserve
    keep if row <= 9
    useest rd
    forest, name(a2_a) title("A  Post-implementation discontinuity") xmin(0.08) xmax(3) xlab(0.1 "0·1" 0.25 "0·25" 0.5 "0·5" 1 "1" 2 "2") color("$C_AFTER") textcol xtitle("Odds ratio for in-hospital death") rmargin(30)
    restore
    useest dd
    forest, name(a2_b) title("B  Difference-in-discontinuities") xmin(0.02) xmax(6) xlab(0.05 "0·05" 0.1 "0·1" 0.25 "0·25" 0.5 "0·5" 1 "1" 2 "2" 5 "5") color("$C_MAIN") textcol xtitle("Odds ratio for in-hospital death") rmargin(30)
    graph combine a2_a a2_b, cols(2) xsize(10.5) ysize(4.4) iscale(*1.1) graphregion(color(white)) name(A2, replace)
    gexp "FA2_spec"
}

if strpos(" $WHICH ", " A3 ") {
    use "$WIN", clear
    gen byte y_d23  = (death == 1 & h_death <= 72)
    gen byte y_d47  = (death == 1 & h_death > 72 & h_death <= 168)
    gen byte y_d830 = (death == 1 & h_death > 168 & h_death <= 720)
    gen byte y_d30  = (death == 1 & h_death > 720 & h_death < .)
    gen byte y_all  = (death == 1 & h_death < .)
    tempfile a3
    save `a3'
    local o1 y_d23
    local o2 y_d47
    local o3 y_d830
    local o4 y_d30
    local o5 y_all
    local l1 "Death on days 2–3"
    local l2 "Death on days 4–7"
    local l3 "Death on days 8–30"
    local l4 "Death after day 30"
    local l5 "All in-hospital death"
    local c7  "age >= 75"
    local c8  "age < 75"
    local c9  "male == 1"
    local c10 "male == 0"
    local c11 "cci >= 2"
    local c12 "cci < 2"
    local c13 "sepsis == 1"
    local c14 "sepsis == 0"
    local c15 "rv_news >= 5 & rv_news < ."
    local c16 "!(rv_news >= 5 & rv_news < .)"
    local l7  "Age ≥75 years"
    local l8  "Age <75 years"
    local l9  "Male"
    local l10 "Female"
    local l11 "Charlson index ≥2"
    local l12 "Charlson index <2"
    local l13 "Sepsis"
    local l14 "No sepsis"
    local l15 "NEWS ≥5"
    local l16 "NEWS <5"
    tempname TS
    postfile `TS' int row str40 lab double(ev evp n dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_p) using "$OUT/figdata_A3$TAG.dta", replace
    forvalues i = 1/5 {
        trio, h(5) outcome(`o`i'') minside(2)
        post `TS' (`i') ("`l`i''") (r(ev)) (r(evp)) (r(n)) (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p)) (r(pre_or)) (r(pre_p))
    }
    post `TS' (6) ("") (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
    forvalues i = 7/16 {
        use `a3', clear
        qui keep if `c`i''
        trio, h(5) minside(2)
        post `TS' (`i') ("`l`i''") (r(ev)) (r(evp)) (r(n)) (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p)) (r(pre_or)) (r(pre_p))
    }
    postclose `TS'
    use "$OUT/figdata_A3$TAG.dta", clear
    xsheet "A3_timing_subgroup"
    list row lab ev dd_or dd_lo dd_hi dd_p rd_or rd_p, noobs clean
    useest dd
    forest, name(a3_a) title("A  Difference-in-discontinuities") xmin(0.01) xmax(30) xlab(0.01 "0·01" 0.1 "0·1" 1 "1" 10 "10") color("$C_MAIN") textcol xtitle("Odds ratio for in-hospital death") rmargin(30)
    useest rd
    forest, name(a3_b) title("B  Post-implementation discontinuity") xmin(0.01) xmax(30) xlab(0.01 "0·01" 0.1 "0·1" 1 "1" 10 "10") color("$C_AFTER") textcol xtitle("Odds ratio for in-hospital death") rmargin(30)
    graph combine a3_a a3_b, cols(2) xsize(10.5) ysize(6.2) iscale(*1.0) graphregion(color(white)) name(A3, replace)
    gexp "FA3_timing_subgroup"
}

if strpos(" $WHICH ", " A4 ") {
    use "$WIN", clear
    tempname HO
    postfile `HO' int row str30 lab double(ev dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p pre_or pre_lo pre_hi pre_p) using "$OUT/figdata_A4$TAG.dta", replace
    forvalues k = 1/3 {
        preserve
        qui keep if hosp_num == `k'
        trio, h(5) nosat
        foreach e in dd rd pre {
            local `e'_b`k' = r(`e'_b)
            local `e'_s`k' = r(`e'_se)
        }
        post `HO' (`k') ("Hospital `: word `k' of A B C'") (r(ev)) (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p)) (r(pre_or)) (r(pre_lo)) (r(pre_hi)) (r(pre_p))
        restore
    }
    foreach e in dd rd pre {
        poolre ``e'_b1' ``e'_s1' ``e'_b2' ``e'_s2' ``e'_b3' ``e'_s3'
        local `e'_bF = r(bF)
        local `e'_sF = r(seF)
        local `e'_bR = r(bR)
        local `e'_sR = r(seR)
    }
    foreach m in F R {
        local rw = cond("`m'" == "F", 4, 5)
        local nm = cond("`m'" == "F", "Pooled, fixed effect", "Pooled, random effects")
        post `HO' (`rw') ("`nm'") (.) ///
            (exp(`dd_b`m'')) (exp(`dd_b`m'' - 1.96*`dd_s`m'')) (exp(`dd_b`m'' + 1.96*`dd_s`m'')) (2*normal(-abs(`dd_b`m''/`dd_s`m''))) ///
            (exp(`rd_b`m'')) (exp(`rd_b`m'' - 1.96*`rd_s`m'')) (exp(`rd_b`m'' + 1.96*`rd_s`m'')) (2*normal(-abs(`rd_b`m''/`rd_s`m''))) ///
            (exp(`pre_b`m'')) (exp(`pre_b`m'' - 1.96*`pre_s`m'')) (exp(`pre_b`m'' + 1.96*`pre_s`m'')) (2*normal(-abs(`pre_b`m''/`pre_s`m'')))
    }
    forvalues k = 1/3 {
        preserve
        qui drop if hosp_num == `k'
        trio, h(5) minside(2)
        post `HO' (`=5 + `k'') ("Omitting hospital `: word `k' of A B C'") (r(ev)) (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p)) (r(pre_or)) (r(pre_lo)) (r(pre_hi)) (r(pre_p))
        restore
    }
    postclose `HO'
    use "$OUT/figdata_A4$TAG.dta", clear
    xsheet "A4_hospital"
    list row lab ev dd_or dd_p rd_or rd_p pre_or pre_p, noobs clean
    useest dd
    forest, name(a4_a) title("A  Difference-in-discontinuities") xmin(0.005) xmax(20) xlab(0.01 "0·01" 0.1 "0·1" 1 "1" 10 "10") color("$C_MAIN") xtitle("Odds ratio")
    useest rd
    forest, name(a4_b) title("B  Post-implementation") xmin(0.05) xmax(5) xlab(0.1 "0·1" 0.5 "0·5" 2 "2") color("$C_AFTER") xtitle("Odds ratio") noylab
    useest pre
    forest, name(a4_c) title("C  Pre-implementation") xmin(0.05) xmax(60) xlab(0.1 "0·1" 1 "1" 10 "10") color("$C_BEFORE") xtitle("Odds ratio") noylab
    graph combine a4_a a4_b a4_c, cols(3) xsize(8.6) ysize(3.5) iscale(*1.25) graphregion(color(white)) name(A4, replace)
    gexp "FA4_hospital"
}

if strpos(" $WHICH ", " A5 ") {
    use "$WIN", clear
    local S1 "0"
    local S2 "abs(s) < 5 & post == 1 & D == 1 & dnr_case == 1 & death == 0"
    local S3 "abs(s) < 5 & dnr_case == 1 & death == 0"
    local S4 "tip_rand12 == 1"
    local L1 "Observed data"
    local L2 "DNR order during the stay, above the threshold, after implementation"
    local L3 "DNR order during the stay, all four groups"
    local L4 "Survivors at random, above the threshold, after implementation"
    tempname TP
    postfile `TP' int row str90 lab double(n_added dd_or dd_lo dd_hi dd_p rd_or rd_lo rd_hi rd_p) using "$OUT/figdata_A5$TAG.dta", replace
    forvalues i = 1/4 {
        preserve
        qui count if `S`i''
        local nadd = r(N)
        qui replace death = 1 if `S`i''
        trio, h(5)
        post `TP' (`i') ("`L`i''  (+`nadd')") (`nadd') (r(dd_or)) (r(dd_lo)) (r(dd_hi)) (r(dd_p)) (r(rd_or)) (r(rd_lo)) (r(rd_hi)) (r(rd_p))
        restore
    }
    postclose `TP'
    use "$OUT/figdata_A5$TAG.dta", clear
    xsheet "A5_tipping"
    list row lab dd_or dd_lo dd_hi dd_p rd_or rd_p, noobs clean
    gen double r_dd = row - 0.15
    gen double r_rd = row + 0.17
    gen double xt = 1.7
    gen str12 ptxt = cond(dd_p < 0.0001, "p<0.0001", "p=" + cond(dd_p >= 0.0995, string(dd_p, "%4.2f"), cond(dd_p >= 0.00995, string(dd_p, "%5.3f"), string(dd_p, "%6.4f"))))
    gen str80 txt = subinstr(string(dd_or, "%9.2f") + " (" + string(dd_lo, "%9.2f") + "–" + string(dd_hi, "%9.2f") + ")   " + ptxt, ".", "·", .)
    foreach v in dd rd {
        gen double `v'_loc = max(`v'_lo, 0.05)
        gen double `v'_hic = min(`v'_hi, 1.7)
    }
    local yl
    forvalues i = 1/`=_N' {
        local l = lab[`i']
        local yl `yl' `i' `"`l'"'
    }
    twoway (rspike dd_loc dd_hic r_dd, horizontal lcolor("$C_MAIN") lwidth(medthin)) (scatter r_dd dd_or, mcolor("$C_MAIN") msymbol(O) msize(small)) ///
           (rspike rd_loc rd_hic r_rd, horizontal lcolor("$C_AFTER") lwidth(thin)) (scatter r_rd rd_or, msymbol(S) mfcolor(white) mlcolor("$C_AFTER") msize(small)) ///
           (scatter r_dd xt, msymbol(none) mlabel(txt) mlabposition(3) mlabgap(*4) mlabsize(small) mlabcolor(black)), ///
           xline(1, lpattern(dash) lcolor(gs7) lwidth(thin)) xscale(log range(0.05 1.7)) xlabel(0.1 "0·1" 0.25 "0·25" 0.5 "0·5" 1 "1", labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
           yscale(reverse range(0 4.6) noline) ylabel(`yl', angle(0) labsize(small) noticks nogrid) ytitle("") xtitle("Odds ratio for in-hospital death", size(small)) ///
           text(0.1 1.7 "{bf:Difference-in-discontinuities:}" "{bf:OR (95% CI)}", placement(e) justification(left) size(small) margin(l=2)) ///
           legend(order(2 "Difference-in-discontinuities" 4 "Post-implementation discontinuity") rows(1) position(6) region(lstyle(none)) size(small)) ///
           graphregion(color(white) margin(r=40)) plotregion(lstyle(none)) xsize(9.4) ysize(3.5) name(A5, replace)
    gexp "FA5_tipping"
}

if strpos(" $WHICH ", " A6 ") {
    use "$DATA", clear
    gen double rv_w24 = rv_dec
    gen byte post = treat
    keep if rv_dec < .
    keep adm_id hosp_num hB hC clust treat post death h_death h_stay rv_w6 rv_w12 rv_w24 rv_w48 $COVS
    tempfile full
    save `full'
    tempname GR
    postfile `GR' double(W L Hdays h n ev dd_or dd_p rd_or rd_p pre_or pre_p) using "$OUT/figdata_A6$TAG.dta", replace
    foreach W in 6 12 24 48 {
        foreach L in 24 48 72 120 {
            if `L' >= `W' {
                use `full', clear
                qui keep if rv_w`W' < . & h_stay > `L' & h_stay < .
                gen double s = rv_w`W' - $CUT
                qui keep if abs(s) < 8
                gen byte D = (s >= 0)
                foreach Hd in 7 14 30 1000000 {
                    cap drop yy
                    gen byte yy = (death == 1 & h_death <= `L' + `Hd'*24)
                    foreach h in 5 8 {
                        trio, h(`h') outcome(yy) minside(2)
                        post `GR' (`W') (`L') (`Hd') (`h') (r(n)) (r(ev)) (r(dd_or)) (r(dd_p)) (r(rd_or)) (r(rd_p)) (r(pre_or)) (r(pre_p))
                    }
                }
            }
        }
    }
    postclose `GR'
    use "$OUT/figdata_A6$TAG.dta", clear
    gen byte main = (W == 24 & L == 24 & Hdays == 1000000)
    xsheet "A6_grid"
    foreach h in 5 8 {
        qui count if h == `h' & dd_or < .
        local nn = r(N)
        qui count if h == `h' & dd_or < 1
        local n1 = r(N)
        qui count if h == `h' & dd_p < 0.05
        local n2 = r(N)
        qui count if h == `h' & rd_p < 0.05
        local n3 = r(N)
        qui count if h == `h' & pre_or > 1 & pre_or < .
        local n4 = r(N)
        qui count if h == `h' & pre_p < 0.05
    }
    keep if dd_or < .
    bysort h (dd_or): gen int rank = _n
    local lo_dd 0.0025
    local hi_dd 1.6
    local yl_dd 0.01 "0·01" 0.1 "0·1" 0.25 "0·25" 0.5 "0·5" 1 "1"
    local nm_dd "Difference-in-discontinuities"
    local lo_rd 0.02
    local hi_rd 2.6
    local yl_rd 0.05 "0·05" 0.25 "0·25" 0.5 "0·5" 1 "1" 2 "2"
    local nm_rd "Post-implementation discontinuity"
    local lo_pre 0.45
    local hi_pre 13
    local yl_pre 0.5 "0·5" 1 "1" 2 "2" 5 "5" 10 "10"
    local nm_pre "Pre-implementation discontinuity"
    local k = 0
    local glist
    foreach e in dd rd pre {
        foreach h in 5 8 {
            local ++k
            local let : word `k' of A B C D E F
            local yt = cond(`h' == 5, "Odds ratio for in-hospital death", "")
            local xt = cond("`e'" == "pre", "Specification, ranked by the difference-in-discontinuities estimate", "")
            qui su rank if h == `h'
            local xmx = r(max) + 1
            twoway (scatter `e'_or rank if h == `h' & !(`e'_p < 0.05) & main == 0 & inrange(`e'_or, `lo_`e'', `hi_`e''), msymbol(O) mfcolor(white) mlcolor("$C_BEFORE") msize(small)) ///
                   (scatter `e'_or rank if h == `h' & `e'_p < 0.05 & main == 0 & inrange(`e'_or, `lo_`e'', `hi_`e''), msymbol(O) mcolor("$C_AFTER") msize(small)) ///
                   (scatter `e'_or rank if h == `h' & main == 1 & inrange(`e'_or, `lo_`e'', `hi_`e''), msymbol(O) mcolor("$C_MAIN") msize(medlarge)), ///
                   yline(1, lpattern(dash) lcolor(gs7) lwidth(thin)) yscale(log range(`lo_`e'' `hi_`e'')) ylabel(`yl_`e'', angle(0) labsize(small) grid glcolor(gs15) glwidth(vthin)) ///
                   xscale(range(0 `xmx')) xlabel(0(10)60, labsize(small)) xtitle("`xt'", size(small)) ytitle("`yt'", size(small)) legend(off) ///
                   title("{bf:`let'  `nm_`e'', window ±`h'}", position(11) size(small) color(black)) graphregion(color(white)) plotregion(lstyle(none)) name(a6_`k', replace) nodraw
            local glist `glist' a6_`k'
        }
    }
    twoway (scatteri 1 1, msymbol(O) mcolor("$C_AFTER") msize(0.6)) (scatteri 1 1, msymbol(O) mfcolor(white) mlcolor("$C_BEFORE") msize(0.6)) (scatteri 1 1, msymbol(O) mcolor("$C_MAIN") msize(1.0)), ///
           legend(order(1 "p<0·05" 2 "p≥0·05" 3 "Specification reported in the main text") rows(1) ring(0) position(0) region(lstyle(none) fcolor(white)) size(1.9)) ///
           xscale(off) yscale(off) xtitle("") ytitle("") graphregion(color(white)) plotregion(lstyle(none)) fysize(5) name(a6_leg, replace) nodraw
    graph combine `glist', cols(2) iscale(*0.85) imargin(small) graphregion(color(white)) name(a6_main, replace) nodraw
    graph combine a6_main a6_leg, cols(1) xsize(7.2) ysize(8.4) imargin(zero) graphregion(color(white)) name(A6, replace)
    gexp "FA6_grid"
}

log close
