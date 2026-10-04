* Table 1
di as txt _n "===== 1. Table 1 ====="
program define npc, rclass
    args v cond
    qui count if `cond' & `v' < .
    local n = r(N)
    qui count if `cond' & `v' == 1
    local k = r(N)
    return local s = subinstr(string(`k') + " (" + string(100*`k'/`n', "%4.1f") + "%)", ".", "·", .)
end
program define miqr, rclass
    args v cond dp
    preserve
    qui keep if `cond' & `v' < .
    sort `v'
    local n = _N
    foreach q in 25 50 75 {
        local pos = (`n' - 1)*`q'/100 + 1
        local lo = floor(`pos')
        local hi = ceil(`pos')
        local q`q' = `v'[`lo'] + (`pos' - `lo')*(`v'[`hi'] - `v'[`lo'])
    }
    restore
    local f = cond(`dp' == 0, "%3.0f", "%4.1f")
    return local s = subinstr(string(`q50', "`f'") + " (" + string(`q25', "`f'") + "–" + string(`q75', "`f'") + ")", ".", "·", .)
end
program define msd, rclass
    args v cond
    qui su `v' if `cond'
    return local s = subinstr(string(r(mean), "%9.2f") + " (" + string(r(sd), "%9.2f") + ")", ".", "·", .)
end
tempname T1
postfile `T1' str60 row str24 (c1 c2 c3 c4 c5 c6 c7 c8) using "replicate_v12\res_T1.dta", replace
local G1 "post==0 & abs(s)<5 & D==0"
local G2 "post==0 & abs(s)<5 & D==1"
local G3 "post==0 & D==0"
local G4 "post==0 & D==1"
local G5 "post==1 & abs(s)<5 & D==0"
local G6 "post==1 & abs(s)<5 & D==1"
local G7 "post==1 & D==0"
local G8 "post==1 & D==1"
local cells
forvalues g = 1/8 {
    qui count if `G`g''
    local cells `cells' ("`=string(r(N))'")
}
post `T1' ("Admissions") `cells'
local R1 "DeepCARS value, median (IQR)|m1|rv_dec"
local R2 "Hospital A|b|hA"
local R3 "Hospital B|b|hB"
local R4 "Hospital C|b|hC"
local R5 "Age, years, median (IQR)|m0|age"
local R6 "Male|b|male"
local R7 "Charlson index, mean (SD)|sd|cci"
local R8 "Hypertension|b|htn"
local R9 "SOFA score, mean (SD)|sd|sofa"
local R10 "Surgical service|b|surgical_dept"
local R11 "Sepsis|b|sepsis"
local R12 "DNR order at admission|b|dnr_at_admit"
local R13 "Length of stay, days, median (IQR)|m1|total_los"
local R14 "DNR order during stay|b|dnr_case"
local R15 "In-hospital death|b|death"
local R16 "Ward cardiac arrest|b|arr_gw"
local R17 "Unplanned transfer to intensive care after 24 h|b|icu_uit_24h"
cap drop hA
gen byte hA = (hosp_num == 1)
forvalues r = 1/17 {
    tokenize "`R`r''", parse("|")
    local lab "`1'"
    local kind "`3'"
    local var "`5'"
    local cells
    forvalues g = 1/8 {
        if "`kind'" == "b" {
            npc `var' "`G`g''"
            local cells `cells' ("`r(s)'")
        }
        else if "`kind'" == "sd" {
            msd `var' "`G`g''"
            local cells `cells' ("`r(s)'")
        }
        else {
            local dp = cond("`kind'" == "m0", 0, 1)
            miqr `var' "`G`g''" `dp'
            local cells `cells' ("`r(s)'")
        }
    }
    post `T1' ("`lab'") `cells'
}
postclose `T1'
xout "replicate_v12\res_T1.dta" "T1"

