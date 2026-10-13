# DeepCARS: in-hospital mortality at the AI alert threshold

Stata 15 code that reproduces every table and figure of the paper. Data are not included; see the data sharing statement in the paper.

Requires `rdrobust` (eTable 2) and `rddensity` version 3.0 (eTable 4; install from github.com/rdpackages/rddensity, which includes `rddensity_fun.do`). Put the data files in this folder and run `do deepcars_tables.do` and then, in the same session, `do deepcars_figures.do` (to rerun one table, run part 1 of `deepcars_tables.do` and then that part).

- `paper_data_v9.dta`: analysis data, one row per admission
- `process_flags_v1.dta`: treatments within 72 h of admission and in the 24 h before admission (eTables 9 and 12, Table 1)
- `process_t90_v2.dta`: process-of-care markers at the alert threshold, including the treatments on days 1-2 of Table 1 (Table 1, Figure 3, eTables 12 and 13)
- `process_t80_v1.dta` ... `process_t62_v1.dta`: the same markers at the placebo thresholds (eTable 13, eFigure 6)
- `rules_t90_v1.dta`: markers under the reference-time rules R0-R9 (eTable 14)
- `alt_assign_v1.dta`: alternative assignment variables (eTable 10)

| File | Output |
|---|---|
| `deepcars_tables.do` | part 1: landmark cohort, estimation and formatting programs; parts 2-3: Table 1, Table 2; parts 4-17: eTable 1 ... eTable 14 |
| `deepcars_figures.do` | part 1: Figure 3, eFigure 6; part 2: Figure 1, Figure 2, eFigures 1-5 and 7 (starts with `clear all`; self-contained) |

Each table part writes one sheet of `results/JAMA_tables.xlsx` with the cells as printed in the paper (table notes are not included); figures are written to `results/figures/`. Intermediate data sets are kept in `results/work/`.

`data_preparation/` holds the Python scripts that build the marker and treatment files from the hospital records (treatment and vital-sign records with time stamps, AI score records, admission times); definitions follow the eMethods.

Estimates in the paper were computed in Python 3.13; this Stata code reproduces them. Values that depend on the rdrobust version (robust standard errors of eTable 2) may differ in the last digit.
