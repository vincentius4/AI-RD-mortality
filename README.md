# DeepCARS: in-hospital mortality at the AI alert threshold

Stata 15 code that reproduces the tables, figures and text estimates of the paper. Data are not included; see the data sharing statement in the paper.

Requires `rdrobust` and `rddensity`. Put the data files in this folder and run `do 00_master.do`:

- `paper_data_v9.dta`: analysis data, one row per admission (all files)
- `process_flags_v1.dta`: treatment indicators within 72 h of admission and in the 24 h before admission (`10_additions.do`)
- `process_t90_v1.dta`, `process_t80_v1.dta` ... `process_t62_v1.dta`, `rules_t90_v1.dta`, `alt_assign_v1.dta`: process-of-care markers at the alert threshold and at the placebo thresholds, the markers under the alternative reference-time rules, and the alternative assignment variables (`11_process_markers.do`)

`data_preparation/` holds the Python scripts that build `process_flags_v1.dta` and the marker files from the hospital records (treatment and vital-sign records with time stamps, AI score records, admission times); the definitions follow the eMethods. `paper_data_v9.dta` was assembled from the same records in Python (variables defined in the eMethods).

`10_jama_additions.do` adds the model-based risk differences of Table 2 (with and without covariates), the process-of-care estimates within 72 h of admission, and in-hospital mortality adjusted for treatments recorded before admission. `11_process_markers.do` reproduces Figure 3, the alternative assignment variables (eTable 10), the vasopressor windows at hospitals B and C (eTable 12), the marker specifications, linear probability models, hospital-level estimates and placebo thresholds at all 3 hospitals and at hospitals B and C (eTable 13, eFigure 6), and the reference-time rules (eTable 14). Results are written to `replicate_v12/`.

Estimates in the paper were computed in Python 3.13; this Stata code reproduces them. Values that depend on random numbers or on the rdrobust version may differ slightly. The matched reference times of rule R8 use a fixed random seed.
