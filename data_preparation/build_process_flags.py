# -*- coding: utf-8 -*-
"""Build process_flags_v1.dta used by 10_jama_additions.do (DeepCARS, JAMA manuscript): treatment indicators per admission.

Inputs (not distributed): paper_data_v9.dta, new_analysis_data.parquet (admit_dt), 1008/antibiotic.parquet, 1008/blood_culture.parquet,
1008/lactate.parquet, 1008/oxygen.parquet, 1008/vasopressor.parquet (admission_id_global, event_dt, time_resolution, and the fields used below).
Only records with a clock time (time_resolution == "minute") are used.
  a_<k>    at least one record in (0, 72] h after admission        b_<k>    at least one record in [-24, 0] h before admission
  k: abx (intravenous antibiotic), bc (blood culture), lac (lactate), o2 (oxygen escalation), vaso (vasopressor or inotrope outside surgery),
     press (vasopressor outside surgery: norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine)
  a_any, b_any: any of abx, bc, lac, o2, vaso        newpress: a_press == 1 and b_press == 0
Usage: python build_process_flags.py <root of the records> <folder with paper_data_v9.dta and new_analysis_data.parquet> <output folder>
"""
import os
import sys

import numpy as np
import pandas as pd

RAW, DATA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
d = pd.read_stata(os.path.join(DATA, "paper_data_v9.dta"), columns=["adm_id"])
ad = pd.read_parquet(os.path.join(DATA, "new_analysis_data.parquet"), columns=["admission_id_global", "admit_dt"]).rename(columns={"admission_id_global": "adm_id"})
d = d.merge(ad, on="adm_id", how="left")
assert d.admit_dt.notna().all()
B = os.path.join(RAW, "1008")


def load(f, filt=None):
    e = pd.read_parquet(os.path.join(B, f + ".parquet"))
    if filt is not None:
        e = filt(e)
    e = e[e.time_resolution == "minute"].rename(columns={"admission_id_global": "adm_id"})
    e = e[e.adm_id.isin(d.adm_id)].merge(d, on="adm_id")
    e["hr"] = (e.event_dt - e.admit_dt).dt.total_seconds() / 3600
    return e[["adm_id", "hr"]]


PRESS = ["norepinephrine", "epinephrine", "phenylephrine", "vasopressin", "dopamine"]
EV = {"abx": load("antibiotic", lambda e: e[e.route == "iv"]), "bc": load("blood_culture"), "lac": load("lactate"),
      "o2": load("oxygen", lambda e: e[e.o2_escalation == 1]), "vaso": load("vasopressor", lambda e: e[e.during_surgery == 0]),
      "press": load("vasopressor", lambda e: e[(e.during_surgery == 0) & e.agent_class.isin(PRESS)])}
out = d[["adm_id"]].copy()
for k, e in EV.items():
    out[f"a_{k}"] = out.adm_id.isin(set(e.adm_id[(e.hr > 0) & (e.hr <= 72)])).astype(np.int8)
    out[f"b_{k}"] = out.adm_id.isin(set(e.adm_id[(e.hr >= -24) & (e.hr <= 0)])).astype(np.int8)
out["a_any"] = out[[f"a_{k}" for k in ("abx", "bc", "lac", "o2", "vaso")]].max(axis=1).astype(np.int8)
out["b_any"] = out[[f"b_{k}" for k in ("abx", "bc", "lac", "o2", "vaso")]].max(axis=1).astype(np.int8)
out["newpress"] = ((out.a_press == 1) & (out.b_press == 0)).astype(np.int8)
out.to_stata(os.path.join(OUT, "process_flags_v1.dta"), write_index=False, version=117)
print("saved process_flags_v1.dta", out.shape)
