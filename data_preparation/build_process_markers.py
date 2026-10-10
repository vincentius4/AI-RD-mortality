# -*- coding: utf-8 -*-
"""Build the process-of-care marker data files used by 11_process_markers.do (DeepCARS, JAMA manuscript).

Inputs (hospital records compiled for the companion study; not distributed):
  paper_data_v9.dta                         analysis data, one row per admission (adm_id, rv_dec, h_stay, treat, hB, hC, ...)
  new_analysis_data.parquet                 admission time per admission (admission_id_global, admit_dt)
  20260828_pooled/dcars_scores_add_predtrans.parquet   every AI score record (adm, rec_time, pred_trans; internal value = 100 * pred_trans)
  1008/vasopressor.parquet, 1008/oxygen.parquet, 1008/blood_culture.parquet    treatment records
  1010/vital.parquet                        vital-sign records
  Record tables carry admission_id_global, event_dt and time_resolution ("minute" when the clock time is recorded).

Outputs (one row per admission; indicators 0/1):
  process_t90_v1.dta   admissions within 8 points of the alert threshold (internal value 89.5): t0hr (reference time, hours after
                       admission) and the markers of Figure 3, eTable 12 panel B and eTable 13
  process_t<c>_v1.dta  the same markers at the placebo thresholds (displayed scores c = 80, 78, ..., 62; internal c - 0.5), with the
                       reference time applied to that threshold (eTable 13 panel B, eFigure 6)
  rules_t90_v1.dta     markers under the reference-time rules R0-R9 for patients below the threshold (eTable 14; admissions within 5 points)
  alt_assign_v1.dta    alternative assignment variables (eTable 10)

Definitions (eMethods). Reference time: above the threshold, the first record with an internal value at or above it within -1 to 24 h of
admission (the alert time); below it, the time of the 24-hour maximum. Marker names: measure_timing_definition, where
A<h>h = within h hours after the reference time, Bd<days> = calendar days after the admission day (day 0), any = at least one record,
new = at least one record and none in the 24 h before the reference time (calendar days: none on the preceding days from day -1).
Alert-timed markers use records with a clock time; calendar-day markers use all records.

Usage: python build_process_markers.py <root of the records> <folder with paper_data_v9.dta and new_analysis_data.parquet> <output folder>
"""
import os
import sys

import numpy as np
import pandas as pd

RAW, DATA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
CUT = 89.5
PLACEBO = [80, 78, 76, 74, 72, 70, 68, 66, 64, 62]
LO, HI = -26.0, 130.0          # hours after admission kept for the records
HL = [0.25, 0.5, 1, 2, 3, 6, 12, 24, 48, 72]

base = pd.read_stata(os.path.join(DATA, "paper_data_v9.dta"), columns=["adm_id", "rv_dec", "h_stay", "treat", "death", "hB", "hC"])
base = base[base.rv_dec.notna() & (base.h_stay > 24)].reset_index(drop=True)      # landmark cohort
base["hosp"] = np.where(base.hB == 1, "kangdong", np.where(base.hC == 1, "naeun", "sihwa"))
adm = pd.read_parquet(os.path.join(DATA, "new_analysis_data.parquet"), columns=["admission_id_global", "admit_dt"]).rename(columns={"admission_id_global": "adm_id"})
base = base.merge(adm, on="adm_id", how="left")
assert base.admit_dt.notna().all()
sc = pd.read_parquet(os.path.join(RAW, "20260828_pooled", "dcars_scores_add_predtrans.parquet"), columns=["adm", "rec_time", "pred_trans"]).rename(columns={"adm": "adm_id"})
sc = sc[sc.adm_id.isin(base.adm_id)].merge(base[["adm_id", "admit_dt"]], on="adm_id", how="inner")
sc["vv"] = 100 * sc.pred_trans.astype(float)
sc["hr"] = (sc.rec_time - sc.admit_dt).dt.total_seconds() / 3600
sc = sc[sc.hr >= -1].sort_values(["adm_id", "rec_time"]).reset_index(drop=True)
sc24 = sc[sc.hr <= 24]

RAW_EV = {"vaso_all": pd.read_parquet(os.path.join(RAW, "1008", "vasopressor.parquet")),
          "o2_escalation": pd.read_parquet(os.path.join(RAW, "1008", "oxygen.parquet")),
          "blood_culture": pd.read_parquet(os.path.join(RAW, "1008", "blood_culture.parquet")),
          "vital_any": pd.read_parquet(os.path.join(RAW, "1010", "vital.parquet"))}
RAW_EV["vaso_all"] = RAW_EV["vaso_all"][RAW_EV["vaso_all"].during_surgery == 0]          # outside surgery
RAW_EV["o2_escalation"] = RAW_EV["o2_escalation"][RAW_EV["o2_escalation"].o2_escalation == 1]
for k in RAW_EV:
    RAW_EV[k] = RAW_EV[k][["admission_id_global", "event_dt"] + (["time_resolution"] if "time_resolution" in RAW_EV[k].columns else [])].rename(columns={"admission_id_global": "adm_id"})


def events(name, P):
    """Records of the admissions in P: hours and calendar days after admission, clock-time flag."""
    df = RAW_EV[name]
    df = df[df.adm_id.isin(P.adm_id)].merge(P[["adm_id", "admit_dt"]], on="adm_id")
    df["h"] = (df.event_dt - df.admit_dt).dt.total_seconds() / 3600
    df["day"] = (df.event_dt.dt.normalize() - df.admit_dt.dt.normalize()).dt.days
    df["minute"] = (df.time_resolution == "minute") if "time_resolution" in df.columns else True
    df = df[(df.day >= -2) & (df.day <= 5) & (df.h >= LO) & (df.h <= HI)]
    return df[["adm_id", "h", "day", "minute"]]


def cohort(cutx, h):
    """Admissions within h points of the internal threshold cutx, with the reference time t0 (hours after admission)."""
    P = base.copy()
    P["sx"] = P.rv_dec - cutx
    P = P[np.abs(P.sx) < h].reset_index(drop=True)
    x = sc24[sc24.adm_id.isin(P.adm_id)]
    t_cross = x[x.vv >= cutx].groupby("adm_id").rec_time.first()
    t_max = x.loc[x.groupby("adm_id").vv.idxmax(), ["adm_id", "rec_time"]].set_index("adm_id").rec_time
    P["t0"] = pd.to_datetime(np.where(P.sx >= 0, P.adm_id.map(t_cross), P.adm_id.map(t_max)))
    P = P[P.t0.notna()].reset_index(drop=True)
    P["t0hr"] = (P.t0 - P.admit_dt).dt.total_seconds() / 3600
    return P


class Counter:
    """Counts of records in (lo, hi] hours relative to a per-admission reference time (sorted-key trick)."""

    def __init__(self, P, ev, t0hr, clock_only=True):
        self.n = len(P)
        self.S, self.OFF = 4e5, 2e5
        code = pd.Series(np.arange(len(P)), index=P.adm_id)
        e = ev[ev.minute] if clock_only else ev
        c = code.reindex(e.adm_id).values.astype(float)
        ok = ~np.isnan(c)
        dt = e.h.values[ok] - t0hr[c[ok].astype(int)]
        self.keys = np.sort(c[ok] * self.S + self.OFF + dt)
        self.base = np.arange(self.n) * self.S + self.OFF

    def __call__(self, lo, hi):
        return np.searchsorted(self.keys, self.base + hi, side="right") - np.searchsorted(self.keys, self.base + lo, side="right")


def day_counts(P, ev):
    code = pd.Series(np.arange(len(P)), index=P.adm_id)
    c = code.reindex(ev.adm_id).values.astype(float)
    ok = ~np.isnan(c)
    g = ev[ok].assign(c=c[ok].astype(int)).groupby(["c", "day"]).size().unstack(fill_value=0).reindex(index=np.arange(len(P)), columns=[-2, -1, 0, 1, 2, 3, 4, 5], fill_value=0)
    return g


def markers(P, which):
    """which: list of (measure, timing, definition). Returns a DataFrame of 0/1 indicators."""
    out = {}
    cache = {}
    for m, tim, dfn in which:
        if m not in cache:
            cache[m] = events(m, P)
        ev = cache[m]
        if tim.startswith("A"):
            H = float(tim[1:-1])
            key = ("A", m)
            if key not in cache:
                cache[key] = Counter(P, ev, P.t0hr.values)
            cnt = cache[key]
            after = cnt(0, H)
            col = (after > 0) if dfn == "any" else ((after > 0) & (cnt(-24, 0) == 0))
        else:
            days = [int(ch) for ch in tim[2:]]
            key = ("B", m)
            if key not in cache:
                cache[key] = day_counts(P, ev)
            g = cache[key]
            a = g[days].sum(axis=1).values
            prev = [d for d in range(-1, min(days))]
            pv = g[prev].sum(axis=1).values
            col = (a > 0) if dfn == "any" else ((a > 0) & (pv == 0))
        out[f"{m}_{tim.replace('.', '_')}_{dfn}"] = col.astype(np.int8)
    return pd.DataFrame(out)


MAIN = [("vital_any", "A1h", "any"), ("vital_any", "A2h", "any"), ("vital_any", "A3h", "any"),
        ("vaso_all", "Bd0", "any"), ("vaso_all", "Bd1", "any"), ("vaso_all", "Bd2", "any"), ("vaso_all", "Bd3", "any"), ("vaso_all", "Bd12", "any"),
        ("vaso_all", "A24h", "any"), ("vaso_all", "A48h", "any"), ("vaso_all", "A72h", "any"),
        ("o2_escalation", "A1h", "new"), ("o2_escalation", "A24h", "new"), ("o2_escalation", "A72h", "new"), ("blood_culture", "Bd01", "any")]
NINE = [("vital_any", "A1h", "any"), ("vital_any", "A2h", "any"), ("vital_any", "A3h", "any"), ("vaso_all", "Bd1", "any"), ("vaso_all", "Bd12", "any"),
        ("o2_escalation", "A1h", "new"), ("o2_escalation", "A24h", "new"), ("o2_escalation", "A72h", "new"), ("blood_culture", "Bd01", "any")]


def save(df, name):
    p = os.path.join(OUT, name)
    df.to_stata(p, write_index=False, version=117)
    print("saved", name, df.shape, flush=True)


# 1. alert threshold
P90 = cohort(CUT, 8)
F = markers(P90, MAIN)
save(pd.concat([P90[["adm_id", "t0hr"]], F], axis=1), "process_t90_v1.dta")

# 2. placebo thresholds (the nine markers of eTable 13)
for c in PLACEBO:
    Pc = cohort(c - 0.5, 8)
    Fc = markers(Pc, NINE)
    save(pd.concat([Pc[["adm_id", "t0hr"]], Fc], axis=1), f"process_t{c}_v1.dta")

# 3. reference-time rules R0-R9 (admissions within 5 points of the alert threshold)
P = base.copy()
P["sx"] = P.rv_dec - CUT
P = P[np.abs(P.sx) < 5].reset_index(drop=True)
x = sc24[sc24.adm_id.isin(P.adm_id)].copy()
vmax = x.groupby("adm_id").vv.max()
x["vmax"] = x.adm_id.map(vmax)
t_cross = x[x.vv >= CUT].groupby("adm_id").rec_time.first()
t_max = x.loc[x.groupby("adm_id").vv.idxmax(), ["adm_id", "rec_time"]].set_index("adm_id").rec_time
first_ge = lambda thr: x[x.vv >= thr].groupby("adm_id").rec_time.first()
first_near = lambda d: x[x.vv >= x.vmax - d].groupby("adm_id").rec_time.first()
below = P.sx < 0
above_t0 = P.adm_id.map(t_cross)
RULES = {0: (above_t0, P.adm_id.map(t_max)), 1: (above_t0, P.adm_id.map(first_ge(84.5))), 2: (above_t0, P.adm_id.map(first_ge(79.5))),
         3: (above_t0, P.adm_id.map(first_near(1.0))), 4: (above_t0, P.adm_id.map(first_near(2.5))), 5: (above_t0, P.adm_id.map(first_near(5.0))),
         6: (P.adm_id.map(t_max), P.adm_id.map(t_max)), 7: (P.adm_id.map(first_near(1.0)), P.adm_id.map(first_near(1.0)))}
# R8 matched above-threshold patient (same hospital and period, nearest admission time of day; ties broken at random, fixed seed); R9 hospital-period median
P["t0hr_above"] = (pd.to_datetime(above_t0) - P.admit_dt).dt.total_seconds() / 3600
P["adm_hod"] = P.admit_dt.dt.hour + P.admit_dt.dt.minute / 60
rng = np.random.default_rng(20261010)
m8 = pd.Series(np.nan, index=P.index)
m9 = pd.Series(np.nan, index=P.index)
for (hp, tr), g in P.groupby(["hosp", "treat"]):
    ab = g[(g.sx >= 0) & g.t0hr_above.notna()]
    bl = g[g.sx < 0]
    if len(ab) == 0:
        continue
    m9.loc[bl.index] = ab.t0hr_above.median()
    ah, at = ab.adm_hod.values, ab.t0hr_above.values
    for i, hod in zip(bl.index, bl.adm_hod.values):
        dd = np.abs(ah - hod)
        dd = np.minimum(dd, 24 - dd)
        cand = np.flatnonzero(dd == dd.min())
        m8.loc[i] = at[rng.choice(cand)]
RULES[8] = (above_t0, P.admit_dt + pd.to_timedelta(m8, unit="h"))
RULES[9] = (above_t0, P.admit_dt + pd.to_timedelta(m9, unit="h"))
RULE_MARK = [("vital_any", "A1h", "any"), ("vital_any", "A2h", "any"), ("o2_escalation", "A1h", "any"), ("o2_escalation", "A1h", "new")]
ev_rule = {m: events(m, P) for m in ("vital_any", "o2_escalation")}
F90 = pd.concat([P90[["adm_id"]], F[["vaso_all_Bd12_any", "vaso_all_Bd1_any"]]], axis=1)
frames = []
for r, (ta, tb) in RULES.items():
    Q = P.copy()
    Q["t0"] = pd.to_datetime(np.where(below, pd.to_datetime(tb), pd.to_datetime(ta)))
    Q = Q[Q.t0.notna()].reset_index(drop=True)
    Q["t0hr"] = (Q.t0 - Q.admit_dt).dt.total_seconds() / 3600
    cols = {}
    for m, tim, dfn in RULE_MARK:
        cnt = Counter(Q, ev_rule[m], Q.t0hr.values)
        H = float(tim[1:-1])
        after = cnt(0, H)
        cols[f"{m}_{tim}_{dfn}"] = ((after > 0) if dfn == "any" else ((after > 0) & (cnt(-24, 0) == 0))).astype(np.int8)
    Q = pd.concat([Q[["adm_id", "t0hr"]], pd.DataFrame(cols)], axis=1).merge(F90, on="adm_id", how="left")
    Q.insert(2, "rule", np.int8(r))
    frames.append(Q)
save(pd.concat(frames, ignore_index=True), "rules_t90_v1.dta")

# 4. alternative assignment variables: first value at or above the threshold within 24 h (otherwise the 24-h maximum); whole stay
fc = sc24[sc24.vv >= CUT].groupby("adm_id").vv.first()
vmax_all = sc.groupby("adm_id").vv.max()
first_val = sc[sc.vv >= CUT].groupby("adm_id").vv.first()
A = base[["adm_id", "rv_dec"]].copy()
A["rv_first24"] = A.rv_dec
A.loc[A.adm_id.isin(fc.index), "rv_first24"] = A.loc[A.adm_id.isin(fc.index), "adm_id"].map(fc)
A["rv_whole"] = A.adm_id.map(first_val)
A["rv_whole"] = A.rv_whole.where(A.rv_whole.notna(), A.adm_id.map(vmax_all))
save(A[["adm_id", "rv_first24", "rv_whole"]], "alt_assign_v1.dta")
print("done")
