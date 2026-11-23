#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MR-coloc discordance dissection.

For MTHFR (the only gene with a strict cis instrument AND a coloc result),
plot the variant-level Wald ratio against position.  If one variant drives both
traits, the ratio is flat; if two distinct variants sit in partial LD, the ratio
tracks the distance to the blood-pressure lead.

Output: 05_results/tables/mthfr_wald_by_position.csv
"""
import csv
import importlib.util
import os

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
spec = importlib.util.spec_from_file_location(
    "mr", os.path.join(ROOT, "04_code", "13_cis_mr.py"))
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)

HTN_LEAD = 11823674   # rs12567136, nearest gene CLCN6
TISSUES = ["Colon_Sigmoid", "Colon_Transverse",
           "Small_Intestine_Terminal_Ileum", "Whole_Blood"]

recs = []
for t in TISSUES:
    rows = mr.load_eqtl("MTHFR", t)
    out = mr.load_out("MTHFR", "I9_HYPTENS")
    for r in rows:
        if r["p"] >= 5e-8:
            continue
        o = out.get(r["pos"])
        if not o:
            continue
        h = mr.harmonise(r, o)
        if h is None:
            continue
        f = (r["beta"] / r["se"]) ** 2
        if f < 10:
            continue
        b = h[0] / r["beta"]
        se = h[1] / abs(r["beta"])
        recs.append(dict(tissue=t, rsid=r["rsid"], pos=r["pos"],
                         p_eqtl="%.2e" % r["p"], p_out="%.2e" % o["p"],
                         wald=round(b, 4), wald_se=round(se, 4),
                         dist_htn_kb=round(abs(r["pos"] - HTN_LEAD) / 1000.0, 2)))

p = os.path.join(ROOT, "05_results", "tables", "mthfr_wald_by_position.csv")
with open(p, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(recs[0].keys()))
    w.writeheader()
    w.writerows(recs)
print("written", p, len(recs), "variant-tissue rows")

# distance-stratified summary
import statistics
def band(k):
    if k < 20: return "<20 kb"
    if k < 50: return "20-50 kb"
    return ">=50 kb"
agg = {}
for r in recs:
    agg.setdefault(band(r["dist_htn_kb"]), []).append(r["wald"])
print()
print("%-10s %5s %10s %10s %10s" % ("band", "n", "median", "min", "max"))
for k in ["<20 kb", "20-50 kb", ">=50 kb"]:
    v = agg.get(k, [])
    if v:
        print("%-10s %5d %10.3f %10.3f %10.3f" % (
            k, len(v), statistics.median(v), min(v), max(v)))
