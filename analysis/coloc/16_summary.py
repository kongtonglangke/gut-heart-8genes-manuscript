#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Master summary: merge the pieces into one gene-level evidence table.

Inputs
  05_results/tables/coloc_<pheno>.csv         (gene x tissue x pheno)
  05_results/tables/matrix_raw.csv            (detectability)
  05_results/tables/mr_results.csv
  05_results/tables/expression_tau.csv

Outputs
  05_results/tables/coloc_master.csv
  05_results/tables/evidence_master.md
"""
import csv
import os
from collections import defaultdict

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T = os.path.join(ROOT, "05_results", "tables")

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]
PHENOS = [("I9_HYPTENS", "Hypertension"),
          ("I9_HEARTFAIL", "Heart failure"),
          ("I9_AF", "Atrial fibrillation"),
          ("N14_CHRONKIDNEYDIS", "Chronic kidney disease"),
          ("T2D", "Type 2 diabetes")]
GUT = ["Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum"]


def read(p):
    return list(csv.DictReader(open(p, encoding="utf-8"))) if os.path.exists(p) else []


# ---- coloc master ----
master = []
for ph, label in PHENOS:
    rows = read(os.path.join(T, "coloc_%s.csv" % ph))
    best = {}
    for r in rows:
        g = r["gene"]
        v = float(r["PPH4"])
        if g not in best or v > float(best[g]["PPH4"]):
            best[g] = r
    for g in GENES:
        r = best.get(g)
        if not r:
            continue
        master.append(dict(
            pheno=ph, pheno_label=label, gene=g, tissue=r["tissue"],
            n_snp=r["n_snp"], minP_eqtl=r["minP_eqtl"], minP_out=r["minP_out"],
            PPH3=r["PPH3"], PPH4=r["PPH4"],
            gut_best="yes" if r["tissue"] in GUT else "no"))

with open(os.path.join(T, "coloc_master.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(master[0].keys()))
    w.writeheader()
    w.writerows(master)

# ---- evidence table ----
tau = {r["gene"]: r for r in read(os.path.join(T, "expression_tau.csv"))}
mr = defaultdict(dict)
for r in read(os.path.join(T, "mr_results.csv")):
    if r["iv_set"] != "strict_5e-8" or r.get("status") != "ok":
        continue
    mr[r["gene"]][r["pheno"]] = r

mat = defaultdict(dict)
for r in read(os.path.join(T, "matrix_raw.csv")):
    mat[r["gene"]][r["tissue"]] = r

lines = []
lines.append("| gene | Tau | gut eGene (of 3) | best gut eQTL minP | "
             "HTN coloc PPH4 | AF coloc PPH4 | HF coloc PPH4 | cis-MR HTN OR (95% CI) |")
lines.append("|---|---|---|---|---|---|---|---|")

for g in GENES:
    t = tau.get(g, {})
    n_eg = sum(1 for x in GUT if g in mat and x in mat[g]
               and float(mat[g][x].get("minP", 1)) < 5e-8)
    gp = min((float(mat[g][x]["minP"]) for x in GUT if g in mat and x in mat[g]),
             default=float("nan"))
    def ph4(ph):
        v = [m for m in master if m["gene"] == g and m["pheno"] == ph]
        return "%.3f" % float(v[0]["PPH4"]) if v else "-"
    m = mr.get(g, {}).get("I9_HYPTENS")
    mrs = ("%.3f (%.3f-%.3f)" % (float(m["OR"]), float(m["ci_lo"]), float(m["ci_hi"]))
           if m else "no IV")
    lines.append("| %s | %s | %d/3 | %.1e | %s | %s | %s | %s |" % (
        g, t.get("tau", "-"), n_eg, gp, ph4("I9_HYPTENS"), ph4("I9_AF"),
        ph4("I9_HEARTFAIL"), mrs))

md = "\n".join(lines) + "\n"
open(os.path.join(T, "evidence_master.md"), "w", encoding="utf-8").write(md)
print(md)

# ---- how many positives ----
tot = len(master)
pos = [m for m in master if float(m["PPH4"]) >= 0.8]
print("coloc cells (gene x pheno, best tissue) = %d ; PPH4>=0.8 = %d" % (tot, len(pos)))
