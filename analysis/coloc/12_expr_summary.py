#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M1 summary: tissue specificity (Tau), gut vs blood vs heart, and top tissues.

Tau = sum(1 - x_i/max(x)) / (n - 1),  x_i = median TPM in tissue i
      (Yanai et al. 2005). 0 = ubiquitous, 1 = single-tissue specific.

Output:
  05_results/tables/expression_tau.csv
  05_results/tables/expression_key_tissues.csv
  05_results/tables/expression_top_tissues.md
"""
import csv
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "03_data_raw", "expression_median_tpm.csv")
OUT = os.path.join(ROOT, "05_results", "tables")
os.makedirs(OUT, exist_ok=True)

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]

GUT = ["Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum"]
BLOOD = ["Whole_Blood"]
HEART = ["Heart_Left_Ventricle", "Heart_Atrial_Appendage"]
KIDNEY = ["Kidney_Cortex", "Kidney_Medulla"]

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
by_gene = defaultdict(dict)
for r in rows:
    by_gene[r["gene"]][r["tissue"]] = float(r["median_tpm"])


def tau(vals):
    m = max(vals)
    if m <= 0:
        return float("nan")
    n = len(vals)
    return sum(1.0 - v / m for v in vals) / (n - 1)


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else float("nan")


# ---- Tau ----
p1 = os.path.join(OUT, "expression_tau.csv")
with open(p1, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene", "n_tissues", "max_tpm", "max_tissue", "median_tpm",
                "tau", "tau_class", "gut_mean", "blood", "heart_mean",
                "kidney_mean", "gut_vs_blood", "gut_vs_heart"])
    for g in GENES:
        d = by_gene[g]
        vals = list(d.values())
        mx = max(vals)
        mx_t = max(d, key=d.get)
        med = sorted(vals)[len(vals) // 2]
        t = tau(vals)
        cls = "tissue-specific" if t >= 0.85 else ("intermediate" if t >= 0.6 else "broadly expressed")
        gm = mean([d.get(x) for x in GUT])
        bm = d.get("Whole_Blood", float("nan"))
        hm = mean([d.get(x) for x in HEART])
        km = mean([d.get(x) for x in KIDNEY])
        w.writerow([g, len(vals), round(mx, 2), mx_t, round(med, 2),
                    round(t, 4), cls, round(gm, 2), round(bm, 2),
                    round(hm, 2), round(km, 2),
                    round(gm / bm, 3) if bm and bm > 0 else "",
                    round(gm / hm, 3) if hm and hm > 0 else ""])
print("written", p1)

# ---- key tissues ----
KEY = GUT + BLOOD + HEART + KIDNEY + ["Liver", "Adrenal_Gland", "Lung", "Muscle_Skeletal"]
p2 = os.path.join(OUT, "expression_key_tissues.csv")
with open(p2, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene"] + KEY)
    for g in GENES:
        w.writerow([g] + [round(by_gene[g].get(t, float("nan")), 2) for t in KEY])
print("written", p2)

# ---- top tissues markdown ----
lines = ["| gene | max TPM | tissue | rank-2 tissue (TPM) | rank-3 tissue (TPM) | Tau | class |",
         "|---|---|---|---|---|---|---|"]
for g in GENES:
    d = by_gene[g]
    top = sorted(d.items(), key=lambda kv: -kv[1])[:3]
    t = tau(list(d.values()))
    cls = "tissue-specific" if t >= 0.85 else ("intermediate" if t >= 0.6 else "broadly expressed")
    lines.append("| %s | %.1f | %s | %s (%.1f) | %s (%.1f) | %.3f | %s |" % (
        g, top[0][1], top[0][0], top[1][0], top[1][1], top[2][0], top[2][1], t, cls))
p3 = os.path.join(OUT, "expression_top_tissues.md")
open(p3, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("written", p3)
print()
print("\n".join(lines))
