#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M5 (complete): cell-type localisation from the FULL HPA single-cell matrix.

The 16.3 MB `rna_single_cell_type.tsv` was downloaded and filtered to the 8
Ensembl ids, giving a complete 8 x 154 matrix (no reporting threshold), which
removes the caveat of the earlier per-gene JSON approach.

Outputs
  05_results/tables/celltype_gut_focus.csv     gene x focused cell types
  05_results/tables/celltype_full_summary.csv  per-gene maxima by compartment
"""
import csv
import os
from collections import defaultdict

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T = os.path.join(ROOT, "05_results", "tables")

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]

GUT = ["colonocytes", "enterocytes", "goblet cells", "paneth cells",
       "tuft cells", "enteric stem cells", "enteric transient amplifying cells"]
CARDIO = ["cardiomyocytes", "vascular smooth muscle cells",
          "vascular endothelial cells", "pericytes", "fibroblasts",
          "fibro-adipogenic progenitors", "epicardial cells"]
KIDNEY = ["proximal tubule cells", "podocytes",
          "loop of henle epithelial cells",
          "renal collecting duct principal cells",
          "renal collecting duct intercalated cells",
          "distal convoluted tubule cells"]
IMMUNE = ["macrophages", "monocytes", "neutrophils", "t-cells", "b-cells",
          "plasma cells", "nk-cells", "mast cells", "cdc", "pdcs",
          "innate lymphoid cells", "kupffer cells"]

rows = list(csv.DictReader(open(os.path.join(T, "celltype_full_matrix.csv"),
                                encoding="utf-8")))
val = defaultdict(dict)
for r in rows:
    val[r["gene"]][r["cell_type"]] = float(r["ncpm"])

# ---- focused table ----
p1 = os.path.join(T, "celltype_gut_focus.csv")
with open(p1, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene", "compartment", "cell_type", "ncpm"])
    for g in GENES:
        for comp, lst in (("gut_epithelial", GUT), ("cardiovascular", CARDIO),
                          ("kidney", KIDNEY), ("immune", IMMUNE)):
            for ct in lst:
                w.writerow([g, comp, ct, round(val[g].get(ct, 0.0), 1)])
print("written", p1)

# ---- per-gene summary ----
p2 = os.path.join(T, "celltype_full_summary.csv")
summary = []
with open(p2, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene", "max_cell", "max_ncpm", "gut_top_cell", "gut_top_ncpm",
                "cardio_top_cell", "cardio_top_ncpm", "immune_top_cell",
                "immune_top_ncpm", "kidney_top_cell", "kidney_top_ncpm",
                "gut_share_of_max"])
    for g in GENES:
        d = val[g]
        mx_ct = max(d, key=d.get)
        mx = d[mx_ct]
        def top(lst):
            sub = {k: d.get(k, 0.0) for k in lst}
            k = max(sub, key=sub.get)
            return k, sub[k]
        gc, gv = top(GUT)
        cc, cv = top(CARDIO)
        ic, iv = top(IMMUNE)
        kc, kv = top(KIDNEY)
        row = [g, mx_ct, round(mx, 1), gc, round(gv, 1), cc, round(cv, 1),
               ic, round(iv, 1), kc, round(kv, 1),
               "%.3f" % (gv / mx if mx else 0)]
        w.writerow(row)
        summary.append(row)
print("written", p2)
print()
print("%-7s %-32s %8s %-30s %8s %8s" % ("gene", "top cell type", "nCPM",
                                        "top gut cell", "nCPM", "gut/max"))
for r in summary:
    print("%-7s %-32s %8.0f %-30s %8.0f %8s" % (r[0], r[1][:32], r[2],
                                                r[3][:30], r[4], r[11]))
