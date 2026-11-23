#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M5: cell-type localisation of the 8 candidate genes.

Source: Human Protein Atlas per-gene JSON (consensus single-cell dataset,
        Siletti/Karlsson/... integrated; 80+ cell types).
Why not the bulk matrix: the 16 MB tsv.zip download is throttled on this
network; the per-gene JSON carries the same nCPM values for one gene and is
two orders of magnitude smaller.

Fields used
  RNA single cell type specific nCPM   -> per immune/stromal/epithelial cell type
  RNA blood cell specific nTPM         -> per circulating blood cell type
  RNA tissue specific nTPM             -> per tissue (cross-check of GTEx)

Output
  05_results/tables/celltype_ncpm.csv
  05_results/tables/celltype_summary.md
"""
import csv
import json
import os
import time
import urllib.request

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
RAW = os.path.join(ROOT, "03_data_raw", "hpa")
OUT = os.path.join(ROOT, "05_results", "tables")
os.makedirs(RAW, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

ENSEMBL = {
    "LRP1": "ENSG00000123384", "CTSD": "ENSG00000117984",
    "MTHFR": "ENSG00000177000", "AUTS2": "ENSG00000158321",
    "FURIN": "ENSG00000140564", "CST3": "ENSG00000101439",
    "FCER1G": "ENSG00000158869", "TBXAS1": "ENSG00000059377",
}

GUT_CELLS = ["Enterocytes", "Paneth cells", "Goblet cells",
             "Intestinal endocrine cells", "Tuft cells",
             "Enteric glial cells", "Smooth muscle cells"]
IMMUNE = ["Macrophages", "Monocytes", "Classical monocytes", "Kupffer cells",
          "Dendritic cells", "cDC", "pDC", "Neutrophils", "NK cells",
          "T-cells", "B-cells", "Plasma cells", "Mast cells"]
UA = {"User-Agent": "Mozilla/5.0"}


def fetch(gene):
    p = os.path.join(RAW, "%s.json" % gene)
    if os.path.exists(p) and os.path.getsize(p) > 2000:
        return json.load(open(p, encoding="utf-8"))
    url = "https://www.proteinatlas.org/%s.json" % ENSEMBL[gene]
    for i in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            data = urllib.request.urlopen(req, timeout=90).read()
            open(p, "wb").write(data)
            time.sleep(0.4)
            return json.loads(data)
        except Exception as e:
            print("  %s retry %d (%s)" % (gene, i, type(e).__name__))
            time.sleep(3)
    return None


def as_dict(v):
    if isinstance(v, dict):
        return v
    if isinstance(v, str):
        try:
            return json.loads(v.replace("'", '"'))
        except Exception:
            return {}
    return {}


def main():
    recs, summary = [], []
    for g in ENSEMBL:
        d = fetch(g)
        if not d:
            print("%-7s FAILED" % g)
            continue
        sc = as_dict(d.get("RNA single cell type specific nCPM"))
        blood = as_dict(d.get("RNA blood cell specific nTPM"))
        tis = as_dict(d.get("RNA tissue specific nTPM"))
        print("%-7s cell types=%3d  blood types=%2d  tissues=%2d  spec=%s" % (
            g, len(sc), len(blood), len(tis),
            d.get("RNA single cell type specificity")))
        for ct, val in sc.items():
            try:
                recs.append(dict(gene=g, compartment="single_cell_type",
                                 cell_type=ct, value=float(val), unit="nCPM"))
            except ValueError:
                pass
        for ct, val in blood.items():
            try:
                recs.append(dict(gene=g, compartment="blood_cell",
                                 cell_type=ct, value=float(val), unit="nTPM"))
            except ValueError:
                pass

        # gut epithelial + immune ranking
        gutvals = {k: float(v) for k, v in sc.items()
                   if k in GUT_CELLS or "nterocyte" in k or "oblet" in k}
        immvals = {k: float(v) for k, v in sc.items()
                   if k in IMMUNE or "cell" in k and ("phage" in k or "cyte" in k)}
        top = sorted(((k, float(v)) for k, v in sc.items()), key=lambda x: -x[1])[:3]
        summary.append(dict(
            gene=g,
            specificity=d.get("RNA single cell type specificity", ""),
            distribution=d.get("RNA single cell type distribution", ""),
            top1="%s (%.0f)" % top[0] if top else "",
            top2="%s (%.0f)" % top[1] if len(top) > 1 else "",
            top3="%s (%.0f)" % top[2] if len(top) > 2 else "",
            max_gut=round(max(gutvals.values()), 1) if gutvals else float("nan"),
            max_gut_cell=max(gutvals, key=gutvals.get) if gutvals else "",
            max_immune=round(max(immvals.values()), 1) if immvals else float("nan"),
            max_immune_cell=max(immvals, key=immvals.get) if immvals else ""))

    p1 = os.path.join(OUT, "celltype_ncpm.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "compartment", "cell_type",
                                          "value", "unit"])
        w.writeheader()
        w.writerows(recs)
    print("written", p1, len(recs), "rows")

    p2 = os.path.join(OUT, "celltype_summary.csv")
    keys = ["gene", "specificity", "distribution", "top1", "top2", "top3",
            "max_gut", "max_gut_cell", "max_immune", "max_immune_cell"]
    with open(p2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(summary)
    print("written", p2)

    lines = ["| gene | HPA specificity | top cell type (nCPM) | max gut epithelial (cell) | max immune (cell) |",
             "|---|---|---|---|---|"]
    for r in summary:
        lines.append("| %s | %s | %s | %.0f (%s) | %.0f (%s) |" % (
            r["gene"], r["specificity"], r["top1"], r["max_gut"],
            r["max_gut_cell"], r["max_immune"], r["max_immune_cell"]))
    p3 = os.path.join(OUT, "celltype_summary.md")
    open(p3, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("written", p3)
    print()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
