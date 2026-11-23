# -*- coding: utf-8 -*-
"""Step 2 - extract protein-level (SomaScan v4) evidence for the eight candidates
plus the neighbouring gene FES
from the deCODE / Ferkingstad 2021 supplementary workbook.

ST01 lists all 5,284 aptamers and whether each entered the analysis, which settles
the coverage question outright.  ST02 holds the 28,191 sentinel and secondary pQTL
associations with effect size, P value, cis/trans class and LD-class membership.

Outputs
  pqtl_panel_coverage.csv   one row per candidate: on array? analysed? SeqId, UniProt
  pqtl_associations.csv     every pQTL row for the four proteins
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xlsx_stream import Workbook

WB = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/_decode_source_data/41588_2021_978_MOESM4_ESM.xlsx"
GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3",
         "FCER1G", "TBXAS1", "FES"]
# ST02 columns we keep
KEEP = {2: "gene", 9: "seqid", 10: "variant", 11: "ld_class_size", 13: "chr",
        14: "pos", 17: "amin", 18: "amaj", 20: "maf", 21: "cis_trans",
        22: "dist_tss", 23: "closest_genes", 25: "beta_adj",
        27: "minuslog10p_adj", 38: "rank_pqtl_per_protein",
        42: "reported", 43: "reported_rsid"}


def main():
    wb = Workbook(WB)

    # ---------------- ST01: panel coverage ----------------
    cov = []
    it = wb.rows("ST01")
    hdr = None
    for r in it:
        if hdr is None:
            hdr = r
            continue
        if len(r) < 5:
            continue
        gene = r[3].upper()
        if gene in GENES:
            cov.append({"gene": r[3], "seqid": r[0], "protein": r[1],
                        "uniprot": r[6] if len(r) > 6 else "",
                        "included": r[5] if len(r) > 5 else "",
                        "ncand": "yes"})
    seen = {c["gene"].upper() for c in cov}
    for g in GENES:
        if g not in seen:
            cov.append({"gene": g, "seqid": "", "protein": "", "uniprot": "",
                        "included": "NOT ON ARRAY", "ncand": "no"})
    cov.sort(key=lambda c: GENES.index(c["gene"].upper()))
    with open(os.path.join(HERE, "pqtl_panel_coverage.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "seqid", "protein", "uniprot",
                                          "included", "ncand"])
        w.writeheader()
        w.writerows(cov)
    print("=== SomaScan v4 面板覆盖 ===")
    for c in cov:
        print("  %-6s %-12s %-10s %-12s %s"
              % (c["gene"], c["seqid"], c["uniprot"], c["included"], c["protein"]))

    # ---------------- ST02: pQTL associations ----------------
    out = []
    it = wb.rows("ST02")
    next(it)                      # title line
    next(it)                      # header line
    n = 0
    for r in it:
        n += 1
        if len(r) < 11:
            continue
        if r[2].upper() not in GENES:
            continue
        rec = {v: (r[k] if k < len(r) else "") for k, v in KEEP.items()}
        out.append(rec)
    print("\n扫描 ST02 行数: %d，命中 %d 条" % (n, len(out)))

    with open(os.path.join(HERE, "pqtl_associations.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(KEEP.values()))
        w.writeheader()
        w.writerows(out)

    print("\n=== 各基因 cis 方向 pQTL 数 / 全部 ===")
    for g in GENES:
        sub = [o for o in out if o["gene"].upper() == g]
        cis = [o for o in sub if o["cis_trans"].lower() == "cis"]
        print("  %-6s 全部 %3d，cis %3d" % (g, len(sub), len(cis)))
        for o in sorted(cis, key=lambda x: -float(x["minuslog10p_adj"] or 0))[:4]:
            print("      %-12s pos=%-11s MAF=%-6s beta=%+8s  -log10P=%s"
                  % (o["variant"], o["pos"], o["maf"][:5],
                     (o["beta_adj"] or "")[:8], o["minuslog10p_adj"][:7]))


if __name__ == "__main__":
    main()
