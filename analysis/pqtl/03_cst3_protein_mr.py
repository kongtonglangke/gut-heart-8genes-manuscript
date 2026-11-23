# -*- coding: utf-8 -*-
"""Step 3 - protein-level test for CST3, and the LD classes needed to judge overlap
with the hypertension signal.

For each CST3 cis-pQTL we keep the full LD-class variant list reported by deCODE,
then pull the same window out of FinnGen R12 hypertension and compute a Wald ratio
(circulating protein -> hypertension) plus an LD-overlap check.
"""
import csv
import os
import sys
from math import erf, sqrt, exp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, r"F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr")

from xlsx_stream import Workbook
from _tabix import Tabix

WB = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/_decode_source_data/41588_2021_978_MOESM4_ESM.xlsx"
HTN_URL = ("https://storage.googleapis.com/finngen-public-data-r12/"
           "summary_stats/release/finngen_R12_I9_HYPTENS.gz")
HTN_TBI = r"F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr/data/finngen_htn.tbi"
# FinnGen R12 column order (same as 30_region_scan.py)
F = dict(pos=1, ref=2, alt=3, rsid=4, gene=5, p=6, mlogp=7, beta=8, se=9, af=10)

CST3_WIN = ("20", 23526706, 23738556)      # gene body +/- 100 kb, GRCh38
FURIN_WIN = ("15", 90768588, 90983564)     # gene body +/- 100 kb, GRCh38


def decode_cst3():
    wb = Workbook(WB)
    it = wb.rows("ST02")
    next(it)
    next(it)
    out = []
    for r in it:
        if len(r) < 13 or r[2].upper() != "CST3":
            continue
        out.append({
            "variant": r[10], "chr": r[13], "pos": r[14],
            "amin": r[17], "amaj": r[18], "maf": r[20],
            "cis_trans": r[21], "beta": r[25], "mlogp": r[27],
            "ld_size": r[11], "ld_class": r[12],
        })
    return out


def finngen(win, tag):
    tb = Tabix(HTN_URL, tbi_path=HTN_TBI)
    chrom, beg, end = win
    rows = tb.query(chrom, beg, end)
    recs, by_rs = {}, {}
    for ln in rows:
        f = ln.split("\t")
        if len(f) < 11:
            continue
        try:
            rec = {
                "rsid": f[F["rsid"]], "pos": int(f[F["pos"]]),
                "ref": f[F["ref"]], "alt": f[F["alt"]],
                "p": float(f[F["p"]]), "beta": float(f[F["beta"]]),
                "se": float(f[F["se"]]), "af": float(f[F["af"]]),
            }
        except (ValueError, IndexError):
            continue
        recs[rec["pos"]] = rec
        if rec["rsid"] and rec["rsid"] != ".":
            by_rs[rec["rsid"]] = rec
    print("%s 窗口 %s:%d-%d -> FinnGen 变异 %d 个（可配 rsID %d 个）"
          % (tag, chrom, beg, end, len(recs), len(by_rs)))
    return recs, by_rs


def main():
    cst3 = decode_cst3()
    print("=== CST3 蛋白层 pQTL（deCODE, n=35,559）===")
    for c in cst3:
        print("  %-12s %-5s pos=%-10s MAF=%-6s beta=%+7.4f -log10P=%-8s LDclass=%s"
              % (c["variant"], c["cis_trans"], c["pos"], c["maf"][:5],
                 float(c["beta"]), c["mlogp"][:7], c["ld_size"]))

    htn, htn_rs = finngen(CST3_WIN, "CST3")
    htn15, _ = finngen(FURIN_WIN, "FURIN/FES")

    print("\n=== CST3 cis-pQTL 变异在 FinnGen R12 高血压中的表现（Wald ratio）===")
    out = []
    for c in cst3:
        if c["cis_trans"].lower() != "cis":
            continue
        hit = htn_rs.get(c["variant"])
        if not hit:
            print("  %-12s FinnGen 窗口内无此 rsID（未通过质控或未定型）" % c["variant"])
            continue
        bx, sy = float(c["beta"]), hit["se"]
        ratio = hit["beta"] / bx
        se_r = abs(sy / bx)
        z = ratio / se_r
        p_r = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))
        print("  %-12s pQTL beta=%+.4f | HTN beta=%+.4f (P=%.2e) | OR=%.3f (%.3f-%.3f) P=%.3f"
              % (c["variant"], bx, hit["beta"], hit["p"], exp(ratio),
                 exp(ratio - 1.96 * se_r), exp(ratio + 1.96 * se_r), p_r))
        out.append({"variant": c["variant"], "pos": c["pos"], "pqtl_beta": bx,
                    "mlogp": c["mlogp"], "htn_beta": hit["beta"],
                    "htn_p": hit["p"], "or": exp(ratio), "p_wald": p_r,
                    "ld_size": c["ld_size"], "ld_class": c["ld_class"]})

    if not out:
        raise SystemExit("CST3: 没有可配对的 pQTL 变异，无法计算 Wald ratio")

    with open(os.path.join(HERE, "cst3_protein_mr.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    # ---- LD overlap between the lead pQTL and the hypertension signal ----
    lead = max(out, key=lambda r: float(r["mlogp"]))
    ld = set(x for x in lead["ld_class"].split(",") if x)
    print("\n=== 最强 CST3 pQTL（%s）LD class 与高血压信号的交接 ===" % lead["variant"])
    print("  LD class 变异数: %d" % len(ld))
    sig = [r for r in htn.values() if r["p"] < 5e-8]
    print("  高血压窗口内全基因组显著变异: %d" % len(sig))
    overlap = [r for r in sig if r["rsid"] in ld]
    print("  与 pQTL LD class 相交的显著高血压变异: %d" % len(overlap))
    for r in overlap[:8]:
        print("     %-12s pos=%d P=%.2e" % (r["rsid"], r["pos"], r["p"]))
    if sig:
        top = min(sig, key=lambda r: r["p"])
        print("  高血压窗口最强: %-12s pos=%d P=%.2e" % (top["rsid"], top["pos"], top["p"]))
        print("  它是否在 pQTL LD class 内: %s" % (top["rsid"] in ld))

    # ---- FURIN window: confirm the hypertension signal and the absent protein arm ----
    print("\n=== FURIN/FES 窗口高血压信号（供对照）===")
    sig15 = sorted([r for r in htn15.values() if r["p"] < 5e-8],
                   key=lambda r: r["p"])[:5]
    for r in sig15:
        print("  %-12s pos=%d P=%.2e" % (r["rsid"], r["pos"], r["p"]))


if __name__ == "__main__":
    main()
