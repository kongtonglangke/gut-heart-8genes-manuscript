"""
M3-1  暴露侧窗口 -> 结局侧同窗口取数（FinnGen R12，GRCh38）

为什么用 tabix + HTTP Range 而不是 FinnGen API：
  实测 r12.finngen.fi/api/variant/{chr}-{pos}-{ref}-{alt} 对每个表型只保留
  约 8% 的变异（按结局显著性筛选），MAF 0.37-0.45 的常见变异同样为 null，
  不能用于 coloc（coloc 要求区域内变异两侧齐全）。改走原始 sumstats：
  本地已有与远端字节数完全一致的 .tbi（1,777,032 B），用 P0_mr/_tabix.py
  做 HTTP Range 区域检索，只取命中的压缩块。

FinnGen R12 列序（实测 + 与 P0_mr/30_region_scan.py 一致）：
  0 chrom  1 pos  2 ref  3 alt  4 rsids  5 nearest_genes
  6 pval   7 mlogp  8 beta  9 sebeta  10 af_alt ...

用法：
  python 05_fetch_outcome.py I9_HYPTENS
"""
import csv
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr")
from _tabix import Tabix  # noqa: E402

GTEX = os.path.join(ROOT, "03_data_raw", "gtex_full")
OUTD = os.path.join(ROOT, "03_data_raw", "outcome")
TBI_DIR = os.path.join(ROOT, "03_data_raw", "finngen_tbi")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(TBI_DIR, exist_ok=True)

GCS = ("https://storage.googleapis.com/finngen-public-data-r12/"
       "summary_stats/release/")
# 本地已有的 .tbi（P0_mr/data），字节数与远端一致，直接复用，避免重复下载
LOCAL_TBI = {"I9_HYPTENS": "F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr/data/finngen_htn.tbi"}

CHROM = {"LRP1": "12", "CTSD": "11", "MTHFR": "1", "AUTS2": "7",
         "FURIN": "15", "CST3": "20", "FCER1G": "1", "TBXAS1": "7"}

PAD = 50000  # 两侧各外扩 50 kb，覆盖调控区
IDX = dict(pos=1, ref=2, alt=3, rsid=4, gene=5, p=6, beta=8, se=9, af=10)


def gtex_window(gene):
    """取该基因在结肠乙状部的变异范围作为窗口基准。"""
    path = os.path.join(GTEX, "%s_Colon_Sigmoid.csv" % gene)
    pos = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            pos.append(int(r["pos"]))
    return min(pos), max(pos)


def tbi_path(pheno):
    if pheno in LOCAL_TBI and os.path.exists(LOCAL_TBI[pheno]):
        return LOCAL_TBI[pheno]
    return os.path.join(TBI_DIR, "finngen_R12_%s.gz.tbi" % pheno)


def parse(lines):
    rows = []
    for ln in lines:
        f = ln.split("\t")
        try:
            rows.append({
                "chrom": f[0], "pos": int(f[IDX["pos"]]),
                "ref": f[IDX["ref"]], "alt": f[IDX["alt"]],
                "rsid": f[IDX["rsid"]], "gene": f[IDX["gene"]],
                "p": float(f[IDX["p"]]), "beta": float(f[IDX["beta"]]),
                "se": float(f[IDX["se"]]), "af": float(f[IDX["af"]]),
            })
        except (ValueError, IndexError):
            continue
    return rows


def main():
    pheno = sys.argv[1] if len(sys.argv) > 1 else "I9_HYPTENS"
    tb_path = tbi_path(pheno)
    if not os.path.exists(tb_path):
        raise SystemExit("缺索引: %s（需先下载远端 .tbi）" % tb_path)
    tb = Tabix(GCS + "finngen_R12_%s.gz" % pheno, tbi_path=tb_path)
    print("[%s] 索引载入，refs=%d" % (pheno, len(tb.names)), flush=True)

    for gene, chrom in CHROM.items():
        a, b = gtex_window(gene)
        beg, end = a - PAD, b + PAD
        t0 = time.time()
        try:
            rows = parse(tb.query(chrom, beg, end))
        except Exception as e:  # noqa: BLE001
            print("%-7s FAIL %s %s" % (gene, type(e).__name__, str(e)[:80]), flush=True)
            continue
        out = os.path.join(OUTD, "%s_%s.csv" % (gene, pheno))
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["chrom", "pos", "ref", "alt",
                                               "rsid", "gene", "p", "beta", "se", "af"])
            w.writeheader()
            w.writerows(rows)
        minp = min((r["p"] for r in rows), default=float("nan"))
        print("%-7s chr%-2s %d-%d  n=%4d  minP=%.2e  %.1fs"
              % (gene, chrom, beg, end, len(rows), minp, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
