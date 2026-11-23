"""
PoC：8 基因在肠道组织的全量 cis-eQTL 抽取可行性验证
—— 只做一件事：证明"数据真的能取到"，不产出任何结论。

复用 P0_mr/_tabix.py（纯 Python tabix）与已缓存的 .tbi 索引。
窗口规则：基因体 ± 100 kb（GRCh38，保留基因体内所有变异）。
"""
import os, sys, csv, time

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/"
P0MR = r"F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr"
TBI = os.path.join(P0MR, "data", "gtex_tbi")
OUT = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes/03_data_raw/gtex_full"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, P0MR)
from _tabix import Tabix

# 基因体坐标来自 Ensembl REST（2026-11-01 实测），窗口 = 体 ±100kb
GENES = {
    "FURIN":  ("chr15", 90768588, 90883564, "ENSG00000140564"),
    "CST3":   ("chr20", 23626706, 23638556, "ENSG00000101439"),
}
PAD = 100000


def run(tissue, genes):
    tbi = os.path.join(TBI, "%s.tbi" % tissue)
    assert os.path.exists(tbi), "缺索引 %s" % tbi
    tb = Tabix(BASE + "%s.tsv.gz" % tissue, tbi_path=tbi)
    for g, (chrom, s, e, ensg) in genes.items():
        t0 = time.time()
        rows = tb.query(chrom, s - PAD, e + PAD)
        recs = []
        for r in rows:
            f = r.split("\t")
            if len(f) < 17 or f[4] != ensg:
                continue
            try:
                recs.append(dict(gene=g, tissue=tissue, rsid=f[17] or f[0],
                                 pos=int(f[13]), ref=f[14], alt=f[15],
                                 p=float(f[2]), beta=float(f[8]), se=float(f[9]),
                                 median_tpm=float(f[7]), maf=float(f[5]), an=int(f[10])))
            except (ValueError, IndexError):
                continue
        p = os.path.join(OUT, "%s_%s_poc.csv" % (g, tissue))
        if recs:
            with open(p, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
                w.writeheader()
                w.writerows(recs)
        if recs:
            mp = min(x["p"] for x in recs)
            best = min(recs, key=lambda x: x["p"])
            print("%-8s %-32s n=%5d  minP=%.3e  at %s  TPM=%.2f  beta=%.3f  [%.1fs]"
                  % (g, tissue, len(recs), mp, best["rsid"], recs[0]["median_tpm"],
                     best["beta"], time.time() - t0))
        else:
            print("%-8s %-32s n=0（无记录）[%.1fs]" % (g, tissue, time.time() - t0))
        sys.stdout.flush()


if __name__ == "__main__":
    t = sys.argv[1] if len(sys.argv) > 1 else "Colon_Sigmoid"
    run(t, GENES)
