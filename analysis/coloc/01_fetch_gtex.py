"""
01_fetch_gtex.py —— 8 基因 × 目标组织的全量 cis-eQTL 抽取
复用 P0_mr/_tabix.py 与已缓存的 .tbi；只补下缺失的组织索引。

窗口：基因体 ± 100 kb（GRCh38）。AUTS2（1.2 Mb）用体 ± 50 kb 以免区间过大，
      并在日志里显式说明（AUTS2 是全基因组最长基因之一，窗口策略必须交代）。

用法：
  python 01_fetch_gtex.py Colon_Sigmoid
  python 01_fetch_gtex.py Colon_Sigmoid Colon_Transverse Small_Intestine_Terminal_Ileum
"""
import os, sys, csv, time

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/"
P0MR = r"F:/longxiawenjian/gaoxueyaqiguanzhou/00_Shared_Resources/P0_mr"
TBI_SRC = os.path.join(P0MR, "data", "gtex_tbi")
ROOT = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
TBI = os.path.join(ROOT, "03_data_raw", "gtex_tbi")
OUT = os.path.join(ROOT, "03_data_raw", "gtex_full")
os.makedirs(TBI, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, P0MR)
from _tabix import Tabix

# Ensembl REST 2026-07 实测坐标（GRCh38）
GENES = {
    "LRP1":   ("chr12", 57128476, 57213394, "ENSG00000123384"),
    "CTSD":   ("chr11", 1752752, 1764573, "ENSG00000117984"),
    "MTHFR":  ("chr1", 11785723, 11806455, "ENSG00000177000"),
    "AUTS2":  ("chr7", 69598296, 70793506, "ENSG00000158321"),
    "FURIN":  ("chr15", 90868588, 90883564, "ENSG00000140564"),
    "CST3":   ("chr20", 23626706, 23638556, "ENSG00000101439"),
    "FCER1G": ("chr1", 161215175, 161220709, "ENSG00000158869"),
    "TBXAS1": ("chr7", 139777051, 140020325, "ENSG00000059377"),
}
PAD = {"AUTS2": 50000}          # 超长基因收窄窗口
DEFAULT_PAD = 100000


def get_tbi(tissue):
    for d in (TBI, TBI_SRC):
        p = os.path.join(d, "%s.tbi" % tissue)
        if os.path.exists(p) and os.path.getsize(p) > 1_000_000:
            return p
    import urllib.request
    p = os.path.join(TBI, "%s.tbi" % tissue)
    for k in range(4):
        try:
            req = urllib.request.Request(BASE + "%s.tsv.gz.tbi" % tissue,
                                         headers={"User-Agent": "Mozilla/5.0"})
            d = urllib.request.urlopen(req, timeout=600).read()
            if len(d) < 1_000_000:
                raise IOError("short %d" % len(d))
            open(p + ".part", "wb").write(d)
            os.replace(p + ".part", p)
            return p
        except Exception:
            if k == 3:
                raise
            time.sleep(5)


def fetch(tissue):
    tb = Tabix(BASE + "%s.tsv.gz" % tissue, tbi_path=get_tbi(tissue))
    log = []
    for g, (chrom, s, e, ensg) in GENES.items():
        pad = PAD.get(g, DEFAULT_PAD)
        t0 = time.time()
        rows = tb.query(chrom, s - pad, e + pad)
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
        dst = os.path.join(OUT, "%s_%s.csv" % (g, tissue))
        if recs:
            with open(dst, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
                w.writeheader()
                w.writerows(recs)
            best = min(recs, key=lambda x: x["p"])
            line = "%-8s %-32s n=%5d minP=%.3e %-14s TPM=%8.2f beta=%+.3f %4.0fs" % (
                g, tissue, len(recs), best["p"], best["rsid"],
                recs[0]["median_tpm"], best["beta"], time.time() - t0)
        else:
            line = "%-8s %-32s n=0 %4.0fs" % (g, tissue, time.time() - t0)
        print(line); sys.stdout.flush()
        log.append(line)
    with open(os.path.join(ROOT, "05_results", "tables", "fetch_log.txt"),
              "a", encoding="utf-8") as fh:
        fh.write("\n## %s  %s\n" % (tissue, time.strftime("%Y-%m-%d %H:%M")))
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    for t in (sys.argv[1:] or ["Colon_Sigmoid"]):
        fetch(t)
