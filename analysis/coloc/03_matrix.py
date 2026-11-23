"""
03_matrix.py —— 从 gtex_full/*.csv 汇总 8 基因 x 6 组织的表达-调控矩阵

输出：
  05_results/tables/matrix_raw.csv   长表（gene, tissue, n, minP, best_snp, tpm, beta）
  05_results/tables/matrix_wide.md   宽表（minP 与 TPM 各一张，Markdown）
不产出任何解释性文字——解释写在报告里，避免数字与文字分家。
"""
import os, csv, glob, collections

ROOT = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
SRC = os.path.join(ROOT, "03_data_raw", "gtex_full")
OUT = os.path.join(ROOT, "05_results", "tables")
os.makedirs(OUT, exist_ok=True)

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]
TISSUES = ["Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum",
           "Heart_Left_Ventricle", "Kidney_Cortex", "Whole_Blood"]
SHORT = {"Colon_Sigmoid": "结肠乙状", "Colon_Transverse": "结肠横",
         "Small_Intestine_Terminal_Ileum": "回肠末端",
         "Heart_Left_Ventricle": "左心室", "Kidney_Cortex": "肾皮质",
         "Whole_Blood": "全血"}
N = {"Colon_Sigmoid": 318, "Colon_Transverse": 368, "Small_Intestine_Terminal_Ileum": 174,
     "Heart_Left_Ventricle": 386, "Kidney_Cortex": 73, "Whole_Blood": 670}

rec = {}
for f in glob.glob(os.path.join(SRC, "*.csv")):
    base = os.path.basename(f)[:-4]
    for t in TISSUES:
        if base.endswith("_" + t):
            g = base[:-(len(t) + 1)]
            break
    else:
        continue
    if g not in GENES:
        continue
    rows = list(csv.DictReader(open(f, encoding="utf-8")))
    if not rows:
        continue
    best = min(rows, key=lambda x: float(x["p"]))
    rec[(g, t)] = dict(n=len(rows), minp=float(best["p"]), snp=best["rsid"],
                       tpm=float(best["median_tpm"]), beta=float(best["beta"]))

# 长表
with open(os.path.join(OUT, "matrix_raw.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "tissue", "eqtl_n", "nsnp", "minP", "best_snp", "median_tpm", "beta"])
    for g in GENES:
        for t in TISSUES:
            r = rec.get((g, t))
            if r:
                w.writerow([g, t, N[t], r["n"], "%.3e" % r["minp"], r["snp"],
                            "%.2f" % r["tpm"], "%.3f" % r["beta"]])

def fmt_p(p):
    if p < 5e-8:
        return "**%.2e**" % p          # eGene
    if p < 1e-5:
        return "%.2e" % p
    return "%.1e" % p

L = []
L.append("### 表 A　各组织内最小 cis-eQTL P 值（**加粗 = eGene，P<5×10⁻⁸**）\n")
L.append("| 基因 | " + " | ".join("%s (n=%d)" % (SHORT[t], N[t]) for t in TISSUES) + " |")
L.append("|---|" + "---|" * len(TISSUES))
for g in GENES:
    cells = []
    for t in TISSUES:
        r = rec.get((g, t))
        cells.append(fmt_p(r["minp"]) if r else "—")
    L.append("| **%s** | %s |" % (g, " | ".join(cells)))

L.append("\n### 表 B　中位表达量（TPM）\n")
L.append("| 基因 | " + " | ".join(SHORT[t] for t in TISSUES) + " |")
L.append("|---|" + "---|" * len(TISSUES))
for g in GENES:
    cells = []
    for t in TISSUES:
        r = rec.get((g, t))
        cells.append("%.1f" % r["tpm"] if r else "—")
    L.append("| **%s** | %s |" % (g, " | ".join(cells)))

L.append("\n### 表 C　eGene 计数（P<5×10⁻⁸）\n")
L.append("| 基因 | " + " | ".join(SHORT[t] for t in TISSUES) + " | 合计 |")
L.append("|---|" + "---|" * len(TISSUES) + "---|")
tot = collections.Counter()
for g in GENES:
    cells, c = [], 0
    for t in TISSUES:
        r = rec.get((g, t))
        k = bool(r and r["minp"] < 5e-8)
        cells.append("✓" if k else "·")
        c += k
    tot[g] = c
    L.append("| **%s** | %s | %d |" % (g, " | ".join(cells), c))

open(os.path.join(OUT, "matrix_wide.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
print("\n已写入", OUT)
