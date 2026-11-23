# -*- coding: utf-8 -*-
"""Rebuild the 240 colocalization tests from raw GTEx / FinnGen files.

The original analysis (04_code/06_coloc.R) harmonised the two datasets by
position, flipped the eQTL effect when the alleles were swapped, dropped
variants with MAF outside (0.005, 0.995) and de-duplicated by position. This
script reproduces that exactly so that the simulation runs on the same designs
as the published tests, and writes the variant-level arrays the simulator needs.

Outputs (data/):
    tests_240.csv     one row per test with m, sample sizes and PPH0-PPH4
    variants_240.csv  long format, one row per variant per test
"""
import os

import numpy as np
import pandas as pd

ROOT = r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
GTEX = os.path.join(ROOT, "03_data_raw", "gtex_full")
OUTD = os.path.join(ROOT, "03_data_raw", "outcome")
DEST = os.path.join(ROOT, "09_simulation", "data")

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]
TISSUES = ["Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum",
           "Heart_Left_Ventricle", "Kidney_Cortex", "Whole_Blood"]
PHENOS = ["I9_HYPTENS", "I9_HEARTFAIL", "I9_AF", "N14_CHRONKIDNEYDIS", "T2D"]

P1 = P2 = 1e-4
P12 = 1e-5
SD_PRIOR_O = 0.2    # coloc: cc -> sd.prior fixed at 0.2

# coloc does NOT use sdY = 1 for a quantitative trait: when sdY is not supplied
# it estimates it with coloc:::sdY.est, i.e. the through-origin slope of
# 2*N*maf*(1-maf) on 1/varbeta, then takes the square root. The prior SD on beta
# is then 0.15 * sdY. Skipping this step changes every posterior, so it is
# reproduced here exactly.


def sdy_est(vbeta, maf, n):
    x = 1.0 / vbeta
    y = 2.0 * n * maf * (1.0 - maf)
    return float(np.sqrt(np.sum(x * y) / np.sum(x * x)))


def logsum(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def labf(z, V, sd_prior):
    r = sd_prior ** 2 / (sd_prior ** 2 + V)
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def logdiff(a, b):
    """log(exp(a) - exp(b)), stable; b <= a by construction."""
    d = b - a
    with np.errstate(divide="ignore", invalid="ignore"):
        return a + np.log1p(-np.exp(np.minimum(d, 0.0)))


def pph(l1, l2):
    """Posterior probabilities under the five coloc.abf hypotheses."""
    s1, s2 = logsum(l1), logsum(l2)
    ls = logsum(l1 + l2)
    l3 = np.log(P1) + np.log(P2) + logdiff(s1 + s2, ls)
    lh = np.array([0.0,
                   np.log(P1) + s1,
                   np.log(P2) + s2,
                   l3,
                   np.log(P12) + ls])
    t = logsum(lh)
    out = np.exp(lh - t)
    return np.nan_to_num(out, nan=0.0)


def read_eqtl(gene, tissue):
    f = os.path.join(GTEX, "%s_%s.csv" % (gene, tissue))
    if not os.path.exists(f):
        return None
    d = pd.read_csv(f)
    if len(d) < 50:
        return None
    d = d[np.isfinite(d["se"]) & (d["se"] > 0) & np.isfinite(d["beta"]) & np.isfinite(d["p"])]
    d = d.drop_duplicates(subset=["pos", "ref", "alt"])
    d["N"] = d["an"] / 2.0
    return d


def read_outcome(gene, pheno):
    f = os.path.join(OUTD, "%s_%s.csv" % (gene, pheno))
    if not os.path.exists(f):
        return None
    d = pd.read_csv(f)
    d = d[np.isfinite(d["se"]) & (d["se"] > 0) & np.isfinite(d["beta"]) & np.isfinite(d["p"])]
    d = d.drop_duplicates(subset=["pos", "ref", "alt"])
    return d


def harmonise(e, o):
    m = e.merge(o, on="pos", suffixes=(".e", ".o"))
    if not len(m):
        return None
    same = (m["ref.e"] == m["ref.o"]) & (m["alt.e"] == m["alt.o"])
    flip = (m["ref.e"] == m["alt.o"]) & (m["alt.e"] == m["ref.o"])
    m = m[same | flip].copy()
    if not len(m):
        return None
    f = (m["ref.e"] == m["alt.o"]) & (m["alt.e"] == m["ref.o"])
    m.loc[f, "beta.e"] = -m.loc[f, "beta.e"]
    m["maf"] = np.minimum(m["af"], 1 - m["af"])
    m = m[(m["maf"] > 0.005) & (m["maf"] < 0.995)]
    m = m.drop_duplicates(subset=["pos"])
    return m


def main():
    os.makedirs(DEST, exist_ok=True)
    tests, variants = [], []
    tid = 0
    for gene in GENES:
        for pheno in PHENOS:
            o = read_outcome(gene, pheno)
            if o is None or len(o) < 50:
                continue
            for tissue in TISSUES:
                e = read_eqtl(gene, tissue)
                if e is None:
                    continue
                m = harmonise(e, o)
                if m is None or len(m) < 50:
                    continue
                V_e = m["se.e"].to_numpy(float) ** 2
                V_o = m["se.o"].to_numpy(float) ** 2
                z_e = m["beta.e"].to_numpy(float) / m["se.e"].to_numpy(float)
                z_o = m["beta.o"].to_numpy(float) / m["se.o"].to_numpy(float)
                maf = m["maf"].to_numpy(float)
                n_e = float(np.median(m["N"]))
                sd_prior_e = 0.15 * sdy_est(V_e, maf, n_e)
                L1 = labf(z_e, V_e, sd_prior_e)
                L2 = labf(z_o, V_o, SD_PRIOR_O)
                pp = pph(L1, L2)
                tid += 1
                tests.append(dict(
                    test_id=tid, gene=gene, tissue=tissue, pheno=pheno, m=len(m),
                    n_eqtl=n_e, sdY=sd_prior_e / 0.15,
                    minP_e=float(m["p.e"].min()), minP_o=float(m["p.o"].min()),
                    zmax_e=float(np.abs(z_e).max()), zmax_o=float(np.abs(z_o).max()),
                    PPH0=pp[0], PPH1=pp[1], PPH2=pp[2], PPH3=pp[3], PPH4=pp[4]))
                variants.append(pd.DataFrame(dict(
                    test_id=tid, pos=m["pos"].to_numpy(int), maf=maf,
                    V_e=V_e, V_o=V_o, z_e=z_e, z_o=z_o)))
    t = pd.DataFrame(tests)
    v = pd.concat(variants, ignore_index=True)
    t.to_csv(os.path.join(DEST, "tests_240.csv"), index=False)
    v.to_csv(os.path.join(DEST, "variants_240.csv"), index=False)
    print("tests:", len(t), " variants:", len(v))
    print(t[["m", "n_eqtl"]].describe().loc[["min", "50%", "max"]].to_string())
    print("\nPPH4 > 0.8 :", int((t["PPH4"] > 0.8).sum()))
    print("both z>5.45:", int(((t["zmax_e"] > 5.45) & (t["zmax_o"] > 5.45)).sum()))


if __name__ == "__main__":
    main()
