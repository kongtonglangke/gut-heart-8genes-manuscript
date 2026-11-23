# -*- coding: utf-8 -*-
"""Posterior predictive check: is the observed PPH4 consistent with H3?

The anchor simulation fixes the causal variant at genome-wide significance
(z = 5.45). That is a fair operating point for a typical test, but the observed
signal strengths differ by more than an order of magnitude between tests, and a
stronger signal raises PPH4 under every hypothesis. Here the causal effect is
therefore set to the signal actually observed in that test, separately for each
side, so that the simulated H3 distribution can be compared with the observed
posterior directly.

For each of the 240 designs the script reports
    p_H3   fraction of H3 replicates reaching at least the observed PPH4
    p_H4   fraction of H4 replicates reaching at least the observed PPH4
A high p_H3 means the observed posterior is unremarkable when the two traits have
distinct causal variants, i.e. the data do not require a shared variant.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
P1 = P2 = 1e-4
P12 = 1e-5
SD_O = 0.2
B = 2000
RNG = np.random.default_rng(20261116)


def logsum(x):
    m = np.max(x, axis=-1, keepdims=True)
    return (m + np.log(np.sum(np.exp(x - m), axis=-1, keepdims=True)))[..., 0]


def labf(z, V, sd_prior):
    r = sd_prior ** 2 / (sd_prior ** 2 + V)
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def pph4(L1, L2):
    s1, s2 = logsum(L1), logsum(L2)
    ls = logsum(L1 + L2)
    with np.errstate(divide="ignore", invalid="ignore"):
        l3 = np.log(P1) + np.log(P2) + (s1 + s2) + np.log1p(-np.exp(np.minimum(ls - (s1 + s2), 0.0)))
    lh = np.stack([np.zeros_like(s1), np.log(P1) + s1, np.log(P2) + s2,
                   l3, np.log(P12) + ls], axis=-1)
    t = logsum(lh)
    return np.exp(lh[..., 4] - t)


def sdy_est(V, maf, n):
    x, y = 1.0 / V, 2.0 * n * maf * (1.0 - maf)
    return float(np.sqrt(np.sum(x * y) / np.sum(x * x)))


def run(Ve, Vo, maf, n_e, ze, zo, scen, rho=0.0):
    """Place the causal effects at the observed lead variants.

    Strength matters as much as position: the effect given to the causal variant
    is the strongest signal actually observed in that trait, so the simulated
    posterior is on the same scale as the observed one. Randomising the position
    instead would drop both effects on a null variant and drive PPH4 to zero for
    reasons that have nothing to do with the hypothesis.
    """
    m = len(Ve)
    sd_e = 0.15 * sdy_est(Ve, maf, n_e)
    Z = RNG.standard_normal((B, m))
    W = RNG.standard_normal((B, m))
    ar = np.arange(B)
    ie = int(np.argmax(np.abs(ze)))
    io = int(np.argmax(np.abs(zo)))
    se_, so_ = np.sign(ze[ie]), np.sign(zo[io])
    ze_max, zo_max = abs(ze[ie]), abs(zo[io])

    if scen == "H3":
        if io == ie:
            return None            # the two leads coincide: H3 is not reachable
        Z[ar, ie] = ze_max * se_
        W[ar, io] = zo_max * so_
        if rho:                    # tag signal at the other lead
            Z[ar, io] += rho * zo_max * so_
            W[ar, ie] += rho * ze_max * se_
    elif scen == "H4":
        Z[ar, ie] = ze_max * se_
        W[ar, ie] = zo_max * se_
    else:
        raise ValueError(scen)
    return pph4(labf(Z, Ve[None, :], sd_e), labf(W, Vo[None, :], SD_O))


def main():
    tests = pd.read_csv(os.path.join(DATA, "r_coloc_240.csv"))
    var = pd.read_csv(os.path.join(DATA, "r_variants_240.csv"))
    blocks = {t: (g["V_e"].to_numpy(float), g["V_o"].to_numpy(float),
                  g["maf"].to_numpy(float), g["z_e"].to_numpy(float), g["z_o"].to_numpy(float))
              for t, g in var.groupby("test_id")}
    rows = []
    for _, t in tests.iterrows():
        tid = int(t["test_id"])
        Ve, Vo, maf, ze, zo = blocks[tid]
        f3 = run(Ve, Vo, maf, t["n_eqtl"], ze, zo, "H3", rho=0.0)
        f4 = run(Ve, Vo, maf, t["n_eqtl"], ze, zo, "H4")
        obs = float(t["PPH4"])
        if f3 is None:
            continue
        rows.append(dict(test_id=tid, gene=t["gene"], tissue=t["tissue"], pheno=t["pheno"],
                         PPH4_obs=obs, PPH4_H3_med=float(np.median(f3)),
                         PPH4_H4_med=float(np.median(f4)),
                         p_H3=float(np.mean(f3 >= obs)), p_H4=float(np.mean(f4 >= obs)),
                         zmax_e=float(np.abs(ze).max()), zmax_o=float(np.abs(zo).max())))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(DATA, "sim_matched.csv"), index=False)

    inf = d[(d.zmax_e >= 5.45) & (d.zmax_o >= 5.45)]
    print("informative designs:", len(inf))
    print(inf[["gene", "tissue", "pheno", "PPH4_obs", "PPH4_H3_med",
               "PPH4_H4_med", "p_H3", "p_H4"]].round(3).to_string(index=False))
    print("\nall 240: median p_H3 = %.3f, median p_H4 = %.3f"
          % (d.p_H3.median(), d.p_H4.median()))
    print("tests whose observed PPH4 is above the 95th centile of H3: %d/240"
          % int((d.p_H3 < 0.05).sum()))


if __name__ == "__main__":
    main()
