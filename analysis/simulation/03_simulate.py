# -*- coding: utf-8 -*-
"""Operating characteristics of the two-control colocalization design.

Every simulated replicate uses the real design of one of the 240 tests: the same
number of shared variants and the same per-variant sampling variances, which are
what set the power of the test. Only the ground truth is imposed.

Scenarios
    H0   neither trait has a causal variant in the window
    H1   eQTL only          H2  outcome only
    H3   two distinct causal variants (the coincidence case)
    H4   one shared causal variant

For H3 the two causal variants may be in linkage disequilibrium: a pair in LD
rho gives the eQTL a marginal signal rho*Z at the outcome's causal variant and
vice versa, which is what a tag variant actually looks like.

Two rules are compared, as prespecified in the manuscript:
    uncalibrated  call "shared" whenever PPH4 >= 0.8, and call everything else
                  negative
    calibrated    call "shared" when PPH4 >= 0.8; call "negative" only when the
                  design has demonstrated that it can reach that threshold (a
                  positive control in the same tissue does); otherwise call the
                  test uninformative

The Python posteriors reproduce coloc 5.2.3 to 3.7e-15 given the same input rows
(see 02_validate.R).
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

P1 = P2 = 1e-4
P12 = 1e-5
SD_O = 0.2
THRESH = 0.8
GWAS_Z = 5.45
Z_GRID = [4.0, 5.45, 7.0]
RHO_GRID = [0.0, 0.5, 0.8]
B_MAIN = 1000
B_GRID = 400
RNG = np.random.default_rng(20261115)


def logsum(x):
    m = np.max(x, axis=-1, keepdims=True)
    return (m + np.log(np.sum(np.exp(x - m), axis=-1, keepdims=True)))[..., 0]


def labf(z, V, sd_prior):
    r = sd_prior ** 2 / (sd_prior ** 2 + V)
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def pph4(L1, L2):
    """PPH4 for a batch: L1, L2 are (B, m) arrays, returns (B,)."""
    s1, s2 = logsum(L1), logsum(L2)
    ls = logsum(L1 + L2)
    d = ls - (s1 + s2)
    with np.errstate(divide="ignore", invalid="ignore"):
        l3 = np.log(P1) + np.log(P2) + (s1 + s2) + np.log1p(-np.exp(np.minimum(d, 0.0)))
    lh = np.stack([np.zeros_like(s1),
                   np.log(P1) + s1,
                   np.log(P2) + s2,
                   l3,
                   np.log(P12) + ls], axis=-1)
    t = logsum(lh)
    return np.exp(lh[..., 4] - t)


def sdy_est(V, maf, n):
    x = 1.0 / V
    y = 2.0 * n * maf * (1.0 - maf)
    return float(np.sqrt(np.sum(x * y) / np.sum(x * x)))


def simulate_design(Ve, Vo, maf, n_e, B, scen, Z=GWAS_Z, rho=0.0):
    """Return PPH4 replicates (B,) for one design under one scenario."""
    m = len(Ve)
    sd_e = 0.15 * sdy_est(Ve, maf, n_e)
    ze = RNG.standard_normal((B, m))
    zo = RNG.standard_normal((B, m))

    if scen == "H0":
        pass
    elif scen == "H1":
        j = RNG.integers(0, m)
        ze[:, j] = Z * RNG.choice([-1.0, 1.0], size=B)
    elif scen == "H2":
        j = RNG.integers(0, m)
        zo[:, j] = Z * RNG.choice([-1.0, 1.0], size=B)
    elif scen == "H4":
        j = RNG.integers(0, m)
        s = RNG.choice([-1.0, 1.0], size=B)
        ze[:, j] = Z * s
        zo[:, j] = Z * s
    elif scen == "H3":
        # two distinct causal variants, separated so that they are never the same
        idx = RNG.integers(0, m, size=(B, 2))
        same = idx[:, 0] == idx[:, 1]
        idx[same, 1] = (idx[same, 1] + 1) % m
        s1 = RNG.choice([-1.0, 1.0], size=B)
        s2 = RNG.choice([-1.0, 1.0], size=B)
        ze[np.arange(B), idx[:, 0]] = Z * s1
        zo[np.arange(B), idx[:, 1]] = Z * s2
        if rho > 0:
            # tag signal: each causal variant leaves rho*Z at the other one
            ze[np.arange(B), idx[:, 1]] += rho * Z * s2
            zo[np.arange(B), idx[:, 0]] += rho * Z * s1
    else:
        raise ValueError(scen)

    L1 = labf(ze, Ve[None, :], sd_e)
    L2 = labf(zo, Vo[None, :], SD_O)
    return pph4(L1, L2)


def main():
    tests = pd.read_csv(os.path.join(DATA, "r_coloc_240.csv"))
    var = pd.read_csv(os.path.join(DATA, "r_variants_240.csv"))
    blocks = {tid: (g["V_e"].to_numpy(float), g["V_o"].to_numpy(float),
                    g["maf"].to_numpy(float), g["z_e"].to_numpy(float),
                    g["z_o"].to_numpy(float))
              for tid, g in var.groupby("test_id")}
    print("designs:", len(tests), "variants:", len(var))

    rows = []
    for _, t in tests.iterrows():
        tid = int(t["test_id"])
        Ve, Vo, maf, ze_obs, zo_obs = blocks[tid]
        rec = dict(test_id=tid, gene=t["gene"], tissue=t["tissue"],
                   pheno=t["pheno"], m=len(Ve), n_eqtl=t["n_eqtl"],
                   PPH4_obs=t["PPH4"],
                   zmax_e=float(np.abs(ze_obs).max()),
                   zmax_o=float(np.abs(zo_obs).max()))
        rec["informative_obs"] = (rec["zmax_e"] >= GWAS_Z) and (rec["zmax_o"] >= GWAS_Z)

        for scen in ("H0", "H1", "H2", "H4"):
            f = simulate_design(Ve, Vo, maf, t["n_eqtl"], B_MAIN, scen)
            rec["PPH4_%s_med" % scen] = float(np.median(f))
            rec["call_%s" % scen] = float(np.mean(f >= THRESH))
        for rho in RHO_GRID:
            f = simulate_design(Ve, Vo, maf, t["n_eqtl"], B_MAIN, "H3", rho=rho)
            rec["call_H3_rho%.1f" % rho] = float(np.mean(f >= THRESH))
            rec["PPH4_H3_rho%.1f_med" % rho] = float(np.median(f))
        rows.append(rec)

    d = pd.DataFrame(rows)
    d["power"] = d["call_H4"]          # P(detect a true shared variant)
    d["fpr_ld0"] = d["call_H3_rho0.0"]
    d["fpr_ld5"] = d["call_H3_rho0.5"]
    d["fpr_ld8"] = d["call_H3_rho0.8"]
    d.to_csv(os.path.join(DATA, "sim_by_design.csv"), index=False)

    # ---- per-tissue summary -------------------------------------------------
    g = d.groupby("tissue").agg(
        n=("test_id", "size"), n_eqtl=("n_eqtl", "median"), m=("m", "median"),
        power=("power", "median"), frac_power80=("power", lambda x: (x >= 0.8).mean()),
        fpr_ld0=("fpr_ld0", "mean"), fpr_ld5=("fpr_ld5", "mean"),
        fpr_ld8=("fpr_ld8", "mean")).reset_index()
    g.to_csv(os.path.join(DATA, "sim_by_tissue.csv"), index=False)

    # ---- effect-size / LD grid ---------------------------------------------
    grid = []
    for Z in Z_GRID:
        for rho in RHO_GRID:
            f4, f3 = [], []
            for _, t in tests.iterrows():
                Ve, Vo, maf, _, _ = blocks[int(t["test_id"])]
                f4.append(np.mean(simulate_design(Ve, Vo, maf, t["n_eqtl"], B_GRID, "H4", Z=Z) >= THRESH))
                f3.append(np.mean(simulate_design(Ve, Vo, maf, t["n_eqtl"], B_GRID, "H3", Z=Z, rho=rho) >= THRESH))
            grid.append(dict(Z=Z, rho=rho, power_H4=float(np.mean(f4)),
                             fpr_H3=float(np.mean(f3))))
            print("  Z=%.2f rho=%.1f  power=%.3f  fpr=%.3f"
                  % (Z, rho, grid[-1]["power_H4"], grid[-1]["fpr_H3"]))
    pd.DataFrame(grid).to_csv(os.path.join(DATA, "sim_grid.csv"), index=False)

    # ---- headline numbers ---------------------------------------------------
    print("\n== headline ==")
    print("designs with >=80%% power to detect a shared variant: %d/240 (%.0f%%)"
          % ((d.power >= 0.8).sum(), 100 * (d.power >= 0.8).mean()))
    print("median power across designs: %.3f" % d.power.median())
    print("FPR under H3, LD r=0   : %.4f" % d.fpr_ld0.mean())
    print("FPR under H3, LD r=0.5 : %.4f" % d.fpr_ld5.mean())
    print("FPR under H3, LD r=0.8 : %.4f" % d.fpr_ld8.mean())
    print("real tests informative  : %d/240" % d.informative_obs.sum())
    print("\nby tissue:")
    print(g.round(3).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
