# -*- coding: utf-8 -*-
"""Figure 5: operating characteristics of the calibrated design.

Panel a  probability of detecting a shared causal variant, by tissue
Panel b  false-positive rate when the two causal variants are in linkage
         disequilibrium
Panel c  posterior obtained per design under each scenario, against the two
         control values (ACE 0.730, FES 0.985)

Springer requirements met here: 174 mm width, height <= 234 mm, RGB, >= 300 dpi,
no title or legend burnt into the image, panel letters in the upper left corner.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "fig")

MM = 1 / 25.4
W, H = 174 * MM, 150 * MM
DPI = 400
INK = "#1a1a1a"
GRID = "#d9d9d9"
C_BAR = "#7a8fa6"
C_HI = "#b5482f"
C_L1, C_L2, C_L3 = "#9db4c8", "#5c7f9e", "#2f4f6f"

# Guidelines: Helvetica or Arial, 8-12 pt, no titles inside the artwork.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
})


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def main():
    d = pd.read_csv(os.path.join(DATA, "sim_by_design.csv"))
    grid = pd.read_csv(os.path.join(DATA, "sim_grid.csv"))

    fig = plt.figure(figsize=(W, H))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.08],
                  hspace=0.55, wspace=0.34,
                  left=0.245, right=0.975, top=0.935, bottom=0.115)

    # ---------------- a: power by tissue ------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    order = (d.groupby("tissue")["power"].mean().sort_values().index.tolist())
    SHORT = {"Whole_Blood": "Whole blood",
             "Heart_Left_Ventricle": "Heart left ventricle",
             "Colon_Sigmoid": "Colon sigmoid",
             "Small_Intestine_Terminal_Ileum": "Terminal ileum",
             "Colon_Transverse": "Colon transverse",
             "Kidney_Cortex": "Kidney cortex"}
    lbl = [SHORT[t] for t in order]
    y = np.arange(len(order))
    for i, t in enumerate(order):
        sub = d[d.tissue == t]
        lo, hi = wilson(int((sub.power >= 0.8).sum()), len(sub))
        v = (sub.power >= 0.8).mean()
        col = C_HI if t == "Kidney_Cortex" else C_BAR
        ax.barh(i, v, height=0.62, color=col, edgecolor="white", linewidth=0.4, zorder=3)
        ax.plot([lo, hi], [i, i], color=INK, lw=0.8, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(lbl, fontsize=8)
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Designs reaching 80% power", fontsize=8)
    ax.axvline(0.8, color=INK, lw=0.6, ls=(0, (2, 2)), zorder=2)
    ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.30, 1.05, "a", transform=ax.transAxes, fontsize=11, fontweight="bold")

    # ---------------- b: false positives under LD ---------------------------
    ax = fig.add_subplot(gs[0, 1])
    rhos = [0.0, 0.5, 0.8]
    x = np.arange(3)
    series = [("All designs", d, C_BAR), ("Kidney cortex", d[d.tissue == "Kidney_Cortex"], C_HI)]
    for k, (name, sub, col) in enumerate(series):
        vals = [sub["fpr_ld%d" % round(r*10)].mean() for r in rhos]
        off = (k - 0.5) * 0.22
        ax.bar(x + off, vals, width=0.20, color=col, edgecolor="white",
               linewidth=0.4, label=name, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["0", "0.5", "0.8"])
    ax.set_xlabel("LD between causal variants (r)", fontsize=8)
    ax.set_ylabel("False-positive rate", fontsize=8)
    ax.set_ylim(0, 0.85)
    ax.legend(frameon=False, fontsize=8, loc="upper left", handlelength=1.0)
    ax.grid(axis="y", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.30, 1.05, "b", transform=ax.transAxes, fontsize=11, fontweight="bold")

    # ---------------- c: posteriors per design, all scenarios ---------------
    ax = fig.add_subplot(gs[1, :])
    keys = [("PPH4_H3_rho0.0_med", "Distinct\nr = 0", C_L1),
            ("PPH4_H3_rho0.5_med", "Distinct\nr = 0.5", C_L2),
            ("PPH4_H3_rho0.8_med", "Distinct\nr = 0.8", C_L3),
            ("PPH4_H4_med", "Shared", "#4a7c59")]
    data = [d[k].to_numpy() for k, _, _ in keys]
    parts = ax.violinplot(data, positions=np.arange(4), widths=0.78,
                          showextrema=False, showmedians=False)
    for body, (_, _, col) in zip(parts["bodies"], keys):
        body.set_facecolor(col)
        body.set_edgecolor("white")
        body.set_linewidth(0.4)
        body.set_alpha(0.85)
    for i, arr in enumerate(data):
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        ax.vlines(i, q1, q3, color=INK, lw=3.0, zorder=4)
        ax.plot([i], [med], marker="_", color="white", ms=7, mew=1.2, zorder=5)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels([lbl for _, lbl, _ in keys], fontsize=8)
    ax.set_ylabel("Posterior for a shared variant (PPH4)", fontsize=8)
    ax.set_xlim(-0.62, 4.62)
    ax.set_ylim(-0.03, 1.06)
    for v, t, ls in ((0.985, "positive control 0.985", (0, (5, 2))),
                     (0.800, "threshold 0.80", (0, (2, 2))),
                     (0.730, "negative control 0.730", (0, (5, 2)))):
        ax.axhline(v, color=INK, lw=0.6, ls=ls, zorder=2)
        ax.text(-0.55, v + 0.028, t, fontsize=8, va="bottom", ha="left", color="#333333",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5), zorder=6)
    ax.grid(axis="y", color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.075, 1.06, "c", transform=ax.transAxes, fontsize=11, fontweight="bold")
    ax.set_xlabel("Ground truth across the 240 designs", fontsize=8, labelpad=4)

    os.makedirs(FIG, exist_ok=True)
    png = os.path.join(FIG, "Figure5_calibration.png")
    fig.savefig(png, dpi=DPI, facecolor="white")
    plt.close(fig)

    im = Image.open(png).convert("RGB")
    w, h = im.size
    dpi = max(300, w / (174 / 25.4), h / (234 / 25.4))
    out = os.path.join(FIG, "Figure5_calibration.tif")
    im.save(out, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    print("size px:", w, h, " width mm: %.1f  height mm: %.1f"
          % (w / dpi * 25.4, h / dpi * 25.4), " dpi:", round(dpi))


if __name__ == "__main__":
    main()
