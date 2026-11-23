#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M4: cis-MR of gut eQTL -> cardiometabolic outcomes (single-instrument design).

Why single instrument
  No LD reference panel is reachable from this network, so multi-IV IVW on
  cis-eQTL would stack correlated variants and inflate precision.  We therefore
  use ONE instrument per gene x tissue: the lowest-P cis variant that is
  present and allele-harmonisable in both datasets.  Precision is lower, but
  the estimate is valid.

Instrument tiers
  strict  : P < 5e-8 ;  relaxed : P < 1e-5  (reported separately, never pooled)

Output
  05_results/tables/mr_results.csv
"""
import csv
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GTEX = os.path.join(ROOT, "03_data_raw", "gtex_full")
OUTC = os.path.join(ROOT, "03_data_raw", "outcome")
OUT = os.path.join(ROOT, "05_results", "tables")
os.makedirs(OUT, exist_ok=True)

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]
GUT = ["Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum"]
P_STRICT, P_RELAX = 5e-8, 1e-5
F_MIN = 10.0


def load_eqtl(gene, tissue):
    p = os.path.join(GTEX, "%s_%s.csv" % (gene, tissue))
    if not os.path.exists(p):
        return []
    out = []
    for r in csv.DictReader(open(p, encoding="utf-8")):
        try:
            out.append(dict(rsid=r["rsid"], pos=int(r["pos"]), ref=r["ref"].upper(),
                            alt=r["alt"].upper(), p=float(r["p"]),
                            beta=float(r["beta"]), se=float(r["se"])))
        except (ValueError, KeyError):
            pass
    return out


def load_out(gene, pheno):
    p = os.path.join(OUTC, "%s_%s.csv" % (gene, pheno))
    if not os.path.exists(p):
        return {}
    d = {}
    for r in csv.DictReader(open(p, encoding="utf-8")):
        try:
            d[int(r["pos"])] = dict(rsid=r["rsid"], ref=r["ref"].upper(),
                                    alt=r["alt"].upper(), p=float(r["p"]),
                                    beta=float(r["beta"]), se=float(r["se"]))
        except (ValueError, KeyError):
            pass
    return d


def harmonise(e, o):
    if e["ref"] == o["ref"] and e["alt"] == o["alt"]:
        return o["beta"], o["se"], "same"
    if e["ref"] == o["alt"] and e["alt"] == o["ref"]:
        return -o["beta"], o["se"], "flipped"
    return None


def pick_iv(rows, out, thr):
    """Lowest-P variant under thr that exists in outcome and harmonises, F>10."""
    cand = sorted([r for r in rows if r["p"] < thr], key=lambda r: r["p"])
    for rank, iv in enumerate(cand):
        f = (iv["beta"] / iv["se"]) ** 2
        if f < F_MIN:
            continue
        o = out.get(iv["pos"])
        if not o:
            continue
        h = harmonise(iv, o)
        if h is None:
            continue
        return iv, o, h, f, rank
    return None


def wald(b_exp, se_exp, b_out, se_out):
    b = b_out / b_exp
    se = se_out / abs(b_exp)
    z = b / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return b, se, p


def two_side(p):
    return "%.3e" % p


def main():
    phenos = sys.argv[1:] or ["I9_HYPTENS"]
    recs = []
    for pheno in phenos:
        for g in GENES:
            best = None
            for t in GUT:
                rows = load_eqtl(g, t)
                if not rows:
                    continue
                m = min(rows, key=lambda r: r["p"])
                if best is None or m["p"] < best[2]["p"]:
                    best = (t, rows, m)
            if best is None:
                recs.append(dict(pheno=pheno, gene=g, status="no eQTL data"))
                continue
            tissue, rows, m = best
            out = load_out(g, pheno)
            if not out:
                recs.append(dict(pheno=pheno, gene=g, tissue=tissue,
                                 status="no outcome data"))
                continue

            for thr, tag in ((P_STRICT, "strict_5e-8"), (P_RELAX, "relaxed_1e-5")):
                got = pick_iv(rows, out, thr)
                if got is None:
                    recs.append(dict(pheno=pheno, gene=g, tissue=tissue,
                                     iv_set=tag, n_iv=0,
                                     minP_eqtl="%.2e" % m["p"],
                                     status="no valid instrument"))
                    continue
                iv, o, (b_out, se_out, flip), f, rank = got
                b, se, p = wald(iv["beta"], iv["se"], b_out, se_out)
                src = "lead" if rank == 0 else "fallback_rank%d" % (rank + 1)
                recs.append(dict(
                    pheno=pheno, gene=g, tissue=tissue, iv_set=tag, n_iv=1,
                    rsid=iv["rsid"], pos=iv["pos"],
                    allele="%s>%s" % (iv["ref"], iv["alt"]),
                    b_exp=round(iv["beta"], 4), p_exp="%.2e" % iv["p"],
                    F=round(f, 1), b_out=round(b_out, 4),
                    p_out="%.2e" % o["p"], harmonise=flip,
                    iv_source=src,
                    beta=round(b, 5), se=round(se, 5),
                    OR=round(math.exp(b), 4),
                    ci_lo=round(math.exp(b - 1.96 * se), 4),
                    ci_hi=round(math.exp(b + 1.96 * se), 4),
                    pval=two_side(p), status="ok"))

    p = os.path.join(OUT, "mr_results.csv")
    keys = ["pheno", "gene", "tissue", "iv_set", "n_iv", "rsid", "pos", "allele",
            "b_exp", "p_exp", "F", "b_out", "p_out", "harmonise", "iv_source",
            "beta", "se", "OR", "ci_lo", "ci_hi", "pval", "status"]
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            w.writerow(r)
    print("written", p, "rows=", len(recs))
    print()
    print("%-7s %-30s %-14s %-12s %8s %20s %11s %7s %s" % (
        "gene", "tissue", "iv_set", "IV", "OR", "95%CI", "P", "F", "source"))
    for r in recs:
        if r.get("status") != "ok":
            continue
        print("%-7s %-30s %-14s %-12s %8.3f %8.3f-%-9.3f %11s %7.1f %s" % (
            r["gene"], r["tissue"], r["iv_set"], r["rsid"], r["OR"],
            r["ci_lo"], r["ci_hi"], r["pval"], r["F"], r["iv_source"]))


if __name__ == "__main__":
    main()
