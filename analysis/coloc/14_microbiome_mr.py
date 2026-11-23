#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M7: gut microbiome taxon -> hypertension, two-sample MR.

Exposure : MiBioGen (211 taxa, GRCh37) genome-wide significant hits,
           locally mirrored in P0_mr/data/mbg_allhits.txt
Outcome  : FinnGen R12 I9_HYPTENS (GRCh38), read by remote tabix

Build mismatch
  MiBioGen is GRCh37, FinnGen R12 is GRCh38.  We lift over by rsID using the
  Ensembl REST variation endpoint (batch POST, 200 ids per call) instead of a
  chain file, so no local liftOver binary is needed.

Instrument selection
  Per taxon: sort by P, greedy distance clumping at 1 Mb, keep at most 5.
  Distance clumping is a conservative substitute for LD clumping, which is not
  possible without an LD reference panel on this network.

Output
  05_results/tables/mr_microbiome_htn.csv
"""
import csv
import json
import math
import os
import sys
import time
import urllib.request
from collections import defaultdict

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou"
PROJ = os.path.join(ROOT, "Gut_Heart_RAS_8genes")
MBG = os.path.join(ROOT, "P0_mr", "data", "mbg_allhits.txt")
RES = os.path.join(PROJ, "05_results", "tables")
os.makedirs(RES, exist_ok=True)

sys.path.insert(0, os.path.join(ROOT, "P0_mr"))
from _tabix import Tabix  # noqa: E402

FINNGEN = ("https://storage.googleapis.com/finngen-public-data-r12/"
           "summary_stats/release/finngen_R12_%s.gz")
TBI_DIR = os.path.join(PROJ, "03_data_raw", "finngen_tbi")
P_TBI_LOCAL = os.path.join(ROOT, "P0_mr", "data", "finngen_htn.tbi")

PHENO = "I9_HYPTENS"
N_CASE, N_CONT = 154630, 345634
CLUMP_BP = 1_000_000
MAX_IV = 5
UA = {"User-Agent": "Mozilla/5.0"}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_mbg():
    rows = []
    with open(MBG, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                rows.append(dict(
                    taxon=r["bac"].strip('"'), chrom=r["chr"].strip('"'),
                    bp37=int(r["bp"]), rsid=r["rsID"].strip('"'),
                    ref=r["ref.allele"].strip('"').upper(),
                    eff=r["eff.allele"].strip('"').upper(),
                    beta=float(r["beta"]), se=float(r["SE"]),
                    p=float(r["P.weightedSumZ"]),
                    n=int(r["N"])))
            except (ValueError, KeyError):
                pass
    return rows


def select_ivs(rows):
    by = defaultdict(list)
    for r in rows:
        by[r["taxon"]].append(r)
    sel = []
    for taxon, xs in by.items():
        xs.sort(key=lambda x: x["p"])
        keep = []
        for x in xs:
            if all(x["chrom"] != k["chrom"] or abs(x["bp37"] - k["bp37"]) > CLUMP_BP
                   for k in keep):
                keep.append(x)
            if len(keep) >= MAX_IV:
                break
        sel.extend(keep)
    return sel


def ensembl_map(rsids, batch=180):
    out = {}
    ids = list(rsids)
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        body = json.dumps({"ids": chunk}).encode()
        for attempt in range(5):
            try:
                req = urllib.request.Request(
                    "https://rest.ensembl.org/variation/homo_sapiens",
                    data=body,
                    headers={"Content-Type": "application/json",
                             "Accept": "application/json"})
                d = json.loads(urllib.request.urlopen(req, timeout=120).read())
                out.update(d)
                break
            except Exception as e:
                print("  ensembl batch %d fail (%s), retry" % (i, type(e).__name__))
                time.sleep(4 + 4 * attempt)
        print("  mapped %d/%d" % (min(i + batch, len(ids)), len(ids)), flush=True)
        time.sleep(0.5)
    return out


def pick_mapping(v):
    """Choose the GRCh38 primary-assembly mapping, prefer SNV."""
    maps = [m for m in v.get("mappings", []) if m.get("assembly_name") == "GRCh38"]
    if not maps:
        maps = v.get("mappings", [])
    if not maps:
        return None
    snv = [m for m in maps if m.get("coord_system") == "chromosome"]
    m = (snv or maps)[0]
    al = m.get("allele_string", "")
    if "/" not in al or len(al.split("/")) != 2:
        return None
    a, b = al.split("/")
    if len(a) != 1 or len(b) != 1:
        return None
    return dict(chrom=m["seq_region_name"], pos=int(m["start"]),
                a1=a.upper(), a2=b.upper())


def allele_beta(row_ref, row_alt, row_beta, e_ref, e_eff):
    """Return beta of the MiBioGen effect allele, or None."""
    r, a = row_ref.upper(), row_alt.upper()
    if (r, a) == (e_ref, e_eff):
        return row_beta
    if (r, a) == (e_eff, e_ref):
        return -row_beta
    cr = COMP.get(r)
    if cr is None:
        return None
    if (cr, COMP.get(a)) == (e_ref, e_eff):
        return row_beta
    if (cr, COMP.get(a)) == (e_eff, e_ref):
        return -row_beta
    return None


def wald(b_exp, b_out, se_out):
    b = b_out / b_exp
    se = se_out / abs(b_exp)
    z = b / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return b, se, p


def ivw(pairs):
    w = [1.0 / (s ** 2) for _, s in pairs]
    b = sum(wi * bi for wi, (bi, _) in zip(w, pairs)) / sum(w)
    se = math.sqrt(1.0 / sum(w))
    z = b / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return b, se, p


def bh_fdr(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * n / rank)
        q[i] = val
        prev = val
    return q


def main():
    rows = load_mbg()
    print("MiBioGen rows=%d  taxa=%d" % (len(rows), len({r['taxon'] for r in rows})))
    ivs = select_ivs(rows)
    print("selected instruments=%d  unique rsID=%d" % (len(ivs), len({r['rsid'] for r in ivs})))

    cache = os.path.join(PROJ, "03_data_raw", "hpa", "..", "ensembl_map.json")
    cache = os.path.abspath(cache)
    if os.path.exists(cache):
        mapped = json.load(open(cache, encoding="utf-8"))
        print("loaded Ensembl cache", len(mapped))
    else:
        mapped = ensembl_map({r["rsid"] for r in ivs})
        json.dump(mapped, open(cache, "w", encoding="utf-8"))
        print("cached", cache)

    # ---- outcome lookup ----
    tb = Tabix(FINNGEN % PHENO,
               tbi_path=P_TBI_LOCAL if os.path.exists(P_TBI_LOCAL)
               else os.path.join(TBI_DIR, "finngen_R12_%s.gz.tbi" % PHENO))
    print("tabix ready, refs=%d" % len(tb.names), flush=True)

    positions = {}
    for r in ivs:
        v = mapped.get(r["rsid"])
        if not v:
            continue
        m = pick_mapping(v)
        if not m:
            continue
        positions[r["rsid"]] = m
    uniq = sorted({(m["chrom"], m["pos"]) for m in positions.values()})
    print("querying %d unique positions" % len(uniq), flush=True)

    got = {}
    t0 = time.time()
    for i, (c, p) in enumerate(uniq):
        try:
            lines = tb.query(c, p, p)
        except Exception as e:
            print("  query fail %s:%d %s" % (c, p, type(e).__name__))
            continue
        for ln in lines:
            f = ln.split("\t")
            if len(f) < 11:
                continue
            try:
                if int(f[1]) != p:
                    continue
                got[(c, p)] = dict(ref=f[2], alt=f[3], rsid=f[4],
                                   p=float(f[6]), beta=float(f[8]),
                                   se=float(f[9]))
            except (ValueError, IndexError):
                pass
        if (i + 1) % 100 == 0:
            print("  %d/%d  hit=%d  %.0fs" % (i + 1, len(uniq), len(got),
                                              time.time() - t0), flush=True)
    print("outcome lookup done: %d/%d positions with data, %.0fs" % (
        len(got), len(uniq), time.time() - t0))

    # ---- per-taxon MR ----
    by_taxon = defaultdict(list)
    for r in ivs:
        m = positions.get(r["rsid"])
        if not m:
            continue
        o = got.get((m["chrom"], m["pos"]))
        if not o:
            continue
        b_out = allele_beta(o["ref"], o["alt"], o["beta"], r["eff"], r["ref"])
        if b_out is None:
            continue
        f_stat = (r["beta"] / r["se"]) ** 2
        if f_stat < 10:
            continue
        b, se, p = wald(r["beta"], b_out, o["se"])
        by_taxon[r["taxon"]].append(dict(rsid=r["rsid"], b=b, se=se, p=p,
                                         f=f_stat, p_exp=r["p"], p_out=o["p"]))

    recs = []
    for taxon, ps in by_taxon.items():
        if not ps:
            continue
        if len(ps) == 1:
            b, se, p = ps[0]["b"], ps[0]["se"], ps[0]["p"]
        else:
            b, se, p = ivw([(x["b"], x["se"]) for x in ps])
        recs.append(dict(taxon=taxon, n_iv=len(ps), beta=round(b, 5),
                         se=round(se, 5), OR=round(math.exp(b), 4),
                         ci_lo=round(math.exp(b - 1.96 * se), 4),
                         ci_hi=round(math.exp(b + 1.96 * se), 4),
                         pval=p, mean_F=round(sum(x["f"] for x in ps) / len(ps), 1),
                         min_p_exp="%.1e" % min(x["p_exp"] for x in ps),
                         min_p_out="%.1e" % min(x["p_out"] for x in ps),
                         ivs=";".join(x["rsid"] for x in ps)))

    qs = bh_fdr([r["pval"] for r in recs])
    for r, q in zip(recs, qs):
        r["fdr"] = q
        r["pval"] = "%.3e" % r["pval"]
        r["fdr"] = "%.3e" % q
    recs.sort(key=lambda r: float(r["pval"]))

    out = os.path.join(RES, "mr_microbiome_htn.csv")
    keys = ["taxon", "n_iv", "beta", "se", "OR", "ci_lo", "ci_hi", "pval",
            "fdr", "mean_F", "min_p_exp", "min_p_out", "ivs"]
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(recs)
    print("written", out, "taxa with MR =", len(recs))
    print()
    print("%-45s %3s %8s %18s %10s %10s" % ("taxon", "k", "OR", "95%CI", "P", "FDR"))
    for r in recs[:15]:
        print("%-45s %3d %8.3f %8.3f-%-9.3f %10s %10s" % (
            r["taxon"][:45], r["n_iv"], r["OR"], r["ci_lo"], r["ci_hi"],
            r["pval"], r["fdr"]))


if __name__ == "__main__":
    main()
