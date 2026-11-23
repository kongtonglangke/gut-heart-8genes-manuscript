#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M7 (fast path): gut microbiota -> hypertension, single-instrument cis-MR.

Differences from 14_microbiome_mr.py
  * ONE instrument per taxon (lowest-P variant, F>10).  With 1,055 instruments
    the tabix layer had to fetch 417 distinct 64-128 kb blocks; with 211 leads
    it is roughly half that, and single-instrument Wald ratios are what we can
    defend anyway without an LD panel.
  * hard per-query timeout: a variant living in an ultra-dense region (the LCT
    peak has thousands of FinnGen records) can otherwise stall for many minutes.
  * progress printed every 20 variants so the run can be monitored.

Output: 05_results/tables/mr_microbiome_htn.csv
"""
import csv
import json
import math
import os
import socket
import sys
import time
import urllib.request
from collections import defaultdict

socket.setdefaulttimeout(60)

ROOT = "F:/longxiawenjian/gaoxueyaqiguanzhou"
PROJ = os.path.join(ROOT, "Gut_Heart_RAS_8genes")
MBG = os.path.join(ROOT, "P0_mr", "data", "mbg_allhits.txt")
RES = os.path.join(PROJ, "05_results", "tables")
CACHE = os.path.join(PROJ, "03_data_raw", "ensembl_map.json")
RESULT_JSON = os.path.join(PROJ, "03_data_raw", "mbg_outcome_cache.json")

sys.path.insert(0, os.path.join(ROOT, "P0_mr"))
from _tabix import Tabix  # noqa: E402

FINNGEN = ("https://storage.googleapis.com/finngen-public-data-r12/"
           "summary_stats/release/finngen_R12_I9_HYPTENS.gz")
TBI = os.path.join(ROOT, "P0_mr", "data", "finngen_htn.tbi")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
QUERY_TIMEOUT = 90


def load_leads():
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
                    p=float(r["P.weightedSumZ"])))
            except (ValueError, KeyError):
                pass
    by = defaultdict(list)
    for r in rows:
        by[r["taxon"]].append(r)
    leads = []
    for taxon, xs in by.items():
        xs = [x for x in xs if (x["beta"] / x["se"]) ** 2 > 10]
        if not xs:
            continue
        leads.append(min(xs, key=lambda x: x["p"]))
    return leads


def g38(rs):
    v = MAPPED.get(rs)
    if not v:
        return None
    maps = [m for m in v.get("mappings", []) if m.get("assembly_name") == "GRCh38"]
    if not maps:
        return None
    m = maps[0]
    al = m.get("allele_string", "")
    if "/" not in al:
        return None
    parts = al.split("/")
    if len(parts) != 2 or len(parts[0]) != 1 or len(parts[1]) != 1:
        return None
    return dict(chrom=m["seq_region_name"], pos=int(m["start"]),
                a1=parts[0].upper(), a2=parts[1].upper())


def allele_beta(ref, alt, beta, e_ref, e_eff):
    r, a = ref.upper(), alt.upper()
    if (r, a) == (e_ref, e_eff):
        return beta
    if (r, a) == (e_eff, e_ref):
        return -beta
    cr = COMP.get(r)
    if cr is None:
        return None
    if (cr, COMP.get(a)) == (e_ref, e_eff):
        return beta
    if (cr, COMP.get(a)) == (e_eff, e_ref):
        return -beta
    return None


def main():
    global MAPPED
    MAPPED = json.load(open(CACHE, encoding="utf-8"))
    leads = load_leads()
    print("taxa with F>10 lead: %d" % len(leads), flush=True)

    pos = {}
    for r in leads:
        m = g38(r["rsid"])
        if m:
            pos[r["rsid"]] = m
    uniq = sorted({(m["chrom"], m["pos"]) for m in pos.values()})
    print("unique positions: %d" % len(uniq), flush=True)

    cache = {}
    if os.path.exists(RESULT_JSON):
        cache = json.load(open(RESULT_JSON, encoding="utf-8"))
        print("outcome cache: %d" % len(cache), flush=True)

    tb = Tabix(FINNGEN, tbi_path=TBI)
    t0 = time.time()
    slow = 0
    for i, (c, p) in enumerate(uniq):
        key = "%s:%d" % (c, p)
        if key in cache:
            continue
        ts = time.time()
        try:
            lines = tb.query(c, p, p)
            hit = None
            for ln in lines:
                f = ln.split("\t")
                if len(f) < 11:
                    continue
                try:
                    if int(f[1]) != p:
                        continue
                    hit = dict(ref=f[2], alt=f[3], rsid=f[4], p=float(f[6]),
                               beta=float(f[8]), se=float(f[9]))
                    break
                except (ValueError, IndexError):
                    pass
            cache[key] = hit
        except Exception:
            cache[key] = None
        dt = time.time() - ts
        if dt > 20:
            slow += 1
        if (i + 1) % 20 == 0:
            print("  %d/%d  %.0fs  slow=%d" % (i + 1, len(uniq), time.time() - t0, slow),
                  flush=True)
            json.dump(cache, open(RESULT_JSON, "w", encoding="utf-8"))
    json.dump(cache, open(RESULT_JSON, "w", encoding="utf-8"))
    print("lookup done: %d/%d with data" % (
        sum(1 for v in cache.values() if v), len(uniq)), flush=True)

    recs = []
    for r in leads:
        m = pos.get(r["rsid"])
        if not m:
            continue
        o = cache.get("%s:%d" % (m["chrom"], m["pos"]))
        if not o:
            continue
        b_out = allele_beta(o["ref"], o["alt"], o["beta"], r["eff"], r["ref"])
        if b_out is None:
            continue
        b = b_out / r["beta"]
        se = o["se"] / abs(r["beta"])
        z = b / se
        pv = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
        recs.append(dict(taxon=r["taxon"], n_iv=1, rsid=r["rsid"], beta=round(b, 5),
                         se=round(se, 5), OR=round(math.exp(b), 4),
                         ci_lo=round(math.exp(b - 1.96 * se), 4),
                         ci_hi=round(math.exp(b + 1.96 * se), 4),
                         pval=pv, F=round((r["beta"] / r["se"]) ** 2, 1),
                         p_exp="%.1e" % r["p"], p_out="%.1e" % o["p"]))

    ps = [r["pval"] for r in recs]
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    q = [0.0] * len(ps)
    prev = 1.0
    for rank in range(len(ps), 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * len(ps) / rank)
        q[i] = prev
    for r, qq in zip(recs, q):
        r["fdr"] = qq
    recs.sort(key=lambda r: r["pval"])

    out = os.path.join(RES, "mr_microbiome_htn.csv")
    keys = ["taxon", "n_iv", "rsid", "beta", "se", "OR", "ci_lo", "ci_hi",
            "pval", "fdr", "F", "p_exp", "p_out"]
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            rr = dict(r)
            rr["pval"] = "%.3e" % r["pval"]
            rr["fdr"] = "%.3e" % r["fdr"]
            w.writerow(rr)
    print("written", out, "taxa =", len(recs), flush=True)
    print("%-45s %8s %18s %10s %10s" % ("taxon", "OR", "95%CI", "P", "FDR"))
    for r in recs[:12]:
        print("%-45s %8.3f %8.3f-%-9.3f %10.2e %10.2e" % (
            r["taxon"][:45], r["OR"], r["ci_lo"], r["ci_hi"], r["pval"], r["fdr"]))


if __name__ == "__main__":
    main()
