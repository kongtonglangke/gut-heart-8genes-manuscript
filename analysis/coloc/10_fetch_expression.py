#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
M1: GTEx v8 54-tissue median TPM for the 8 candidate genes.

Pipeline:
  1. Ensembl REST  -> Ensembl gene id (unversioned)
  2. GTEx reference/gene -> GTEx gencodeId (versioned, v26/GRCh38)
  3. GTEx medianGeneExpression -> 54 tissues median TPM

Output:
  03_data_raw/expression_median_tpm.csv   (long format)
  03_data_raw/tissue_sample_n.csv         (tissue -> n samples)
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

GENES = ["LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1"]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "03_data_raw")
os.makedirs(OUTDIR, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0"}


def get_json(url, timeout=60, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            last = e
            time.sleep(2 + 2 * i)
    raise last


# Verified 2026-11-01 via Ensembl REST (GRCh38). Hard-coded because the
# batch lookup endpoint returns sporadic 500/503 and is not needed at runtime.
ENSEMBL = {
    "LRP1": "ENSG00000123384",
    "CTSD": "ENSG00000117984",
    "MTHFR": "ENSG00000177000",
    "AUTS2": "ENSG00000158321",
    "FURIN": "ENSG00000140564",
    "CST3": "ENSG00000101439",
    "FCER1G": "ENSG00000158869",
    "TBXAS1": "ENSG00000059377",
}


def ensembl_ids():
    return dict(ENSEMBL)


def gtex_gencode(eid):
    url = ("https://gtexportal.org/api/v2/reference/gene?geneId=%s"
           "&genomeBuild=GRCh38/hg38" % urllib.parse.quote(eid))
    d = get_json(url)
    rows = d.get("data", [])
    return rows[0]["gencodeId"] if rows else None


def median_tpm(gid):
    url = ("https://gtexportal.org/api/v2/expression/medianGeneExpression"
           "?gencodeId=%s&datasetId=gtex_v8" % urllib.parse.quote(gid))
    d = get_json(url, timeout=90)
    return d.get("data", [])


def main():
    eids = ensembl_ids()
    print("Ensembl ids:", eids)

    recs = []
    for g in GENES:
        eid = eids.get(g)
        if not eid:
            print("%-7s no Ensembl id, skip" % g)
            continue
        gid = gtex_gencode(eid)
        if not gid:
            print("%-7s no GTEx gencodeId, skip" % g)
            continue
        rows = median_tpm(gid)
        print("%-7s %-20s tissues=%d" % (g, gid, len(rows)))
        for r in rows:
            recs.append(dict(gene=g, ensembl=eid, gencode=gid,
                             tissue=r["tissueSiteDetailId"],
                             median_tpm=round(float(r["median"]), 4),
                             unit=r.get("unit", "TPM")))
        time.sleep(0.4)

    p1 = os.path.join(OUTDIR, "expression_median_tpm.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "ensembl", "gencode", "tissue",
                                          "median_tpm", "unit"])
        w.writeheader()
        w.writerows(recs)
    print("written %s  rows=%d" % (p1, len(recs)))

    # tissue sample sizes (for power statements)
    t = get_json("https://gtexportal.org/api/v2/dataset/tissueSiteDetail?datasetId=gtex_v8")
    p2 = os.path.join(OUTDIR, "tissue_sample_n.csv")
    with open(p2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["tissue", "n_samples", "colorHex"])
        w.writeheader()
        for r in t.get("data", []):
            for col in ("rnaSeqSampleCount", "totalSampleCount"):
                pass
            w.writerow(dict(tissue=r["tissueSiteDetailId"],
                            n_samples=r.get("rnaSeqSampleCount") or r.get("totalSampleCount"),
                            colorHex=r.get("colorHex", "")))
    print("written %s" % p2)


if __name__ == "__main__":
    sys.exit(main())
