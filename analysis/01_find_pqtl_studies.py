# -*- coding: utf-8 -*-
"""Step 1 - find plasma pQTL studies on the GWAS Catalog for the four candidates.

The GWAS Catalog keeps one study per protein per cohort, so the trait name has to
be probed in several spellings.  For every hit we record the accession, sample
size, ancestry and whether full summary statistics (and therefore a possible
tabix index) are deposited.

Usage:  python 01_find_pqtl_studies.py
"""
import json
import os
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = os.path.join(HERE, "_tmp")
os.makedirs(TMP, exist_ok=True)
API = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"

# protein -> candidate trait spellings on the GWAS Catalog
QUERIES = {
    "FURIN": ["Furin levels", "FURIN levels", "furin measurement"],
    "FES": ["Fes levels", "FES levels", "Tyrosine-protein kinase Fes/Fps levels",
            "Fes/Fps levels"],
    "CST3": ["Cystatin-C levels", "Cystatin C levels", "Cystatin-3 levels",
             "cystatin C measurement"],
    "MTHFR": ["MTHFR levels", "Methylenetetrahydrofolate reductase levels",
              "MTHFR measurement"],
}
FIELDS = ("accession_id,full_summary_stats_available,full_summary_stats,"
          "initial_sample_size,discovery_ancestry,disease_trait,pubmed_id")


def get(url, tries=3, timeout=45):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            raw = urllib.request.urlopen(req, timeout=timeout).read().decode()
            return json.loads(raw)
        except Exception as e:
            if k == tries - 1:
                print("    ! %s %s" % (type(e).__name__, str(e)[:60]))
                return None
            time.sleep(2)


def main():
    report = {}
    for gene, names in QUERIES.items():
        print("=" * 70)
        print(gene)
        seen = {}
        for name in names:
            url = "%s?diseaseTrait=%s&size=50" % (API, urllib.parse.quote(name))
            d = get(url)
            tot = (d or {}).get("page", {}).get("totalElements", 0)
            print("  %-46s hits=%s" % (name, tot))
            for s in ((d or {}).get("_embedded", {}) or {}).get("studies", []):
                acc = s.get("accession_id")
                if not acc or acc in seen:
                    continue
                seen[acc] = {
                    "accession": acc,
                    "trait": s.get("disease_trait"),
                    "n": s.get("initial_sample_size"),
                    "ancestry": ", ".join(s.get("discovery_ancestry") or []),
                    "pubmed": s.get("pubmed_id"),
                    "has_sumstats": s.get("full_summary_stats_available"),
                    "ftp_dir": s.get("full_summary_stats"),
                }
            time.sleep(0.5)
        rows = sorted(seen.values(), key=lambda r: (not r["has_sumstats"],
                                                    r["accession"]))
        report[gene] = rows
        for r in rows:
            print("    %-14s %-32s %s" % (r["accession"], (r["n"] or "")[:32],
                                          r["ancestry"][:40]))
    out = os.path.join(HERE, "pqtl_studies.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    print("\nwritten:", out)
    print("\n=== 覆盖汇总 ===")
    for g, rows in report.items():
        ok = [r for r in rows if r["has_sumstats"]]
        print("  %-6s studies=%d  with-sumstats=%d" % (g, len(rows), len(ok)))


if __name__ == "__main__":
    main()
