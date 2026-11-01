# -*- coding: utf-8 -*-
"""Look up APC for candidate journals on DOAJ (authoritative for fully-OA titles).

Hybrid/subscription journals are not indexed by DOAJ; for those the free route is
the non-OA (subscription) publishing option, which is checked separately.

Usage:  python 09_journals.py
"""
import json
import urllib.parse
import urllib.request

TITLES = [
    "Scientific Reports", "PLOS ONE", "BMC Medical Genomics", "Frontiers in Genetics",
    "Journal of Human Genetics", "Human Genetics", "Molecular Genetics & Genomic Medicine",
    "BMC Cardiovascular Disorders", "Egyptian Journal of Medical Human Genetics",
    "Molecular Vision", "Global Medical Genetics", "Journal of Human Hypertension",
    "Hypertension Research", "BMC Medicine", "Genome Medicine", "Human Genomics",
    "Journal of Negative Results in BioMedicine", "BBA Advances",
    "Clinical and Translational Discovery", "Journal of Circulating Biomarkers",
]


def q(title):
    url = ("https://doaj.org/api/search/journals/"
           + urllib.parse.quote('bibjson.title:"%s"' % title)
           + "?pageSize=3")
    try:
        with urllib.request.urlopen(url, timeout=30) as f:
            return json.loads(f.read().decode())
    except Exception as e:
        return {"_err": str(e)[:40]}


rows = []
for t in TITLES:
    d = q(t)
    res = (d.get("results") or []) if isinstance(d, dict) else []
    if not res:
        rows.append((t, "-", "not indexed by DOAJ (hybrid / subscription)"))
        continue
    b = res[0]["bibjson"]
    name = b.get("title", t)
    apc = b.get("apc", {}) or {}
    if not apc.get("has_apc"):
        rows.append((name, "0", "no APC"))
    else:
        mx = (apc.get("max") or [{}])[0]
        rows.append((name, "%s %s" % (mx.get("price"), mx.get("currency")),
                     "APC charged"))

w = max(len(r[0]) for r in rows) + 2
print("%-*s %-14s %s" % (w, "journal", "APC", "note"))
print("-" * (w + 40))
for r in rows:
    print("%-*s %-14s %s" % (w, r[0][:w - 1], r[1], r[2]))
