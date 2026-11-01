# -*- coding: utf-8 -*-
"""Append the two independent FURIN references and renumber by first appearance.

Also drops the single stray PMID so that all entries carry the same metadata set
(DOI), which was one of the reviewer's consistency points.

Usage:  python 08_add_functional_refs.py
"""
import io
import re

MD = (r"F:\longxiawenjian\gaoxueyaqiguanzhou\Manuscript_GutHeart_8genes"
      r"\Manuscript_GutHeart_8genes_v1.md")

NEW = {
    13: ("Yang X, Yang W, McVey DG, et al. FURIN expression in vascular endothelial "
         "cells is modulated by a coronary artery disease-associated genetic variant and "
         "influences monocyte transendothelial migration. *J Am Heart Assoc*. "
         "2020;9(4):e014333. doi:10.1161/JAHA.119.014333."),
    14: ("Chen R, Shen D, Mahieux A, et al. FURIN induces vascular endothelial cell "
         "death and efferocytosis through enzymatic cleavage of thrombospondin-1. "
         "*J Physiol*. 2026;604(14):5874\u20135891. doi:10.1113/JP290390."),
}

t = io.open(MD, encoding="utf-8").read()
body, sep, refs = t.partition("# References")
assert sep

# ---- collect existing entries
entries = {}
for m in re.finditer(r"^(\d+)\. (.+?)(?=\n\d+\. |\n\n## |\Z)", refs, re.M | re.S):
    entries[int(m.group(1))] = " ".join(m.group(2).split())
assert len(entries) == 12, "expected 12 existing entries, found %d" % len(entries)
entries.update(NEW)

# ---- drop the lone PMID so every entry carries the same metadata
entries[6] = entries[6].replace(" PMID: 27019110.", "")

# ---- renumber by first appearance
order = []
for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", body):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in order:
            order.append(n)
assert sorted(order) == sorted(entries), "citation set != reference set: %s" % (
    set(order) ^ set(entries))
mapping = {old: i + 1 for i, old in enumerate(order)}


def remap(m):
    return "[" + ",".join(str(n) for n in
                          sorted(mapping[int(x)] for x in m.group(1).split(","))) + "]"


body = re.sub(r"\[(\d+(?:\s*,\s*\d+)*)\]", remap, body)

new_lines = ["%d. %s" % (i, entries[old]) for i, old in enumerate(order, 1)]
refs_new = "\n\n" + "\n".join(new_lines) + "\n"
io.open(MD, "w", encoding="utf-8", newline="\n").write(body + sep + refs_new)

# ---- verify
txt = io.open(MD, encoding="utf-8").read()
b2, _, r2 = txt.partition("# References")
o2 = []
for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", b2):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in o2:
            o2.append(n)
l2 = [int(m.group(1)) for m in re.finditer(r"^(\d+)\. ", r2, re.M)]
print("首次出现顺序:", o2)
print("文献条目:", l2)
assert o2 == list(range(1, len(l2) + 1))
assert l2 == list(range(1, len(l2) + 1))
print("OK: %d 条文献，编号连续且按首次引用排序" % len(l2))
print("\n=== 受影响条目 ===")
for i, old in enumerate(order, 1):
    if old in NEW or old == 6:
        print("  [%2d] <- 原 [%2d]  %s" % (i, old, entries[old][:78]))
