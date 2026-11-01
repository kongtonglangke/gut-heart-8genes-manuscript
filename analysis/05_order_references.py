# -*- coding: utf-8 -*-
"""Renumber the citations so that they run 1..n in order of first appearance.

The source manuscript numbered its references by theme rather than by citation
order, so the Introduction cited entry 11 before the Methods cited entry 2.  A
clean permutation fixes this without touching any wording.

Usage:  python 05_order_references.py
"""
import io
import re

MD = r"F:\longxiawenjian\gaoxueyaqiguanzhou\Manuscript_GutHeart_8genes\Manuscript_GutHeart_8genes_v1.md"

text = io.open(MD, encoding="utf-8").read()
body, sep, refs = text.partition("# References")
assert sep

# ---- 1. first-appearance order
order = []
for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", body):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in order:
            order.append(n)
mapping = {old: i + 1 for i, old in enumerate(order)}
print("原编号首次出现顺序:", order)
if order == list(range(1, len(order) + 1)):
    raise SystemExit("already in citation order; nothing to do")

# ---- 2. rewrite the body citations
def remap(m):
    nums = sorted(mapping[int(x)] for x in m.group(1).split(","))
    return "[" + ",".join(str(n) for n in nums) + "]"

body = re.sub(r"\[(\d+(?:\s*,\s*\d+)*)\]", remap, body)

# ---- 3. reorder the reference list
entries = {}
for m in re.finditer(r"^(\d+)\. (.+?)(?=\n\d+\. |\n\n|\Z)", refs, re.M | re.S):
    entries[int(m.group(1))] = " ".join(m.group(2).split())
assert len(entries) == len(order), "条目数 %d 与引用数 %d 不符" % (len(entries), len(order))

new_lines = []
for i, old in enumerate(order, 1):
    new_lines.append("%d. %s" % (i, entries[old]))

new_refs = "\n\n" + "\n".join(new_lines) + "\n"
io.open(MD, "w", encoding="utf-8", newline="\n").write(body + sep + new_refs)

# ---- 4. verification
txt = io.open(MD, encoding="utf-8").read()
b2, _, r2 = txt.partition("# References")
o2 = []
for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", b2):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in o2:
            o2.append(n)
l2 = [int(m.group(1)) for m in re.finditer(r"^(\d+)\. ", r2, re.M)]
print("新首次出现顺序:", o2)
print("文献条目:", l2)
assert o2 == list(range(1, len(l2) + 1)), "重排后仍非递增"
assert l2 == list(range(1, len(l2) + 1))
print("OK: 引用按首次出现顺序 1..%d，条目一一对应" % len(l2))
for i, old in enumerate(order, 1):
    print("  [%2d] <- 原 [%2d]  %s" % (i, old, entries[old][:62]))
