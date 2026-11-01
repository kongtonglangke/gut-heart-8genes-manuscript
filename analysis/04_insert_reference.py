# -*- coding: utf-8 -*-
"""Insert Ferkingstad 2021 as reference [8] and shift the old 8-11 to 9-12.

Renumbering is a single regex pass with a mapping function, so a shifted number can
never be shifted twice.  The newly added [8] is held behind a sentinel during the
pass.  Afterwards the citation sequence is verified: no gaps, no orphans.

Usage:  python 04_insert_reference.py
"""
import io
import re

MD = r"F:\longxiawenjian\gaoxueyaqiguanzhou\Manuscript_GutHeart_8genes\Manuscript_GutHeart_8genes_v1.md"
SENT = "\x01DECODE\x01"

NEW_REF = ("8. Ferkingstad E, Sulem P, Atlason BA, et al. Large-scale integration of "
           "the plasma proteome with genetics and disease. *Nat Genet*. "
           "2021;53(12):1712\u20131721. doi:10.1038/s41588-021-00978-w.")

text = io.open(MD, encoding="utf-8").read()
head, sep, refs = text.partition("# References")
assert sep, "References heading not found"

if "Ferkingstad" in refs:
    raise SystemExit("already applied; nothing to do")

# 1. protect the newly added deCODE citation
anchor = "study [8], in which 4,907"
assert head.count(anchor) == 1, "deCODE citation anchor not unique (%d)" % head.count(anchor)
head = head.replace(anchor, "study " + SENT + ", in which 4,907")

# 2. shift every citation in [8, 11] by +1 in one pass
def shift(m):
    nums = [int(x) for x in m.group(1).split(",")]
    return "[" + ",".join(str(n + 1 if 8 <= n <= 11 else n) for n in nums) + "]"

head = re.sub(r"\[(\d+(?:\s*,\s*\d+)*)\]", shift, head)

# 3. restore the sentinel
head = head.replace(SENT, "[8]")

# 4. renumber the reference list (descending, so entries never collide), then insert
lines = refs.split("\n")
for old in (11, 10, 9, 8):
    for i, ln in enumerate(lines):
        if re.match(r"^%d\. " % old, ln):
            lines[i] = "%d. " % (old + 1) + ln.split(". ", 1)[1]
idx = next(i for i, ln in enumerate(lines) if re.match(r"^9\. ", ln))
lines.insert(idx, NEW_REF)
refs = "\n".join(lines)

io.open(MD, "w", encoding="utf-8", newline="\n").write(head + sep + refs)

# ---------------- verification ----------------
txt = io.open(MD, encoding="utf-8").read()
body, _, rl = txt.partition("# References")
cited = set()
for m in re.finditer(r"\[(\d+(?:\s*,\s*\d+)*)\]", body):
    for k in m.group(1).split(","):
        cited.add(int(k.strip()))
listed = [int(m.group(1)) for m in re.finditer(r"^(\d+)\. ", rl, re.M)]
print("正文引用编号:", sorted(cited))
print("文献条目编号:", sorted(listed))
print("引用数 %d / 条目数 %d" % (len(cited), len(listed)))
assert sorted(cited) == list(range(1, len(listed) + 1)), "引用编号有缺口"
assert listed == list(range(1, len(listed) + 1)), "文献条目编号不连续"
order = []
for m in re.finditer(r"\[(\d+)\]", body):
    order.append(int(m.group(1)))
first = []
for n in order:
    if n not in first:
        first.append(n)
print("首次出现顺序:", first)
print("OK" if first == sorted(first) else "WARN 首次引用顺序非递增")
