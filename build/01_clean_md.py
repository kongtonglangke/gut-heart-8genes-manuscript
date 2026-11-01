# -*- coding: utf-8 -*-
"""Step 1 - strip every Chinese-language auxiliary note from the manuscript.

The source markdown is machine-assembled and carries Chinese working notes
(working title, submission to-dos, audit trails, citation to-do lists).  The
user asked for all of them to be removed so that the docx contains only the
submission-grade English text.  This script is read-only with respect to the
source file: it writes a cleaned copy next to it.

Usage:  python 01_clean_md.py
"""
import io
import os
import re
import sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "Manuscript_GutHeart_8genes_v1.md")
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "Manuscript_GutHeart_8genes_v1.clean.md")

CJK = re.compile(r"[\u4e00-\u9fff\u3000-\u303f\uff00-\uffef]")

# inline Chinese fragments that must disappear while the English around them stays
INLINE_KILL = ["（通讯作者）", "（白静）", "（可选，部分期刊需要）"]


def clean(lines):
    out = []
    i = 0
    n = len(lines)
    dropped = []
    while i < n:
        raw = lines[i]
        s = raw.strip()

        # 1) generated-file comment
        if s.startswith("<!--"):
            dropped.append((i + 1, "comment"))
            i += 1
            continue

        # 2) Chinese working title
        if s.startswith("Chinese working title"):
            dropped.append((i + 1, "working title"))
            i += 1
            continue

        # 3) the "citations still to add" block (heading + its bullet list)
        if s.startswith("**待补充引用"):
            while i < n and lines[i].strip():
                dropped.append((i + 1, "citation to-do"))
                i += 1
            continue

        # 4) submission to-do notes (single-line Chinese italics block)
        if s.startswith("*") and CJK.search(s) and s.endswith("*"):
            dropped.append((i + 1, "submission note"))
            i += 1
            continue

        # 5) free-standing Chinese parenthetical notes
        if s.startswith("（") and CJK.search(s):
            dropped.append((i + 1, "parenthetical note"))
            i += 1
            continue

        # 6) the citation-number audit table section
        if s.startswith("## 引编号使用自查"):
            while i < n and not lines[i].startswith("# "):
                dropped.append((i + 1, "citation audit"))
                i += 1
            continue

        for frag in INLINE_KILL:
            raw = raw.replace(frag, "")
        out.append(raw.rstrip())
        i += 1
    return out, dropped


def main():
    src = io.open(SRC, encoding="utf-8").read()
    lines = src.split("\n")
    before = sum(1 for l in lines if CJK.search(l))

    kept, dropped = clean(lines)

    # collapse runs of >1 blank line created by the deletions
    final = []
    for l in kept:
        if not l.strip() and final and not final[-1].strip():
            continue
        final.append(l)

    text = "\n".join(final).rstrip() + "\n"
    remaining = [l for l in text.split("\n") if CJK.search(l)]
    if remaining:
        sys.stderr.write("CJK still present after cleaning:\n")
        for l in remaining[:20]:
            sys.stderr.write("  " + l[:120] + "\n")
        raise SystemExit("cleaning incomplete")

    with io.open(DST, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)

    print("source lines : %d" % len(lines))
    print("CJK lines    : %d -> 0" % before)
    print("dropped      : %d" % len(dropped))
    for ln, why in dropped:
        print("   -%5d  %s" % (ln, why))
    print("output       : %s" % DST)
    print("words        : %d" % len(text.split()))


if __name__ == "__main__":
    main()
