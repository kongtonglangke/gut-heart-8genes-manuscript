# -*- coding: utf-8 -*-
"""Step 3 - verify the generated docx against the skill's 13-item checklist.

Usage:  python 03_verify.py
"""
import io
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, r"C:\Users\user\.workbuddy\skills\docx-academic-format\scripts")

from docx import Document                                    # noqa: E402
from docx.shared import Pt                                   # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH                # noqa: E402
from docx.oxml.ns import qn                                  # noqa: E402
from docx.table import Table                                 # noqa: E402
from docx.text.paragraph import Paragraph                    # noqa: E402

from tri_line_table import verify_table_borders              # noqa: E402
import importlib
b = importlib.import_module("02_md2docx")

SRC = os.path.join(SCRIPT_DIR, "Manuscript_GutHeart_8genes_v1.clean.md")
DST = (sys.argv[1] if len(sys.argv) > 1
       else os.path.join(SCRIPT_DIR, "Manuscript_GutHeart_8genes_v1.docx"))

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
STATS = {"P", "F", "n", "\u03b2"}
UPRIGHT_ABBR = ("OR", "CI", "RR", "HR", "SD", "SE", "AUC", "FDR", "TPM")
CJK = re.compile(r"[\u4e00-\u9fff\u3000-\u303f\uff00-\uffef]")

fails, warns, oks = [], [], []


def check(cond, msg):
    (oks if cond else fails).append(msg)
    return cond


# ---------------------------------------------------------------- helpers

def body_seq(doc):
    seq = []
    for c in doc.element.body.iterchildren():
        if c.tag == W + "p":
            seq.append(("p", Paragraph(c, doc)))
        elif c.tag == W + "tbl":
            seq.append(("t", Table(c, doc)))
    return seq


def canon(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


def cell_texts(tbl):
    out = []
    for row in tbl.rows:
        for c in row.cells:
            out.append(c.text)
    return out


def desup(s):
    """Mirror the renderer: unicode superscript glyphs come out as ASCII."""
    return re.sub("[" + "".join(b.SUP_MAP) + "]+",
                  lambda m: "".join(b.SUP_MAP[c] for c in m.group(0)), s)


def expected_stream():
    md = io.open(SRC, encoding="utf-8").read()
    out = []
    for kind, payload in b.parse(md):
        if kind == "tbl":
            out.extend(desup(cell) for row in payload for cell in row)
        else:
            out.append(desup(payload))
    return out


# ------------------------------------------------------------------ 0. load
doc = Document(DST)
seq = body_seq(doc)
paras = [o for k, o in seq if k == "p"]
tables = [o for k, o in seq if k == "t"]

# 1. no empty paragraphs (a paragraph holding only a figure is not a blank line)
def has_fig(p):
    return bool(p._element.findall(".//" + qn("w:drawing")))


fig_paras = [p for p in paras if has_fig(p)]
emp = [i for i, p in enumerate(paras)
       if not p.text.strip() and not has_fig(p)]
check(not emp, "1. empty paragraphs == 0 (found %s)" % (emp[:5] or 0))
check(len(fig_paras) == 10,
      "1b. embedded figures: %d / 10" % len(fig_paras))

# 2. no bold anywhere in running text
bold = [r.text for p in paras for r in p.runs if r.font.bold]
check(not bold, "2. bold runs in running text == 0 (found %s)" % bold[:5])

# 2b. table bold is allowed, but only for the cells the source actually marks
def md_table_bold():
    md_txt = io.open(SRC, encoding="utf-8").read()
    found = []
    for line in md_txt.split("\n"):
        if line.strip().startswith("|"):
            for m in re.finditer(r"\*\*(.+?)\*\*", line):
                found.append(canon(desup(b.normalise(m.group(1)))))
    return sorted(found)


doc_bold_cells = []
for t in tables:
    for row in t.rows:
        for c in row.cells:
            bl = "".join(r.text for r in c.paragraphs[0].runs if r.font.bold)
            if bl.strip():
                doc_bold_cells.append(canon(bl))
exp_bold = md_table_bold()
check(sorted(doc_bold_cells) == exp_bold,
      "2b. table bold: %d cells, matches the %d source-marked cells"
      % (len(doc_bold_cells), len(exp_bold)))
if sorted(doc_bold_cells) != exp_bold:
    warns.append("bold mismatch: docx=%s | md=%s"
                 % (sorted(doc_bold_cells), exp_bold))

# 2c. in Table 2 the bold marks must equal the eGene threshold exactly, cell by cell
PCELL = re.compile(r"([\d.]+)\s*\u00d7\s*10\u2212(\d+)")
if len(tables) > 1:
    t2 = tables[1]
    wrong = []
    n_sig = 0
    for row in t2.rows[1:]:
        gene = row.cells[0].text.strip()
        for c in row.cells[1:]:
            m = PCELL.search(c.text)
            if not m:
                continue
            pval = float(m.group(1)) * 10 ** (-int(m.group(2)))
            is_bold = any(r.font.bold for r in c.paragraphs[0].runs
                          if r.text.strip())
            sig = pval < 5e-8
            n_sig += sig
            if sig != is_bold:
                wrong.append("%s %s bold=%s P=%.1e" % (gene, c.text, is_bold, pval))
    check(not wrong,
          "2c. Table 2 bold == eGene status cell by cell (%d bold, %d with P<5e-8) %s"
          % (len(doc_bold_cells) - 5, n_sig, wrong[:3]))

# 3. level-1 headings
h1 = [p for p in paras if p.runs and p.runs[0].font.size == Pt(15)]
h1txt = {p.text for p in h1}
bad = []
for p in h1:
    rf = p.runs[0]._element.rPr.rFonts
    if rf.get(qn("w:eastAsia")) != "黑体":
        bad.append(("font", p.text[:30]))
    if (p.paragraph_format.space_before, p.paragraph_format.space_after) \
            != (Pt(18), Pt(18)):
        bad.append(("spacing", p.text[:30]))
check(len(h1) > 0 and not bad,
      "3. level-1 headings: n=%d, 15pt/黑体/18+18pt  %s" % (len(h1), bad[:3]))

# 4. level-2 headings
h2 = [p for p in paras if p.runs and p.runs[0].font.size == Pt(14)
      and p.text not in h1txt]
bad = [p.text[:30] for p in h2
       if p.paragraph_format.space_before != Pt(18)
       or p.paragraph_format.space_after != Pt(0)]
check(len(h2) > 0 and not bad,
      "4. level-2 headings: n=%d, 14pt/before 18pt  %s" % (len(h2), bad[:3]))

# 5. body: 12 pt, SimSun eastAsia, 1.5 spacing
body = [p for p in paras if p.runs and p.runs[0].font.size == Pt(12)]
bad = []
for p in body:
    rf = p.runs[0]._element.rPr.rFonts
    if rf.get(qn("w:eastAsia")) != "宋体":
        bad.append(("font", p.text[:30]))
    if p.paragraph_format.line_spacing != 1.5:
        bad.append(("line", p.text[:30]))
check(len(body) > 0 and not bad,
      "5. body paragraphs: n=%d, 12pt/宋体/1.5  %s" % (len(body), bad[:3]))

# 6. every run uses Times New Roman
bad = []
for p in paras:
    for r in p.runs:
        if r._element.findall(".//" + qn("w:drawing")):
            continue          # the run that carries a figure has no text rPr
        rf = r._element.rPr.rFonts if r._element.rPr is not None else None
        if rf is None or rf.get(qn("w:ascii")) != "Times New Roman" \
                or rf.get(qn("w:hAnsi")) != "Times New Roman":
            bad.append(r.text[:20])
for t in tables:
    for row in t.rows:
        for c in row.cells:
            for r in c.paragraphs[0].runs:
                rf = r._element.rPr.rFonts
                if rf.get(qn("w:ascii")) != "Times New Roman":
                    bad.append(r.text[:20])
check(not bad, "6. all runs ascii/hAnsi == Times New Roman (%s)" % bad[:3])

# 7. italics only for statistical symbols or source-marked spans
md = io.open(SRC, encoding="utf-8").read()
md_spans = set()
for m in re.finditer(r"\*+([^*\n]+?)\*+", md):
    md_spans.add(m.group(1).strip())
ital = []
for p in paras:
    ital += [r.text.strip() for r in p.runs if r.font.italic and r.text.strip()]
for t in tables:
    for row in t.rows:
        for c in row.cells:
            ital += [r.text.strip() for r in c.paragraphs[0].runs
                     if r.font.italic and r.text.strip()]
bad = [t for t in ital if t not in STATS and t not in md_spans]
check(len(ital) > 0 and not bad,
      "7. italic runs: n=%d, all are statistics or source-italics  %s"
      % (len(ital), bad[:5]))
bad = [t for t in ital if t.split(" ")[0] in UPRIGHT_ABBR]
check(not bad, "7b. upright abbreviations never italicised (%s)" % bad[:5])
check(sum(1 for t in ital if t in STATS) > 0,
      "7c. statistical symbols italicised: %d"
      % sum(1 for t in ital if t in STATS))

# 8. three-line tables
bad = []
for i, t in enumerate(tables, 1):
    good, prob = verify_table_borders(t)
    if not good:
        bad.append((i, prob[:2]))
check(len(tables) > 0 and not bad,
      "8. three-line tables: %d tables  %s" % (len(tables), bad))

# 9/10. table captions and cell fonts
cap = [p for p in paras if p.runs and p.runs[0].font.size == Pt(10.5)
       and (p.text.startswith("Table ") or p.text.startswith("Figure ")
            or p.text.startswith("Supplementary Figure"))]
tbl_caps = [p for p in cap if p.text.startswith("Table ")]
bad = [p.text[:40] for p in tbl_caps
       if p.alignment != WD_ALIGN_PARAGRAPH.CENTER]
check(len(tbl_caps) == 6 and not bad,
      "9. table captions: n=%d, 10.5pt centred  %s" % (len(tbl_caps), bad))
# captions must not end with punctuation
bad = [p.text[:40] for p in cap if p.text.rstrip().endswith((".", ";", ":"))]
check(not bad, "9b. captions carry no terminal punctuation (%s)" % bad)

bad = []
for t in tables:
    for row in t.rows:
        for c in row.cells:
            for r in c.paragraphs[0].runs:
                rf = r._element.rPr.rFonts
                if rf.get(qn("w:eastAsia")) != "宋体" or r.font.size != Pt(10.5):
                    bad.append(r.text[:15])
check(not bad, "10. table runs: 宋体 + 10.5pt (%s)" % bad[:3])

# 11. body justified
check(any(p.alignment == WD_ALIGN_PARAGRAPH.JUSTIFY for p in body),
      "11. body paragraphs justified (%d/%d)"
      % (sum(1 for p in body
             if p.alignment == WD_ALIGN_PARAGRAPH.JUSTIFY), len(body)))

# 12/13. spacing conventions (English text: no stray % or CJK)
alltext = "\n".join(p.text for p in paras) + "\n" + \
    "\n".join(x for t in tables for x in cell_texts(t))
check(not CJK.search(alltext),
      "12. zero CJK characters in the document")
# `95% CI` is the standard medical-journal form; only free-standing
# percentages would need the skill's "33 %" spacing rule
stray_pct = []
for m in re.finditer(r"%", alltext):
    ctx = alltext[max(0, m.start() - 14):m.start() + 12]
    if not re.search(r"\d{2}%\s*(CI|confidence)", ctx):
        stray_pct.append(ctx)
check(not stray_pct,
      "13. no free-standing percent signs needing a space (%s)" % stray_pct[:4])
check(not re.search(r"(?<!\*)\*(?!\*)", alltext),
      "13b. no unpaired asterisk survived")

# fidelity: no lost or duplicated words
exp = "".join(canon(x) for x in expected_stream())
act = "".join(canon(x) for k, o in seq
              for x in (cell_texts(o) if k == "t"
                        else ([] if has_fig(o) else [o.text])))
check(exp == act, "14. text fidelity: marker-stripped streams identical")
if exp != act:
    for i in range(min(len(exp), len(act))):
        if exp[i] != act[i]:
            warns.append("first divergence at %d: md ...%s | docx ...%s"
                         % (i, exp[max(0, i - 40):i + 40],
                            act[max(0, i - 40):i + 40]))
            break
    warns.append("length md=%d docx=%d" % (len(exp), len(act)))

# extras
sup = sum(1 for p in paras for r in p.runs
          if r._element.rPr is not None
          and r._element.rPr.find(qn("w:vertAlign")) is not None)
sup += sum(1 for t in tables for row in t.rows for c in row.cells
           for r in c.paragraphs[0].runs
           if r._element.rPr.find(qn("w:vertAlign")) is not None)
check(sup > 0, "15. genuine superscript runs: %d" % sup)
check("^" not in alltext, "16. no raw caret exponents left (%s)"
      % alltext.count("^"))
check("**" not in alltext and "|" not in alltext,
      "17. no markdown markers leaked into the document")
check("`" not in alltext, "17b. no backtick code spans survived")
check("+/-" not in alltext and "\u00b1" in alltext,
      "17c. '+/-' unified to '±' (%d occurrences)"
      % alltext.count("\u00b1"))
check("Wald ratios are single-instrument estimates." not in alltext,
      "17d. Table 5 note not duplicated by a trailing sentence")

# ------------------------------------------------------------------ report
print("=" * 68)
print("VERIFICATION  %s" % os.path.basename(DST))
print("=" * 68)
for m in oks:
    print("  PASS  %s" % m)
for m in warns:
    print("  WARN  %s" % m)
if fails:
    for m in fails:
        print("  FAIL  %s" % m)
print("-" * 68)
print("pass %d / fail %d / warn %d"
      % (len(oks), len(fails), len(warns)))
print("paragraphs %d | tables %d | italic runs %d"
      % (len(paras), len(tables), len(ital)))
sys.exit(1 if fails else 0)
