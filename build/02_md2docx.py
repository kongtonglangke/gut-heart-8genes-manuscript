# -*- coding: utf-8 -*-
"""Step 2 - render the cleaned manuscript markdown into a submission-grade docx.

Generated with python-docx 1.2.0, following the `docx-academic-format` skill:
  * three-line tables (1.5 / 1.0 / 1.5 lb, everything else suppressed)
  * Times New Roman for every Latin run, SimSun / SimHei in the eastAsia slot
  * 12 pt body, 1.5 line spacing, justified; 15 pt level-1 / 14 pt level-2
  * "blank line" implemented as space_before / space_after, never as an
    empty paragraph
  * no bold anywhere; italics only for statistical symbols and gene names

Deviations from the stock skill scripts are deliberate and are listed in
`03_verify.py`'s report:
  * `fmt_lib.add_text_with_italics` / `fix_italic.py` italicise `OR` and `CI`,
    which contradicts the skill's own authoritative table (OR / CI / RR / HR /
    SD / SE / AUC must stay upright).  Italics are therefore resolved inline
    here instead of by the post-processing passes.
  * `fix_spacing2.py` is a Chinese-typesetting pass; the cleaned text is
    English throughout, so spacing is inherited verbatim from the source.

Usage:  python 02_md2docx.py
"""
import io
import os
import re
import struct
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = r"C:\Users\user\.workbuddy\skills\docx-academic-format\scripts"
sys.path.insert(0, SKILL_DIR)

from docx import Document                                    # noqa: E402
from docx.shared import Pt, Cm                               # noqa: E402
from docx.enum.text import (WD_ALIGN_PARAGRAPH,              # noqa: E402
                            WD_LINE_SPACING, WD_BREAK)
from docx.enum.table import WD_TABLE_ALIGNMENT               # noqa: E402
from docx.enum.section import WD_SECTION                     # noqa: E402
from docx.oxml.ns import qn                                  # noqa: E402
from docx.oxml import OxmlElement                            # noqa: E402

from tri_line_table import set_table_borders                 # noqa: E402

SRC = os.path.join(SCRIPT_DIR, "Manuscript_GutHeart_8genes_v1.clean.md")
DST = os.path.join(SCRIPT_DIR, "Manuscript_GutHeart_8genes_v1.docx")

# --------------------------------------------------------------- type scale
SZ_H1 = Pt(15)        # 小三
SZ_H2 = Pt(14)        # 四号
SZ_BODY = Pt(12)      # 小四
SZ_SMALL = Pt(10.5)   # 五号
CN_BODY = "宋体"
CN_HEAD = "黑体"
EN = "Times New Roman"
GAP = Pt(18)          # 12 pt x 1.5 -> one blank line
ZERO = Pt(0)

CONTENT_CM = 16.0     # A4 minus 2 x 2.5 cm
FIGDIR = r"F:\longxiawenjian\gaoxueyaqiguanzhou\Gut_Heart_RAS_8genes\05_results\figures_final"
FIG_WIDTH_CM = 15.5
FIG_MAX_H_CM = 19.0
MISSING_FIGS = []
SRC_RE = re.compile(r"Source:\s*(\S+\.png)\.?\s*$")


def png_px(path):
    with open(path, "rb") as f:
        head = f.read(33)
    return struct.unpack(">II", head[16:24])


def add_figure(doc, name):
    """Place the rendered panel directly under its legend."""
    path = os.path.join(FIGDIR, name)
    if not os.path.exists(path):
        MISSING_FIGS.append(name)
        return
    px_w, px_h = png_px(path)
    w_cm = FIG_WIDTH_CM
    h_cm = w_cm * px_h / float(px_w)
    if h_cm > FIG_MAX_H_CM:
        h_cm = FIG_MAX_H_CM
        w_cm = h_cm * px_w / float(px_h)
    doc.add_picture(path, width=Cm(w_cm), height=Cm(h_cm))
    fig = doc.paragraphs[-1]
    pf = fig.paragraph_format
    pf.space_before = Pt(6)
    pf.space_after = Pt(14)
    pf.line_spacing = 1.0
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.keep_together = True
    fig.alignment = WD_ALIGN_PARAGRAPH.CENTER


SUP_MAP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3",
           "\u2074": "4", "\u2075": "5", "\u2076": "6", "\u2077": "7",
           "\u2078": "8", "\u2079": "9", "\u207b": "\u2212",
           "\u207a": "+", "\u207f": "n"}
SUP_RE = re.compile("[" + "".join(SUP_MAP) + "]+")

# statistical / mathematical variables that take italics; multi-letter
# abbreviations (OR, CI, SD, SE, AUC, FDR, TPM) deliberately excluded
STAT_ITALIC = re.compile(
    r"(?<![A-Za-z0-9])(P)(?=\s*[=<>\u2264\u2265]\s*[\d.])"
    r"|(?<![A-Za-z0-9])(F)(?=\s*[=<>]\s*[\d.])"
    r"|(?<![A-Za-z0-9])(n)(?=\s*=\s*\d)"
    r"|(\u03b2)")

MD_ITALIC = re.compile(r"(\*[^*\n]+?\*)")
TABLE_POW = re.compile(r"x\s+10\^(-?\d+)")
BARE_POW = re.compile(r"(?<![\d.])10\^(-?\d+)")


def to_super(digits):
    inv = {v: k for k, v in SUP_MAP.items()}
    out = []
    for ch in digits:
        out.append(inv.get(ch, inv.get("\u2212") if ch == "-" else ch))
    return "".join(out)


def normalise(text):
    """`x 10^-8` / `10^-8` -> real superscript glyphs, then handled below.

    Markdown code spans (file names in the figure legends) lose their backticks,
    and `+/-` is unified to `±` because the Methods already use the latter.
    """
    text = text.replace("`", "").replace("+/-", "\u00b1")
    text = TABLE_POW.sub(lambda m: "\u00d7 10" + to_super(m.group(1)), text)
    text = BARE_POW.sub(lambda m: "10" + to_super(m.group(1)), text)
    return text


# ------------------------------------------------------------- run emission

def set_font(run, size, cn_font, italic=False, sup=False, bold=False):
    run.font.size = size
    run.font.bold = bool(bold)
    run.font.italic = italic
    rPr = run._element.get_or_add_rPr()
    rf = rPr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts")
        rPr.insert(0, rf)
    rf.set(qn("w:ascii"), EN)
    rf.set(qn("w:hAnsi"), EN)
    rf.set(qn("w:eastAsia"), cn_font)
    rf.set(qn("w:cs"), EN)
    if sup:
        va = rPr.find(qn("w:vertAlign"))
        if va is None:
            va = OxmlElement("w:vertAlign")
            rPr.append(va)
        va.set(qn("w:val"), "superscript")


STRAY_STARS = [0]


def strip_bold(text):
    """Drop `**...**` emphasis markers before italics are resolved.

    Order matters: the italic splitter matches the *innermost* asterisk pair, so
    a leftover `**` makes it swallow the wrong span (e.g. `**Yan Wang**` would
    pair the second `*` with the first `*` of the next name).  Bold is never
    rendered, so the markers are simply removed.
    """
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    # source typo `**gene*` -> a proper single-asterisk italic span
    text = re.sub(r"\*\*([^*\n]+?)\*", r"*\1*", text)
    return text


def emit(p, text, size, cn_font=CN_BODY, md_italic=False, allow_bold=False):
    """Write one markdown inline string as runs.

    Bold is never rendered in running text.  It is rendered *inside table cells*
    when `allow_bold` is set, because there the source markers carry meaning the
    table note depends on (eGene status in Table 2, the summary row in Table 3)
    rather than acting as decorative emphasis.
    """
    if allow_bold:
        for part in re.split(r"(\*\*.+?\*\*)", text):
            if not part:
                continue
            if part.startswith("**") and part.endswith("**") and len(part) > 4:
                _emit_italics(p, part[2:-2], size, cn_font, md_italic, True)
            else:
                _emit_italics(p, part, size, cn_font, md_italic, False)
        return
    _emit_italics(p, strip_bold(text), size, cn_font, md_italic, False)


def _emit_italics(p, text, size, cn_font, md_italic, bold):
    for seg in MD_ITALIC.split(text):
        if not seg:
            continue
        forced = md_italic
        if seg.startswith("*") and seg.endswith("*") and len(seg) > 2:
            forced, seg = True, seg[1:-1]
        elif "*" in seg:
            # an unpaired marker survived; drop it rather than print it
            STRAY_STARS[0] += seg.count("*")
            seg = seg.replace("*", "")
        pos = 0
        for m in STAT_ITALIC.finditer(seg):
            piece = seg[pos:m.start()]
            _emit_plain(p, piece, size, cn_font, forced, bold)
            r = p.add_run(m.group(0))
            set_font(r, size, cn_font, italic=True, sup=False, bold=bold)
            pos = m.end()
        _emit_plain(p, seg[pos:], size, cn_font, forced, bold)


def _emit_plain(p, seg, size, cn_font, italic, bold=False):
    pos = 0
    for m in SUP_RE.finditer(seg):
        if m.start() > pos:
            r = p.add_run(seg[pos:m.start()])
            set_font(r, size, cn_font, italic=italic, bold=bold)
        r = p.add_run("".join(SUP_MAP[c] for c in m.group(0)))
        set_font(r, size, cn_font, italic=italic, sup=True, bold=bold)
        pos = m.end()
    if pos < len(seg):
        r = p.add_run(seg[pos:])
        set_font(r, size, cn_font, italic=italic, bold=bold)


# ------------------------------------------------------------------ blocks

def parse(md):
    lines = md.split("\n")
    blocks, i, in_refs = [], 0, False
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        if s == "---":
            i += 1
            continue

        m = re.match(r"^(#{1,3})\s+(.*)$", s)
        if m:
            lvl, text = len(m.group(1)), m.group(2)
            if lvl <= 2:
                in_refs = (text.strip() == "References")
            blocks.append(("h%d" % lvl, text))
            i += 1
            continue

        if s.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in
                         lines[i].strip().strip("|").split("|")]
                if not set("".join(cells)) <= set("-: "):
                    rows.append([normalise(c) for c in cells])
                i += 1
            blocks.append(("tbl", rows))
            continue

        if in_refs and re.match(r"^\d+\.\s", s):
            blocks.append(("ref", normalise(s)))
            i += 1
            continue

        if re.match(r"^[-*]\s+", s) and not s.startswith("**"):
            buf = [re.sub(r"^[-*]\s+", "", s)]
            i += 1
            while i < len(lines):
                nxt = lines[i].strip()
                if (not nxt or nxt.startswith("|") or nxt.startswith("#")
                        or re.match(r"^[-*]\s+", nxt)):
                    break
                buf.append(nxt)
                i += 1
            blocks.append(("li", " ".join(buf)))
            continue

        buf = [s]
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if (not nxt or nxt == "---" or nxt.startswith("|")
                    or nxt.startswith("#") or nxt.startswith("![")
                    or re.match(r"^[-*]\s+", nxt)
                    or (in_refs and re.match(r"^\d+\.\s", nxt))):
                break
            buf.append(nxt)
            i += 1
        blocks.append(("p", normalise(" ".join(buf))))
    return blocks


# ------------------------------------------------------------------ tables

def add_table(doc, rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]

    plain = [[re.sub(r"[*]", "", c) for c in r] for r in rows]
    weight = [min(max(max(len(r[c]) for r in plain), 8), 46)
              for c in range(ncol)]
    tot = float(sum(weight))
    total = int(CONTENT_CM / 2.54 * 1440)
    widths = [int(total * w / tot) for w in weight]
    widths[-1] += total - sum(widths)

    t = doc.add_table(rows=len(rows), cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    tblPr = t._tbl.tblPr
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)
    mar = OxmlElement("w:tblCellMar")
    for tag, w in (("w:top", 40), ("w:left", 80),
                   ("w:bottom", 40), ("w:right", 80)):
        e = OxmlElement(tag)
        e.set(qn("w:w"), str(w))
        e.set(qn("w:type"), "dxa")
        mar.append(e)
    tblPr.append(mar)

    for ri, row in enumerate(rows):
        hdr = (ri == 0)
        for ci, cell in enumerate(row):
            c = t.cell(ri, ci)
            c.width = __import__("docx").shared.Emu(
                int(widths[ci] / 1440.0 * 914400))
            p = c.paragraphs[0]
            for r in list(p.runs):
                r._element.getparent().remove(r._element)
            pf = p.paragraph_format
            pf.space_before = Pt(2)
            pf.space_after = Pt(2)
            pf.line_spacing = 1.0
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            short = len(plain[ri][ci]) <= 12
            p.alignment = (WD_ALIGN_PARAGRAPH.CENTER if (hdr or short)
                           else WD_ALIGN_PARAGRAPH.LEFT)
            if cell == "":
                r = p.add_run("")
                set_font(r, SZ_SMALL, CN_BODY)
            else:
                emit(p, cell, SZ_SMALL, CN_BODY, allow_bold=True)
            tcPr = c._tc.get_or_add_tcPr()
            va = OxmlElement("w:vAlign")
            va.set(qn("w:val"), "center")
            tcPr.append(va)

    trPr = t.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader")
    th.set(qn("w:val"), "true")
    trPr.append(th)

    set_table_borders(t)
    return t


# ------------------------------------------------------------- page number

def add_page_number(par):
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = par.add_run()
    set_font(r, SZ_SMALL, CN_BODY)
    for kind, txt in (("begin", None), (None, " PAGE "), ("end", None)):
        if kind:
            e = OxmlElement("w:fldChar")
            e.set(qn("w:fldCharType"), kind)
            r._r.append(e)
        else:
            e = OxmlElement("w:instrText")
            e.set(qn("xml:space"), "preserve")
            e.text = txt
            r._r.append(e)


CAPTION_RE = re.compile(r"^(Figure|Table)\s|\ASupplementary (Figure|Table)\s")


def main():
    md = io.open(SRC, encoding="utf-8").read()
    doc = Document()

    st = doc.styles["Normal"]
    st.font.name = EN
    st.font.size = SZ_BODY
    st.element.rPr.rFonts.set(qn("w:eastAsia"), CN_BODY)

    sec = doc.sections[0]
    sec.page_width = Cm(21.0)
    sec.page_height = Cm(29.7)
    sec.top_margin = Cm(2.5)
    sec.bottom_margin = Cm(2.5)
    sec.left_margin = Cm(2.5)
    sec.right_margin = Cm(2.5)
    sec.footer_distance = Cm(1.5)
    add_page_number(sec.footer.paragraphs[0])

    for kind, payload in parse(md):
        if kind == "h1":
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.space_before, pf.space_after = GAP, GAP
            pf.line_spacing = 1.0
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            pf.keep_with_next = True
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            emit(p, payload, SZ_H1, CN_HEAD)
        elif kind == "h2":
            if CAPTION_RE.match(payload):
                p = doc.add_paragraph()
                pf = p.paragraph_format
                pf.space_before, pf.space_after = Pt(12), Pt(4)
                pf.line_spacing = 1.0
                pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
                pf.keep_with_next = True
                p.alignment = (WD_ALIGN_PARAGRAPH.CENTER
                               if payload.startswith("Table")
                               else WD_ALIGN_PARAGRAPH.LEFT)
                emit(p, payload, SZ_SMALL, CN_HEAD)
            else:
                p = doc.add_paragraph()
                pf = p.paragraph_format
                pf.space_before, pf.space_after = GAP, ZERO
                pf.line_spacing = 1.0
                pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
                pf.keep_with_next = True
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
                emit(p, payload, SZ_H2, CN_HEAD)
        elif kind == "h3":
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.space_before, pf.space_after = ZERO, ZERO
            pf.line_spacing = 1.5
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            emit(p, payload, SZ_BODY, CN_HEAD)
        elif kind == "tbl":
            add_table(doc, payload)
        elif kind == "ref":
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.space_before, pf.space_after = ZERO, Pt(4)
            pf.line_spacing = 1.0
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            pf.left_indent = Cm(0.75)
            pf.first_line_indent = Cm(-0.75)
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            emit(p, payload, SZ_SMALL, CN_BODY)
        elif kind == "li":
            p = doc.add_paragraph(style="List Bullet")
            pf = p.paragraph_format
            pf.space_before, pf.space_after = ZERO, ZERO
            pf.line_spacing = 1.5
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            emit(p, payload, SZ_BODY, CN_BODY)
        else:
            p = doc.add_paragraph()
            pf = p.paragraph_format
            pf.space_before, pf.space_after = ZERO, ZERO
            pf.line_spacing = 1.5
            pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            emit(p, payload, SZ_BODY, CN_BODY)
            m = SRC_RE.search(payload)
            if m:
                add_figure(doc, m.group(1))

    if MISSING_FIGS:
        raise SystemExit("missing figure files: %s" % MISSING_FIGS)
    # the target may be open in a viewer, which locks it on Windows
    FALLBACK = DST.replace(".docx", "_new.docx")
    OUT = DST
    try:
        doc.save(DST)
    except PermissionError:
        doc.save(FALLBACK)
        OUT = FALLBACK
        print("LOCKED  %s is open in a viewer; wrote %s instead"
              % (os.path.basename(DST), os.path.basename(FALLBACK)))
    print("stray  : %d unpaired asterisk(s) removed" % STRAY_STARS[0])
    print("written: %s" % OUT)
    print("figures: %d" % len(doc.inline_shapes))
    print("tables : %d" % len(doc.tables))
    print("paras  : %d" % len(doc.paragraphs))
    print("size   : %.1f KB" % (os.path.getsize(OUT) / 1024.0))


if __name__ == "__main__":
    main()
