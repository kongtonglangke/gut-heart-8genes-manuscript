# -*- coding: utf-8 -*-
"""Minimal streaming reader for large .xlsx workbooks.

No openpyxl / pandas on this machine, so the OOXML package is read directly.
Rows are yielded as lists of strings; iteration can stop early, which matters
because the master pQTL sheet alone holds >1 M cells.
"""
import re
import zipfile

from lxml import etree

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _col_index(ref):
    """'BC12' -> 54 (0-based column index)."""
    letters = re.match(r"([A-Z]+)", ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


class Workbook:
    def __init__(self, path):
        self.zip = zipfile.ZipFile(path)
        self.shared = self._shared_strings()
        self.sheets = self._sheet_map()

    def _shared_strings(self):
        try:
            raw = self.zip.read("xl/sharedStrings.xml")
        except KeyError:
            return []
        out = []
        for _, el in etree.iterparse(__import__("io").BytesIO(raw),
                                     tag=NS + "si", huge_tree=True):
            out.append("".join(t.text or "" for t in el.iter(NS + "t")))
            el.clear()
        return out

    def _sheet_map(self):
        wb = etree.fromstring(self.zip.read("xl/workbook.xml"))
        rels = etree.fromstring(self.zip.read("xl/_rels/workbook.xml.rels"))
        tgt = {r.get("Id"): r.get("Target") for r in rels}
        out = {}
        for s in wb.iter(NS + "sheet"):
            rid = s.get(RELNS + "id")
            t = tgt.get(rid, "")
            if not t:
                continue
            out[s.get("name")] = "xl/" + t.lstrip("/").replace("xl/", "", 1)
        return out

    def rows(self, name, limit=None):
        """Yield rows as lists of strings."""
        path = self.sheets[name]
        with self.zip.open(path) as fh:
            n = 0
            for _, row in etree.iterparse(fh, tag=NS + "row", huge_tree=True):
                cells = {}
                for c in row.iter(NS + "c"):
                    ref = c.get("r")
                    if not ref:
                        continue
                    v = c.find(NS + "v")
                    isel = c.find(NS + "is")
                    if c.get("t") == "s" and v is not None:
                        try:
                            val = self.shared[int(v.text)]
                        except (ValueError, IndexError):
                            val = ""
                    elif isel is not None:
                        val = "".join(t.text or "" for t in isel.iter(NS + "t"))
                    else:
                        val = v.text if v is not None else ""
                    cells[_col_index(ref)] = (val or "").strip()
                width = max(cells) + 1 if cells else 0
                yield [cells.get(i, "") for i in range(width)]
                row.clear()
                n += 1
                if limit and n >= limit:
                    return
