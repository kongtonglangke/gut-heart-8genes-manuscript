# gut-heart 8-gene triage — manuscript working copy

Private working repository for the manuscript **"Genetic triage of 'gut–heart axis'
candidates: colocalization and Mendelian randomization of eight hypertension
transcriptomic markers"** (Yan Wang, Jing Bai).

## Layout

| Path | Contents |
|---|---|
| `manuscript/` | `*_v1.md` is the editable source of truth; `*_v1.src.md` is the untouched original (read-only); `*_v1.clean.md` is the Chinese-free intermediate; `*_v1d.docx` is the current rendered submission file |
| `build/` | `01_clean_md.py` → `02_md2docx.py` → `03_verify.py`: the three-stage build with a 27-item automated check |
| `analysis/` | pQTL extraction and review-application scripts, the 1000G LD extraction, and the faceted Figure 1 redraw |
| `figures/` | 300–400 dpi panels as embedded in the manuscript |
| `docs/` | revision log with the per-comment verdict for each review round |

## Rebuilding the Word file

```
python build/01_clean_md.py      # strip Chinese working notes -> *.clean.md
python build/02_md2docx.py       # render -> *.docx, embedding the 10 figures
python build/03_verify.py        # 27 assertions; must print "pass 27 / fail 0"
```

`02_md2docx.py` expects the figure directory at
`Gut_Heart_RAS_8genes/05_results/figures_final/`. If the target `.docx` is open in a
viewer on Windows the save falls back to `*_new.docx` and says so.

## Data sources (all public, no controlled access)

- GTEx v8 cis-eQTL, as uniformly reprocessed by the eQTL Catalogue
- FinnGen release 12 (hypertension and four other outcomes)
- Human Protein Atlas single-cell type expression
- MiBioGen microbiome GWAS
- deCODE plasma pQTL supplementary tables (Ferkingstad 2021)
- 1000 Genomes high-coverage panel (GRCh38), for LD

## Notes

- The manuscript is English throughout; the Chinese in `docs/` is commentary.
- Every numeric value in the figures and tables is produced by the scripts, not typed.
- `*_v1.src.md` is kept read-only as an immutable record of the received draft.
