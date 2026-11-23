# Calibrated colocalization triage of eight gut–heart candidate genes

Code and derived result tables for the manuscript:

> **Two internal controls calibrate colocalization: a three-layer genetic triage of eight gut–heart candidate genes for hypertension**
> Yan Wang, Jing Bai. Submitted to *Annals of Human Genetics*.

This repository contains the analysis scripts and the derived tables behind every
number in the paper. No individual-level or identifiable data are included;
all inputs are public summary statistics, listed under **Data sources** below.

---

## What is here

```
analysis/
  coloc/         phenome-wide cis-eQTL × outcome colocalization (coloc.abf)
  simulation/    operating-characteristic simulation on the 240 real tests
  ld/            1000G EUR reference panel, fine-mapping, conditional analysis
  pqtl/          deCODE plasma pQTL extraction and cis-pQTL tests
results/
  coloc/         per-outcome colocalization tables, expression and cell-type summaries
  simulation/    per-design power, false-positive grid, matched analysis
  pqtl/          pQTL association table and panel coverage
docs/
  data_sources.md      accession numbers and download routes for every input
  file_manifest.md     one line per file: what it is and which table it feeds
```

## Data sources

| Input | Version | Access |
|---|---|---|
| GTEx cis-eQTL (via eQTL Catalogue) | GTEx v8 | https://www.ebi.ac.uk/eqtl/ |
| FinnGen | release 12 | https://r12.finngen.fi/ |
| Human Protein Atlas single-cell | v2 | https://www.proteinatlas.org |
| MiBioGen microbiome GWAS | 2021 | https://mibiogen.org/ |
| deCODE plasma proteome | Ferkingstad 2021 | https://www.decode.com/summarydata/ |
| 1000 Genomes high-coverage panel | GRCh38 | https://www.internationalgenome.org |

Six eQTL tissues were queried: sigmoid colon (n = 318), transverse colon (n = 368),
terminal ileum (n = 174), left ventricle (n = 386), kidney cortex (n = 73) and
whole blood (n = 670).

## Requirements

- R 4.3.3 with `coloc` 5.2.3, `susieR` 0.14.2, `data.table`, `ggplot2`
- Python 3.13 with `pandas`, `numpy`, `scipy`
- `tabix`/`bgzip` for remote range queries
- About 2 GB of disk if the full simulation inputs are regenerated

## Reproducing

Download the inputs listed in `docs/data_sources.md` into `data_raw/`, then run in order:

```bash
# 1 colocalization, 240 tests
python analysis/coloc/03_matrix.py
Rscript analysis/coloc/06_coloc.R

# 2 linkage-disequilibrium panel and fine-mapping
python analysis/ld/01_extract_ld.py
Rscript analysis/ld/02_ld_susie.R
Rscript analysis/ld/04_coloc_ld.R

# 3 calibration simulation (rebuilds the 240 designs, then simulates)
python analysis/simulation/01_build_tests.py
Rscript analysis/simulation/02_validate.R
python analysis/simulation/03_simulate.py
python analysis/simulation/04_matched.py
python analysis/simulation/05_figure5.py
```

`analysis/simulation/01_build_tests.py` regenerates the two large per-variant
files (`variants_240.csv`, `r_variants_240.csv`, about 40 MB together) that are
deliberately **not** committed here. The smaller derived tables are committed,
so every number in the manuscript can be checked without re-running anything.

## Outputs that map to the paper

| File | Feeds |
|---|---|
| `results/coloc/coloc_master.csv` | Table 1, Table S1 |
| `results/coloc/coloc_positive_control.csv` | Table S4 (FES 0.985, ACE 0.730) |
| `results/coloc/expression_tau.csv` | Figure 1, Table S2 |
| `results/coloc/celltype_gut_focus.csv` | Figure S3 |
| `results/simulation/sim_by_tissue.csv` | Figure 5a, Table S9 |
| `results/simulation/sim_grid.csv` | Figure 5b |
| `results/simulation/sim_by_design.csv` | Figure 5c, Table S10 |
| `results/pqtl/pqtl_associations.csv` | Table S6 |

## Citation

If you use this code, please cite the manuscript above. A machine-readable record
is provided in `CITATION.cff`.

## Licence

Code: MIT (`LICENSE`). Derived tables: CC BY 4.0 (`LICENSE-DATA`).
