# Data sources and accessions

All inputs are public summary statistics. Nothing individual-level was accessed.

## 1 cis-eQTL — GTEx v8 via the eQTL Catalogue

| Item | Value |
|---|---|
| Resource | eQTL Catalogue, uniformly reprocessed GTEx v8 |
| URL | https://www.ebi.ac.uk/eqtl/ |
| FTP pattern | `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/<tissue>.tsv.gz` |
| Index | same path with `.tbi` suffix, queried by `tabix` |
| Build | GRCh38 |
| Column order | variant, r2, pvalue, mt_obj, gene, maf, gene_id, median_tpm, beta, se, an, ac, chrom, pos, ref, alt, type, rsid |

Tissues and sample sizes used: colon sigmoid (318), colon transverse (368),
small intestine terminal ileum (174), heart left ventricle (386),
kidney cortex (73), whole blood (670).

## 2 Outcomes — FinnGen release 12

| Item | Value |
|---|---|
| URL | https://r12.finngen.fi/ |
| Build | GRCh38 |
| Endpoints | `I9_HYPTENS` 154,630 cases / 345,634 controls; `I9_HEARTFAIL` 37,653/462,695; `I9_AF` 63,532/527,935; `N14_CHRONKIDNEYDIS` 5,712/482,219; `T2D` 33,043/308,872 |

Note: FinnGen R12 does not include chromosome X, so *ACE2* cannot be tested in
this framework.

## 3 Single-cell expression — Human Protein Atlas

| Item | Value |
|---|---|
| URL | https://www.proteinatlas.org |
| File | `rna_single_cell_type.tsv.zip` |
| Unit | normalised counts per million (nCPM) |
| Cell types | 154 profiled; 32 in the cardiovascular, intestinal epithelial, immune and kidney compartments displayed |

## 4 Microbiome — MiBioGen

| Item | Value |
|---|---|
| URL | https://mibiogen.org/ (the old `mibiogen.gcc.rug.nl` host is retired) |
| File | `MiBioGen_QmbQTL_summary_<rank>.zip` |
| Taxa | 211 |
| Instrument rule | single lead variant per taxon, F > 10 |

## 5 Plasma proteome — deCODE

| Item | Value |
|---|---|
| URL | https://www.decode.com/summarydata/ |
| Study | Ferkingstad et al. 2021, *Nature Genetics* 53:1712–1721 |
| Platform | SomaScan v4, 4,907 aptamers, 35,559 participants |
| cis definition | P < 5 × 10⁻⁸ within 1 Mb of the transcription start site |

## 6 Reference panel — 1000 Genomes

| Item | Value |
|---|---|
| URL | https://www.internationalgenome.org |
| Release | high-coverage phased panel, GRCh38 |
| Samples | 3,202 total; 633 European-ancestry samples used |
| Retained variants | 884 at 1p36 (*MTHFR*), 710 at 15q26.1 (*FES*/*FURIN*) |

## Harmonisation

Alleles were matched on position with explicit strand checks. Three counts are
kept distinct throughout: **tested** (all variants in the 100 kb window with
summary statistics), **analysed** (harmonised, ambiguous and monomorphic sites
removed), **fine-mapped** (retained after intersecting with the reference panel).
For *MTHFR* these were 2,803 / 909 / 884; for the *ACE* control, 351 / 270 / 242.
