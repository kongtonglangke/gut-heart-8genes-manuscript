## M3  colocalisation: GTEx v8 cis-eQTL  x  FinnGen R12 hypertension
## Input : 03_data_raw/gtex_full/<gene>_<tissue>.csv
##         03_data_raw/outcome/<gene>_<pheno>.csv
## Output: 05_results/tables/coloc_<pheno>.csv
## Note  : source kept ASCII-only on purpose (C locale + Chinese source =
##         parse errors on this machine).

suppressMessages({
  library(coloc)
})

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
GTEX <- file.path(ROOT, "03_data_raw", "gtex_full")
OUTD <- file.path(ROOT, "03_data_raw", "outcome")
RESD <- file.path(ROOT, "05_results", "tables")
dir.create(RESD, showWarnings = FALSE, recursive = TRUE)

CHROM <- c(LRP1 = "12", CTSD = "11", MTHFR = "1", AUTS2 = "7",
           FURIN = "15", CST3 = "20", FCER1G = "1", TBXAS1 = "7")
TISSUES <- c("Colon_Sigmoid", "Colon_Transverse",
             "Small_Intestine_Terminal_Ileum", "Heart_Left_Ventricle",
             "Kidney_Cortex", "Whole_Blood")

PHENO <- commandArgs(trailingOnly = TRUE)
if (length(PHENO) == 0) PHENO <- "I9_HYPTENS"

## Sample sizes read from the FinnGen R12 variant endpoint (metadata only;
## that endpoint is NOT used for statistics).
PHENO_N <- list(
  I9_HYPTENS         = c(154630, 345634),
  I9_HEARTFAIL       = c(37653, 462695),
  I9_AF              = c(63532, 252810),
  N14_CHRONKIDNEYDIS = c(12787, 480448),
  T2D                = c(82878, 403489),
  I9_CARDMYO         = c(7649, 375343)
)
if (!PHENO %in% names(PHENO_N)) stop("no sample size recorded for ", PHENO)
NCASE <- PHENO_N[[PHENO]][1]
NCONT <- PHENO_N[[PHENO]][2]
NOUT <- NCASE + NCONT
SFRAC <- NCASE / NOUT

P12_MAIN <- 1e-5
P12_GRID <- c(1e-6, 1e-5, 1e-4)

read_eqtl <- function(gene, tissue) {
  f <- file.path(GTEX, sprintf("%s_%s.csv", gene, tissue))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) < 50) return(NULL)
  d$chrom <- CHROM[[gene]]
  d$N <- d$an / 2
  d <- d[d$se > 0 & is.finite(d$se) & is.finite(d$beta) & is.finite(d$p), ]
  d <- d[!duplicated(paste(d$pos, d$ref, d$alt)), ]
  d
}

read_out <- function(gene, pheno) {
  f <- file.path(OUTD, sprintf("%s_%s.csv", gene, pheno))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$se > 0 & is.finite(d$se) & is.finite(d$beta) & is.finite(d$p), ]
  d <- d[!duplicated(paste(d$pos, d$ref, d$alt)), ]
  d
}

## merge on position, harmonise effect alleles, return aligned data.frame
harmonise <- function(e, o) {
  m <- merge(e, o, by = "pos", suffixes = c(".e", ".o"))
  if (nrow(m) == 0) return(NULL)
  same <- m$ref.e == m$ref.o & m$alt.e == m$alt.o
  flip <- m$ref.e == m$alt.o & m$alt.e == m$ref.o
  m <- m[same | flip, , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  fl <- m$ref.e == m$alt.o & m$alt.e == m$ref.o
  m$beta.e[fl] <- -m$beta.e[fl]
  m$maf <- pmin(m$af, 1 - m$af)
  m <- m[m$maf > 0.005 & m$maf < 0.995, , drop = FALSE]
  m <- m[!duplicated(m$pos), ]
  m
}

rows <- list()
for (gene in names(CHROM)) {
  o <- read_out(gene, PHENO)
  if (is.null(o) || nrow(o) < 50) {
    cat(sprintf("%-7s outcome window missing or too small\n", gene))
    next
  }
  for (tissue in TISSUES) {
    e <- read_eqtl(gene, tissue)
    if (is.null(e)) next
    m <- harmonise(e, o)
    if (is.null(m) || nrow(m) < 50) {
      cat(sprintf("%-7s %-30s shared SNPs < 50 (n=%s)\n", gene, tissue,
                  ifelse(is.null(m), 0, nrow(m))))
      next
    }
    D1 <- list(snp = paste0("v", m$pos), beta = m$beta.e, varbeta = m$se.e^2,
               N = median(m$N), MAF = m$maf, type = "quant")
    D2 <- list(snp = paste0("v", m$pos), beta = m$beta.o, varbeta = m$se.o^2,
               N = NOUT, MAF = m$maf, type = "cc", s = SFRAC)
    res <- tryCatch(suppressWarnings(coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4,
                                               p12 = P12_MAIN)),
                    error = function(err) NULL)
    if (is.null(res)) {
      cat(sprintf("%-7s %-30s coloc.abf FAILED\n", gene, tissue))
      next
    }
    pp <- res$summary
    ## sensitivity across prior p12
    sens <- sapply(P12_GRID, function(p12) {
      r2 <- tryCatch(suppressWarnings(coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4,
                                                p12 = p12)),
                     error = function(err) NULL)
      if (is.null(r2)) NA_real_ else as.numeric(r2$summary["PP.H4.abf"])
    })
    top <- res$results[which.max(res$results$SNP.PP.H4), ]
    rows[[length(rows) + 1]] <- data.frame(
      gene = gene, tissue = tissue, pheno = PHENO, n_snp = nrow(m),
      minP_eqtl = min(m$p.e), minP_out = min(m$p.o),
      PPH0 = as.numeric(pp[["PP.H0.abf"]]), PPH1 = as.numeric(pp[["PP.H1.abf"]]),
      PPH2 = as.numeric(pp[["PP.H2.abf"]]), PPH3 = as.numeric(pp[["PP.H3.abf"]]),
      PPH4 = as.numeric(pp[["PP.H4.abf"]]),
      H4_H3 = as.numeric(pp[["PP.H4.abf"]]) /
        max(as.numeric(pp[["PP.H3.abf"]]), 1e-12),
      top_snp = top$snp, top_snp_pp = top$SNP.PP.H4,
      H4_p12_1e6 = sens[1], H4_p12_1e5 = sens[2], H4_p12_1e4 = sens[3],
      stringsAsFactors = FALSE)
    cat(sprintf("%-7s %-30s n=%4d  PPH3=%.3f PPH4=%.3f  H4/H3=%6.2f\n",
                gene, tissue, nrow(m), as.numeric(pp[["PP.H3.abf"]]),
                as.numeric(pp[["PP.H4.abf"]]),
                as.numeric(pp[["PP.H4.abf"]]) /
                  max(as.numeric(pp[["PP.H3.abf"]]), 1e-12)))
  }
}

if (length(rows) == 0) stop("no results")
out <- do.call(rbind, rows)
out <- out[order(-out$PPH4), ]
f <- file.path(RESD, sprintf("coloc_%s.csv", PHENO))
write.csv(out, f, row.names = FALSE)
cat("\nSaved:", f, "\n")
cat(sprintf("positive (PPH4>=0.8 & H4/H3>=3): %d / %d\n",
            sum(out$PPH4 >= 0.8 & out$H4_H3 >= 3), nrow(out)))
