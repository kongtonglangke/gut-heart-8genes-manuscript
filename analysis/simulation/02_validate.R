## Validate the Python re-implementation of coloc.abf against the R package.
## Runs the original pipeline end to end (harmonise by position, flip alleles,
## MAF filter) for all 8 genes x 6 tissues x 5 outcomes and writes the posteriors
## to data/r_coloc_240.csv. Source kept ASCII-only: this machine is under a C
## locale and Chinese characters in R source cause parse errors.

suppressMessages(library(coloc))

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
GTEX <- file.path(ROOT, "03_data_raw", "gtex_full")
OUTD <- file.path(ROOT, "03_data_raw", "outcome")
DEST <- file.path(ROOT, "09_simulation", "data")

GENES <- c("LRP1", "CTSD", "MTHFR", "AUTS2", "FURIN", "CST3", "FCER1G", "TBXAS1")
TISSUES <- c("Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum",
             "Heart_Left_Ventricle", "Kidney_Cortex", "Whole_Blood")
PHENOS <- c("I9_HYPTENS", "I9_HEARTFAIL", "I9_AF", "N14_CHRONKIDNEYDIS", "T2D")
PHENO_N <- list(I9_HYPTENS = c(154630, 345634), I9_HEARTFAIL = c(37653, 462695),
                I9_AF = c(63532, 252810), N14_CHRONKIDNEYDIS = c(12787, 480448),
                T2D = c(82878, 403489))

read_eqtl <- function(gene, tissue) {
  f <- file.path(GTEX, sprintf("%s_%s.csv", gene, tissue))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) < 50) return(NULL)
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

out <- list()
var <- list()
tid <- 0
for (gene in GENES) {
  for (pheno in PHENOS) {
    o <- read_out(gene, pheno)
    if (is.null(o) || nrow(o) < 50) next
    nn <- PHENO_N[[pheno]]
    nout <- sum(nn); sfrac <- nn[1] / nout
    for (tissue in TISSUES) {
      e <- read_eqtl(gene, tissue)
      if (is.null(e)) next
      m <- harmonise(e, o)
      if (is.null(m) || nrow(m) < 50) next
      D1 <- list(snp = paste0("v", m$pos), beta = m$beta.e, varbeta = m$se.e^2,
                 N = median(m$N), MAF = m$maf, type = "quant")
      D2 <- list(snp = paste0("v", m$pos), beta = m$beta.o, varbeta = m$se.o^2,
                 N = nout, MAF = m$maf, type = "cc", s = sfrac)
      r <- tryCatch(suppressWarnings(coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)),
                    error = function(err) NULL)
      if (is.null(r)) next
      tid <- tid + 1
      s <- r$summary
      out[[length(out) + 1]] <- data.frame(
        test_id = tid, gene = gene, tissue = tissue, pheno = pheno,
        m_r = as.integer(s[["nsnps"]]), n_eqtl = median(m$N),
        PPH0 = s[["PP.H0.abf"]], PPH1 = s[["PP.H1.abf"]], PPH2 = s[["PP.H2.abf"]],
        PPH3 = s[["PP.H3.abf"]], PPH4 = s[["PP.H4.abf"]])
      var[[length(var) + 1]] <- data.frame(
        test_id = tid, pos = m$pos, maf = m$maf,
        V_e = m$se.e^2, V_o = m$se.o^2, z_e = m$beta.e / m$se.e, z_o = m$beta.o / m$se.o)
    }
  }
}
res <- do.call(rbind, out)
vres <- do.call(rbind, var)
write.csv(res, file.path(DEST, "r_coloc_240.csv"), row.names = FALSE)
write.csv(vres, file.path(DEST, "r_variants_240.csv"), row.names = FALSE)
cat("R tests:", nrow(res), " variants:", nrow(vres), "\n")
