## Reusable pairwise coloc: one eQTL csv  x  one outcome csv
## Used for the positive control (ACE kidney cortex vs FinnGen hypertension)
## and for any single-gene follow-up.
##   Rscript 07_coloc_pair.R <eqtl.csv> <outcome.csv> <chrom> <label>

suppressMessages(library(coloc))
a <- commandArgs(trailingOnly = TRUE)
eqtl_f <- a[1]; out_f <- a[2]; chrom <- a[3]; label <- a[4]
NCASE <- 154630; NCONT <- 345634
NOUT <- NCASE + NCONT; SFRAC <- NCASE / NOUT

e <- read.csv(eqtl_f, stringsAsFactors = FALSE)
o <- read.csv(out_f, stringsAsFactors = FALSE)
e$N <- e$an / 2
e <- e[e$se > 0 & is.finite(e$se) & is.finite(e$beta), ]
o <- o[o$se > 0 & is.finite(o$se) & is.finite(o$beta), ]
e <- e[!duplicated(paste(e$pos, e$ref, e$alt)), ]
o <- o[!duplicated(paste(o$pos, o$ref, o$alt)), ]

m <- merge(e, o, by = "pos", suffixes = c(".e", ".o"))
same <- m$ref.e == m$ref.o & m$alt.e == m$alt.o
flip <- m$ref.e == m$alt.o & m$alt.e == m$ref.o
m <- m[same | flip, , drop = FALSE]
fl <- m$ref.e == m$alt.o & m$alt.e == m$ref.o
m$beta.e[fl] <- -m$beta.e[fl]
m$maf <- pmin(m$af, 1 - m$af)
m <- m[m$maf > 0.005 & m$maf < 0.995, ]
m <- m[!duplicated(m$pos), ]

cat(sprintf("%s: eQTL n=%d  outcome n=%d  shared n=%d  (flipped %d)\n",
            label, nrow(e), nrow(o), nrow(m), sum(fl)))

D1 <- list(snp = paste0("v", m$pos), beta = m$beta.e, varbeta = m$se.e^2,
           N = median(m$N), MAF = m$maf, type = "quant")
D2 <- list(snp = paste0("v", m$pos), beta = m$beta.o, varbeta = m$se.o^2,
           N = NOUT, MAF = m$maf, type = "cc", s = SFRAC)
for (p12 in c(1e-6, 1e-5, 1e-4)) {
  r <- suppressWarnings(coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = p12))
  s <- r$summary
  cat(sprintf("  p12=%.0e  PPH0=%.4f PPH1=%.4f PPH2=%.4f PPH3=%.4f PPH4=%.4f\n",
              p12, s[["PP.H0.abf"]], s[["PP.H1.abf"]], s[["PP.H2.abf"]],
              s[["PP.H3.abf"]], s[["PP.H4.abf"]]))
}
cat(sprintf("  minP eQTL=%.2e  minP outcome=%.2e\n", min(m$p.e), min(m$p.o)))
