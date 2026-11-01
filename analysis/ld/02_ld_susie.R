# Step 2 - LD from 1000G EUR, susieR credible sets, conditional analysis,
# and coloc.susie.  Source kept pure ASCII on purpose (the R console here
# runs under a C locale and mangles non-ASCII literals).

suppressPackageStartupMessages({
  library(susieR)
  library(coloc)
})

BASE <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Gut_Heart_RAS_8genes"
LD   <- file.path(BASE, "08_ld")
OUT  <- file.path(LD, "out")

dir.create(OUT, showWarnings = FALSE)

say <- function(...) cat(..., "\n", sep = "")

read_gt <- function(gene) {
  f <- file.path(LD, sprintf("gt_%s.tsv", gene))
  tab <- read.table(f, header = FALSE, sep = "\t",
                    colClasses = c("character", "integer", "character",
                                   "character", "numeric", rep("numeric", 633)))
  list(pos = tab[[2]], ref = tab[[3]], alt = tab[[4]], af = tab[[5]],
       G = as.matrix(tab[, 6:ncol(tab)]))
}

# dosage correlation among the requested positions
ld_matrix <- function(gt, want) {
  idx <- match(want, gt$pos)
  keep <- !is.na(idx)
  G <- gt$G[idx[keep], , drop = FALSE]
  # drop monomorphic sites, which make cor() return NaN
  sds <- apply(G, 1, sd)
  ok <- !is.na(sds) & sds > 0
  G <- G[ok, , drop = FALSE]
  R <- cor(t(G))
  list(R = R, pos = want[keep][ok])
}

# ---- analysis for one gene/locus ------------------------------------------
run_locus <- function(gene, chrom, eqtl_n, out_n, label) {
  say("\n===== ", gene, " / ", label, " =====")

  gt <- read_gt(gene)
  eq <- read.csv(file.path(BASE, "03_data_raw/gtex_full",
                           sprintf("%s_Colon_Transverse.csv", gene)),
                 stringsAsFactors = FALSE)
  ou <- read.csv(file.path(BASE, "03_data_raw/outcome",
                           sprintf("%s_I9_HYPTENS.csv", gene)),
                 stringsAsFactors = FALSE)

  # variants present in all three sources
  want <- Reduce(intersect, list(gt$pos, eq$pos, ou$pos))
  want <- sort(want)
  say("shared variants: ", length(want))

  ld <- ld_matrix(gt, want)
  pos <- ld$pos
  R <- ld$R
  say("LD matrix: ", nrow(R), " x ", ncol(R))

  eq <- eq[match(pos, eq$pos), ]
  ou <- ou[match(pos, ou$pos), ]

  # align to the 1000G alt allele; there were no flips at this stage but the
  # check is kept so the script cannot silently mis-sign an effect
  stopifnot(all(eq$ref == gt$ref[match(pos, gt$pos)]))

  z_eq <- eq$beta / eq$se
  z_ou <- ou$beta / ou$se

  # ---- univariate susie on the exposure ----------------------------------
  s_eq <- try(susie_rss(z = z_eq, R = R, n = eqtl_n, L = 10,
                        estimate_residual_variance = FALSE,
                        max_iter = 500, verbose = FALSE), silent = TRUE)
  s_ou <- try(susie_rss(z = z_ou, R = R, n = out_n, L = 10,
                        estimate_residual_variance = FALSE,
                        max_iter = 500, verbose = FALSE), silent = TRUE)

  say("susie eQTL : ", class(s_eq)[1],
      if (!inherits(s_eq, "try-error")) paste0(" | sets=", length(s_eq$sets$cs)) else "")
  say("susie outcome: ", class(s_ou)[1],
      if (!inherits(s_ou, "try-error")) paste0(" | sets=", length(s_ou$sets$cs)) else "")

  if (!inherits(s_eq, "try-error") && length(s_eq$sets$cs)) {
    for (k in seq_along(s_eq$sets$cs)) {
      v <- pos[s_eq$sets$cs[[k]]]
      say(sprintf("  eQTL CS%d: n=%d pip_max=%.3f  rsids=%s",
                  k, length(v), max(s_eq$pip[s_eq$sets$cs[[k]]]),
                  paste(eq$rsid[s_eq$sets$cs[[k]]], collapse = ",")))
    }
  }
  if (!inherits(s_ou, "try-error") && length(s_ou$sets$cs)) {
    for (k in seq_along(s_ou$sets$cs)) {
      v <- pos[s_ou$sets$cs[[k]]]
      say(sprintf("  outcome CS%d: n=%d pip_max=%.3f  rsids=%s",
                  k, length(v), max(s_ou$pip[s_ou$sets$cs[[k]]]),
                  paste(ou$rsid[s_ou$sets$cs[[k]]], collapse = ",")))
    }
  }

  # ---- coloc.susie -------------------------------------------------------
  d1 <- list(beta = eq$beta, varbeta = eq$se^2, snp = eq$rsid,
             position = pos, type = "quant", sdY = 1, N = eqtl_n, LD = R)
  d2 <- list(beta = ou$beta, varbeta = ou$se^2, snp = ou$rsid,
             position = pos, type = "cc", s = 154630 / (154630 + 345634),
             N = 154630 + 345634, LD = R)

  cs <- try(coloc.susie(s_eq, s_ou, d1, d2), silent = TRUE)
  if (inherits(cs, "try-error")) {
    say("coloc.susie FAILED: ", conditionMessage(attr(cs, "condition")))
  } else if (is.null(cs$summary)) {
    say("coloc.susie returned no credible-set pairs")
  } else {
    say("coloc.susie pairs:")
    print(cs$summary[, c("idx1", "idx2", "PP.H0.abf", "PP.H1.abf",
                         "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
    write.csv(cs$summary, file.path(OUT, sprintf("coloc_susie_%s.csv", gene)),
              row.names = FALSE)
  }

  list(pos = pos, R = R, eq = eq, ou = ou, z_eq = z_eq, z_ou = z_ou,
       s_eq = s_eq, s_ou = s_ou)
}

# ---- conditional z-score test ---------------------------------------------
# z_cond for SNP1 given SNP2, using the joint LD matrix (no individual data)
cond_z <- function(z, R, i1, i2) {
  R12 <- R[i1, i2, drop = FALSE]
  R22 <- R[i2, i2, drop = FALSE]
  b <- as.numeric(R12 %*% solve(R22, z[i2]))
  v <- as.numeric(1 - R12 %*% solve(R22, t(R12)))
  list(z = (z[i1] - b) / sqrt(v), shrink = b / z[i1])
}

report_cond <- function(res, gene, a, b, tag) {
  pos <- res$pos
  ia <- which(res$eq$rsid == a)
  ib <- which(res$eq$rsid == b)
  if (!length(ia) || !length(ib)) {
    say(sprintf("  %s: variant missing (a=%d b=%d)", tag, length(ia), length(ib)))
    return(invisible())
  }
  say(sprintf("\n-- %s -- %s (pos %d) vs %s (pos %d), %.1f kb apart",
              tag, a, pos[ia], b, pos[ib], abs(pos[ia] - pos[ib]) / 1000))
  say(sprintf("   marginal z: %s = %+.2f (P=%.2e) | %s = %+.2f (P=%.2e)",
              a, res$z_eq[ia], 2 * pnorm(-abs(res$z_eq[ia])),
              b, res$z_eq[ib], 2 * pnorm(-abs(res$z_eq[ib]))))
  for (pair in list(c(ia, ib), c(ib, ia))) {
    cc <- cond_z(res$z_eq, res$R, pair[1], pair[2])
    say(sprintf("   conditioning on %s: z(%s) = %+.3f (P=%.2e)",
                res$eq$rsid[pair[2]], res$eq$rsid[pair[1]],
                cc$z, 2 * pnorm(-abs(cc$z))))
  }
}

# ---- main -----------------------------------------------------------------
mthfr <- run_locus("MTHFR", 1, 368, 500264, "transverse colon vs FinnGen R12")
report_cond(mthfr, "MTHFR", "rs3818762", "rs12567136", "MTHFR two-variant test")
report_cond(mthfr, "MTHFR", "rs12404124", "rs198397", "MTHFR sign-flip pair")

furin <- run_locus("FURIN", 15, 368, 500264, "transverse colon vs FinnGen R12")
report_cond(furin, "FURIN", "rs6227", "rs1894400", "FURIN/FES lead vs BP lead")

saveRDS(list(mthfr = mthfr, furin = furin), file.path(OUT, "ld_susie.rds"))
say("\ndone -> ", OUT)
