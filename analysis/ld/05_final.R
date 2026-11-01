# Step 5 - final LD-based analysis.
# Fine-mapping, formal conditional analysis and coloc.susie for MTHFR
# (the negative) and FES/FURIN (the positive control), using 1000G EUR LD.

suppressPackageStartupMessages({
  library(susieR)
  library(coloc)
})

BASE <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Gut_Heart_RAS_8genes"
LD <- file.path(BASE, "08_ld")
OUT <- file.path(LD, "out")
dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(..., "\n", sep = "")

N_EQTL <- 368
N_CASE <- 154630
N_CTRL <- 345634
N_OUT <- N_CASE + N_CTRL
S_OUT <- N_CASE / N_OUT

cache <- new.env()
get_gt <- function(ldfile) {
  if (!is.null(cache[[ldfile]])) return(cache[[ldfile]])
  tab <- read.table(file.path(LD, ldfile), header = FALSE, sep = "\t",
                    colClasses = c("character", "integer", "character",
                                   "character", "numeric", rep("numeric", 633)))
  g <- list(pos = tab[[2]], ref = tab[[3]], alt = tab[[4]],
            G = as.matrix(tab[, 6:ncol(tab)]))
  cache[[ldfile]] <- g
  g
}

build <- function(gene, ldfile, outcome_rs) {
  gt <- get_gt(ldfile)
  eq <- read.csv(file.path(BASE, "03_data_raw/gtex_full",
                           sprintf("%s_Colon_Transverse.csv", gene)),
                 stringsAsFactors = FALSE)
  ou <- read.csv(file.path(BASE, "03_data_raw/outcome", outcome_rs),
                 stringsAsFactors = FALSE)
  want <- sort(Reduce(intersect, list(gt$pos, eq$pos, ou$pos)))
  idx <- match(want, gt$pos)
  G <- gt$G[idx, , drop = FALSE]
  sds <- apply(G, 1, sd)
  ok <- !is.na(sds) & sds > 0
  G <- G[ok, , drop = FALSE]
  pos <- want[ok]
  R <- cor(t(G))
  eq <- eq[match(pos, eq$pos), ]
  ou <- ou[match(pos, ou$pos), ]
  stopifnot(all(eq$ref == gt$ref[match(pos, gt$pos)]))
  dimnames(R) <- list(eq$rsid, eq$rsid)
  list(pos = pos, R = R, eq = eq, ou = ou,
       z1 = eq$beta / eq$se, z2 = ou$beta / ou$se)
}

# conditional z for a set of index variants given conditioning set
cond_z <- function(z, R, i1, i2) {
  R12 <- R[i1, i2, drop = FALSE]
  R22 <- R[i2, i2, drop = FALSE]
  b <- R12 %*% solve(R22, z[i2])
  v <- 1 - R12 %*% solve(R22, t(R12))
  as.numeric((z[i1] - b) / sqrt(v))
}

fmt_p <- function(z) sprintf("%.2e", 2 * pnorm(-abs(z)))

analyse <- function(gene, ldfile, outcome_rs, a, b) {
  say("\n================ ", gene, " ================")
  d <- build(gene, ldfile, outcome_rs)
  say("joint variants: ", length(d$pos))

  ia <- which(d$eq$rsid == a); ib <- which(d$eq$rsid == b)
  say(sprintf("r(%s, %s) = %+.3f", a, b, d$R[ia, ib]))

  say("\n-- marginal and conditional --")
  say(sprintf("  %-11s eQTL z=%+.2f (P=%s) | HTN z=%+.2f (P=%s)",
              a, d$z1[ia], fmt_p(d$z1[ia]), d$z2[ia], fmt_p(d$z2[ia])))
  say(sprintf("  %-11s eQTL z=%+.2f (P=%s) | HTN z=%+.2f (P=%s)",
              b, d$z1[ib], fmt_p(d$z1[ib]), d$z2[ib], fmt_p(d$z2[ib])))
  say(sprintf("  eQTL  z(%s | %s) = %+.3f (P=%s)", b, a,
              cond_z(d$z1, d$R, ib, ia), fmt_p(cond_z(d$z1, d$R, ib, ia))))
  say(sprintf("  eQTL  z(%s | %s) = %+.3f (P=%s)", a, b,
              cond_z(d$z1, d$R, ia, ib), fmt_p(cond_z(d$z1, d$R, ia, ib))))
  say(sprintf("  HTN   z(%s | %s) = %+.3f (P=%s)", b, a,
              cond_z(d$z2, d$R, ib, ia), fmt_p(cond_z(d$z2, d$R, ib, ia))))
  say(sprintf("  HTN   z(%s | %s) = %+.3f (P=%s)", a, b,
              cond_z(d$z2, d$R, ia, ib), fmt_p(cond_z(d$z2, d$R, ia, ib))))

  # what a single causal variant at a would predict at b
  say(sprintf("\n  single-variant prediction: if a is causal, expected eQTL z(b) = %+.2f (obs %+.2f)",
              d$R[ia, ib] * d$z1[ia], d$z1[ib]))
  say(sprintf("                             if a is causal, expected HTN  z(b) = %+.2f (obs %+.2f)",
              d$R[ia, ib] * d$z2[ia], d$z2[ib]))

  say("\n-- susie_rss fine-mapping --")
  s1 <- susie_rss(z = d$z1, R = d$R, n = N_EQTL, L = 10,
                  estimate_residual_variance = FALSE, max_iter = 2000, verbose = FALSE)
  s2 <- susie_rss(z = d$z2, R = d$R, n = N_OUT, L = 10,
                  estimate_residual_variance = FALSE, max_iter = 2000, verbose = FALSE)
  show <- function(s, dat, tag) {
    if (!length(s$sets$cs)) { say("  ", tag, ": no credible set"); return(invisible()) }
    for (k in seq_along(s$sets$cs)) {
      ii <- s$sets$cs[[k]]
      lead <- ii[which.max(s$pip[ii])]
      say(sprintf("  %s CS%d n=%2d maxPIP=%.3f minCorr=%.2f lead=%s (%s)",
                  tag, k, length(ii), max(s$pip[ii]),
                  s$sets$purity$min.abs.corr[k], dat$rsid[lead], dat$pos[lead]))
    }
  }
  show(s1, d$eq, "eQTL   ")
  show(s2, d$ou, "HTN    ")

  say("\n-- coloc.abf (reproduction) --")
  d1 <- list(beta = d$eq$beta, varbeta = d$eq$se^2, snp = d$eq$rsid,
             position = d$pos, type = "quant", sdY = 1, N = N_EQTL, LD = d$R)
  # both datasets must carry the same SNP identifiers as the LD matrix dimnames
  d2 <- list(beta = d$ou$beta, varbeta = d$ou$se^2, snp = d$eq$rsid,
             position = d$pos, type = "cc", s = S_OUT, N = N_OUT, LD = d$R)
  ab <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  p <- ab$summary
  say(sprintf("  nSNP=%d PPH3=%.4g PPH4=%.4g", p[["nsnps"]], p[["PP.H3.abf"]], p[["PP.H4.abf"]]))
  top <- ab$results[order(ab$results[["SNP.PP.H4"]], decreasing = TRUE)[1:3], ]
  say("  top SNP.PP.H4:")
  for (i in seq_len(nrow(top))) {
    say(sprintf("    %-12s level1=%s level2=%s H4=%s", top$snp[i],
                format(top[["pvalues1"]][i], digits = 3),
                format(top[["pvalues2"]][i], digits = 3),
                format(top[["SNP.PP.H4"]][i], digits = 4)))
  }

  say("\n-- coloc.susie --")
  cs <- try(coloc.susie(s1, s2, d1, d2), silent = TRUE)
  if (inherits(cs, "try-error")) {
    say("  FAILED: ", conditionMessage(attr(cs, "condition")))
  } else if (is.null(cs$summary)) {
    say("  no credible-set pairs")
  } else {
    for (i in seq_len(nrow(cs$summary))) {
      r <- cs$summary[i, ]
      say(sprintf("  CS%d x CS%d : H3=%.4f H4=%.4f", r$idx1, r$idx2,
                  r[["PP.H3.abf"]], r[["PP.H4.abf"]]))
    }
    write.csv(cs$summary, file.path(OUT, sprintf("coloc_susie_%s.csv", gene)),
              row.names = FALSE)
  }
  invisible(d)
}

analyse("FES", "gt_FURIN.tsv", "FURIN_I9_HYPTENS.csv", "rs6227", "rs1894400")
say("\nwritten -> ", OUT)
