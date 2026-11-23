# Step 4 - formal colocalization with the 1000G EUR LD matrix.
# coloc.abf is rerun on exactly the same variants so the two approaches are
# comparable, then coloc.susie is fitted with a named LD matrix.

suppressPackageStartupMessages({
  library(susieR)
  library(coloc)
})

BASE <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
LD <- file.path(BASE, "08_ld")
OUT <- file.path(LD, "out")
dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(..., "\n", sep = "")

N_EQTL <- 368
N_CASE <- 154630
N_CTRL <- 345634
N_OUT <- N_CASE + N_CTRL

build <- function(gene) {
  tab <- read.table(file.path(LD, sprintf("gt_%s.tsv", gene)), header = FALSE,
                    sep = "\t",
                    colClasses = c("character", "integer", "character",
                                   "character", "numeric", rep("numeric", 633)))
  gt <- list(pos = tab[[2]], ref = tab[[3]], alt = tab[[4]],
             G = as.matrix(tab[, 6:ncol(tab)]))
  eq <- read.csv(file.path(BASE, "03_data_raw/gtex_full",
                           sprintf("%s_Colon_Transverse.csv", gene)),
                 stringsAsFactors = FALSE)
  ou <- read.csv(file.path(BASE, "03_data_raw/outcome",
                           sprintf("%s_I9_HYPTENS.csv", gene)),
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
  dn <- eq$rsid
  dimnames(R) <- list(dn, dn)
  list(pos = pos, R = R, eq = eq, ou = ou)
}

for (gene in c("MTHFR", "FURIN")) {
  say("\n########## ", gene, " ##########")
  b <- build(gene)
  say("variants: ", length(b$pos))

  d1 <- list(beta = b$eq$beta, varbeta = b$eq$se^2, snp = b$eq$rsid,
             position = b$pos, type = "quant", sdY = 1, N = N_EQTL, LD = b$R)
  d2 <- list(beta = b$ou$beta, varbeta = b$ou$se^2, snp = b$ou$rsid,
             position = b$pos, type = "cc", s = N_CASE / N_OUT,
             N = N_OUT, LD = b$R)

  say("\n-- coloc.abf, priors p1=p2=1e-4, p12=1e-5 --")
  ab <- coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  pp <- ab$summary
  say(sprintf("nSNP=%d  PPH0=%.3g PPH1=%.3g PPH2=%.3g PPH3=%.3g PPH4=%.3g",
              pp[["nsnps"]], pp[["PP.H0.abf"]], pp[["PP.H1.abf"]], pp[["PP.H2.abf"]],
              pp[["PP.H3.abf"]], pp[["PP.H4.abf"]]))
  say("top 5 by SNP.PP.H4:")
  o <- order(ab$results$SNP.PP.H4, decreasing = TRUE)[1:5]
  print(ab$results[o, c("snp", "pvalues1", "pvalues2", "SNP.PP.H4")])

  say("\n-- susie_rss --")
  s1 <- susie_rss(z = b$eq$beta / b$eq$se, R = b$R, n = N_EQTL, L = 10,
                  estimate_residual_variance = FALSE, max_iter = 1000,
                  verbose = FALSE)
  s2 <- susie_rss(z = b$ou$beta / b$ou$se, R = b$R, n = N_OUT, L = 10,
                  estimate_residual_variance = FALSE, max_iter = 1000,
                  verbose = FALSE)
  say("eQTL credible sets: ", length(s1$sets$cs), " | outcome: ", length(s2$sets$cs))

  summarise_sets <- function(s, dat, tag) {
    if (!length(s$sets$cs)) { say("  ", tag, ": none"); return(invisible()) }
    for (k in seq_along(s$sets$cs)) {
      ii <- s$sets$cs[[k]]
      say(sprintf("  %s CS%d: n=%2d maxpip=%.3f cov=%.2f purity=%.2f | %s",
                  tag, k, length(ii), max(s$pip[ii]),
                  s$sets$coverage[k], s$sets$purity$min.abs.corr[k],
                  paste(head(dat$rsid[ii], 6), collapse = ",")))
    }
  }
  summarise_sets(s1, b$eq, "eQTL")
  summarise_sets(s2, b$ou, "outcome")

  # where does the eQTL lead sit on the outcome side?
  lead <- which.max(b$eq$beta / b$eq$se)
  say(sprintf("\neQTL lead %s -> outcome z = %+.2f (P = %.2e)",
              b$eq$rsid[lead], b$ou$beta[lead] / b$ou$se[lead],
              2 * pnorm(-abs(b$ou$beta[lead] / b$ou$se[lead]))))

  say("\n-- coloc.susie --")
  cs <- try(coloc.susie(s1, s2, d1, d2), silent = TRUE)
  if (inherits(cs, "try-error")) {
    say("FAILED: ", conditionMessage(attr(cs, "condition")))
  } else if (is.null(cs$summary)) {
    say("no credible-set pairs returned")
  } else {
    print(cs$summary)
    write.csv(cs$summary, file.path(OUT, sprintf("coloc_susie_%s.csv", gene)),
              row.names = FALSE)
  }

  # ---- joint conditional model for the disputed MTHFR pair ---------------
  if (gene == "MTHFR") {
    say("\n-- joint fit conditioning eQTL on rs3818762 --")
    j <- which(b$eq$rsid == "rs3818762")
    s1c <- susie_rss(z = b$eq$beta / b$eq$se, R = b$R, n = N_EQTL, L = 5,
                     estimate_residual_variance = FALSE, max_iter = 1000,
                     verbose = FALSE, prior_condition_number = 1e6)
    # explicit conditioning: regress out the lead from the z vector
    zc <- (b$eq$beta / b$eq$se)[-j] - b$R[-j, j] * (b$eq$beta[j] / b$eq$se[j])
    zc <- zc / sqrt(pmax(1 - b$R[-j, j]^2, 1e-6))
    k <- which(b$eq$rsid[-j] == "rs12567136")
    say(sprintf("  rs12567136 conditioned on rs3818762: z = %+.3f (P = %.2e)",
                zc[k], 2 * pnorm(-abs(zc[k]))))
    say(sprintf("  largest remaining |z|: %s z=%+.2f",
                b$eq$rsid[-j][which.max(abs(zc))], zc[which.max(abs(zc))]))
  }
}
