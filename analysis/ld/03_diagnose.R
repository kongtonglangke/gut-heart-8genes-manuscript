# Step 3 - diagnostic check of the conditional analysis.
# Prints the actual LD coefficients so the sign-flip claim can be judged
# against the correlation structure rather than against intuition.

suppressPackageStartupMessages(library(susieR))

BASE <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
LD <- file.path(BASE, "08_ld")
say <- function(...) cat(..., "\n", sep = "")

read_gt <- function(gene) {
  tab <- read.table(file.path(LD, sprintf("gt_%s.tsv", gene)), header = FALSE,
                    sep = "\t",
                    colClasses = c("character", "integer", "character",
                                   "character", "numeric", rep("numeric", 633)))
  list(pos = tab[[2]], ref = tab[[3]], alt = tab[[4]], G = as.matrix(tab[, 6:ncol(tab)]))
}

for (gene in c("MTHFR", "FURIN")) {
  say("\n########## ", gene, " ##########")
  gt <- read_gt(gene)
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

  say("variants in joint set: ", length(pos))

  # allele frequencies in 1000G EUR, to show the LD is not driven by rarities
  af <- rowMeans(G) / 2
  say("MAF range: ", sprintf("%.3f - %.3f", min(pmin(af, 1 - af)), max(pmin(af, 1 - af))))

  show_pair <- function(a, b) {
    ia <- which(eq$rsid == a); ib <- which(eq$rsid == b)
    if (!length(ia) || !length(ib)) { say("  missing: ", a, " / ", b); return(invisible()) }
    r <- R[ia, ib]
    say(sprintf("\n%s (%.0f kb) vs %s : r = %+.3f  r2 = %.3f", a,
                pos[ia] / 1000, b, r, r^2))
    say(sprintf("   eQTL betas %+.4f / %+.4f   z %+.2f / %+.2f",
                eq$beta[ia], eq$beta[ib], eq$beta[ia] / eq$se[ia], eq$beta[ib] / eq$se[ib]))
    say(sprintf("   outcome betas %+.4f / %+.4f   P %.2e / %.2e",
                ou$beta[ia], ou$beta[ib], ou$p[ia], ou$p[ib]))
    # if a single causal variant sits at a, the marginal z at b is r * z at a
    say(sprintf("   expected z(b) if a is the sole causal variant: %+.2f (observed %+.2f)",
                r * eq$beta[ia] / eq$se[ia], eq$beta[ib] / eq$se[ib]))
  }

  if (gene == "MTHFR") {
    show_pair("rs3818762", "rs12567136")
    show_pair("rs12404124", "rs198397")
    show_pair("rs3818762", "rs1023252")
    show_pair("rs3818762", "rs12404124")
    show_pair("rs3818762", "rs198397")
    # strongest eQTL variant in the window
    k <- which.min(eq$p)
    say(sprintf("\nstrongest eQTL: %s pos=%.0f beta=%+.4f P=%.3e", eq$rsid[k], pos[k],
                eq$beta[k], eq$p[k]))
    say(sprintf("strongest HTN  : %s pos=%.0f beta=%+.4f P=%.3e",
                ou$rsid[which.min(ou$p)], pos[which.min(ou$p)],
                ou$beta[which.min(ou$p)], min(ou$p)))
  } else {
    show_pair("rs6227", "rs1894400")
    k <- which.min(eq$p)
    say(sprintf("\nstrongest eQTL: %s pos=%.0f beta=%+.4f P=%.3e", eq$rsid[k], pos[k],
                eq$beta[k], eq$p[k]))
    say(sprintf("strongest HTN  : %s pos=%.0f beta=%+.4f P=%.3e",
                ou$rsid[which.min(ou$p)], pos[which.min(ou$p)],
                ou$beta[which.min(ou$p)], min(ou$p)))
    say(sprintf("rs6227 eQTL P=%.3e  (window rank %d of %d)",
                eq$p[which(eq$rsid == "rs6227")],
                rank(eq$p)[which(eq$rsid == "rs6227")], nrow(eq)))
  }
}
