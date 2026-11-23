## Fig 9: cis-MR forest (strict instrument tier only).
suppressMessages({library(ggplot2); library(data.table)})
ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T <- file.path(ROOT, "05_results", "tables")
FIG <- file.path(ROOT, "05_results", "figures")

phs <- c(I9_HYPTENS = "Hypertension", I9_HEARTFAIL = "Heart failure",
         I9_AF = "Atrial fibrillation", N14_CHRONKIDNEYDIS = "Chronic kidney disease",
         T2D = "Type 2 diabetes")

d <- fread(file.path(T, "mr_results.csv"))
d <- d[iv_set == "strict_5e-8" & status == "ok"]
d[, OR := as.numeric(OR)][, ci_lo := as.numeric(ci_lo)][, ci_hi := as.numeric(ci_hi)]
d[, p := as.numeric(pval)]
d[, pheno := factor(phs[pheno], levels = rev(phs))]
d[, gene := factor(gene, levels = c("CST3", "MTHFR", "LRP1"))]
d[, sig := p < 0.05]

p <- ggplot(d, aes(OR, pheno, colour = sig)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22, linewidth = 0.5) +
  geom_point(size = 2.1) +
  facet_wrap(~ gene, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = c(`TRUE` = "#cb181d", `FALSE` = "grey35"),
                      labels = c(`TRUE` = "P<0.05", `FALSE` = "not significant"),
                      name = NULL) +
  labs(x = "OR of outcome per SD increase in intestinal expression (95% CI)",
       y = NULL,
       title = "cis-Mendelian randomization with strict instruments (3 of 8 genes have one)",
       subtitle = "one instrument per gene-tissue pair, so the estimate is a single-instrument Wald ratio") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        plot.title.position = "plot")
FIGF <- file.path(ROOT, "05_results", "figures_final")
ggsave(file.path(FIGF, "Supplementary_Figure_S3_cisMR_forest.png"), p,
       width = 7.4, height = 3.3, dpi = 400)
cat("ok\n")
