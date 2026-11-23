## Fig 8: where the colocalisation test is actually powered.
suppressMessages({library(ggplot2); library(data.table)})
ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T <- file.path(ROOT, "05_results", "tables")
FIG <- file.path(ROOT, "05_results", "figures")

phs <- c(I9_HYPTENS = "Hypertension", I9_HEARTFAIL = "Heart failure",
         I9_AF = "Atrial fibrillation", N14_CHRONKIDNEYDIS = "Chronic kidney disease",
         T2D = "Type 2 diabetes")
d <- rbindlist(lapply(names(phs), function(p) {
  x <- fread(file.path(T, sprintf("coloc_%s.csv", p)))
  x[, pheno := phs[[p]]]; x
}))
d[, e := -log10(as.numeric(minP_eqtl))]
d[, o := -log10(as.numeric(minP_out))]
d[, PPH4 := as.numeric(PPH4)]
d[, both := e > -log10(5e-8) & o > -log10(5e-8)]

lab <- d[both == TRUE, .(gene, tissue, e, o)][order(-o)][1]

p <- ggplot(d, aes(e, o)) +
  annotate("rect", xmin = -log10(5e-8), xmax = Inf, ymin = -log10(5e-8), ymax = Inf,
           fill = "#fdd0a2", alpha = 0.35) +
  geom_hline(yintercept = -log10(5e-8), linetype = 3, colour = "grey35") +
  geom_vline(xintercept = -log10(5e-8), linetype = 3, colour = "grey35") +
  geom_point(aes(colour = PPH4), size = 1.8, alpha = 0.9) +
  scale_colour_gradientn(colours = c("#4575b4", "#fee090", "#d73027"),
                         limits = c(0, 0.5), oob = scales::squish,
                         name = "PPH4") +
  annotate("text", x = 13, y = 26,
           label = "both sides significant: 4 tests,\nall MTHFR x hypertension,\nall PPH4 = 0.000-0.001",
           size = 2.5, hjust = 0, colour = "grey20") +
  annotate("segment", x = 12.8, y = 30, xend = 10.5, yend = 42.5,
           arrow = arrow(length = unit(0.12, "cm")), colour = "grey30", linewidth = 0.3) +
  labs(x = expression(-log[10]*"(minimum eQTL "*italic(P)*" in the cis window)"),
       y = expression(-log[10]*"(minimum outcome "*italic(P)*" in the same window)"),
       title = "Signal availability across 240 cis-eQTL x outcome tests",
       subtitle = "the test is only interpretable in the shaded quadrant; everywhere else at least one trait is null") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), plot.title.position = "plot")
ggsave(file.path(FIG, "fig8_signal_quadrants.png"), p, width = 6.4, height = 4.3, dpi = 300)
cat("ok\n")
