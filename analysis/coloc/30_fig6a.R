## Fig 6a: regional association at the MTHFR locus.
## eQTL panel on top, hypertension panel below, stacked explicitly so the
## order cannot depend on factor levels.

suppressMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
FIG <- file.path(ROOT, "05_results", "figures_final")

eq <- fread(file.path(ROOT, "03_data_raw", "gtex_full", "MTHFR_Colon_Transverse.csv"))
ou <- fread(file.path(ROOT, "03_data_raw", "outcome", "MTHFR_I9_HYPTENS.csv"))
eq[, `:=`(mb = pos / 1e6, y = -log10(p))]
ou[, `:=`(mb = pos / 1e6, y = -log10(p))]

LEAD_E <- 11790946 / 1e6
LEAD_O <- 11823674 / 1e6

base_theme <- theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8, hjust = 0))

mk <- function(d, col, ymax, show_x) {
  p <- ggplot(d, aes(mb, y)) +
    geom_point(size = 0.7, alpha = 0.85, colour = col) +
    geom_vline(xintercept = c(LEAD_E, LEAD_O), linetype = "dashed",
               linewidth = 0.35, colour = "grey25") +
    annotate("text", x = LEAD_E, y = ymax * 0.97, label = "rs3818762",
             size = 2.4, hjust = -0.05, vjust = 0.4) +
    annotate("text", x = LEAD_O, y = ymax * 0.86, label = "rs12567136",
             size = 2.4, hjust = -0.05, vjust = 0.4) +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
    base_theme
  if (show_x) p + labs(x = "chr1 position (Mb, GRCh38)")
  else p + labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

p1 <- mk(eq, "#2b7bba", max(eq$y), FALSE) +
  labs(y = expression(-log[10]*italic(P)),
       title = sprintf("MTHFR, colon transverse (%s variants)",
                       format(nrow(eq), big.mark = ",")))

p2 <- mk(ou, "#c0392b", max(ou$y), TRUE) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dotted",
             linewidth = 0.35, colour = "grey35") +
  labs(y = expression(-log[10]*italic(P)),
       title = "FinnGen R12 hypertension")

p <- (p1 / p2) +
  plot_annotation(
    title = "MTHFR: expression lead and blood-pressure lead lie 32.7 kb apart",
    theme = theme(plot.title = element_text(size = 9)))

out <- file.path(FIG, "Figure_6a_MTHFR_region.png")
ggsave(out, p, width = 7.0, height = 5.4, dpi = 400)
cat("written:", out, "\n")
cat("eQTL max:", round(max(eq$y), 2), "| outcome max:", round(max(ou$y), 2),
    "| variants:", nrow(eq), "/", nrow(ou), "\n")
