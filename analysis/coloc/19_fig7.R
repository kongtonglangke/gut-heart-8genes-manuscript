## Fig 7: MR-coloc discordance dissection for MTHFR.
suppressMessages({library(ggplot2); library(data.table)})
ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T <- file.path(ROOT, "05_results", "tables")
FIG <- file.path(ROOT, "05_results", "figures")

d <- fread(file.path(T, "mthfr_wald_by_position.csv"))
d[, pos_mb := pos / 1e6]
d[, tissue := factor(tissue, levels = c("Colon_Transverse", "Colon_Sigmoid",
                                        "Small_Intestine_Terminal_Ileum", "Whole_Blood"))]

p <- ggplot(d, aes(pos_mb, wald)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(colour = tissue), size = 0.9, alpha = 0.7) +
  geom_vline(xintercept = 11.790946, linetype = 2, colour = "#2c7fb8") +
  geom_vline(xintercept = 11.823674, linetype = 2, colour = "#cb181d") +
  annotate("text", x = 11.790946, y = -0.33, label = "eQTL lead\nrs3818762",
           size = 2.3, hjust = 1.1, colour = "#2c7fb8") +
  annotate("text", x = 11.823674, y = -0.33, label = "BP lead\nrs12567136",
           size = 2.3, hjust = -0.1, colour = "#cb181d") +
  scale_colour_brewer(palette = "Set2", name = NULL) +
  labs(x = "chr1 position (Mb, GRCh38)",
       y = "Wald ratio (log OR per SD of expression)",
       title = "MTHFR: cis-MR is significant but the Wald ratio is not constant across the locus",
       subtitle = "a single shared variant would give a flat line; the median drifts with distance from the BP lead") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        plot.title.position = "plot",
        plot.margin = margin(6, 8, 6, 26))
FIGF <- file.path(ROOT, "05_results", "figures_final")
ggsave(file.path(FIGF, "Figure_6b_MTHFR_wald_gradient.png"), p,
       width = 7.0, height = 4.0, dpi = 400)
cat("ok\n")
