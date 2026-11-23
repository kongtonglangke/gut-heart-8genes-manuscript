## Fig 6 (complete matrix version): cell-type localisation.
suppressMessages({library(ggplot2); library(data.table)})
ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T <- file.path(ROOT, "05_results", "tables")
FIG <- file.path(ROOT, "05_results", "figures")

GENES <- c("LRP1","CTSD","MTHFR","AUTS2","FURIN","CST3","FCER1G","TBXAS1")
d <- fread(file.path(T, "celltype_gut_focus.csv"))
d[, gene := factor(gene, levels = GENES)]

## order cell types by compartment then by mean nCPM
ord <- d[, .(m = mean(ncpm)), by = .(compartment, cell_type)][order(compartment, -m)]
d[, cell_type := factor(cell_type, levels = ord$cell_type)]
d[, gene := factor(gene, levels = rev(GENES))]

p <- ggplot(d, aes(cell_type, gene, fill = log10(ncpm + 1))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
                       name = expression(log[10]*"(nCPM+1)")) +
  labs(x = NULL, y = NULL,
       title = "Cell-type localisation of the eight candidates across 32 cell types of four gut-heart compartments",
       subtitle = "complete matrix, not thresholded: 0 means at or below the HPA reporting floor (0.1 nCPM)") +
  theme_bw(base_size = 8.5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.4),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 7),
        legend.key.height = unit(1.1, "cm"),
        plot.title.position = "plot")
FIGF <- file.path(ROOT, "05_results", "figures_final")
ggsave(file.path(FIGF, "Supplementary_Figure_S1_celltypes.png"), p,
       width = 9.6, height = 3.6, dpi = 400)
cat("ok\n")
