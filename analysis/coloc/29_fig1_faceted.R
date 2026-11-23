## Figure 1, faceted redraw.
## The single-panel original packs 54 tissue labels into one axis, which is
## illegible at print width.  Grouping is taken verbatim from 17_figures.R so that
## the compartments, the ordering rule and the colour ramp all stay identical.
## Source kept ASCII-only (C locale + UTF-8 source = parse errors here).
suppressMessages({library(ggplot2); library(data.table)})

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T    <- file.path(ROOT, "05_results", "tables")
FIGD <- file.path(ROOT, "05_results", "figures_final")

GENES <- c("LRP1","CTSD","MTHFR","AUTS2","FURIN","CST3","FCER1G","TBXAS1")

grp <- function(t) {
  ifelse(grepl("Colon|Intestine|Stomach|Esophagus", t), "Gut",
  ifelse(grepl("Heart|Artery", t), "Cardiovascular",
  ifelse(grepl("Kidney|Liver|Adrenal|Pancreas|Salivary|Thyroid", t), "Metabolic",
  ifelse(grepl("Whole_Blood|Spleen|Cells_", t), "Immune",
  ifelse(grepl("Brain|Nerve", t), "Neuro",
  ifelse(grepl("Skin|Adipose|Muscle|Fibroblast|Ovary|Uterus|Vagina|Breast|Prostate|Testis|Bladder|Cervix|Fallopian", t), "Structural",
         "Other"))))))
}

e <- fread(file.path(ROOT, "03_data_raw", "expression_median_tpm.csv"))
e[, log2tpm := log2(median_tpm + 1)]
e[, group := grp(tissue)]
e[, gene := factor(gene, levels = GENES)]

## order tissues by decreasing LRP1 within compartment, exactly as the original
ord <- e[gene == "LRP1"][order(group, -log2tpm)][, tissue]
e[, tissue := factor(tissue, levels = ord)]
e[, group := factor(group, levels = c("Gut", "Cardiovascular", "Metabolic",
                                      "Immune", "Neuro", "Structural", "Other"))]
e <- e[!is.na(group)]

p <- ggplot(e, aes(gene, tissue, fill = log2tpm)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#f7f7f7", "#fdd49e", "#fdbb84", "#e34a33", "#7f0000"),
    name = expression(log[2]*"(median TPM + 1)"),
    breaks = c(0, 5, 10), limits = c(0, 11)) +
  guides(fill = guide_colourbar(barwidth = unit(5.0, "cm"),
                              barheight = unit(0.28, "cm"),
                              title.position = "top", title.hjust = 0.5,
                              ticks.colour = "black", frame.colour = "black")) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size = 7.5, angle = 0),
        axis.text.y = element_text(size = 6.2),
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.28, "cm"),
        legend.box.spacing = unit(0.15, "cm"),
        plot.margin = margin(t = 4, r = 10, b = 4, l = 4),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.2),
        legend.margin = margin(t = 2, r = 6, b = 2, l = 6),
        panel.spacing.y = unit(0.12, "cm"))

out <- file.path(FIGD, "Figure_1_expression_atlas.png")
if (!file.exists(file.path(FIGD, "Figure_1_expression_atlas_ORIGINAL.png"))) {
  file.copy(out, file.path(FIGD, "Figure_1_expression_atlas_ORIGINAL.png"))
}
ggsave(out, p, width = 6.6, height = 7.9, dpi = 400)
cat("written:", out, "\n")
cat("n tissues:", length(unique(e$tissue)), "| n genes:", length(unique(e$gene)), "\n")
