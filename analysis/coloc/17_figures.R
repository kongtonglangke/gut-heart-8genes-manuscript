## Figures for the gut-heart 8-gene project.
## Source kept ASCII-only (C locale + UTF-8 source = parse errors here).
suppressMessages({
  library(ggplot2); library(data.table); library(patchwork); library(scales)
})

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T    <- file.path(ROOT, "05_results", "tables")
FIG  <- file.path(ROOT, "05_results", "figures")
GEX  <- file.path(ROOT, "03_data_raw", "gtex_full")
OUTD <- file.path(ROOT, "03_data_raw", "outcome")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

GENES <- c("LRP1","CTSD","MTHFR","AUTS2","FURIN","CST3","FCER1G","TBXAS1")
GUT   <- c("Colon_Sigmoid","Colon_Transverse","Small_Intestine_Terminal_Ileum")
theme_set(theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(),
                  strip.background = element_rect(fill = "grey92", colour = NA)))

## ------------------------------------------------------------------ Fig 1
e <- fread(file.path(T, "..", "..", "03_data_raw", "expression_median_tpm.csv"))
e[, log2tpm := log2(median_tpm + 1)]

grp <- function(t) {
  ifelse(grepl("Colon|Intestine|Stomach|Esophagus", t), "Gut",
  ifelse(grepl("Heart|Artery", t), "Cardiovascular",
  ifelse(grepl("Kidney|Liver|Adrenal|Pancreas|Salivary|Thyroid", t), "Metabolic",
  ifelse(grepl("Whole_Blood|Spleen|Cells_", t), "Immune",
  ifelse(grepl("Brain|Nerve", t), "Neuro",
  ifelse(grepl("Skin|Adipose|Muscle|Fibroblast|Ovary|Uterus|Vagina|Breast|Prostate|Testis|Bladder|Cervix|Fallopian", t), "Structural",
         "Other"))))))
}
e[, group := grp(tissue)]
e[, gene := factor(gene, levels = GENES)]
ord <- e[gene == "LRP1"][order(group, -log2tpm)][, tissue]
e[, tissue := factor(tissue, levels = ord)]

p1 <- ggplot(e, aes(tissue, gene, fill = log2tpm)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(colours = c("#f7f7f7", "#fdd49e", "#fdbb84", "#e34a33", "#7f0000"),
                       name = expression(log[2]*"(TPM+1)")) +
  labs(x = NULL, y = NULL,
       title = "Tissue expression atlas of the eight candidates (GTEx v8, 54 tissues, median TPM)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "right", legend.key.height = unit(0.9, "cm"))
ggsave(file.path(FIG, "fig1_expression_atlas.png"), p1, width = 9.2, height = 3.4, dpi = 300)

## ------------------------------------------------------------------ Fig 2
m <- fread(file.path(T, "matrix_raw.csv"))
m[, minP := as.numeric(minP)]
m[, neglog := -log10(minP)]
m[, eGene := minP < 5e-8]
tis_order <- c("Colon_Sigmoid","Colon_Transverse","Small_Intestine_Terminal_Ileum",
               "Heart_Left_Ventricle","Kidney_Cortex","Whole_Blood")
m[, tissue := factor(tissue, levels = tis_order)]
m[, gene := factor(gene, levels = rev(GENES))]
m[, lab := ifelse(eGene, sprintf("%.0f", neglog), sprintf("%.1f", neglog))]

p2 <- ggplot(m, aes(tissue, gene, fill = neglog)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 2.6, colour = "grey15") +
  geom_tile(data = m[eGene == TRUE], colour = "black", fill = NA, linewidth = 0.8) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",
                       midpoint = 4, name = expression(-log[10]*"(min "*italic(P)*")")) +
  labs(x = NULL, y = NULL,
       title = "Detectability of cis-eQTL signal (eight genes x six tissues)",
       subtitle = "black box = eGene (P<5e-8); numbers are -log10 min P") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7))
ggsave(file.path(FIG, "fig2_detectability.png"), p2, width = 5.6, height = 3.2, dpi = 300)

## ------------------------------------------------------------------ Fig 3
phs <- c(I9_HYPTENS = "Hypertension", I9_HEARTFAIL = "Heart failure",
         I9_AF = "Atrial fibrillation", N14_CHRONKIDNEYDIS = "Chronic kidney disease",
         T2D = "Type 2 diabetes")
cl <- rbindlist(lapply(names(phs), function(p) {
  f <- file.path(T, sprintf("coloc_%s.csv", p))
  if (!file.exists(f)) return(NULL)
  d <- fread(f); d[, pheno := phs[[p]]]; d
}))
cl[, PPH3 := as.numeric(PPH3)][, PPH4 := as.numeric(PPH4)]
cl[, inGut := tissue %in% GUT]

p3 <- ggplot(cl, aes(PPH3, PPH4)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = 0.8, linetype = 2, colour = "grey45") +
  geom_point(aes(colour = pheno, shape = inGut), size = 1.7, alpha = 0.85) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "gut tissue", `FALSE` = "non-gut tissue"),
                     name = NULL) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  annotate("text", x = 0.55, y = 0.95, label = "shared", size = 2.6, colour = "grey30") +
  annotate("text", x = 0.55, y = 0.03, label = "distinct", size = 2.6, colour = "grey30") +
  labs(x = expression("PPH3  (two distinct causal variants)"),
       y = expression("PPH4  (shared causal variant)"),
       title = "Colocalization of gut and cardiac cis-eQTL with five cardiometabolic outcomes",
       subtitle = sprintf("%d gene x tissue x outcome tests; 0 reach PPH4>=0.8", nrow(cl))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
ggsave(file.path(FIG, "fig3_coloc_scatter.png"), p3, width = 6.4, height = 4.2, dpi = 300)

## ------------------------------------------------------------------ Fig 4/5
region_plot <- function(gene, chrom, tissue, pheno, mark_rs, title, fout,
                        mark_rs2 = NULL, extra = NULL) {
  ev <- fread(file.path(GEX, sprintf("%s_%s.csv", gene, tissue)))
  ov <- fread(file.path(OUTD, sprintf("%s_%s.csv", gene, pheno)))
  ev[, neglog := -log10(as.numeric(p))]
  ov[, neglog := -log10(as.numeric(p))]
  a <- ggplot(ev, aes(pos / 1e6, neglog)) +
    geom_point(size = 0.6, colour = "#2c7fb8", alpha = 0.6) +
    labs(x = NULL, y = expression(-log[10]*italic(P)~"(eQTL)"),
         subtitle = sprintf("%s, %s (n=%d variants)", gene, gsub("_", " ", tissue), nrow(ev)))
  b <- ggplot(ov, aes(pos / 1e6, neglog)) +
    geom_point(size = 0.6, colour = "#cb181d", alpha = 0.6) +
    labs(x = paste0("chr", chrom, " position (Mb, GRCh38)"),
         y = expression(-log[10]*italic(P)~"(FinnGen)")) +
    geom_hline(yintercept = -log10(5e-8), linetype = 3, colour = "grey30")
  for (rs in c(mark_rs, mark_rs2)) {
    if (is.null(rs)) next
    if (rs %in% ev$rsid) {
      x <- ev[rsid == rs, pos] / 1e6; y <- ev[rsid == rs, neglog]
      a <- a + geom_vline(xintercept = x, linetype = 2, colour = "grey20") +
        annotate("text", x = x, y = Inf, label = rs, vjust = 1.6, hjust = -0.1,
                 size = 2.4, colour = "grey20")
    }
    if (rs %in% ov$rsid) {
      x <- ov[rsid == rs, pos] / 1e6; y <- ov[rsid == rs, neglog]
      b <- b + geom_vline(xintercept = x, linetype = 2, colour = "grey20") +
        annotate("text", x = x, y = Inf, label = rs, vjust = 1.6, hjust = -0.1,
                 size = 2.4, colour = "grey20")
    }
  }
  p <- a / b + plot_annotation(title = title)
  ggsave(fout, p, width = 6.4, height = 4.6, dpi = 300)
  invisible(p)
}

region_plot("FURIN", "15", "Colon_Sigmoid", "I9_HYPTENS", "rs6227",
            "FURIN: the strongest hypertension variant lies in the 3'UTR and has no local eQTL effect",
            file.path(FIG, "fig4_furin_region.png"))

region_plot("MTHFR", "1", "Colon_Transverse", "I9_HYPTENS", "rs3818762",
            "MTHFR: expression lead and blood-pressure lead lie 32.7 kb apart",
            file.path(FIG, "fig5_mthfr_region.png"), mark_rs2 = "rs12567136")

## ------------------------------------------------------------------ Fig 6
ct <- fread(file.path(T, "celltype_ncpm.csv"))
ct[, gene := factor(gene, levels = GENES)]
keep <- ct[, .(maxv = max(value)), by = .(gene, cell_type)][maxv >= 100, unique(cell_type)]
ct2 <- ct[cell_type %in% keep]
ct2[, val2 := pmin(value, 4000)]
p6 <- ggplot(ct2, aes(cell_type, gene, fill = log10(val2 + 1))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#f7fbff", "#9ecae1", "#2171b5", "#08306b"),
                       name = expression(log[10]*"(nCPM+1)")) +
  labs(x = NULL, y = NULL,
       title = "Cell-type localisation (Human Protein Atlas single-cell)",
       subtitle = "cell types with >=100 nCPM in at least one gene; blank = below HPA reporting threshold") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5))
ggsave(file.path(FIG, "fig6_celltype.png"), p6, width = 8.4, height = 3.4, dpi = 300)

cat("figures written to", FIG, "\n")
