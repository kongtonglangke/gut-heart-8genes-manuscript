## Fig 10: gut microbiota -> hypertension, two-sample MR (211 taxa).
suppressMessages({library(ggplot2); library(data.table); library(ggrepel)})
ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
T <- file.path(ROOT, "05_results", "tables")
FIG <- file.path(ROOT, "05_results", "figures")
f <- file.path(T, "mr_microbiome_htn.csv")
if (!file.exists(f)) { cat("no data yet\n"); quit(status = 0) }

d <- fread(f)
d[, OR := as.numeric(OR)][, p := as.numeric(pval)][, fdr := as.numeric(fdr)]
d[, log10p := -log10(p)]
d[, sig := fdr < 0.05]
d[, short := sub("^(class|order|family|genus|species)\\.", "", taxon)]
d[, short := sub("\\.id\\.[0-9]+$", "", short)]

p <- ggplot(d, aes(OR, log10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, colour = "grey45") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_point(aes(colour = sig), size = 1.9, alpha = 0.85) +
  geom_text_repel(data = d[sig == TRUE][order(-log10p)][1:min(8, .N)],
                  aes(label = short), size = 2.3, max.overlaps = 20,
                  min.segment.length = 0, segment.size = 0.2) +
  scale_colour_manual(values = c(`TRUE` = "#cb181d", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "FDR<0.05", `FALSE` = "not significant"),
                      name = NULL) +
  scale_x_log10() +
  labs(x = "OR of hypertension per unit increase in taxon abundance (log scale)",
       y = expression(-log[10]*italic(P)),
       title = "J  Gut microbiota -> hypertension: two-sample MR across 211 MiBioGen taxa",
       subtitle = sprintf("%d taxa estimable; %d reach FDR<0.05", nrow(d), sum(d$sig))) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title.position = "plot")
ggsave(file.path(FIG, "fig10_microbiome_mr.png"), p, width = 6.6, height = 4.4, dpi = 300)
cat("ok\n")
