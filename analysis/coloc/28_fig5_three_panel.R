## Figure 5: three-panel view of the FURIN/FES locus (chr15, GRCh38).
## Same variants, same window, same outcome GWAS: the signal belongs to FES.
suppressMessages({
  library(ggplot2); library(data.table); library(patchwork)
})

ROOT <- "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes"
GEX  <- file.path(ROOT, "03_data_raw", "gtex_full")
OUTD <- file.path(ROOT, "03_data_raw", "outcome")
FIGF <- file.path(ROOT, "05_results", "figures_final")

theme_set(theme_bw(base_size = 8.5) +
            theme(panel.grid = element_blank(),
                  plot.title.position = "plot",
                  title = element_text(size = 9)))

furin <- fread(file.path(GEX, "FURIN_Colon_Transverse.csv"))
fes   <- fread(file.path(GEX, "FES_Colon_Transverse.csv"))
out   <- fread(file.path(OUTD, "FURIN_I9_HYPTENS.csv"))

prep <- function(d) d[, .(pos, rsid, p = as.numeric(p))][!is.na(p)][, neglog := -log10(p)]

a <- prep(furin)
b <- prep(fes)
c <- prep(out)

panel <- function(d, ylab, subtitle, colr) {
  ggplot(d, aes(pos / 1e6, neglog)) +
    geom_point(size = 0.6, colour = colr, alpha = 0.65) +
    labs(x = NULL, y = ylab, subtitle = subtitle) +
    theme(axis.title.y = element_text(size = 7.5))
}

pa <- panel(a, expression(-log[10]*italic(P)~"(eQTL)"),
            "a  FURIN cis-eQTL, transverse colon", "#2c7fb8")
pb <- panel(b, expression(-log[10]*italic(P)~"(eQTL)"),
            "b  FES cis-eQTL, transverse colon", "#0571b0")
pc <- panel(c, expression(-log[10]*italic(P)~"(hypertension)"),
            "c  FinnGen R12 hypertension", "#cb181d") +
  scale_x_continuous(name = "chr15 position (Mb, GRCh38)") +
  geom_hline(yintercept = -log10(5e-8), linetype = 3, colour = "grey30")

## lead variant marker
lead <- "rs6227"
ybot <- floor(min(a$neglog, b$neglog, c$neglog))
for (p in list(pa, pb, pc)) NULL
xtop <- unique(c(a[rsid == lead, pos], c[rsid == lead, pos])) / 1e6

pa <- pa + geom_vline(xintercept = xtop, linetype = 2, colour = "grey20") +
  annotate("text", x = xtop, y = Inf, label = lead, vjust = 1.9, hjust = -0.1,
           size = 2.3, colour = "grey20")
pb <- pb + geom_vline(xintercept = xtop, linetype = 2, colour = "grey20") +
  annotate("text", x = xtop, y = Inf, label = lead, vjust = 1.9, hjust = -0.1,
           size = 2.3, colour = "grey20")
pc <- pc + geom_vline(xintercept = xtop, linetype = 2, colour = "grey20") +
  annotate("text", x = xtop, y = Inf, label = lead, vjust = 1.9, hjust = -0.1,
           size = 2.3, colour = "grey20")

p <- (pa / pb / pc) +
  plot_annotation(
    title = "The 15q26.1 hypertension signal is a FES eQTL and not a FURIN eQTL",
    subtitle = sprintf("identical cis window and variant set in all three panels; %s:\nP = 3.2e-34 (hypertension), 4.2e-9 (FES), 0.35 (FURIN)", lead),
    theme = theme(plot.subtitle = element_text(size = 7.5, lineheight = 1.1)))

ggsave(file.path(FIGF, "Figure_5_FURIN_FES_region.png"), p, width = 8.0, height = 7.5, dpi = 300)
cat("ok\n")
