# Fig 10: microbiome -> hypertension single-instrument MR volcano
# ASCII-only source (project rule: avoid C-locale encoding issues)
library(ggplot2)

d <- read.csv("F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes/05_results/tables/mr_microbiome_htn.csv",
              stringsAsFactors = FALSE)
d$logor <- log(d$OR)
d$mlogp <- -log10(d$pval)
d$grp <- ifelse(d$fdr < 0.05 & d$OR > 1, "FDR<0.05, risk",
         ifelse(d$fdr < 0.05 & d$OR < 1, "FDR<0.05, protective", "ns"))

lab <- d[d$taxon %in% c("genus.Bifidobacterium.id.436",
                        "genus.Romboutsia.id.11347",
                        "genus.Victivallis.id.2256"), ]
lab$name <- sub("\\.id\\.\\d+$", "", lab$taxon)
lab$name <- paste0(lab$name, "\n(", lab$rsid, ")")

p <- ggplot(d, aes(logor, mlogp, colour = grp)) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_text(data = lab, aes(label = name), colour = "black",
            size = 3, vjust = -0.9, lineheight = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  scale_colour_manual(values = c("FDR<0.05, risk" = "#C0392B",
                                 "FDR<0.05, protective" = "#2471A3",
                                 "ns" = "grey75")) +
  labs(title = "Gut microbiota -> hypertension: single-instrument cis-MR (107 taxa)",
       subtitle = "MiBioGen lead variants (F>10) x FinnGen R12 I9_HYPTENS; dashed line P=0.05",
       caption = "Actinobacteria/Bifidobacterium clade shares one instrument: rs182549 (LCT/MCM6 locus)",
       x = "log(OR) per SD taxon abundance", y = "-log10(P)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Gut_Heart_RAS_8genes/05_results/figures/fig10_microbiome_mr.png",
       p, width = 7.5, height = 5.5, dpi = 200)
cat("fig10 saved; n =", nrow(d), "; FDR<0.05 =", sum(d$fdr < 0.05), "\n")
