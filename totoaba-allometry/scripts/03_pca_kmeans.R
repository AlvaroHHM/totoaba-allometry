# ============================================================
# 03_pca_kmeans.R — PCA + K-means (Fig. 1b)
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(factoextra); library(fpc); library(cowplot)
})
set.seed(123)

dat <- read_csv("data/medidas_alometria.csv", show_col_types = FALSE) %>%
  filter(complete.cases(SL, Eye, Head, Trunk, Notochord))
cat(sprintf("N filas tras filtrado: %d\n", nrow(dat)))

bac <- dat %>% dplyr::select(Notochord, Trunk, Head, Eye, SL)
Tab <- scale(bac)
respca <- princomp(Tab)

var_exp <- round(100 * (respca$sdev^2 / sum(respca$sdev^2))[1:2], 1)
cat(sprintf("Varianza explicada: PC1 = %.1f%%, PC2 = %.1f%%\n",
            var_exp[1], var_exp[2]))

contrib <- sweep(respca$loadings[, 1:2]^2, 2,
                 colSums(respca$loadings[, 1:2]^2), "/") * 100
print(round(contrib, 1))

k1 <- kmeansruns(Tab, k = 3, runs = 1000)
cat("\nTamano de clusters:\n"); print(table(k1$cluster))

palette_aaas <- c("#008B45CC","#3B4992CC","#EE0000CC","#631879CC",
                  "#008280CC","#BB0021CC","#5F559BCC","#A20056CC","#808180CC")

p <- fviz_cluster(k1, data = Tab, geom = "point", shape = 19,
                  ggtheme = theme_cowplot(), palette = palette_aaas,
                  repel = TRUE, main = "Morphometric PCA",
                  ellipse.type = "t")
ggsave("output/figures/fig1b_pca.png", p, width = 6, height = 5, dpi = 300)

scores <- as_tibble(respca$scores[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(DPH, Temperature, SL)) %>%
  mutate(Cluster = factor(k1$cluster))
write_csv(scores, "output/tables/pca_scores.csv")
write_csv(as_tibble(contrib, rownames = "variable"),
          "output/tables/pca_contributions.csv")
message("03_pca_kmeans.R completado.")
