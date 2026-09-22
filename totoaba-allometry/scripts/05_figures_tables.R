# ============================================================
# 05_figures_tables.R — Fig. 1a + Fig. S1
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(cowplot); library(corrplot)
})

if (file.exists("data/flexion_noto.csv")) {
  flex <- read_csv("data/flexion_noto.csv", show_col_types = FALSE)
  cat("Columnas flexion_noto:", paste(names(flex), collapse = ", "), "\n")

  flex_26 <- flex %>%
    filter(Diet == "N", as.numeric(Temperature) == 26) %>%
    mutate(Flexion_pct = (Flexion * 100) / 10, DPH = as.integer(DPH)) %>%
    group_by(DPH) %>%
    summarise(mean = mean(Flexion_pct, na.rm = TRUE),
              sd = sd(Flexion_pct, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    arrange(DPH)
  print(flex_26)
  write_csv(flex_26, "output/tables/Fig1a_flexion_summary.csv")

  p1a <- ggplot(flex_26, aes(factor(DPH), mean)) +
    geom_col(fill = "#3B4992CC", alpha = 0.8) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    labs(x = "DPH", y = "Notochord flexion (%)") +
    theme_cowplot(font_size = 12)
  ggsave("output/figures/fig1a_flexion.png", p1a,
         width = 5, height = 4, dpi = 300)
  message("Fig. 1a generada")
} else {
  warning("flexion_noto.csv no encontrado")
}

if (file.exists("data/medidas_alometria.csv")) {
  ma <- read_csv("data/medidas_alometria.csv", show_col_types = FALSE)
  cat("Columnas medidas_alometria:", paste(names(ma), collapse = ", "), "\n")

  vars <- intersect(c("SL","Eye","Head","Trunk","Notochord"), names(ma))
  if (length(vars) == 5) {
    Tab <- scale(ma[, vars])
    cor_matrix <- cor(Tab, method = "spearman")
    p_values <- cor.mtest(Tab, conf.level = 0.95)

    png("output/figures/figS1_correlation.png",
        width = 1800, height = 1800, res = 300)
    corrplot(cor_matrix, order = "hclust", hclust.method = "complete",
             addrect = 2, method = "circle", type = "lower",
             tl.col = "black", addCoef.col = TRUE, addCoefasPercent = TRUE,
             p.mat = p_values$p, insig = "blank", diag = FALSE,
             addgrid = FALSE, number.cex = 0.7)
    dev.off()
    message("Fig. S1 generada")
  } else {
    warning("medidas_alometria.csv no tiene las 5 variables necesarias")
  }
} else {
  warning("medidas_alometria.csv no encontrado")
}

message("05_figures_tables.R completado.")
