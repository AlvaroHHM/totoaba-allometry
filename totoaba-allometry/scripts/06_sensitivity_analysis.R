# ============================================================
# 06_sensitivity_analysis.R — Tabla S2 + FigS3 + FigS4
# Tres metodos: DPH-based, CDD-based, continuo Tm × CDD
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(lme4); library(lmerTest); library(emmeans)
  library(multcomp); library(MASS); library(robustbase)
  library(performance); library(effectsize)
  library(ggplot2); library(cowplot); library(patchwork)
})
set.seed(123)

# --- 1. Cargar datos ---
medidas <- read_csv("data/RMR.csv", show_col_types = FALSE) %>%
  mutate(
    id       = factor(id),
    Tm       = as.numeric(Temperature),
    DPH      = as.integer(DPH),
    O2       = as.numeric(O2),
    Peso     = as.numeric(Weight),
    etapa    = factor(etapa_E, levels = c("E1","E2","E3")),
    nombre_etapa = factor(Stage,
                          levels = c("Preflexion","Flexion",
                                     "Postflexion","Juvenile")),
    Individuo = factor(paste(Temperature, id, sep = "_"))
  )

cat("N total:", nrow(medidas), "\n")
cat("Etapas DPH-based (etapa_E):\n"); print(table(medidas$etapa, medidas$Tm))
cat("\nEtapas CDD-based (Stage):\n"); print(table(medidas$nombre_etapa, medidas$Tm))

# --- 2. Funcion generica para modelo robusto ---
fit_robust_model <- function(data, group_var, group_label = NULL,
                             temp_var = "Tm", response = "O2",
                             weight = "Peso",
                             random = c("(1|id)", "(1|DPH)")) {

  formula_str <- paste0("log(", response, ") ~ log(", weight, ") + ",
                        temp_var, " * ", group_var, " + ",
                        paste(random, collapse = " + "))
  formula_obj   <- as.formula(formula_str)
  group_formula <- as.formula(paste("~", group_var))

  cat("\n========================================\n")
  cat("Modelo con grupo:", group_var, "\n")
  cat("========================================\n")

  modelo_std <- lmer(formula_obj, data = data, REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 100000)))

  residuos <- residuals(modelo_std)
  pesos    <- MASS::psi.huber(residuos / sd(residuos), k = 1.345)

  modelo_rob <- lmer(formula_obj, data = data, weights = pesos, REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 100000)))

  pendientes   <- emtrends(modelo_rob, group_formula, var = temp_var)
  confint_pend <- confint(pendientes, level = 0.95)
  df_pend      <- as.data.frame(confint_pend)

  colnames(df_pend)[1] <- "Group"
  colnames(df_pend)[grep("trend", colnames(df_pend))] <- "m"
  df_pend$Q10 <- exp(10 * df_pend$m)

  comparaciones <- pairs(pendientes, adjust = "tukey")
  letras        <- cld(pendientes, Letters = letters, alpha = 0.05)
  df_pend$Tukey <- letras$.group

  anova_res <- anova(modelo_rob, ddf = "Satterthwaite", type = 3)
  r2        <- r2_nakagawa(modelo_rob)

  df_pend$Method <- ifelse(is.null(group_label), group_var, group_label)
  df_pend <- df_pend[, c("Method", "Group", "m", "SE", "df",
                          "lower.CL", "upper.CL", "Q10", "Tukey")]

  list(table = df_pend, model = modelo_rob, anova = anova_res,
       r2 = r2, comparisons = comparaciones, emtrends = pendientes)
}

# --- 3. Analisis 1: DPH-based (E1, E2, E3) ---
result_dph <- fit_robust_model(data = medidas,
                               group_var = "etapa",
                               group_label = "DPH-based")

cat("\n=== RESULTADOS: DPH-BASED ===\n"); print(result_dph$table)
cat("\n=== ANOVA ===\n"); print(result_dph$anova)
cat("\n=== R2 ===\n");    print(result_dph$r2)

# --- 4. Analisis 2: CDD-based (etapas morfologicas) ---
cat("\n=== MUESTRAS POR ETAPA Y TEMPERATURA (CDD-based) ===\n")
tabla_muestras <- table(medidas$nombre_etapa, medidas$Tm)
print(tabla_muestras)

etapas_validas <- names(which(rowSums(tabla_muestras > 0) >= 2))
cat("\nEtapas validas:", paste(etapas_validas, collapse = ", "), "\n")

medidas_cdd <- medidas %>%
  filter(nombre_etapa %in% etapas_validas) %>%
  droplevels()

result_cdd <- fit_robust_model(data = medidas_cdd,
                               group_var = "nombre_etapa",
                               group_label = "CDD-based")

cat("\n=== RESULTADOS: CDD-BASED ===\n"); print(result_cdd$table)
cat("\n=== ANOVA ===\n"); print(result_cdd$anova)
cat("\n=== R2 ===\n");    print(result_cdd$r2)

# --- 5. Analisis 3: Modelo continuo Tm × CDD ---
medidas$CDD_std <- as.numeric(scale(medidas$CDD_en_etapa))

modelo_cont_std <- lmer(
  log(O2) ~ log(Peso) + Tm * CDD_std + (1 | id) + (1 | DPH),
  data = medidas, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

residuos_cont <- residuals(modelo_cont_std)
pesos_cont    <- MASS::psi.huber(residuos_cont / sd(residuos_cont), k = 1.345)

modelo_cont_rob <- lmer(
  log(O2) ~ log(Peso) + Tm * CDD_std + (1 | id) + (1 | DPH),
  data = medidas, weights = pesos_cont, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

cat("\n=== MODELO CONTINUO: Tm x CDD ===\n")
print(summary(modelo_cont_rob)$coefficients)

anova_cont <- anova(modelo_cont_rob, ddf = "Satterthwaite", type = 3)
cat("\n=== ANOVA CONTINUO ===\n"); print(anova_cont)

r2_cont <- r2_nakagawa(modelo_cont_rob)
cat("\n=== R2 ===\n"); print(r2_cont)

efectos_cont <- effectsize::eta_squared(modelo_cont_rob, partial = TRUE,
                                         ci = 0.95, method = "Satterthwaite")
cat("\n=== TAMANOS DEL EFECTO ===\n"); print(efectos_cont)

coefs   <- fixef(modelo_cont_rob)
b_Tm    <- coefs["Tm"]
b_inter <- coefs["Tm:CDD_std"]

cdd_mean <- attr(scale(medidas$CDD_en_etapa), "scaled:center")
cdd_sd   <- attr(scale(medidas$CDD_en_etapa), "scaled:scale")

cdd_valores     <- c(36, 144, 216, 360)
cdd_std_valores <- (cdd_valores - cdd_mean) / cdd_sd
slopes_por_cdd  <- b_Tm + b_inter * cdd_std_valores
q10_por_cdd     <- exp(10 * slopes_por_cdd)

tabla_continua <- data.frame(
  CDD   = cdd_valores,
  Stage = c("Preflexion", "Flexion", "Postflexion", "Juvenile"),
  m     = round(slopes_por_cdd, 4),
  Q10   = round(q10_por_cdd, 3)
)
cat("\n=== SLOPES DE Tm EN FUNCION DE CDD ===\n"); print(tabla_continua)

# --- 6. Tabla S2 unificada ---
add_n_by_temp <- function(df, data, group_var) {
  n_table <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(n_20 = sum(Tm == 20), n_24 = sum(Tm == 24),
              n_26 = sum(Tm == 26), n_28 = sum(Tm == 28),
              n_total = n(), .groups = "drop")
  colnames(n_table)[1] <- "Group"
  merge(df, n_table, by = "Group", all.x = TRUE)
}

tabla_dph <- add_n_by_temp(result_dph$table, medidas, "etapa")
tabla_cdd <- add_n_by_temp(result_cdd$table, medidas_cdd, "nombre_etapa")

tabla_cont <- data.frame(
  Method = "Continuous (Tm x CDD)",
  Group  = tabla_continua$Stage,
  m      = tabla_continua$m, SE = NA, df = NA,
  lower.CL = NA, upper.CL = NA, Q10 = tabla_continua$Q10, Tukey = NA,
  n_20 = c(3,0,6,0), n_24 = c(3,2,3,0),
  n_26 = c(3,2,1,3), n_28 = c(3,2,1,3),
  n_total = c(12,6,11,6)
)

cols_keep <- c("Method","Group","m","SE","Q10","Tukey",
               "n_20","n_24","n_26","n_28","n_total")

tabla_s2 <- rbind(
  tabla_dph[, cols_keep],
  tabla_cdd[, cols_keep],
  tabla_cont[, cols_keep]
)
tabla_s2$m   <- round(tabla_s2$m, 3)
tabla_s2$SE  <- round(tabla_s2$SE, 3)
tabla_s2$Q10 <- round(tabla_s2$Q10, 3)

cat("\n\n============================================\n")
cat("TABLA S2: SENSIBILIDAD TERMICA POR METODO\n")
cat("============================================\n")
print(tabla_s2)
write_csv(tabla_s2, "output/tables/Tabla_S2_sensitivity.csv")

# --- 7. Figuras ---
plot_dph <- result_dph$table
plot_dph$Method <- "DPH-based (E1-E3)"
plot_dph$Group  <- factor(plot_dph$Group, levels = c("E1","E2","E3"))

plot_cdd <- result_cdd$table
plot_cdd$Method <- "CDD-based (biological stage)"

plot_all <- rbind(
  plot_dph[, c("Method","Group","m","lower.CL","upper.CL","Q10","Tukey")],
  plot_cdd[, c("Method","Group","m","lower.CL","upper.CL","Q10","Tukey")]
)

p_combined <- ggplot(plot_all, aes(x = Group, y = m)) +
  geom_point(size = 4, shape = 21, fill = "black", color = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.8) +
  geom_text(aes(label = Tukey), vjust = -1.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Method, scales = "free_x") +
  labs(x = "Developmental stage",
       y = "Thermal sensitivity coefficient (m)",
       title = "Comparison of thermal sensitivity across grouping strategies") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))
ggsave("output/figures/FigS3_sensitivity_comparison.png", p_combined,
       width = 12, height = 5, dpi = 300)

p_continuous <- ggplot(tabla_continua, aes(x = CDD, y = m)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 4, shape = 21, fill = "steelblue", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  annotate("rect", xmin = 36,  xmax = 144, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "green") +
  annotate("rect", xmin = 144, xmax = 216, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "yellow") +
  annotate("rect", xmin = 216, xmax = 360, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "orange") +
  labs(x = "Accumulated degree-days (CDD)",
       y = "Thermal sensitivity coefficient (m)",
       title = "Ontogenetic change in thermal sensitivity (continuous model)") +
  theme_cowplot()
ggsave("output/figures/FigS4_continuous_sensitivity.png", p_continuous,
       width = 8, height = 5, dpi = 300)

# --- 8. Resumen ---
p_inter <- anova_cont["Tm:CDD_std", "Pr(>F)"]

cat("\n\n==========================================================\n")
cat("  RESUMEN: SENSITIVITY ANALYSIS\n")
cat("==========================================================\n")
cat(sprintf("  Metodo 1 (DPH-based):   %d etapas, %d larvas\n",
            nrow(result_dph$table), nrow(medidas)))
cat(sprintf("  Metodo 2 (CDD-based):   %d etapas, %d larvas\n",
            nrow(result_cdd$table), nrow(medidas_cdd)))
cat(sprintf("  Metodo 3 (Continuo):    1 modelo, %d larvas\n", nrow(medidas)))
cat("----------------------------------------------------------\n")
cat(sprintf("  Interaccion Tm x CDD:   %s (p = %.4f)\n",
            ifelse(b_inter > 0, "positiva", "negativa"), p_inter))
cat(sprintf("  R2 condicional:         %.1f%%\n", r2_cont$R2_conditional * 100))
cat("==========================================================\n")
message("06_sensitivity_analysis.R completado.")
