# ============================================================
# 04_lmm_rmr.R — LMM RMR (Tabla 2 + S2 + Fig. 1c)
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(lme4); library(lmerTest); library(emmeans)
  library(multcomp); library(performance); library(effectsize)
 library(MASS) 
library(ggplot2); library(cowplot); library(rlang)
})
set.seed(123)
rmr <- read_csv("data/RMR.csv", show_col_types = FALSE) %>%
  mutate(Temperature = as.numeric(Temperature), DPH = as.integer(DPH),
         Stage = factor(Stage, levels = c("Preflexion","Flexion",
                                          "Postflexion","Juvenile")),
         log_O2 = log(O2), log_Weight = log(Weight))

# === SECCIÓN A (con pesos Huber) ===
dat_A <- rmr %>% filter(DPH %in% c(7, 16, 24)) %>%
  mutate(Stage_E = factor(case_when(
    DPH == 7 ~ "E1", DPH == 16 ~ "E2", DPH == 24 ~ "E3"),
    levels = c("E1","E2","E3")), DPH_f = factor(DPH))
cat(sprintf("\n[Sección A] N obs = %d\n", nrow(dat_A)))

# Modelo estandar
modelo_A_std <- lmer(log_O2 ~ log_Weight + Temperature * Stage_E +
                       (1|id) + (1|DPH_f), data = dat_A,
                     control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 2e5)))

# Pesos Huber
residuos_A <- residuals(modelo_A_std)
pesos_A    <- MASS::psi.huber(residuos_A / sd(residuos_A), k = 1.345)

# Modelo robusto
modelo_A <- update(modelo_A_std, weights = pesos_A)

print(summary(modelo_A))
tabla2 <- anova(modelo_A, ddf = "Satterthwaite", type = 3)
print(tabla2)
# ... (resto igual)
eta2 <- effectsize::eta_squared(modelo_A, partial = TRUE,
                                ci = 0.95, method = "Satterthwaite")
print(eta2)
r2 <- r2_nakagawa(modelo_A)
cat(sprintf("\nR2 marginal = %.1f %% | R2 condicional = %.1f %%\n",
            100 * r2$R2_marginal, 100 * r2$R2_conditional))

pend_A <- emtrends(modelo_A, ~ Stage_E, var = "Temperature")
conf_A <- as.data.frame(confint(pend_A))
conf_A$Q10 <- exp(10 * conf_A$Temperature.trend)
letras_A <- cld(pend_A, Letters = letters, alpha = 0.05)
conf_A$Grupo <- letras_A$.group
print(conf_A)
print(pairs(pend_A, adjust = "tukey"))

write_csv(as.data.frame(tabla2) %>% tibble::rownames_to_column("Source"),
          "output/tables/Tabla2_LMM_primary.csv")
write_csv(conf_A, "output/tables/Tabla2_Q10_primary.csv")

p1c <- ggplot(conf_A, aes(Stage_E, Temperature.trend)) +
  geom_point(size = 4, shape = 21, fill = "white", color = "black", stroke = 1.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
  geom_text(aes(y = upper.CL + 0.005, label = Grupo),
            size = 5, fontface = "bold", vjust = 0) +
  geom_text(aes(y = lower.CL - 0.005,
                label = sprintf("Q[10]==%.2f", Q10)),
            parse = TRUE, size = 3.5, vjust = 1) +
  labs(x = "Larval stage (DPH)", y = "RMR-T slope") +
  theme_cowplot(font_size = 12)
ggsave("output/figures/fig1c_q10.png", p1c, width = 5, height = 4.5, dpi = 300)

dat_B <- rmr %>% filter(!is.na(Stage))
cat(sprintf("\n[Seccion B] N obs = %d\n", nrow(dat_B)))
print(table(dat_B$Stage))

modelo_B_std <- lmer(log_O2 ~ log_Weight + Temperature * Stage +
                       (1|id) + (1|DPH), data = dat_B,
                     control = lmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 2e5)))
residuos_B <- residuals(modelo_B_std)
pesos_B <- MASS::psi.huber(residuos_B / sd(residuos_B), k = 1.345)
modelo_B <- update(modelo_B_std, weights = pesos_B)
cat("\nSingularidad:", lme4::isSingular(modelo_B), "\n")

pend_B <- emtrends(modelo_B, ~ Stage, var = "Temperature")
conf_B <- as.data.frame(confint(pend_B))
conf_B$Q10 <- exp(10 * conf_B$Temperature.trend)
letras_B <- cld(pend_B, Letters = letters, alpha = 0.05)
conf_B$Grupo <- letras_B$.group
print(conf_B)
print(pairs(pend_B, adjust = "tukey"))
write_csv(conf_B, "output/tables/TablaS2_Q10_morphological.csv")

r2_B <- r2_nakagawa(modelo_B)
cat(sprintf("\n[Seccion B] R2 marginal = %.1f %% | R2 condicional = %.1f %%\n",
            100 * r2_B$R2_marginal, 100 * r2_B$R2_conditional))
message("04_lmm_rmr.R completado.")
