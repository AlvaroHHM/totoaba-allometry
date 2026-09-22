# ============================================================
# 02_breakpoints.R — Breakpoints alometricos + Tabla 3
# Excepciones: 26 Head (BP 1) y 28 Eye (BP 2) usan npsi=2
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(segmented); library(boot); library(ggplot2); library(cowplot)
})
set.seed(123)

datos <- read_csv("data/medidas.csv", show_col_types = FALSE) %>%
  mutate(Temperature = as.numeric(Temperature), DPH = as.integer(DPH)) %>%
  rename(tm = Temperature, Tail = Notochord)

temperaturas <- c(20, 24, 26, 28)
estructuras  <- c("Eye", "Head", "Trunk", "Tail")

# Tabla de excepciones: cuando el manuscrito uso un BP del modelo npsi=2
excepciones_bp <- tibble::tribble(
  ~Temperature, ~Structure, ~bp_idx,
  26,           "Head",      1L,   # 26 Head -> primer BP de npsi=2 (3.25)
  28,           "Eye",       2L    # 28 Eye  -> segundo BP de npsi=2 (4.30)
)

usa_excepcion <- function(temp, est) {
  any(excepciones_bp$Temperature == temp & excepciones_bp$Structure == est)
}

bp_idx_excepcion <- function(temp, est) {
  idx <- excepciones_bp$bp_idx[excepciones_bp$Temperature == temp &
                                excepciones_bp$Structure == est]
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

seleccionar_bp <- function(temp, est, datos) {
  d <- datos %>% filter(tm == temp) %>%
    mutate(log_x = log(SL), log_y = log(.data[[est]])) %>%
    dplyr::select(SL, log_x, log_y, DPH) %>% drop_na()
  if (nrow(d) < 10) return(NULL)

  lm0 <- lm(log_y ~ log_x, data = d)

  seg1 <- tryCatch(segmented(lm0, seg.Z = ~log_x, npsi = 1),
                   error = function(e) NULL)
  if (is.null(seg1)) return(NULL)

  seg2 <- tryCatch(segmented(lm0, seg.Z = ~log_x, npsi = 2),
                   error = function(e) NULL)

  p_val <- if (!is.null(seg2)) {
    tryCatch(anova(seg1, seg2)[["Pr(>F)"]][2], error = function(e) NA)
  } else NA

  if (usa_excepcion(temp, est) && !is.null(seg2)) {
    idx    <- bp_idx_excepcion(temp, est)
    psi_mm <- exp(seg2$psi[, "Est."])[idx]
    coefs  <- coef(seg2)
    s1     <- coefs["log_x"]
    delta1 <- coefs["U1.log_x"]
    delta2 <- coefs["U2.log_x"]
    slopes <- c(s1, s1 + delta1, s1 + delta1 + delta2)
    return(tibble(Temperature = temp, Structure = est, N = 1L,
                  p_anova = p_val, BP1_mm = psi_mm, BP2_mm = NA_real_,
                  slope_S1 = slopes[1], slope_S2 = slopes[2],
                  slope_S3 = NA_real_))
  }

  psi_mm <- exp(seg1$psi[, "Est."])
  coefs  <- coef(seg1)
  s1     <- coefs["log_x"]
  delta  <- coefs["U1.log_x"]
  slopes <- c(s1, s1 + delta)

  tibble(Temperature = temp, Structure = est, N = 1L, p_anova = p_val,
         BP1_mm = psi_mm[1], BP2_mm = NA_real_,
         slope_S1 = slopes[1], slope_S2 = slopes[2], slope_S3 = NA_real_)
}

bootstrap_bp <- function(temp, est, datos, n_bp, B = 200) {
  d <- datos %>% filter(tm == temp) %>%
    mutate(log_x = log(SL), log_y = log(.data[[est]])) %>%
    dplyr::select(SL, log_x, log_y) %>% drop_na()
  if (nrow(d) < 10) return(tibble())

  npsi_uso <- if (usa_excepcion(temp, est)) 2L else 1L
  idx_bp   <- if (usa_excepcion(temp, est)) bp_idx_excepcion(temp, est) else 1L

  f <- function(data, idx_boot) {
    dd <- data[idx_boot, ]
    m0 <- lm(log_y ~ log_x, data = dd)
    s  <- tryCatch(segmented(m0, seg.Z = ~log_x, npsi = npsi_uso),
                   error = function(e) NULL)
    if (is.null(s)) return(NA_real_)
    return(exp(s$psi[, "Est."])[idx_bp])
  }

  br <- tryCatch(boot(d, f, R = B), error = function(e) NULL)
  if (is.null(br)) return(tibble())

  ci1 <- quantile(br$t[, 1], c(0.025, 0.975), na.rm = TRUE)
  tibble(Temperature = temp, Structure = est,
         BP1_CI_low = ci1[1], BP1_CI_high = ci1[2],
         stable_BP1 = diff(ci1) < median(br$t[, 1], na.rm = TRUE) * 0.5)
}

res_anova <- map2_dfr(
  rep(temperaturas, each = length(estructuras)),
  rep(estructuras, times = length(temperaturas)),
  ~ seleccionar_bp(.x, .y, datos))

res_boot <- map2_dfr(
  rep(temperaturas, each = length(estructuras)),
  rep(estructuras, times = length(temperaturas)),
  ~ {
    info <- res_anova %>% filter(Temperature == .x, Structure == .y)
    if (nrow(info) == 0 || is.na(info$N[1])) return(tibble())
    bootstrap_bp(.x, .y, datos, info$N[1], B = 200)
  })

res <- res_anova %>% left_join(res_boot, by = c("Temperature", "Structure"))

decisiones <- tribble(
  ~Temperature, ~Structure, ~N_final,
  20,"Eye",0L, 20,"Head",1L, 20,"Tail",1L, 20,"Trunk",1L,
  24,"Eye",0L, 24,"Head",0L, 24,"Tail",1L, 24,"Trunk",1L,
  26,"Eye",1L, 26,"Head",1L, 26,"Tail",1L, 26,"Trunk",0L,
  28,"Eye",1L, 28,"Head",0L, 28,"Tail",0L, 28,"Trunk",1L)

res <- res %>% left_join(decisiones, by = c("Temperature", "Structure")) %>%
  mutate(N_final = replace_na(N_final, 0L))

eq_by_temp <- tribble(
  ~Temperature, ~a_sl, ~b_sl,
  20, 2.61, 0.0047, 24, 2.33, 0.0037,
  26, 1.67, 0.0049, 28, 2.08, 0.0038)
a_gen <- 2.46; b_gen <- 0.0035
cdd_sp  <- function(SL, a, b) (log(SL) - log(a)) / b
cdd_gen <- function(SL) (log(SL) - log(a_gen)) / b_gen
asignar_etapa <- function(cdd) case_when(
  cdd < 36 ~ "Egg", cdd < 144 ~ "Preflexion",
  cdd < 216 ~ "Flexion", cdd < 360 ~ "Postflexion", TRUE ~ "Juvenile")

tabla3 <- res %>%
  filter(N_final >= 1, !is.na(BP1_mm)) %>%
  left_join(eq_by_temp, by = "Temperature") %>%
  mutate(
    CDD_sp  = cdd_sp(BP1_mm, a_sl, b_sl),
    CDD_gen = cdd_gen(BP1_mm),
    Stage_sp  = asignar_etapa(CDD_sp),
    Stage_gen = asignar_etapa(CDD_gen),
    dCDD    = CDD_sp - CDD_gen,
    Dev_sp  = (CDD_sp / 360) * 100,
    Dev_gen = (CDD_gen / 360) * 100,
    CI_95   = sprintf("(%.2f-%.2f)", BP1_CI_low, BP1_CI_high)) %>%
  dplyr::select(Temperature, Structure, N_final, BP1_mm, CI_95,
         CDD_sp, Stage_sp, CDD_gen, Stage_gen, dCDD, Dev_sp, Dev_gen) %>%
  arrange(Temperature, Structure)

write_csv(tabla3, "output/tables/Tabla3_breakpoints.csv")
print(tabla3)

pS2 <- res %>% filter(!is.na(BP1_mm)) %>%
  ggplot(aes(factor(Temperature), BP1_mm, color = Structure)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = BP1_CI_low, ymax = BP1_CI_high), width = 0.15) +
  facet_wrap(~Structure, scales = "free_y") +
  labs(x = "Temperature (C)", y = "Breakpoint (mm SL)") +
  theme_classic(base_size = 11) + theme(legend.position = "none")

ggsave("output/figures/figS2_breakpoints.png", pS2,
       width = 7, height = 5, dpi = 300)
message("02_breakpoints.R completado.")
