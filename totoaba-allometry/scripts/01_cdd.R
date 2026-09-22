# ============================================================
# 01_cdd.R — CDD, regla 10 C, Fig. 1d
# ============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot)
})

cdd <- read_csv("data/CDDT0.csv", show_col_types = FALSE) %>%
  mutate(Temperature = as.numeric(Temperature),
         T0 = as.numeric(T0), DPH = as.integer(DPH))

fit_by_t0 <- cdd %>%
  group_by(T0) %>%
  group_modify(~ {
    m <- nls(TL_cm ~ a * exp(b * CDD), data = .x,
             start = list(a = 0.2, b = 0.005))
    tibble(a = coef(m)["a"], b = coef(m)["b"],
           R2 = 1 - sum(residuals(m)^2) /
                    sum((.x$TL_cm - mean(.x$TL_cm))^2),
           SEE = sqrt(sum(residuals(m)^2) / (nrow(.x) - 2)),
           RE  = (SEE / mean(.x$TL_cm)) * 100, n = nrow(.x))
  }) %>% ungroup()
print(fit_by_t0)

best <- fit_by_t0 %>% slice_max(R2, n = 1)
cat(sprintf("\n Mejor T0 = %.0f C | R2 = %.3f | a = %.3f | b = %.4f\n",
            best$T0, best$R2, best$a, best$b))

dat10 <- cdd %>% filter(T0 == 10)
m10 <- nls(TL_cm ~ a * exp(b * CDD), data = dat10,
           start = list(a = 0.2, b = 0.005))
a10 <- coef(m10)["a"]; b10 <- coef(m10)["b"]

med <- read_csv("data/medidas.csv", show_col_types = FALSE)
fit_by_temp <- med %>%
  group_by(Temperature) %>%
  group_modify(~ {
    m <- nls(SL ~ a * exp(b * CDD), data = .x, start = list(a = 2, b = 0.004))
    tibble(a = coef(m)["a"], b = coef(m)["b"],
           R2 = 1 - sum(residuals(m)^2) / sum((.x$SL - mean(.x$SL))^2))
  }) %>% ungroup()
print(fit_by_temp)

x_seq <- seq(0, max(dat10$CDD), length.out = 200)
pred <- data.frame(CDD = x_seq, TL_cm = a10 * exp(b10 * x_seq))
p1d <- ggplot(dat10, aes(CDD, TL_cm)) +
  geom_point(aes(color = factor(Temperature)), size = 2.5, alpha = 0.8) +
  geom_line(data = pred, color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = c(36, 144, 216, 360),
             linetype = "solid", color = "gray40") +
  labs(x = "CDD10 (C*d)", y = "Total length (cm)",
       color = "T (C)") +
  theme_cowplot(font_size = 12)
ggsave("output/figures/fig1d_cdd_general.png", p1d,
       width = 5.5, height = 4.5, dpi = 300)

write_csv(fit_by_t0,  "output/tables/T0_validation.csv")
write_csv(fit_by_temp, "output/tables/CDD_by_temperature.csv")
write_csv(tibble(a = a10, b = b10, R2 = best$R2,
                 SEE = best$SEE, RE = best$RE),
          "output/tables/CDD_general.csv")
message("01_cdd.R completado.")
