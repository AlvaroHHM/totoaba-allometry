# ============================================================
# 00_run_all.R — Orquestador del pipeline
# ============================================================
options(stringsAsFactors = FALSE)
scripts <- c("scripts/01_cdd.R",
             "scripts/02_breakpoints.R",
             "scripts/03_pca_kmeans.R",
             "scripts/04_lmm_rmr.R",
             "scripts/05_figures_tables.R",
             "scripts/06_sensitivity_analysis.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
for (s in scripts) {
  if (!file.exists(s)) { warning(sprintf("No existe: %s", s)); next }
  message(sprintf("\n=== %s ===", s))
  t0 <- Sys.time()
  source(s, echo = FALSE)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("=== %.1f s ===\n", dt))
}
message("Pipeline completo.")
