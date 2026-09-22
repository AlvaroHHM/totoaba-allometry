# Totoaba alometría — Reproducible analysis pipeline

**Manuscript:** *Accumulated degree days reveal thermal plasticity in larval growth allometry of Totoaba macdonaldi*

**Authors:** Hernández-Montiel, Á. H.; Galaviz-Espinoza, M. A.; Gisbert, E.; Alcaraz, C.; López, L. M.; Giffard-Mena, I.

**Institutions:** Faculty of Marine Sciences, Universidad Autónoma de Baja California (UABC), Mexico; IRTA Sant Carles de la Ràpita, Spain.

---

## Overview

This repository contains the complete, reproducible R pipeline used to generate the figures and tables of the manuscript. The pipeline consists of seven scripts (00–06) that run end-to-end in approximately 5 minutes.

---

## Repository structure

    totoaba-allometry/
    ├── data/                         # 5 processed CSV files
    │   ├── medidas.csv               # Morphometric measurements
    │   ├── RMR.csv                   # Respirometry (O2, Weight, Stage)
    │   ├── CDDT0.csv                 # Accumulated degree-days (T0)
    │   ├── flexion_noto.csv          # Notochord flexion
    │   └── medidas_alometria.csv     # Morphometric data for PCA
    ├── scripts/                      # 7 R scripts (00-06)
    │   ├── 00_run_all.R              # Pipeline orchestrator
    │   ├── 01_cdd.R                  # CDD + 10 C rule
    │   ├── 02_breakpoints.R          # Allometric breakpoints
    │   ├── 03_pca_kmeans.R           # PCA + K-means
    │   ├── 04_lmm_rmr.R              # Linear mixed models (RMR)
    │   ├── 05_figures_tables.R       # Notochord flexion + correlation
    │   └── 06_sensitivity_analysis.R # Thermal sensitivity
    ├── output/
    │   ├── figures/                  # 8 generated PNG figures
    │   └── tables/                   # 11 generated CSV tables
    ├── README.md
    ├── LICENSE
    └── .gitignore

---

## Requirements

R >= 4.3

Install required packages:

    install.packages(c("readr","dplyr","tidyr","purrr","ggplot2","cowplot",
      "ggsci","patchwork","segmented","boot","lme4","lmerTest","emmeans",
      "multcomp","performance","effectsize","factoextra","FactoMineR",
      "fpc","cluster","corrplot","rstatix","ggpubr","MASS","robustbase","rlang"))

System: Any OS with R and bash (tested on Linux/macOS; Windows via WSL or Git Bash).

---

## Usage

Clone the repository and run the full pipeline:

    git clone https://github.com/YOUR-USERNAME/totoaba-allometry.git
    cd totoaba-allometry
    Rscript scripts/00_run_all.R

Run individual scripts:

    Rscript -e 'source("scripts/01_cdd.R")'
    Rscript -e 'source("scripts/02_breakpoints.R")'
    Rscript -e 'source("scripts/03_pca_kmeans.R")'
    Rscript -e 'source("scripts/04_lmm_rmr.R")'
    Rscript -e 'source("scripts/05_figures_tables.R")'
    Rscript -e 'source("scripts/06_sensitivity_analysis.R")'

Note: The pipeline must be run from the repository root (where data/ and scripts/ are located).

---

## Outputs

### Figures (8 PNGs)

| File | Script | Description |
|---|---|---|
| fig1a_flexion.png | 05 | Notochord flexion at 26 C |
| fig1b_pca.png | 03 | PCA + K-means (3 clusters) |
| fig1c_q10.png | 04 | Q10 by developmental stage |
| fig1d_cdd_general.png | 01 | Exponential TL vs CDD (T0 = 10 C) |
| figS1_correlation.png | 05 | Spearman correlation matrix |
| figS2_breakpoints.png | 02 | Breakpoints by structure x temperature |
| FigS3_sensitivity_comparison.png | 06 | Comparison of grouping methods |
| FigS4_continuous_sensitivity.png | 06 | Continuous thermal sensitivity |

### Tables (11 CSVs)

| File | Script | Content |
|---|---|---|
| T0_validation.csv | 01 | R2, a, b, SEE, RE per T0 |
| CDD_by_temperature.csv | 01 | Temperature-specific coefficients |
| CDD_general.csv | 01 | General model coefficients |
| Tabla3_breakpoints.csv | 02 | Table 3 of the manuscript |
| pca_scores.csv | 03 | PC1/PC2 scores + clusters |
| pca_contributions.csv | 03 | Variable contributions |
| Tabla2_LMM_primary.csv | 04 | Type III ANOVA with Huber weights |
| Tabla2_Q10_primary.csv | 04 | Slopes + Q10 per stage |
| TablaS2_Q10_morphological.csv | 04 | Complementary morphological analysis |
| Fig1a_flexion_summary.csv | 05 | Flexion summary by DPH |
| Tabla_S2_sensitivity.csv | 06 | Table S2 (3 methods) |

---

## Key results

| Parameter | Value |
|---|---|
| Optimal T0 | 10 C |
| R2 at T0 = 10 C | 0.848 |
| General SL equation | SL = 2.46 * exp(0.0035 * CDD) |
| PCA clusters | 329 / 280 / 51 |
| PC1 variance explained | 96.5 % |
| LMM R2 conditional | 99.3 % |
| Q10 (Preflexion / Flexion / Postflexion) | 0.15 / 0.90 / 0.81 |
| Tm x CDD interaction | negative (p = 0.015) |

---

## Data availability

Included in this repository:

- 5 processed CSV files (data/) ready for the pipeline
- 8 PNG figures (output/figures/)
- 11 CSV tables (output/tables/)

Not included:

- Raw source data (available from the corresponding author upon request)
- Manuscript drafts or personal documentation

---

## Note on figures

The figures in `output/figures/` are the unedited, reproducible outputs 
of the R scripts. The published versions were subsequently polished in 
vector graphics software (e.g., Inkscape, Adobe Illustrator) for layout 
and typography. The underlying data and code are identical.

To regenerate the unedited versions:

    Rscript scripts/00_run_all.R

---

## Citation

If you use this code or data, please cite:

Hernández-Montiel, Á. H., Galaviz-Espinoza, M. A., Gisbert, E., Alcaraz, C., López, L. M., & Giffard-Mena, I. (2026). Accumulated degree days reveal thermal plasticity in larval growth allometry of Totoaba macdonaldi. [Journal name].

---

## License

MIT License - see LICENSE for details.

---

## Contact

Corresponding author: Álvaro H. Hernández-Montiel (a334571@uabc.edu.mx)
