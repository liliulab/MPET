# MPET: Multi-level Protein Expression Trajectory Analysis

MPET is an R package for modeling regulatory mechanisms of surface protein expression using multimodal single-cell datasets. It provides tools to identify condition-specific pathways linking mRNA, intermediate regulators (ICTs), and disease conditions across multiple modules.

---

## 📦 Installation

Clone and install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("liliulab/MPET")
```

Ensure all dependencies (e.g., `lme4`, `lmerTest`, `regmed`, `lavaan`, etc.) are installed.

---

## 🧪 Example Workflow

```r
library(MPET)

# Load Seurat object and prepare MPET input
seurat_obj <- readRDS("cov_mrna_adt_sub.RDS")
adt_rna_map <- read.csv("adt_rna_map.csv")
clinical_df <- read.csv("Cov_meta.csv")

mpet_data <- prepare_mpet_data(
  seurat_obj  = seurat_obj,
  adt_rna_map = adt_rna_map,
  clinical_df = clinical_df,
  group_col   = "group",
  sample_col  = "patient_id"
)

# Run Module 1
mod1_res <- run_mpet_module1(
  mpet_data  = mpet_data,
  trios      = readRDS("Trios.RDS"),
  prot       = "ADT_CD69",
  group      = "Severe",
  covariates = c("Age", "Sex")
)

# Run Module 2 (Regularized)
regmed_out <- run_regmed_module2(
  prot         = "ADT_CD69",
  mpet_data    = mpet_data,
  g1           = "Severe",
  g0           = "Healthy",
  g1_mod1_res  = mod1_res,
  g0_mod1_res  = mod1_res
)

# Run Module 2 (Single mediation)
smed_out <- run_smed_module2(
  prot         = "ADT_CD69",
  mpet_data    = mpet_data,
  g1           = "Severe",
  g0           = "Healthy",
  g1_mod1_res  = mod1_res,
  g0_mod1_res  = mod1_res
)
```

---

## 📁 Directory Structure

```
MPET/
├── R/                  # Core function scripts
├── man/                # Documentation
├── DESCRIPTION         # Package metadata
├── NAMESPACE           # Exported functions
├── README.md           # This file
```

---

## 👤 Authors

- **Li Liu** (ASU) – Maintainer | [liliu@asu.edu](mailto:liliu@asu.edu)
- **Rekha Mudappathi** – Contributor

---

## 📄 License

This package is licensed under the MIT License.

---


