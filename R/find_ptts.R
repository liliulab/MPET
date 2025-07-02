#' Identify Putative Transportation Trios (PTTs) Using Mixed Models (Filtered by Group)
#'
#' Runs linear mixed-effects models for each ADT–mRNA–ICT gene trio in a specific group,
#' returning coefficients and p-values.
#'
#' @param trios Data.frame with columns: adt, mrna, all_gene (ICT gene).
#' @param data Data.frame with expression, covariates, group, and sample_col.
#' @param group Character; value of `group` column to filter by (e.g., "Healthy").
#' @param covariates Character vector of covariate names (e.g., c("Sex", "Age")).
#' @param sample_col Column name for random effect (e.g., "patient_id").
#' @param verbose Logical; print progress.
#'
#' @return Data.frame with prot, mrna, gene, coefficients, p-values.
#' @export
find_ptts <- function(trios, data, group, covariates, sample_col = "patient_id", verbose = TRUE) {
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(tidyr)

  extract_lmer_results <- function(mod, prot, mrna, gene, covariates) {
    x <- summary(mod)$coefficients
    x <- as.data.frame(x)
    x <- x[rownames(x) != "(Intercept)", , drop = FALSE]

    coefs <- x[, "Estimate", drop = FALSE]
    pvals <- x[, "Pr(>|t|)", drop = FALSE]
    coefs$cov <- rownames(coefs)
    pvals$cov <- rownames(pvals)

    coef_wide <- pivot_wider(coefs, names_from = cov, values_from = Estimate)
    pval_wide <- pivot_wider(pvals, names_from = cov, values_from = `Pr(>|t|)`)

    all_covs <- c("mrna", "gene", covariates)
    for (cov in setdiff(all_covs, colnames(coef_wide))) {
      coef_wide[[cov]] <- NA
      pval_wide[[cov]] <- NA
    }

    names(coef_wide) <- paste0(names(coef_wide), ".coef")
    names(pval_wide) <- paste0(names(pval_wide), ".pval")

    cbind(
      data.frame(prot = prot, mrna = mrna, gene = gene),
      coef_wide[, paste0(all_covs, ".coef")],
      pval_wide[, paste0(all_covs, ".pval")]
    )
  }

  # Filter by group
  data <- data %>% filter(group == !!group)

  results_all <- list()
  prots <- unique(trios$adt)

  for (prot in prots) {
    if (verbose) message("Processing: ", prot)

    trios_sub <- trios %>% filter(adt == prot & mrna != all_gene)

    res_list <- foreach(i = seq_len(nrow(trios_sub)), .combine = rbind,
                        .packages = c("lme4", "lmerTest", "dplyr", "tidyr"),
                        .errorhandling = "remove") %dopar% {
      row <- trios_sub[i, ]
      mrna <- as.character(row$mrna)
      gene <- as.character(row$all_gene)

      # Select and rename
      cols_needed <- c(sample_col, prot, mrna, gene, covariates)
      dat_sub <- data[, intersect(cols_needed, colnames(data)), drop = FALSE]
      names(dat_sub)[names(dat_sub) == prot] <- "prot"
      names(dat_sub)[names(dat_sub) == mrna] <- "mrna"
      names(dat_sub)[names(dat_sub) == gene] <- "gene"

      dat_sub <- dat_sub %>%
        mutate(across(c(prot, mrna, gene), as.numeric)) %>%
        na.omit()

      if (nrow(dat_sub) > 1) {
        fml <- as.formula(
          paste0("prot ~ mrna + gene + (1|", sample_col, ")",
                 if (length(covariates) > 0) paste0(" + ", paste(covariates, collapse = " + ")) else "")
        )
        model <- lmer(fml, data = dat_sub, REML = FALSE)
        extract_lmer_results(model, prot, mrna, gene, covariates)
      } else {
        NULL
      }
    }

    if (!is.null(res_list)) {
      results_all[[prot]] <- res_list
    }
  }

  bind_rows(results_all) %>%
    filter(rowSums(is.na(.)) < ncol(.))
}
