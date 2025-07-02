#' Run Single Mediation Analysis for a Surface Protein using SEM
#'
#' This function implements an alternative version of Module 2 of MPET using structural equation modeling (SEM).
#' For a single surface protein, it fits individual mediation models for each ICT gene using `lavaan`,
#' based on expression data across two groups.
#'
#' @param prot Character. Name of the surface protein (e.g., "ADT_CD86").
#' @param mpet_data Data frame containing RNA, ADT, ICT expression, and metadata. Must include a `group` column.
#' @param g1 Character. Comparison group 1 (e.g., "Severe").
#' @param g0 Character. Comparison group 2 (e.g., "Healthy").
#' @param g1_mod1_res Data frame of GLM results for group `g1`. Must contain `prot`, `gene`, `gene.Padj`, `gene.coef`.
#' @param g0_mod1_res Data frame of GLM results for group `g0`. Must contain `prot`, `gene`, `gene.Padj`, `gene.coef`.
#' @param pval_cutoff Numeric. FDR-adjusted p-value threshold for selecting significant ICTs (default = 0.05).
#' @param thr_cutoff Numeric. Coefficient threshold for selecting significant ICTs (default = 0.05).
#'
#' @return A data frame of mediation results with alpha, beta, delta, ab estimates, p-values, and mediation classification.
#' @export
run_single_module2 <- function(prot,
                            mpet_data,
                            g1,
                            g0,
                            g1_mod1_res,
                            g0_mod1_res,
                            pval_cutoff = 0.05,
                            thr_cutoff = 0.05) {
  message("\nRunning SEM Module 2 for protein: ", prot)

  # Subset mpet_data and encode binary outcome
  g.xy.df <- mpet_data |>
    dplyr::filter(group %in% c(g1, g0)) |>
    dplyr::mutate(D = ifelse(group == g1, 1, 0)) |>
    dplyr::select(all_of(c("D", prot, "group")), everything())

  # Get significant ICTs from both groups
  x_sig <- g1_mod1_res[g1_mod1_res$prot == prot & g1_mod1_res$gene.Padj < pval_cutoff & abs(g1_mod1_res$gene.coef) > thr_cutoff, ]
  y_sig <- g0_mod1_res[g0_mod1_res$prot == prot & g0_mod1_res$gene.Padj < pval_cutoff & abs(g0_mod1_res$gene.coef) > thr_cutoff, ]
  icts <- unique(na.omit(c(x_sig$gene, y_sig$gene)))
  message("Number of significant ICTs to evaluate via SEM: ", length(icts))

  if (length(icts) == 0) return(NULL)

  # SEM model template
  sem_model <- '
    M ~ a * G
    D ~ b * M
    D ~ d * G
    ab := a * b
    total := d + ab
  '

  # Apply SEM to each ICT
  sem_results <- purrr::map_dfr(icts, function(cur_ict) {
    cur_data <- g.xy.df |>
      dplyr::select(all_of(c("D", prot, cur_ict))) |>
      dplyr::rename(G = !!cur_ict, M = !!prot)

    cur_data$D <- as.ordered(cur_data$D)

    fit <- tryCatch({
      lavaan::sem(sem_model, data = cur_data, estimator = "WLSMV", ordered = "D")
    }, error = function(e) return(NULL))

    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(NULL)

    est <- lavaan::parameterEstimates(fit, standardized = TRUE)

    extract_val <- function(lhs, rhs, op) {
      row <- est[est$lhs == lhs & est$rhs == rhs & est$op == op, ]
      if (nrow(row) > 0) return(row) else return(NULL)
    }

    alpha <- extract_val("M", "G", "~")
    beta  <- extract_val("D", "M", "~")
    delta <- extract_val("D", "G", "~")
    ab    <- est[est$label == "ab", ]

    if (any(lengths(list(alpha, beta, delta, ab)) == 0)) return(NULL)

    mediation <- dplyr::case_when(
      ab$pvalue < 0.05 & delta$pvalue >= 0.05 ~ "full",
      ab$pvalue < 0.05 & delta$pvalue < 0.05  ~ "partial",
      ab$pvalue >= 0.05 & delta$pvalue < 0.05 ~ "null",
      TRUE ~ "none"
    )

    return(data.frame(
      ict = cur_ict,
      alpha = alpha$est,
      alpha_pval = alpha$pvalue,
      beta  = beta$est,
      beta_pval  = beta$pvalue,
      delta = delta$est,
      delta_pval = delta$pvalue,
      ab = ab$est,
      ab_pval = ab$pvalue,
      mediation = mediation,
      stringsAsFactors = FALSE
    ))
  })

  return(sem_results)
}
