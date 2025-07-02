#' Run Regularized Mediation Analysis for a Single Protein
#'
#' This function implements Module 2 of MPET: regularized mediation analysis using `regmed`. 
#' It models mediation effects for a surface protein across two comparison groups using ICTs 
#' filtered by significance and effect size in Module 1.
#'
#' @param prot Character. Name of the surface protein (e.g., "ADT_CD16").
#' @param mpet_data Data frame containing RNA, ADT, ICT expression and metadata. Must include a column named `group`.
#' @param g1 Character. Comparison group 1 (e.g., "Severe").
#' @param g0 Character. Comparison group 2 (e.g., "Healthy").
#' @param g1_mod1_res Data frame of Module 1 GLM results for group `g1`. Must include columns `prot`, `gene`, `gene.Padj`, `gene.coef`.
#' @param g0_mod1_res Data frame of Module 1 GLM results for group `g0`. Must include columns `prot`, `gene`, `gene.Padj`, `gene.coef`.
#' @param pval_cutoff Numeric. FDR-adjusted p-value threshold for selecting significant ICTs (default = 0.05).
#' @param thr_cutoff Numeric. Absolute coefficient threshold for selecting significant ICTs (default = 0.05).
#' @param lambda_grid Numeric vector. Grid of lambda values used for regularization in `mvregmed.grid` (default = seq(0.03, 0.005, -0.005)).
#' @param out_dir (Unused) Character. Output directory. Retained for compatibility.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{regmed_results}{Data frame containing mediation statistics and mediation category for each ICT.}
#'   \item{regmed_model}{Best-fit multivariate regularized mediation model (`mvregmed`) object.}
#' }
#' @export

run_regmed_module2 <- function(prot,
                                mpet_data,
                                g1,
                                g0,
                                g1_mod1_res,
                                g0_mod1_res,
                                pval_cutoff = 0.05,
                                thr_cutoff = 0.05,
                                lambda_grid = seq(0.03, 0.005, -0.005),
                                out_dir) {

  # Subset for groups
  g.xy.df <- mpet_data |>
      filter(group %in% c(g1,g0)) |>
      mutate(D = ifelse(group == g1, 1, 0))

  message("\nProcessing: ", prot)

  x_sig <- g1_mod1_res[g1_mod1_res$gene.Padj < pval_cutoff & abs(g1_mod1_res$gene.coef) > thr_cutoff, ]
  y_sig <- g0_mod1_res[g0_mod1_res$gene.Padj < pval_cutoff & abs(g0_mod1_res$gene.coef) > thr_cutoff, ]

  x.icts <- unique(x_sig$gene)
  y.icts <- unique(y_sig$gene)
  xy.icts <- unique(na.omit(c(x.icts, y.icts)))
  message("Number of initial glm icts: ", length(xy.icts))

  if (length(xy.icts) > 0) {
    tryCatch({
      med <- g.xy.df[[prot]]
      x <- g.xy.df[, xy.icts, drop = FALSE]
      y <- g.xy.df[["D"]]

      fit.grid <- mvregmed.grid(x, med, y, lambda_grid)
      mvfit.best <- mvregmed.grid.bestfit(fit.grid)
      mvfit.edges <- mvregmed.edges(mvfit.best)

      beta.res <- mvfit.edges$beta.df[, c("mediator", "beta")]
      alpha.res <- mvfit.edges$alpha.df[, c("mediator", "x", "alpha")]
      delta.res <- mvfit.edges$delta.df[, c("x", "delta")]

      if (!is.null(alpha.res) | !is.null(beta.res) | !is.null(delta.res)) {
        if (is.null(delta.res)) delta.res <- data.frame(x = character(), delta = numeric())
        if (is.null(beta.res)) beta.res <- data.frame(mediator = character(), beta = numeric())
        if (is.null(alpha.res)) alpha.res <- data.frame(mediator = character(), x = character(), alpha = numeric())

        regmed.res_1 <- merge(alpha.res, delta.res, by = "x", all = TRUE)
        regmed.res_1$mediator[is.na(regmed.res_1$mediator)] <- "nil"
        all_res <- merge(regmed.res_1, beta.res, by = "mediator", all = TRUE)
        all_res$mediator <- prot

        all_res <- all_res |>
          mutate(mediation = case_when(
            is.na(delta) & !is.na(beta) & !is.na(alpha) ~ "full",
            !is.na(delta) & !is.na(beta) & !is.na(alpha) ~ "partial",
            !is.na(delta) & is.na(beta) & is.na(alpha) ~ "nil",
            TRUE ~ "none"
          ))

        all_res_keep <- all_res |> filter(mediation %in% c("partial", "full", "nil"))
        message("Number of retained icts after mediation: ", nrow(all_res_keep))
        return(list(regmed_results = all_res_keep, regmed_model = mvfit.best))      
      }
    }, error = function(e) {
      warning(paste("Error processing:", prot, "-", e$message))
       return(NULL)
    })
  }
}
