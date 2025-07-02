#' Run MPET Module 1 for a Single Protein in a Specific Condition 
#'
#' Fits mixed-effects models for trios associated with one protein using RNA, ADT, and ICT gene expression.
#' Includes optional covariates in the model. Ensures one result row per trio, even if coefficients are dropped.
#'
#' @param mpet_data A data frame with RNA, ADT, and metadata. Must include a "sample"  column.
#' @param trios A data frame of trios for a single protein (columns: adt, mrna, ict_gene).
#' @param prot Character. Name of the surface protein (e.g., "ADT_CD16").
#' @param group Character; Sample group to subset.
#' @param padj_method Method for p-value adjustment (default = "BH").
#' @param covariates Optional character vector of covariate column names to include in the model.
#' @return A data frame with estimated coefficients, p-values, and adjusted p-values for each trio.
#' @export
run_mpet_module1 <- function(mpet_data, 
                             trios,
                             prot,
                             group,
                             padj_method = "BH",
                             covariates = NULL) {


  extract_glm_results <- function(mod, prot, mrna, gene) {
    x <- summary(mod)$coefficients
    x <- as.data.frame(x)
    x <- x[rownames(x) != "(Intercept)", c("Estimate", "Pr(>|t|)")]
    x$cov <- rownames(x)
    x <- x[x$cov %in% c("mrna", "gene"), ]

    # Always return a row with NA-filled values if necessary
    out <- tibble(
      prot       = prot,
      mrna       = mrna,
      gene       = gene,
      mrna.coef  = NA_real_,
      gene.coef  = NA_real_,
      mrna.pval  = NA_real_,
      gene.pval  = NA_real_
    )

    if ("mrna" %in% x$cov) {
      out$mrna.coef <- x$Estimate[x$cov == "mrna"]
      out$mrna.pval <- x$`Pr(>|t|)`[x$cov == "mrna"]
    }
    if ("gene" %in% x$cov) {
      out$gene.coef <- x$Estimate[x$cov == "gene"]
      out$gene.pval <- x$`Pr(>|t|)`[x$cov == "gene"]
    }

    return(out)
  }

  mpet_g_data <- mpet_data[which(mpet_data$group == group), ]
  trios_p <- trios[which(trios$adt == prot), ]
  message("Running MPET Module 1 for protein: ", prot, " | ", nrow(trios_p), " trios")


  results <- foreach(p = 1:nrow(trios_p), .combine = rbind) %do% {
    trio <- trios_p[p, ]
    prot_col <- as.character(trio$adt)
    mrna_col <- as.character(trio$mrna)
    gene_col <- as.character(trio$ict_gene)

    cols <- c("sample", prot_col, mrna_col, gene_col, covariates)
    df <- mpet_g_data[, intersect(cols, colnames(mpet_g_data))]
    names(df)[names(df) == prot_col] <- "prot"
    names(df)[names(df) == mrna_col] <- "mrna"
    names(df)[names(df) == gene_col] <- "gene"
    df <- na.omit(df)

    if (nrow(df) > 1) {
      cov_string <- if (!is.null(covariates) && length(covariates) > 0) paste(covariates, collapse = " + ") else ""
      formula_str <- paste("prot ~ mrna + gene",
                     if (cov_string != "") paste("+", cov_string),
                     "+ (1 | sample)")
      mod_formula <- as.formula(formula_str)

      mod <- tryCatch(
        suppressWarnings(lmer(mod_formula, data = df, REML = FALSE)),
        error = function(e) NULL
      )
      if (!is.null(mod)) {
        return(extract_glm_results(mod, prot_col, mrna_col, gene_col))
      }
    }

    # Return NA row if model fails or insufficient data
    return(data.frame(
      prot       = prot_col,
      mrna       = mrna_col,
      gene       = gene_col,
      mrna.coef  = NA_real_,
      gene.coef  = NA_real_,
      mrna.pval  = NA_real_,
      gene.pval  = NA_real_
    ))
  }

  # Adjust p-values
  if (!is.null(results) && nrow(results) > 0) {
    results <- as.data.frame(results)
    results$mrna.Padj <- p.adjust(results$mrna.pval, method = padj_method)
    results$gene.Padj <- p.adjust(results$gene.pval, method = padj_method)
  }

  return(results)
}
