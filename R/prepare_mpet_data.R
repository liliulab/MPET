#' Prepare MPET Input Data from a Seurat Object
#'
#' Generates RNA expression, ADT expression, metadata, and GLM input for MPET
#' from a single-cell Seurat object that is already normalized, scaled, and subset to a single cell type.
#'
#' @param seurat_obj A preprocessed Seurat object with RNA and ADT assays (subset to a single cell type).
#' @param adt_rna_map A data.frame mapping ADTs to coding genes (columns: adt, gene).
#' @param clinical_df A data.frame of patient-level metadata.
#' @param group_col Character; name of the group/condition column (e.g., "group").
#'        The column name must match the corresponding column in the Seurat metadata.
#' @param sample_col Character; name of the sample or patient ID column (e.g., "patient_id").
#'        The column name must match the corresponding column in the Seurat metadata.
#' @param prefix_adt Logical; if TRUE, add "ADT_" prefix to adt column. Default = TRUE.
#'
#' @return Dataframe of RNA expression, ADT expression, metadata, GLM input, and gene list.
#' @export
prepare_mpet_data <- function(seurat_obj, adt_rna_map, clinical_df,
                              group_col = "group",
                              sample_col = "patient_id",
                              prefix_adt = TRUE) {

  # ---- ADT Expression ----
  DefaultAssay(seurat_obj) <- "IADT"
  adt_data <- as.data.frame(GetAssayData(seurat_obj, layer = "data"))
  adt_data <- adt_data[rownames(adt_data) %in% adt_rna_map$adt, ]
  adt_expr <- t(adt_data)
  if (prefix_adt) {
   colnames(adt_expr) <- paste0("ADT_", colnames(adt_expr))
  }
  adt_expr <- as.data.frame(adt_expr)

  # ---- RNA Expression ----
  DefaultAssay(seurat_obj) <- "RNA"
  rna_data <- as.data.frame(GetAssayData(seurat_obj, layer = "data"))
  rna_data$gene_ens <- rownames(rna_data)
  rna_data <- rna_data %>% separate(gene_ens, c("gene", "ens.id"), sep = ":")
  rna_data <- rna_data[!is.na(rna_data$ens.id), ]
  rna_data <- rna_data[!duplicated(rna_data$gene), ]
  rownames(rna_data) <- rna_data$gene
  rna_data <- rna_data[, !(names(rna_data) %in% c("gene", "ens.id"))]
  rna_expr <- as.data.frame(t(rna_data))

  # ---- Metadata ----
  meta <- seurat_obj@meta.data[, c(group_col, sample_col)]
  meta$cell_id <- rownames(meta)
  meta_merged <- merge(meta, clinical_df, by = c(sample_col, group_col))

  # ---- GLM Input ----
  stopifnot(identical(rownames(rna_expr), rownames(adt_expr)))
  combined <- cbind(rna_expr, adt_expr)
  combined$ID <- gsub("_.*", "", rownames(combined))
  combined$cell_id <- rownames(combined)
  combined_1 <- merge(combined, meta_merged, by = "cell_id")
  names(combined_1)[names(combined_1) == group_col] <- "group"
  names(combined_1)[names(combined_1) == sample_col] <- "sample"
  return(combined_1)

}
