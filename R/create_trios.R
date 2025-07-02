#' Create Surface Protein–Coding Gene–ICT Gene Trios
#'
#' Combines all ADT–coding gene pairs with all ICT genes to produce trios
#' for downstream mixed-effects modeling.
#'
#' @param adt_rna_map Table of ADTs mapped to coding genes.
#' @param ict_genes All ICT genes (vector).
#' @param prefix_adt Logical; if TRUE, add "ADT_" prefix to adt column. Default = TRUE.
#'
#' @return Returns dataframe of ADT–coding gene–ICT gene trios.
#' @export
create_trios <- function(adt_rna_map,
                         ict_genes,
                         prefix_adt = TRUE) {
  if (prefix_adt) {
    adt_rna_map$adt <- paste0("ADT_", adt_rna_map$adt)
  }

  # Repeat ADT–mRNA map to match length of ICT gene list
  adt_rna_gene_map <- adt_rna_map[rep(seq_len(nrow(adt_rna_map)), each = length(ict_genes)), ]
  rownames(adt_rna_gene_map) <- seq_len(nrow(adt_rna_gene_map))

  # Repeat ICT genes to match ADT–mRNA pairs
  adt_rna_gene_map$ict_gene <- rep(ict_genes, times = nrow(adt_rna_map))

  message("Trio table created")

  return(adt_rna_gene_map)
}
