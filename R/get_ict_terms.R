#' Identify ICT Genes Based on GO Annotations (Reusable Ontology Logic)
#'
#' Filters GO annotations from BP, MF, and CC using curated keywords and exclusion terms,
#' returning GO-term–gene associations for intracellular transport genes.
#'
#' @param gene_list Vector of gene symbols.
#' @param orgdb An AnnotationHub::OrgDb object.
#' @param bp_xx, mf_xx, cc_xx Lists of GO child terms from GO.db.
#'
#' @return A data.frame with SYMBOL, GOID, TERM, DEFINITION, and ONTOLOGY for ICT genes.
#' @export
get_ict_terms <- function(gene_list, orgdb,
                          bp_xx = as.list(GOBPOFFSPRING),
                          mf_xx = as.list(GOMFOFFSPRING),
                          cc_xx = as.list(GOCCOFFSPRING)) {
  if (is.null(gene_list)) {
    gene_list <- keys(orgdb, keytype = "SYMBOL")
  }
  
  go_all <- AnnotationDbi::select(orgdb, keys = gene_list, columns = "GO", keytype = "SYMBOL")
  colnames(go_all)[colnames(go_all) == "GO"] <- "GOID"

  # Get GO term metadata for all annotations
  all_terms <- AnnotationDbi::select(GO.db, keys = unique(go_all$GOID),
                                     columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                     keytype = "GOID")

  # Define exclusion patterns
  exclusion_patterns <- c(
    "phagocytosis", "granul", "exocytosis of neurotransmitter", "histamine secretion",
    "serotonin secretion", "protein secretion", "exosomal secretion",
    "cortical granule", "zymogen granule"
  )

  exclusion_terms <- unique(all_terms[grepl(paste(exclusion_patterns, collapse = "|"),
                                            all_terms$TERM, ignore.case = TRUE), ])
  exclusion_gos <- list(
    BP = exclusion_terms$GOID[exclusion_terms$ONTOLOGY == "BP"],
    CC = exclusion_terms$GOID[exclusion_terms$ONTOLOGY == "CC"]
  )

  # Helper: extract relevant GO entries per ontology
  extract_ontology_terms <- function(ontology, keywords, offspring_list, exclude_ids = character(0)) {
    go_onto <- go_all[go_all$ONTOLOGY == ontology, ]
    term_anno <- AnnotationDbi::select(GO.db, keys = unique(go_onto$GOID),
                                       columns = c("TERM", "DEFINITION"), keytype = "GOID")
    matched_terms <- grep(keywords, term_anno$TERM, value = TRUE, ignore.case = TRUE)
    seed_ids <- unique(term_anno$GOID[term_anno$TERM %in% matched_terms])

    # Expand with child terms
    all_ids <- unique(na.omit(unlist(c(seed_ids, offspring_list[seed_ids]))))
    filtered_ids <- setdiff(all_ids, exclude_ids)

    go_filtered <- go_onto[go_onto$GOID %in% filtered_ids, ]
    term_filtered <- term_anno[term_anno$GOID %in% filtered_ids, ] %>%
      mutate(ONTOLOGY = ontology) %>% distinct()

    list(go = go_filtered, terms = term_filtered)
  }

  # Run for BP, MF, CC
  bp_res <- extract_ontology_terms("BP",
    "intracellular transport|intracellular protein transport|vesicle-mediated transport",
    bp_xx, exclusion_gos$BP
  )

  mf_res <- extract_ontology_terms("MF",
    "protein carrier activity|folding chaperon[e]?|chaperone binding",
    mf_xx
  )

  cc_res <- extract_ontology_terms("CC",
    "vacuole|autophago|endoplasmic|lysosome|golgi|endosome|vesicle|vesicular|plasma membrane|intermediate compartment|cytoplasmic microtubule|microtubule bundle",
    cc_xx, exclusion_gos$CC
  )

  # Intersect BP and/or MF with CC
  cc_genes <- unique(cc_res$go$SYMBOL)
  bp_genes <- intersect(cc_genes, unique(bp_res$go$SYMBOL))
  mf_genes <- intersect(cc_genes, unique(mf_res$go$SYMBOL))
  all_ict <- unique(c(bp_genes, mf_genes))

  # Filter and return merged result
  all_go <- bind_rows(
    bp_res$go %>% filter(SYMBOL %in% all_ict),
    mf_res$go %>% filter(SYMBOL %in% all_ict),
    cc_res$go %>% filter(SYMBOL %in% all_ict)
  ) %>% distinct()

  all_terms <- bind_rows(bp_res$terms, mf_res$terms, cc_res$terms) %>% distinct()
  go_all_filtered <- merge(all_go, all_terms, by = c("GOID","ONTOLOGY"))

  return(go_all_filtered)
}


