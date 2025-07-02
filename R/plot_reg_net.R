#' Plot Regulatory Graph from Mediation Results
#'
#' Visualizes a directed graph showing ICT-mRNA-protein-disease relationships 
#' for a specific surface protein based on mediation analysis results.
#'
#' @param prot Character. Surface protein name (e.g., "ADT_HLA-DR").
#' @param mediation_df Data frame of mediation results (output of Module 2).
#' @param mrna_df Data frame containing significant mRNA links.
#' @param glm_res_x Data frame of GLM results for group x.
#' @param glm_res_yz Data frame of GLM results for groups y and z.
#' @param mrna_logfc Data frame of mRNA DEG results (logFC).
#' @param output_path Character. Path to save the output PNG figure.
#' @param offset Numeric. Spacing offset to separate overlapping nodes.
#' @param width,height,res Plot dimensions and resolution for PNG.
#'
#' @return Saves a PNG plot of the graph.
#' @export
plot_regulatory_graph <- function(prot, mediation_df, mrna_df, glm_res_x, glm_res_yz,
                                  mrna_logfc, output_path, offset = 10,
                                  width = 10, height = 10, res = 600) {
  require(igraph)

  mediation_df$mediator <- ifelse(mediation_df$mediator == "med.1", prot, mediation_df$mediator)

  mrnas <- unique(mrna_df$mrna)
  xyz_glm <- rbind(glm_res_x, glm_res_yz)
  p_glm <- xyz_glm[xyz_glm$prot == prot & xyz_glm$mrna %in% mrnas & xyz_glm$mrna.Padj < 0.05, ]
  mrna_means <- p_glm |>
    dplyr::group_by(mrna) |>
    dplyr::summarise(mean_mrna_coef = mean(mrna.coef, na.rm = TRUE)) |>
    dplyr::distinct()

  # Generate links and nodes
  make_links <- function(mediator_type, from_col, to_col) {
    df <- mediation_df[mediation_df$mediation == mediator_type, ]
    data.frame(g = df[[from_col]], l = df[[to_col]])
  }
  f_link <- make_links("full", "x", "mediator")
  p_link <- make_links("partial", "x", "mediator")
  d_link <- mediation_df[!is.na(mediation_df$delta), ] |> dplyr::select(g = x, l = "D")
  m_link <- data.frame(g = prot, l = "D")

  mrna_links <- data.frame(
    g = rep(mrnas, each = 2),
    l = rep(c(prot, "D"), times = length(mrnas))
  )

  all_links <- dplyr::distinct(rbind(f_link, p_link, d_link, m_link, mrna_links))

  # Nodes
  all_nodes <- unique(c(all_links$g, all_links$l))
  node_types <- dplyr::case_when(
    all_nodes == "D" ~ "D",
    all_nodes == prot ~ "M",
    all_nodes %in% mediation_df$x[mediation_df$mediation == "full"] ~ "Full",
    all_nodes %in% mediation_df$x[mediation_df$mediation == "partial"] ~ "Partial",
    all_nodes %in% mediation_df$x[mediation_df$mediation == "nil"] ~ "nil",
    all_nodes %in% mrnas ~ "mRNA",
    TRUE ~ "other"
  )
  node_df <- data.frame(n = all_nodes, t = node_types)

  net <- igraph::graph_from_data_frame(d = all_links, vertices = node_df, directed = TRUE)

  # Layout with offset
  generate_custom_layout <- function(net, offset) {
    layout <- layout_with_fr(net)
    layout_df <- as.data.frame(layout)
    colnames(layout_df) <- c("x", "y")
    layout_df$node <- V(net)$name
    layout_df$type <- V(net)$t

    m_coords <- layout_df[layout_df$type == "M", c("x", "y")]
    d_coords <- layout_df[layout_df$type == "D", c("x", "y")]

    if (nrow(m_coords) == 0 || nrow(d_coords) == 0) return(as.matrix(layout_df[, c("x", "y")]))

    dir_vec <- as.numeric(d_coords - m_coords)
    norm_vec <- dir_vec / sqrt(sum(dir_vec^2))
    perp_vec <- c(-norm_vec[2], norm_vec[1])

    adjust_nodes <- function(df, type) {
      df[df$type == type, c("x", "y")] <- df[df$type == type, c("x", "y")] +
        offset * matrix(rep(perp_vec, sum(df$type == type)), ncol = 2, byrow = TRUE)
      df
    }

    layout_df <- adjust_nodes(layout_df, "Partial")
    layout_df <- adjust_nodes(layout_df, "mRNA")
    return(as.matrix(layout_df[, c("x", "y")]))
  }

  layout <- generate_custom_layout(net, offset)

  # Edge values
  edge_values <- mapply(function(src, tgt) {
    if (src %in% mediation_df$x & tgt == prot) {
      match_row <- mediation_df[mediation_df$x == src & mediation_df$mediator == prot, ]
      return(ifelse(nrow(match_row) > 0, match_row$alpha[1], NA))
    } else if (src == prot & tgt == "D") {
      return(mediation_df$beta[1])
    } else if (tgt == "D" & src %in% mediation_df$x) {
      match_row <- mediation_df[mediation_df$x == src, ]
      return(ifelse(nrow(match_row) > 0, match_row$delta[1], NA))
    } else if (src %in% mrna_means$mrna & tgt == prot) {
      return(mrna_means$mean_mrna_coef[mrna_means$mrna == src])
    } else if (src %in% mrna_means$mrna & tgt == "D") {
      match_row <- mrna_logfc[mrna_logfc$mrna == src & mrna_logfc$prot == prot, ]
      return(ifelse(nrow(match_row) > 0, match_row$mrna.avg_log2FC[1], NA))
    } else return(NA)
  }, all_links$g, all_links$l)

  all_links$edge_val <- edge_values
  all_links$edge_col <- ifelse(edge_values > 0, "firebrick",
                          ifelse(edge_values < 0, "dodgerblue4", "gray60"))
  E(net)$color <- all_links$edge_col

  # Node aesthetics
  colrs <- c("D" = "tomato1", "M" = "magenta2", "Full" = "limegreen", 
             "Partial" = "gold2", "nil" = "dodgerblue1", "mRNA" = "gray30")
  shape_map <- c("D" = "square", "M" = "rectangle", "Full" = "circle",
                 "Partial" = "circle", "nil" = "circle", "mRNA" = "circle")

  V(net)$color <- colrs[V(net)$t]
  V(net)$shape <- shape_map[V(net)$t]
  V(net)$size <- ifelse(V(net)$shape == "rectangle", 25, 10)

  # Save plot
  png(filename = output_path, width = width, height = height, units = "in", res = res)
  plot(net,
       layout = layout,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       vertex.label.font = 2,
       vertex.size = 20,
       edge.arrow.size = 0.4,
       edge.arrow.width = 1.5,
       margin = 0.1)
  dev.off()
}
