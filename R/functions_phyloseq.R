# Plotting functions & palettes -------------------------------------------

# Define pallete list
PALETTES <- list(
  study_group_name = c(
    "CD" = "#dc267f",
    "UC" = "#ffb000",
    "HC"  = "#648fff"
  ),
  IBD_group_name = c(
    "IBD" = "#785ef0",
    "HC"  = "#648fff"
  ),
  continent = c(
    "Asia"          = "#dc267f",
    "Europe"        = "#785ef0",
    "North America" = "#648fff"
  )
)

# ALL
plot_beta_diversity <- function(physeq_object, comparator, color_map = NULL, palettes = PALETTES, auto_color_if_missing = TRUE, show_pairwise = TRUE) {
  # Step 1: Prune taxa
  physeq_object <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_object) > 0, physeq_object)
  
  # Step 2: Ordination with Bray-Curtis distance
  bray <- phyloseq::ordinate(physeq_object, method = "PCoA", distance = "bray")
  
  # Step 3: Extract variance explained for the first two axes
  variance_axis1 <- bray$values$Relative_eig[1] * 100
  variance_axis2 <- bray$values$Relative_eig[2] * 100
  
  axis1_label <- sprintf("PCo1 (%.1f%%)", variance_axis1)
  axis2_label <- sprintf("PCo2 (%.1f%%)", variance_axis2)
  
  # Step 4: PERMANOVA (adonis2)
  arc_dist_matrix <- phyloseq::distance(physeq_object, method = "bray")
  
  adonis_result <- vegan::adonis2(
    as.formula(paste("arc_dist_matrix ~", comparator)),
    data = as(phyloseq::sample_data(physeq_object), "data.frame")
  )
  
  p_value <- adonis_result$`Pr(>F)`[1]
  r2_value <- adonis_result$R2[1]
  title_text <- sprintf("PERMANOVA p = %.3f", p_value)
  
  print(adonis_result)
  
  # Optional: Pairwise ADONIS
  pairwise_result <- NULL
  if (show_pairwise) {
    pairwise_result <- pairwise.adonis(
      x = arc_dist_matrix,
      factors = phyloseq::sample_data(physeq_object)[[comparator]]
    )
    print(pairwise_result)
  }
  
  # Step 5: Base phyloseq ordination plot
  p.s.nmds <- phyloseq::plot_ordination(
    physeq_object,
    bray,
    color = comparator
  )
  
  # Step 6: Build plot
  beta <- p.s.nmds +
    ggplot2::geom_point(size = 2) +
    ggplot2::stat_ellipse(ggplot2::aes(group = !!rlang::sym(comparator)), linetype = 1) +
    cowplot::theme_cowplot(20) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = axis1_label, y = axis2_label, title = title_text) +
    ggplot2::theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.title = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1, fill = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7),
        color = "black",
        linewidth = 0.5
      )
    )
  
  # Step 7: Apply colour mapping
  if (!is.null(color_map)) {
    beta <- beta + ggplot2::scale_color_manual(values = color_map)
  } else if (!is.null(palettes) && comparator %in% names(palettes)) {
    beta <- beta + ggplot2::scale_color_manual(values = palettes[[comparator]])
  } else if (auto_color_if_missing) {
    # auto-generate a named palette from the factor levels present
    levs <- levels(as.factor(phyloseq::sample_data(physeq_object)[[comparator]]))
    auto_cols <- setNames(scales::hue_pal()(length(levs)), levs)
    beta <- beta + ggplot2::scale_color_manual(values = auto_cols)
  }
  
  # Step 8: Add marginal plots
  beta2 <- ggExtra::ggMarginal(beta, type = "boxplot", groupFill = TRUE, size = 10)
  
  # Return
  return(list(
    plot = beta2,
    p_value = p_value,
    r2 = r2_value,
    adonis = adonis_result,
    pairwise = pairwise_result
  ))
}

# SGN
plot_beta_diversity_faceted <- function(physeq_object, comparator, facet_variable, color_map = NULL, palettes = PALETTES, auto_color_if_missing = TRUE, show_pairwise = TRUE, ncol = 3, facet_scales = "fixed") {
  # Step 1: Prune taxa
  physeq_object <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq_object) > 0, physeq_object)
  
  # Step 2: Ordination with Bray-Curtis distance
  bray <- phyloseq::ordinate(physeq_object, method = "PCoA", distance = "bray")
  
  # Step 3: Variance explained
  variance_axis1 <- bray$values$Relative_eig[1] * 100
  variance_axis2 <- bray$values$Relative_eig[2] * 100
  
  axis1_label <- sprintf("PCo1 (%.1f%%)", variance_axis1)
  axis2_label <- sprintf("PCo2 (%.1f%%)", variance_axis2)
  
  # Step 4: PERMANOVA (adonis2)
  arc_dist_matrix <- phyloseq::distance(physeq_object, method = "bray")
  adonis_result <- vegan::adonis2(
    as.formula(paste("arc_dist_matrix ~", comparator)),
    data = as(phyloseq::sample_data(physeq_object), "data.frame")
  )
  print(adonis_result)
  
  # Optional: Pairwise PERMANOVA
  pairwise_result <- NULL
  if (show_pairwise) {
    pairwise_result <- pairwise.adonis(
      arc_dist_matrix,
      as.factor(phyloseq::sample_data(physeq_object)[[comparator]])
    )
    print(pairwise_result)
  }
  
  p_value <- adonis_result$`Pr(>F)`[1]
  title_text <- sprintf("PERMANOVA p = %.3f", p_value)
  
  # Step 5: Ordination plot
  p.s.nmds <- phyloseq::plot_ordination(physeq_object, bray, color = comparator)
  
  # Step 6: Build plot
  facetted_plot <- p.s.nmds +
    ggplot2::geom_point(size = 2) +
    ggplot2::stat_ellipse(ggplot2::aes(group = !!rlang::sym(comparator)), linetype = 1) +
    cowplot::theme_cowplot(20) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = axis1_label,
      y = axis2_label,
      title = title_text
    ) +
    ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_variable)),
                        scales = facet_scales, ncol = ncol) +
    ggplot2::theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7),
        color = "black",
        linewidth = 0.5
      ),
      plot.title = ggplot2::element_text(
        size = 16,
        margin = ggplot2::margin(t = 10, b = 5)
      ),
      panel.border = ggplot2::element_rect(fill = NA, colour = "white", linewidth = 2),
      strip.text = ggplot2::element_text(size = 20),
      strip.background = ggplot2::element_rect(color = "black", fill = "white", linewidth = 1),
      panel.spacing = grid::unit(1, "lines")
    )
  
  # Step 7: Apply colour mapping
  if (!is.null(color_map)) {
    facetted_plot <- facetted_plot + ggplot2::scale_color_manual(values = color_map)
  } else if (!is.null(palettes) && comparator %in% names(palettes)) {
    facetted_plot <- facetted_plot + ggplot2::scale_color_manual(values = palettes[[comparator]])
  } else if (auto_color_if_missing) {
    levs <- levels(as.factor(phyloseq::sample_data(physeq_object)[[comparator]]))
    auto_cols <- setNames(scales::hue_pal()(length(levs)), levs)
    facetted_plot <- facetted_plot + ggplot2::scale_color_manual(values = auto_cols)
  }
  
  return(list(
    plot = facetted_plot,
    p_value = p_value,
    r2 = adonis_result$R2[1],
    adonis = adonis_result,
    pairwise = pairwise_result
  ))
}

adonis_with_subsets <- function(physeq_object, comparator, p_adjust_method = "fdr") {
  # Define the subgroups for subsetting
  groups <- list(
    "Healthy Control" = subset_samples(physeq_object, IBD_group_name == "HC"),
    "CD" = subset_samples(physeq_object, study_group_name == "CD"),
    "UC" = subset_samples(physeq_object, study_group_name == "UC")
  )
  
  # Initialize a list to store results
  adonis_results_list <- list()
  
  # Loop through each group and perform ADONIS
  for (group_name in names(groups)) {
    subset_physeq <- groups[[group_name]]  # Get the subset for the current group
    
    # Calculate Bray-Curtis distance
    dist_matrix <- phyloseq::distance(subset_physeq, method = "bray")
    
    # Perform ADONIS test
    adonis_result <- vegan::adonis2(
      as.formula(paste("dist_matrix ~", comparator)),
      data = as(sample_data(subset_physeq), "data.frame")
    )
    
    # Store the result
    adonis_results_list[[group_name]] <- data.frame(
      Group = group_name,
      Pr_F = adonis_result$`Pr(>F)`[1],
      R2 = adonis_result$R2[1],
      Note = "ADONIS performed successfully"
    )
  }
  
  # Combine all results into a single data frame
  adonis_results_df <- do.call(rbind, adonis_results_list)
  
  # Adjust p-values
  valid_p_values <- !is.na(adonis_results_df$Pr_F)  # Identify valid p-values
  adonis_results_df$Adjusted_P <- NA  # Initialize column for adjusted p-values
  adonis_results_df$Adjusted_P[valid_p_values] <- p.adjust(
    adonis_results_df$Pr_F[valid_p_values], 
    method = p_adjust_method
  )
  
  return(adonis_results_df)
}