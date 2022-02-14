dotplot_beta_PIP <- function(fit,
                             target_names,
                             reorder_targets = target_names,
                             reorder_factors = NULL,
                             exclude_offset = TRUE){
  # Dots will be colored by effect size and sized by PIP value.
  if (!inherits(fit, "gsfa_fit")){
    stop("Input argument \"fit\" should be an object of class ",
         "\"gsfa_fit\", such as an output of fit_gsfa_multivar().")
  }
  beta_pip <- t(fit$posterior_means$Gamma_pm) # factor by target association PIP matrix
  beta_pm <- t(fit$posterior_means$beta_pm) # factor by target association effect matrix
  if (exclude_offset){
    beta_pip <- beta_pip[, -ncol(beta_pip)]
    beta_pm <- beta_pm[, -ncol(beta_pm)]
  }
  colnames(beta_pip) <- target_names
  beta_pip <- beta_pip[, reorder_targets]
  beta_pip_df <- as.data.frame(beta_pip)
  beta_pip_df$Factor <- paste0("Factor ", 1:nrow(beta_pip_df))
  if (!is.null(reorder_factors)){
    beta_pip_df <- beta_pip_df[reorder_factors, ]
  }
  beta_pip_plot_df <- reshape2::melt(beta_pip_df, value.name = "PIP")

  colnames(beta_pm) <- target_names
  beta_pm <- beta_pm[, reorder_targets]
  beta_pm_df <- as.data.frame(beta_pm)
  beta_pm_df$Factor <- paste0("Factor ", 1:nrow(beta_pm_df))
  if (!is.null(reorder_factors)){
    beta_pm_df <- beta_pm_df[reorder_factors, ]
  }
  beta_pm_plot_df <- reshape2::melt(beta_pm_df, id.var = "Factor",
                                    variable.name = "Perturbation",
                                    value.name = "Effect size")
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(PIP = beta_pip_plot_df$PIP,
           Perturbation = factor(Perturbation, levels = reorder_targets))
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(Factor = factor(Factor, levels = beta_pip_df$Factor[1:nrow(beta_pip_df)]))

  plot_out <- ggplot(beta_pm_plot_df) +
    geom_point(aes(x = Factor, y = Perturbation,
                   size = PIP, color = `Effect size`)) +
    scale_color_gradient2(low = "purple3", mid = "grey90", high = "darkorange1") +
    # color scale can be changed externally
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13, hjust = 1),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 12))
  return(plot_out)
}

dotplot_total_effect <- function(fit, gene_indices, gene_names = gene_indices,
                                 reorder_targets = NULL,
                                 plot_max_score = 0.2, plot_min_score = -0.2, plot_by_score = 0.1){
  # Both inputs should be gene by guide/marker matrices,
  # Dots will be colored by effect size and sized by 1-LFSR value.
  lfsr_binning <- function(lfsr){
    if (lfsr <= 0.05){
      return("0 - 0.05")
    } else if (lfsr <= 0.25){
      return("0.05 - 0.25")
    } else {
      return("> 0.25")
    }
  }

  lfsr_matrix <- fit$lfsr
  if (is.null(reorder_targets)){
    reorder_targets <- colnames(lfsr_matrix)
  }
  selected_lfsr_mat <- lfsr_matrix[gene_indices, reorder_targets]
  rownames(selected_lfsr_mat) <- gene_names

  effect_matrix <- fit$posterior_means$W_pm %*%
    t(fit$posterior_means$beta_pm[-nrow(fit$posterior_means$beta_pm), ])
  colnames(effect_matrix) <- colnames(lfsr_matrix)
  selected_effect_mat <- effect_matrix[gene_indices, reorder_targets]
  rownames(selected_effect_mat) <- gene_names

  pip_df <- as.data.frame(selected_lfsr_mat)
  pip_df$gene <- rownames(selected_lfsr_mat)
  pip_plot_df <- reshape2::melt(pip_df, variable.name = "Perturbation",
                                value.name = "LFSR")

  effect_df <- as.data.frame(selected_effect_mat)
  effect_df$gene <- rownames(selected_effect_mat)
  effect_plot_df <- reshape2::melt(effect_df, variable.name = "Perturbation",
                                   value.name = "Effect_size")

  combined_plot_df <- effect_plot_df %>%
    mutate(LFSR = pip_plot_df$LFSR,
           gene = factor(gene, levels = rownames(selected_effect_mat)),
           Perturbation = factor(Perturbation, levels = reorder_targets))
  combined_plot_df$Effect_size[combined_plot_df$Effect_size > max_score] <- max_score
  combined_plot_df$Effect_size[combined_plot_df$Effect_size < min_score] <- min_score
  combined_plot_df <- combined_plot_df %>%
    rowwise() %>%
    mutate(LFSR_bin = lfsr_binning(LFSR)) %>%
    mutate(LFSR_bin = factor(LFSR_bin, levels = c("> 0.25", "0.05 - 0.25", "0 - 0.05")))
  plot_out <- ggplot(combined_plot_df) +
    geom_point(aes(x = Perturbation, y = gene,
                   size = LFSR_bin, color = Effect_size)) +
    scale_color_gradientn(limits = c(min_score, max_score),
                          colours = c("blue3", "blue", "grey90", "red", "red3"),
                          breaks = seq(min_score, max_score, by_score)) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)) +
    labs(color = "Summarized effect", size = "LFSR")
  return(plot_out)
}
