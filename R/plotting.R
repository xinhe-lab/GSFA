#' @title Dotplot of the Effects of Perturbations on Factors
#' @details The size of a dot represents the PIP of association;
#' the color represents the effect size.
#'
#' @param fit An object of class \code{"gsfa_fit"}
#' @param target_names Names of the perturbations, the default uses the
#' row names in \code{fit$posterior_means$beta_pm}
#' @param reorder_targets A subset or rearrangement of \code{target_names}
#' if one wishes to visualize a subset of perturbations or rearrange their order
#' @param reorder_factors A subset or rearrangement of factor indices if
#' one wishes to visualize a subset of factors or rearrange their order
#' @return A \code{ggplot} object.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @export
#'
dotplot_beta_PIP <- function(fit,
                             target_names = NULL,
                             reorder_targets = target_names,
                             reorder_factors = NULL,
                             exclude_offset = TRUE){
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
  if (is.null(target_names)){
    target_names <- colnames(beta_pm)
  } else {
    colnames(beta_pip) <- target_names
    colnames(beta_pm) <- target_names
  }
  if (is.null(reorder_targets)){
    reorder_targets <- target_names
  }

  beta_pip <- beta_pip[, reorder_targets]
  beta_pip_df <- as.data.frame(beta_pip)
  beta_pip_df$Factor <- paste0("Factor ", 1:nrow(beta_pip_df))
  if (!is.null(reorder_factors)){
    beta_pip_df <- beta_pip_df[reorder_factors, ]
  }
  beta_pip_plot_df <- melt(beta_pip_df, value.name = "PIP")

  beta_pm <- beta_pm[, reorder_targets]
  beta_pm_df <- as.data.frame(beta_pm)
  beta_pm_df$Factor <- paste0("Factor ", 1:nrow(beta_pm_df))
  if (!is.null(reorder_factors)){
    beta_pm_df <- beta_pm_df[reorder_factors, ]
  }
  beta_pm_plot_df <- melt(beta_pm_df, id.var = "Factor",
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
          legend.text = element_text(size = 12))
  return(plot_out)
}

#' @title Dotplot of the Total Effects of Perturbations on Genes
#' @details Visualize the effects of perturbations on selected genes.
#' Sizes of the dots represent LFSR bins;
#' colors of the dots represent the summarized effect sizes.
#'
#' @param fit An object of class \code{"gsfa_fit"}
#' @param gene_indices The indices of genes to visualize
#' @param gene_names Names of the genes chosen to visualize, the default uses
#' \code{gene_indices}
#' @param target_names Names of the perturbations, the default uses the column
#' names of \code{fit$lfsr}
#' @param reorder_targets A subset or rearrangement of \code{target_names}
#' if one wishes to visualize a subset of perturbations or rearrange their order
#' @return A \code{ggplot} object.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr rowwise
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#'
#' @export
#'
dotplot_total_effect <- function(fit, gene_indices, gene_names = gene_indices,
                                 target_names = NULL, reorder_targets = NULL,
                                 plot_max_score = NULL){
  # Both inputs should be gene by guide/marker matrices,
  # Dots will be colored by effect size and sized according to LFSR bins.
  require(dplyr)
  require(ggplot2)
  lfsr_binning <- function(lfsr){
    if (lfsr <= 0.05){
      return("0 - 0.05")
    } else if (lfsr <= 0.25){
      return("0.05 - 0.25")
    } else {
      return("> 0.25")
    }
  }
  lfsr_matrix <- fit$lfsr[, -ncol(fit$lfsr)] # remove offset
  if (is.null(target_names)){
    target_names <- colnames(lfsr_matrix)
  } else {
    colnames(lfsr_matrix) <- target_names
  }
  if (is.null(reorder_targets)){
    reorder_targets <- target_names
  }
  selected_lfsr_mat <- lfsr_matrix[gene_indices, reorder_targets]
  rownames(selected_lfsr_mat) <- gene_names

  effect_matrix <- fit$posterior_means$W_pm %*%
    t(fit$posterior_means$beta_pm[-nrow(fit$posterior_means$beta_pm), ])
  colnames(effect_matrix) <- target_names
  selected_effect_mat <- effect_matrix[gene_indices, reorder_targets]
  rownames(selected_effect_mat) <- gene_names

  pip_df <- as.data.frame(selected_lfsr_mat)
  pip_df$gene <- rownames(selected_lfsr_mat)
  pip_plot_df <- melt(pip_df, variable.name = "Perturbation",
                      value.name = "LFSR")

  effect_df <- as.data.frame(selected_effect_mat)
  effect_df$gene <- rownames(selected_effect_mat)
  effect_plot_df <- melt(effect_df, variable.name = "Perturbation",
                         value.name = "Effect_size")

  combined_plot_df <- pip_plot_df %>%
    mutate(Effect_size = effect_plot_df$Effect_size,
           gene = factor(gene, levels = gene_names),
           Perturbation = factor(Perturbation, levels = reorder_targets))
  ## Stratify LSFR values into discrete bins
  combined_plot_df <- combined_plot_df %>%
    rowwise() %>%
    mutate(LFSR_bin = lfsr_binning(LFSR)) %>%
    mutate(LFSR_bin = factor(LFSR_bin, levels = c("> 0.25", "0.05 - 0.25", "0 - 0.05")))

  if (!is.null(plot_max_score)){
    ## Capping effect size values on both ends for more friendly visualization
    ## useful when the data contain only a few extreme values
    plot_min_score <- plot_max_score * (-1)
    plot_by_score <- plot_max_score/2
    combined_plot_df$Effect_size[combined_plot_df$Effect_size > plot_max_score] <- plot_max_score
    combined_plot_df$Effect_size[combined_plot_df$Effect_size < plot_min_score] <- plot_min_score
  }

  plot_out <- ggplot(combined_plot_df) +
    geom_point(aes(x = Perturbation, y = gene,
                   size = LFSR_bin, color = Effect_size)) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)) +
    labs(color = "Summarized effect", size = "LFSR")
  if (!is.null(plot_max_score)){
    plot_out <- plot_out +
      scale_color_gradientn(limits = c(plot_min_score, plot_max_score),
                            colours = c("blue3", "blue", "grey90", "red", "red3"),
                            breaks = seq(plot_min_score, plot_max_score, plot_by_score))
  } else {
    plot_out <- plot_out +
      scale_color_gradient2(low = "blue3", mid = "grey90", high = "red3")
  }
  return(plot_out)
}
