#' @title Deviance Residual Transformation on Gene Count Matrix
#' @description Performs deviance residual transformation (proposed in Townes et al., 2019)
#' and converts given raw gene count matrix into a continuous matrix with quantities
#' analogous to z-scores and approximately follow a normal distribution.
#' @param count_mat A sample by gene numeric matrix that stores gene count values;
#' @return A sample by gene numeric matrix that stores transformed deviance residual values.
#' @export
deviance_residual_transform <- function(count_mat){
  resid_mat <- matrix(nrow = nrow(count_mat), ncol = ncol(count_mat))
  lib_size <- rowSums(count_mat)
  total_size <- sum(lib_size)
  for (j in 1:ncol(count_mat)){
    relative_pi <- sum(count_mat[, j]) / total_size
    for (i in 1:nrow(count_mat)){
      mu_ij <- lib_size[i] * relative_pi
      if (count_mat[i, j] == 0){
        cross_prod <-
          2 * lib_size[i] * log(lib_size[i] / (lib_size[i] - mu_ij))
      } else {
        cross_prod <- 2 * count_mat[i, j] * log(count_mat[i, j] / mu_ij) +
          2 * (lib_size[i] - count_mat[i, j]) * log((lib_size[i] - count_mat[i, j]) / (lib_size[i] - mu_ij))
      }
      resid_mat[i, j] <- sign(count_mat[i, j] - mu_ij) * sqrt(max(cross_prod, 0))
    }
  }
  return(resid_mat)
}

#' @title Pearson Residual Transformation on Gene Count Matrix
#' @description Performs Pearson residual transformation (proposed in Townes et al., 2019)
#' and converts given raw gene count matrix into a continuous matrix with quantities
#' analogous to z-scores and approximately follow a normal distribution.
#' @param count_mat A sample by gene numeric matrix that stores gene count values;
#' @return A sample by gene numeric matrix that stores transformed Pearson residual values.
#' @export
pearson_residual_transform <- function(count_mat){
  resid_mat <- matrix(nrow = nrow(count_mat), ncol = ncol(count_mat))
  lib_size <- rowSums(count_mat)
  total_size <- sum(lib_size)
  for (j in 1:ncol(count_mat)){
    relative_pi <- sum(count_mat[, j]) / total_size
    for (i in 1:nrow(count_mat)){
      mu_ij <- lib_size[i] * relative_pi
      resid_mat[i, j] <- (count_mat[i, j] - mu_ij) / sqrt(mu_ij - mu_ij^2/lib_size[i])
    }
  }
  return(resid_mat)
}

#' @title Feature Selection Based on the Deviance Statistics of Genes
#' @description Select the top genes ranked by deviance statistics.
#' @param dev_resid_mat A sample by gene numeric matrix that stores deviance residuals;
#' @param num_top_genes Number of top genes to keep;
#' @return A numeric vector that stores the indices of top \code{num_top_genes}
#' genes ranked by deviance statistics.
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#'
#' @export
select_top_devres_genes <- function(dev_resid_mat, num_top_genes = 6000){
  dev_stats_df <- data.frame(stats = colSums(dev_resid_mat^2) / nrow(dev_resid_mat),
                             index = 1:ncol(dev_resid_mat))
  keep_index <- dev_stats_df %>% arrange(-stats) %>% slice(1:num_top_genes) %>% pull(index)
  return(keep_index)
  # reduced matrix: dev_resid_mat[, keep_index]
}

#' @title Covariate Correction for Gene Expression Matrix
#' @description Remove the effects of sample covariates from gene expression via linear regression.
#' @param expression_mat A sample by gene numeric matrix that stores gene expression values;
#' @param covariate_df A sample by covariate numeric matrix that stores sample covariates;
#' @return A sample by gene numeric matrix that stores corrected gene expression values.
#' @export
covariate_removal <- function(expression_mat, covariate_df){
  if (nrow(expression_mat) != nrow(covariate_df)){
    stop("Sample size (nrow) of the expression matrix is different from that (nrow) in the covariate matrix!")
  }
  # Linear regression per gene
  corrected_mat <- matrix(nrow = nrow(expression_mat),
                          ncol = ncol(expression_mat))
  regression_df <- covariate_df
  for (i in 1:ncol(expression_mat)) {
    regression_df$expression <- as.numeric(expression_mat[, i])
    lm.summary <- summary(lm(expression ~ ., data = regression_df))
    corrected_mat[, i] <- lm.summary$residuals
  }
  rownames(corrected_mat) <- rownames(expression_mat)
  colnames(corrected_mat) <- colnames(expression_mat)
  return(corrected_mat)
  # scaled_mat <- scale(corrected_mat)
  # used as input Y in GSFA
}
