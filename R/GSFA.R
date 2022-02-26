#' @title Bayesian Guided Sparse Factor Analysis on Perturbed Gene Expression Matrix
#' @description Performs GSFA on given gene expression matrix and matching perturbation
#' information using Gibbs sampling.
#' @details Uses functions implemented in Rcpp from \code{GSFA_gibbs_inference.cpp}.
#' @param Y A sample by gene numeric matrix that stores normalized gene expression values;
#' \code{is.matrix(Y)} should be \code{TRUE};
#' @param G Either a numeric vector or a sample by perturbation numeric matrix
#' that stores sample-level perturbation information;
#' length or nrow of \code{G} should be the same as \code{nrow(Y)};
#' @param K Number of factors to use in the model; only one of \code{K}
#' and \code{fit0} is needed;
#' @param fit0 A list of class 'gsfa_fit' that is obtained from a previous
#' \code{fit_gsfa_multivar} run, so that more iterations of Gibbs sampling
#' can continue from the last updates in it;
#' @param prior_type Type of sparse prior used on gene weights, can be "mixture_normal"
#' or "spike_slab", "mixture_normal" sometimes works better in inducing sparsity;
#' @param init.method Method to initialize the factors, can be one of
#' "svd" (truncated SVD on \code{Y}) or "random";
#' @param niter Total number of Gibbs sampling iterations;
#' @param used_niter Number of iterations (counting from the last iteration)
#' from which the posterior means of parameters are to be computed;
#' @param lfsr_niter Number of iterations (counting from the last iteration)
#' of posterior samples to use for the computation of LFSR;
#' @param return_samples Boolean indicator of whether all posterior samples throughout
#' Gibbs sampling should be returned;
#' @return A list of class 'gsfa_fit' which stores the Gibbs sampling updates
#' and posterior mean estimates, and the prior parameters used during the inference.
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib GSFA
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit0 <- fit_gsfa_multivar(Y, G, 10, init.method = "svd", niter = 500, used_niter = 200)
#' fit1 <- fit_gsfa_multivar(Y, G, fit0 = fit0, niter = 500, used_niter = 200)
#' }
fit_gsfa_multivar <- function(Y, G, K, fit0,
                              prior_type = c("mixture_normal", "spike_slab"),
                              init.method = c("svd", "random"),
                              prior_w_s = 50, prior_w_r = 0.2,
                              prior_beta_s = 20, prior_beta_r = 0.2,
                              niter = 500,
                              used_niter = floor(niter/2),
                              lfsr_niter = used_niter,
                              verbose = TRUE, return_samples = TRUE){
  prior_type <- match.arg(prior_type)
  init.method <- match.arg(init.method)
  stopifnot(is.matrix(Y))
  # Add a column of 1's to G to help infer the offset
  if (is.vector(G) & is.numeric(G)){
    stopifnot(length(G) == nrow(Y))
    G <- cbind(G, rep(1, length(G)))
  } else if (is.matrix(G)){
    stopifnot(nrow(G) == nrow(Y))
    G <- cbind(G, rep(1, nrow(G)))
  } else {
    stop("G should be either a numeric vector or matrix.")
  }
  # Only one of "K" and "fit0" should be provided.
  # Argument "K" will be ignored if "fit0" is given.
  if (!(missing(K) & !missing(fit0) | (!missing(K) & missing(fit0)))){
    stop("Please provide the number of factors, \"K\", ",
         "or a Gibbs initialization, \"fit0\", but not both.")
  }

  if (missing(fit0)){
    fit <- gsfa_gibbs_cpp(Y = Y, G = G, K = K,
                          prior_type = prior_type,
                          initialize = init.method,
                          prior_s = prior_w_s, prior_r = prior_w_r,
                          prior_sb = prior_beta_s, prior_rb = prior_beta_r,
                          niter = niter, ave_niter = used_niter, lfsr_niter = lfsr_niter,
                          verbose = verbose, return_samples = return_samples)
  }
  if (missing(K)){
    # Check input argument "fit0".
    if (!inherits(fit0, "gsfa_fit")){
      stop("Input argument \"fit0\" should be an object of class ",
           "\"gsfa_fit\", such as an output of fit_gsfa_multivar().")
    }
    fit <- restart_gsfa_gibbs_cpp(Y = Y, G = G,
                                  Z = fit0$updates$Z,
                                  F = fit0$updates$F,
                                  W = fit0$updates$W,
                                  Gamma = fit0$updates$Gamma,
                                  beta = fit0$updates$beta,
                                  pi_vec = fit0$updates$pi_vec,
                                  pi_beta = fit0$updates$pi_beta,
                                  psi = fit0$updates$psi,
                                  sigma_w2 = fit0$updates$sigma_w2,
                                  sigma_b2 = fit0$updates$sigma_b2,
                                  c2 = fit0$updates$c2,
                                  prior_params = fit0$prior_params,
                                  prior_type = prior_type,
                                  niter = niter,
                                  ave_niter = used_niter,
                                  lfsr_niter = lfsr_niter,
                                  verbose = verbose,
                                  return_samples = return_samples)
  }
  factor_names <- paste0("Factor_", 1:ncol(fit$posterior_means$Z_pm))
  if (is.null(rownames(Y))){
    sample_names <- 1:nrow(Y)
  } else {
    sample_names <- rownames(Y)
  }
  if (is.null(colnames(Y))){
    gene_names <- 1:ncol(Y)
  } else {
    gene_names <- colnames(Y)
  }
  G <- G[, -ncol(G)]
  if (is.null(colnames(G))){
    perturbation_names <- 1:ncol(G)
  } else {
    perturbation_names <- colnames(G)
  }

  rownames(fit$posterior_means$Z_pm) <- sample_names
  colnames(fit$posterior_means$Z_pm) <- factor_names
  rownames(fit$posterior_means$F_pm) <- gene_names
  colnames(fit$posterior_means$F_pm) <- factor_names
  rownames(fit$posterior_means$W_pm) <- gene_names
  colnames(fit$posterior_means$W_pm) <- factor_names
  rownames(fit$posterior_means$Gamma_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$Gamma_pm) <- factor_names
  rownames(fit$posterior_means$beta_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$beta_pm) <- factor_names
  rownames(fit$lfsr) <- gene_names
  colnames(fit$lfsr) <- c(perturbation_names, "offset")
  class(fit) <- c("gsfa_fit", "list")
  return(fit)
}

#' @title Bayesian Guided Sparse Factor Analysis on Perturbed Gene Expression Matrix
#' @description Performs GSFA on given gene expression matrix and matching perturbation
#' information using Gibbs sampling for samples that come from two groups
#' @details Similar to the function \code{fit_gsfa_multivar()}, but associations
#' between factors and perturbations are estimated for each group of samples separately.
#' @param Y A sample by gene numeric matrix that stores normalized gene expression values;
#' \code{is.matrix(Y)} should be \code{TRUE};
#' @param G Either a numeric vector or a sample by perturbation numeric matrix
#' that stores sample-level perturbation information;
#' length or nrow of \code{G} should be the same as \code{nrow(Y)};
#' @param group a vector of sample size length, with two types of unique values
#' indicating one of the two groups each sample belongs to;
#' @param K Number of factors to use in the model; only one of \code{K}
#' and \code{fit0} is needed;
#' @param fit0 A list of class 'gsfa_fit' that is obtained from a previous
#' \code{fit_gsfa_multivar} run, so that more iterations of Gibbs sampling
#' can continue from the last updates in it;
#' only one of \code{K} and \code{fit0} is needed;
#' @param prior_type Type of sparse prior used on gene weights, can be "mixture_normal"
#' or "spike_slab", "mixture_normal" sometimes works better in inducing sparsity;
#' @param init.method Method to initialize the factors, can be one of
#' "svd" (truncated SVD on \code{Y}) or "random";
#' @param niter Total number of Gibbs sampling iterations;
#' @param used_niter Number of iterations (counting from the last iteration)
#' from which the posterior means of parameters are to be computed;
#' @param lfsr_niter Number of iterations (counting from the last iteration)
#' of posterior samples to use for the computation of LFSR;
#' @param return_samples Boolean indicator of whether all posterior samples throughout
#' Gibbs sampling should be returned;
#' @return A list of class 'gsfa_fit' which stores the Gibbs sampling updates
#' and posterior mean estimates, and the prior parameters used during the inference.
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib GSFA
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit0 <- fit_gsfa_multivar_2groups(Y, G, group, 10, init.method = "svd", niter = 500, used_niter = 200)
#' fit1 <- fit_gsfa_multivar_2groups(Y, G, group, fit0 = fit0, niter = 500, used_niter = 200)
#' }
fit_gsfa_multivar_2groups <- function(Y, G, group, K, fit0,
                                      prior_type = c("mixture_normal", "spike_slab"),
                                      init.method = c("svd", "random"),
                                      prior_w_s = 50, prior_w_r = 0.2,
                                      prior_beta_s = 20, prior_beta_r = 0.2,
                                      niter = 500,
                                      used_niter = floor(niter/2),
                                      lfsr_niter = used_niter,
                                      verbose = TRUE, return_samples = TRUE){
  prior_type <- match.arg(prior_type)
  init.method <- match.arg(init.method)
  stopifnot(is.matrix(Y))
  # Add a column of 1's to G to help infer the offset
  if (is.vector(G) & is.numeric(G)){
    stopifnot(length(G) == nrow(Y))
    G <- cbind(G, rep(1, length(G)))
  } else if (is.matrix(G)){
    stopifnot(nrow(G) == nrow(Y))
    G <- cbind(G, rep(1, nrow(G)))
  } else {
    stop("G should be either a numeric vector or matrix.")
  }

  if (is.vector(group) & length(group) == nrow(Y)){
    if (length(unique(group)) != 2){
      stop(paste0("There should be exactly 2 types of values in \"group\" to",
                  " indicate the 2 sample groups."))
    } else {
      if (!is.factor(group)){
        group <- factor(group)
      }
      print("The number of samples in each group:")
      print(table(group))
      group_levels <- levels(group)
      if (length(group_levels) != 2){
        stop("Please make sure your factor vector, \"group\", has only 2 levels.")
      }
      # Map "group" to a binary vector
      numeric_group <- as.numeric(group) - 1
      print(paste0("Samples of ", group_levels[1], " are assigned to group 0."))
      print(paste0("Samples of ", group_levels[2], " are assigned to group 1."))
    }
  } else {
    stop("\"group\" should be a vector with length equal to nrow(Y).")
  }

  # Only one of "K" and "fit0" should be provided.
  # Argument "K" will be ignored if "fit0" is given.
  if (!(missing(K) & !missing(fit0) | (!missing(K) & missing(fit0)))){
    stop("Please provide the number of factors, \"K\", ",
         "or a Gibbs initialization, \"fit0\", but not both.")
  }

  if (missing(fit0)){
    fit <- gsfa_gibbs_2groups_cpp(Y = Y, G = G,
                                  group = numeric_group, K = K,
                                  prior_type = prior_type,
                                  initialize = init.method,
                                  prior_s = prior_w_s, prior_r = prior_w_r,
                                  prior_sb = prior_beta_s, prior_rb = prior_beta_r,
                                  niter = niter, ave_niter = used_niter,
                                  lfsr_niter = lfsr_niter,
                                  verbose = verbose, return_samples = return_samples)
  }
  if (missing(K)){
    # Check input argument "fit0".
    if (!inherits(fit0, "gsfa_fit")){
      stop("Input argument \"fit0\" should be an object of class ",
           "\"gsfa_fit\", such as an output of fit_gsfa_multivar().")
    }
    fit <- restart_gibbs_2groups_cpp(Y = Y, G = G, group = numeric_group,
                                     Z = fit0$updates$Z,
                                     F = fit0$updates$F,
                                     W = fit0$updates$W,
                                     Gamma0 = fit0$updates$Gamma0, Gamma1 = fit0$updates$Gamma1,
                                     beta0 = fit0$updates$beta0, beta1 = fit0$updates$beta1,
                                     pi_vec = fit0$updates$pi_vec,
                                     pi_beta0 = fit0$updates$pi_beta0, pi_beta1 = fit0$updates$pi_beta1,
                                     psi = fit0$updates$psi,
                                     sigma_w2 = fit0$updates$sigma_w2,
                                     sigma_b20 = fit0$updates$sigma_b20, sigma_b21 = fit0$updates$sigma_b21,
                                     c2 = fit0$updates$c2,
                                     prior_params = fit0$prior_params,
                                     prior_type = prior_type,
                                     niter = niter,
                                     ave_niter = used_niter,
                                     lfsr_niter = lfsr_niter,
                                     verbose = verbose,
                                     return_samples = return_samples)
  }
  factor_names <- paste0("Factor_", 1:ncol(fit$posterior_means$Z_pm))
  if (is.null(rownames(Y))){
    sample_names <- 1:nrow(Y)
  } else {
    sample_names <- rownames(Y)
  }
  if (is.null(colnames(Y))){
    gene_names <- 1:ncol(Y)
  } else {
    gene_names <- colnames(Y)
  }
  G <- G[, -ncol(G)]
  if (is.null(colnames(G))){
    perturbation_names <- 1:ncol(G)
  } else {
    perturbation_names <- colnames(G)
  }

  rownames(fit$posterior_means$Z_pm) <- sample_names
  colnames(fit$posterior_means$Z_pm) <- factor_names
  rownames(fit$posterior_means$F_pm) <- gene_names
  colnames(fit$posterior_means$F_pm) <- factor_names
  rownames(fit$posterior_means$W_pm) <- gene_names
  colnames(fit$posterior_means$W_pm) <- factor_names
  rownames(fit$posterior_means$Gamma0_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$Gamma0_pm) <- factor_names
  rownames(fit$posterior_means$Gamma1_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$Gamma1_pm) <- factor_names
  rownames(fit$posterior_means$beta0_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$beta0_pm) <- factor_names
  rownames(fit$posterior_means$beta1_pm) <- c(perturbation_names, "offset")
  colnames(fit$posterior_means$beta1_pm) <- factor_names
  rownames(fit$lfsr0) <- gene_names
  colnames(fit$lfsr0) <- c(perturbation_names, "offset")
  rownames(fit$lfsr1) <- gene_names
  colnames(fit$lfsr1) <- c(perturbation_names, "offset")
  class(fit) <- c("gsfa_fit", "list")
  return(fit)
}
