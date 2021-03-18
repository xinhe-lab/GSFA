#' @title Bayesian Guided Sparse Factor Analysis on Perturbed Gene Expression Matrix
#' @description Performs GSFA on given gene expression matrix and matching perturbation
#' (guide) information using Gibbs sampling.
#' @details Uses functions implemented in Rcpp from \code{GSFA_gibbs_inference.cpp}.
#' @param Y A numeric matrix that stores sample by gene normalized gene expression values;
#' \code{is.matrix(Y)} should be \code{TRUE};
#' @param G Either a numeric vector or a numeric matrix that stores sample-level guide information;
#' length or nrow of \code{G} should be the same as \code{nrow(Y)};
#' @param K An integer specifying the number of factors in the model; only one of \code{K}
#' and \code{fit0} is needed;
#' @param fit0 A list of class 'gsfa_fit' that is obtained from a previous \code{fit_gsfa_multivar}
#' run, so that more iterations of Gibbs sampling can be carried out on top of it;
#' only one of \code{K} and \code{fit0} is needed;
#' @param init.method The method to initialize the factors, can be one of "svd", "random", "given";
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
#' fit0 <- fit_gsfa_multivar(Y, G, 10, init.method = "svd", niter = 500, average_niter = 200)
#' fit1 <- fit_gsfa_multivar(Y, G, fit0 = fit0, niter = 500, average_niter = 200)
#' }
fit_gsfa_multivar <- function(Y, G, K, fit0,
                              init.method = c("svd", "random", "given"),
                              prior_w_s = 50, prior_w_r = 0.2,
                              prior_beta_s = 10, prior_beta_r = 0.2,
                              niter = 500, average_niter = 200, lfsr_niter = average_niter,
                              verbose = TRUE, return_samples = FALSE){
  init.method <- match.arg(init.method)
  stopifnot(is.matrix(Y))
  # Add an offset of 1
  if (is.vector(G)){
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
                          initialize = init.method,
                          prior_s = prior_w_s, prior_r = prior_w_r,
                          prior_sb = prior_beta_s, prior_rb = prior_beta_r,
                          niter = niter, ave_niter = average_niter, lfsr_niter = lfsr_niter,
                          verbose = verbose, return_samples = return_samples)
  }
  if (missing(K)){
    # Check input argument "fit0".
    if (!inherits(fit0, "gsfa_fit")){
      stop("Input argument \"fit0\" should be an object of class ",
           "\"gsfa_fit\", such as an output of fit_gsfa_multivar().")
    }
    fit <- restart_gibbs_cpp(Y = Y, G = G,
                             Z = fit0$updates$Z, F = fit0$updates$F, W = fit0$updates$W,
                             Gamma = fit0$updates$Gamma, beta = fit0$updates$beta,
                             pi_vec = fit0$updates$pi_vec, pi_beta = fit0$updates$pi_beta,
                             psi = fit0$updates$psi,
                             sigma_w2 = fit0$updates$sigma_w2, sigma_b2 = fit0$updates$sigma_b2,
                             prior_params = fit0$prior_params,
                             niter = niter, ave_niter = average_niter, lfsr_niter = lfsr_niter,
                             verbose = verbose, return_samples = return_samples)
  }
  class(fit) <- c("gsfa_fit", "list")
  return(fit)
}
