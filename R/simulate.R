#' @title Simulate a Continuous Gene Expression Matrix and an Accompanying
#' Perturbation Matrix
#' @description Generate a binary perturbation matrix and a continuous gene
#' expression matrix in a bottom-up fashion according to a hierarchical factor
#' model with normal noise terms.
#'
#' @param N Number of samples to simulate
#' @param P Number of genes to simulate
#' @param K Number of factors to simulate
#' @param M Number of perturbations to simulate
#' @param beta_true A \eqn{M} by \eqn{K} numeric matrix that stores the
#' true effect sizes of perturbation-factor associations; when \code{offset=TRUE},
#' \eqn{M+1} rows should be provided instead.
#' @param pi_true The true density (proportion of nonzero gene loading) of each factor
#' @param G_prob The Bernoulli probability based on which the binary
#' perturbation matrix \code{G} will be generated; determines the frequency
#' of each perturbation in the sample population
#' @param offset Default is FALSE. If TRUE, \code{beta_true} should have
#' \eqn{M+1} rows, with the last row storing the intercept values \eqn{\beta_0}
#'
#' @return A list object with the following elements:
#'
#' \item{Y}{a sample by gene matrix with continuous gene
#' expression values;}
#' \item{G}{a binary sample by perturbation matrix;}
#' \item{Z}{a sample by factor matrix;}
#' \item{F}{a binary gene by factor matrix that indicates whether a gene has
#' non-zero loading in the factor;}
#' \item{U}{a gene by factor matrix with normal effect sizes, and \code{F*U}
#' (element-wise multiplication) gives the loading matrix \code{W}.}
#' @export
#'
#' @examples
#' set.seed(12345)
#' beta_true <- rbind(c(1, 0, 0, 0, 0), c(0, 0.8, 0, 0, 0))
#' sim_data <- normal_data_sim(N = 4000, P = 6000, beta_true = beta_true)
#'
normal_data_sim <- function(N = 400, P = 600,
                            beta_true,
                            K = ncol(beta_true), M = nrow(beta_true),
                            pi_true = rep(0.1, K), sigma_w2_true = rep(0.5, K),
                            psi_true = 1, G_prob = 0.2,
                            offset = FALSE){
  stopifnot(length(sigma_w2_true) == K & length(pi_true) == K & ncol(beta_true) == K)

  G <- matrix(rbinom(n = N * M, size = 1, prob = G_prob),
              nrow = N, ncol = M)
  if (offset){
    ## If TRUE, the last row in beta_true contains all the intercepts
    stopifnot(nrow(beta_true) == (M + 1))
    G <- cbind(G, rep(1, nrow(G)))
  } else {
    stopifnot(nrow(beta_true) == M )
  }

  Phi <- matrix(rnorm(n = N*K, mean = 0, sd = 1),
                nrow = N, ncol = K) # Phi is fixed at ~ N(0,1)
  Z_true <- G %*% beta_true + Phi
  U_true <- matrix(0, nrow = P, ncol = K)
  F_true <- matrix(0, nrow = P, ncol = K)
  for (k in 1:K){
    U_true[, k] <- rnorm(P, mean = 0, sd = sqrt(sigma_w2_true[k]))
    F_true[, k] <- rbinom(P, 1, pi_true[k])
  }
  W_true <- F_true * U_true
  lower_rank <- Z_true %*% t(W_true)
  E <- matrix(rnorm(n = N*P, mean = 0, sd = sqrt(psi_true)),
              nrow = N, ncol = P)
  Y <- lower_rank + E

  if (offset){
    G <- G[, -ncol(G)]
  }

  return(list(G = G, Z = Z_true, U = U_true, F = F_true, Y = Y))
}

#' @title Add a Simulated Gene Count Matrix to a Dataset with Continuous Gene
#' Expression
#' @description Generate gene expression count data from a Poisson model based on
#' the continuous gene expression rates generated from a normal model,
#' incorporating varying scaling factors (library sizes) across cells.
#' @details Count \eqn{c_{ij}} is generated from a Poisson model based on
#' rate \eqn{y_{ij}} and library size (or scaling factor) \eqn{L_i} generated
#' from a normal model:
#' \deqn{L_i \sim N(\text{p_scale_mean}, \text{p_scale_sd}^2)}
#' \deqn{L_i = max(L_i, 0)}
#' \deqn{c_{ij} \sim Pois(L_i exp(y_{ij} + \text{p_offset}))}
#' @param sim_data A dataset of list type generated from \code{normal_data_sim()}
#' @param p_scale_mean The mean of the normal distribution where the library
#' size is drawn from
#' @param p_scale_sd The variance of the normal distribution where the library
#' size is drawn from
#' @param p_offset An offset number added to the input expression rates
#' @return A list that has a sample by gene count matrix added to the original
#' dataset input.
#' @export
#'
#' @examples
#' set.seed(12345)
#' beta_true <- rbind(c(1, 0, 0, 0, 0), c(0, 0.8, 0, 0, 0))
#' sim_data <- normal_data_sim(N = 4000, P = 6000, beta_true = beta_true)
#' sim_data <- poisson_count_sim_var_scale(sim_data)
#'
poisson_count_sim_var_scale <- function(sim_data,
                                        p_scale_mean = 5e5, p_scale_sd = 1e5,
                                        p_offset = 1/5e5){
  N <- nrow(sim_data$Y)
  P <- ncol(sim_data$Y)
  p_scale <- rnorm(N, mean = p_scale_mean, sd = p_scale_sd)
  if (sum(p_scale <= 0) > 0){
    min_pos_val <- min(p_scale[p_scale > 0])
    p_scale[p_scale <= 0] <- min_pos_val
  }
  count_mat <- matrix(0, nrow = N, ncol = P)
  for (i in 1:N){
    count_mat[i, ] <- rpois(P, p_scale[i] * exp(log(p_offset) + sim_data$Y[i, ]))
  }
  sim_data$count <- count_mat
  return(sim_data)
}

#' @description Generate gene expression count data from a Poisson model based on
#' the continuous gene expression rates generated from a normal model,
#' with a fixed scaling factor (library size) across all cells.
#' @param sim_data A dataset of list type generated from \code{normal_data_sim()}
poisson_count_sim <- function(sim_data,
                              p_scale = 5e5, p_offset = 1){
  Y_vec <- as.numeric(sim_data$Y)
  count_vec <- rpois(length(Y_vec), p_scale * exp(log(p_offset / p_scale) + Y_vec))
  sim_data$count <- matrix(count_vec, nrow = nrow(sim_data$Y))
  return(sim_data)
}

#' @description Generate gene expression count data from a multinomial model
#' based on the continuous gene expression rates generated from a normal model,
#' with a fixed library size across all cells.
#' @param sim_data A dataset of list type generated from \code{normal_data_sim()}
multinomial_count_gen <- function(sim_data, lib_size = 1e3){
  N <- nrow(sim_data$Y)
  P <- ncol(sim_data$Y)
  count_data <- matrix(nrow = N, ncol = P)
  for (i in 1:N){
    count_data[i, ] <- rmultinom(1, size = lib_size, prob = exp(sim_data$Y[i, ]))[, 1]
  }
  sim_data$count <- count_data
  return(sim_data)
}
