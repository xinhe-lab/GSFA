#' @title Simulate a dataset with gene expression and perturbation according to a normal model
#' @description Generate a binary perturbation matrix and a continuous gene
#' expression matrix according to a factor model with normal noise terms.
#' @param N Number of samples to simulate
#' @param P Number of genes to simulate
#' @param K Number of factors to simulate
#' @param M Number of perturbations to simulate
#' @param beta_true A perturbation by factor numeric matrix that stores association
#' @param pi_true True sparsity of each factor
#' @param G_prob The Bernoulli probability based on which to generate the binary
#' perturbation matrix G
#' effect sizes;
#' @return A list object with the following elements:
#'
#' \item{Y}{a sample by gene matrix with continuous gene
#' expression values;}
#' \item{G}{a binary sample by perturbation matrix;}
#' \item{Z}{a sample by factor matrix;}
#' \item{F}{a binary gene by factor matrix that indicates whether a gene has
#' non-zero loading in the factor;}
#' \item{U}{a gene by factor matrix with normal effect sizes, and F*U gives W,
#' the loading matrix.}
#' @export
normal_data_sim <- function(N = 400, P = 600, K = 5, M = 2,
                            beta_true = matrix(c(1, 0, 0, 0.8, 0, 0, 0, 0, 0, 0),
                                               ncol = K),
                            pi_true = rep(0.1, K), sigma_w2_true = rep(0.5, K),
                            psi_true = 1, G_prob = 0.05,
                            seed = 46568, offset = FALSE){
  set.seed(seed)
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

#' @title Simulate gene expression count data
#' @description Generate gene expression count data from a Poisson model based on
#' the continuous gene expression rates generated from a normal model.
#' @details Count \eqn{c_{ij}} is generated from a Poisson model based on
#' rate \eqn{y_{ij}} and library size (or scaling factor) \eqn{L_i} generated
#' from a normal model.
#' \deqn{c_{ij} \sim Pois(L_i exp(y_{ij} + \text{p_offset}))}
#' @param sim_data A dataset generated from \code{normal_data_sim()}
#' @param p_scale_mean The mean of the normal distribution where the library
#' size is drawn from
#' @param p_scale_sd The variance of the normal distribution where the library
#' size is drawn from
#' @param p_offset An offset number added to the input expression rates
#' @return A list that has a sample by gene count matrix added to the original
#' dataset input.
#' @export
poisson_count_sim_var_scale <- function(sim_data,
                                        p_scale_mean = 5e5, p_scale_sd = 1e5,
                                        p_offset = 1/5e5){

  p_scale <- rnorm(parameters$N, mean = p_scale_mean, sd = p_scale_sd)
  if (sum(p_scale <= 0) > 0){
    min_pos_val <- min(p_scale[p_scale > 0])
    p_scale[p_scale <= 0] <- min_pos_val
  }
  count_mat <- matrix(0, nrow = parameters$N, ncol = parameters$P)
  for (i in 1:parameters$N){
    count_mat[i, ] <- rpois(parameters$P, p_scale[i] * exp(log(p_offset) + sim_data$Y[i, ]))
  }
  sim_data$count <- count_mat
  return(sim_data)
}

poisson_count_sim <- function(sim_data,
                              p_scale = 5e5, p_offset = 1){
  # sim_data is generated from normal_data_sim()
  Y_vec <- as.numeric(sim_data$Y)
  count_vec <- rpois(length(Y_vec), p_scale * exp(log(p_offset / p_scale) + Y_vec))
  sim_data$count <- matrix(count_vec, nrow = nrow(sim_data$Y))
  return(sim_data)
}

multinomial_count_gen <- function(sim_data, lib_size = 1e3){
  # sim_data is generated from data_gen_multi_markers()
  count_data <- matrix(nrow = parameters$N, ncol = parameters$P)
  for (i in 1:parameters$N){
    count_data[i, ] <- rmultinom(1, size = lib_size, prob = exp(sim_data$Y[i, ]))[, 1]
  }
  sim_data$count <- count_data
  return(sim_data)
}
