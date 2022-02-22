#include <RcppArmadillo.h>
#include "sampling_functions.h"

using namespace Rcpp;

// FUNCTION DEFINITIONS
// ---------------------

//' Perform guided sparse factor analysis (GSFA) on the given expression matrix and guide matrix
//' using Gibbs sampling.
//'
//' @param Y sample by feature numeric matrix
//' @param G sample by covariate numeric matrix
//' @param K number of factors to infer in the model
//' @param prior_type character value indicating which prior to use on gene weights, recommend
//' to use the default "mixture_normal" (mixture-of-normals),
//' but "spkie_slab" (spike-and-slab) is also an option, but is sometimes insufficient to impose sparsity
//' @param initialize character value indicating which initialization method to use, can be one of
//' "svd", "random", or "given"
//' @param niter a numeric value indicating the total number of iterations Gibbs sampling should last
//' @param ave_niter a numeric value indicating how many of the last iterations
//' should be used to compute the posterior means
//' @return a list that stores the gibbs samples in each step, the posterior means,
//' as well as the prior parameters used.
//' @export
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List gsfa_gibbs_cpp(arma::mat Y, arma::mat G, int K,
                    String prior_type="mixture_normal",
                    String initialize="svd",
                    Rcpp::Nullable<Rcpp::NumericMatrix> Z_given=R_NilValue,
                    double prior_s=50, double prior_r=0.2,
                    double prior_sb=20, double prior_rb=0.2,
                    double prior_gp=1, double prior_hp=1,
                    double prior_gb=1, double prior_hb=1,
                    double prior_gw=1, double prior_hw=1,
                    double prior_gc=3, double prior_hc=0.5,
                    int niter=500, int ave_niter=200, int lfsr_niter=200,
                    bool verbose=true, bool return_samples=true){
  // Initialization
  int N = Y.n_rows;
  int P = Y.n_cols;
  int M = G.n_cols;

  List prior_params = List::create(Named("prior_s") = prior_s, Named("prior_r") = prior_r,
                                   Named("prior_sb") = prior_sb, Named("prior_rb") = prior_rb,
                                   Named("prior_gp") = prior_gp, Named("prior_hp") = prior_hp,
                                   Named("prior_gb") = prior_gb, Named("prior_hb") = prior_hb,
                                   Named("prior_gw") = prior_gw, Named("prior_hw") = prior_hw,
                                   Named("prior_gc") = prior_gc, Named("prior_hc") = prior_hc);
  arma::vec prior_pi(2);
  prior_pi(0) = prior_s;
  prior_pi(1) = prior_r;
  arma::vec pi_vec(K);
  pi_vec.fill(0.2);

  arma::vec prior_pibeta(2);
  prior_pibeta(0) = prior_sb;
  prior_pibeta(1) = prior_rb;
  arma::vec pi_beta(M);
  pi_beta.fill(0.2);

  arma::vec prior_psi(2);
  prior_psi(0) = prior_gp;
  prior_psi(1) = prior_hp;
  arma::vec psi(P);
  psi.fill(1.0);

  arma::vec prior_sigma2b(2);
  prior_sigma2b(0) = prior_gb;
  prior_sigma2b(1) = prior_hb;
  arma::vec sigma_b2(M);
  sigma_b2.fill(1.0);

  arma::vec prior_sigma2w(2);
  prior_sigma2w(0) = prior_gw;
  prior_sigma2w(1) = prior_hw;
  arma::vec sigma_w2(K);
  sigma_w2.fill(1.0);

  arma::mat Z,F,W;
  List intial_params;
  if (initialize == "given"){
    if (Z_given.isNull()){
      throw("No initial value of Z provided!");
    } else {
      Rprintf("Initializing with given Z.");
      arma::mat init_Z = Rcpp::as<arma::mat>(Z_given);
      intial_params = initialize_given_Z(K,Y,init_Z);
    }
  } else if (initialize == "svd"){
    if (Z_given.isNotNull()){
      warning("Provided value of Z ignored, initializing using truncated SVD!");
    }
    Rprintf("Initializing with SVD.");
    intial_params = initialize_svd(K,Y);
  } else {
    if (Z_given.isNotNull()){
      warning("Provided value of Z ignored, initializing using random values!");
    }
    Rprintf("Initializing with random values.");
    intial_params = initialize_random(K,Y);
  }

  Z = as<arma::mat>(intial_params["Z"]);
  F = as<arma::mat>(intial_params["F"]);
  W = as<arma::mat>(intial_params["W"]);
  arma::mat beta = arma::zeros<arma::mat>(M,K);
  arma::vec tmp_beta;
  for (int m=0; m<M; m++){
    arma::vec G_col = G.col(m);
    arma::vec tmp_beta = Z.t() * G_col / (sum(G_col%G_col)+1e-4);
    beta.row(m) = arma::conv_to< arma::rowvec >::from(tmp_beta);
  }
  arma::umat init_bool = abs(beta) > Cquantile(vectorise(abs(beta)), 0.5);
  arma::mat Gamma = arma::conv_to< arma::mat >::from(init_bool);
  beta = beta % Gamma;

  arma::cube Z_mtx = arma::zeros<arma::cube>(N,K,niter+1);
  Z_mtx.slice(0) = Z;
  arma::cube F_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  F_mtx.slice(0) = F;
  arma::cube W_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  W_mtx.slice(0) = W;
  arma::cube beta_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta_mtx.slice(0) = beta;
  arma::cube Gamma_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma_mtx.slice(0) = Gamma;
  arma::mat pi_mtx = arma::zeros<arma::mat>(K,niter+1);
  pi_mtx.col(0) = pi_vec;
  arma::mat pi_beta_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta_mtx.col(0) = pi_beta;
  arma::mat sigma_w2_mtx = arma::zeros<arma::mat>(K,niter+1);
  sigma_w2_mtx.col(0) = sigma_w2;

  arma::vec prior_c(2);
  prior_c(0) = prior_gc;
  prior_c(1) = prior_hc;

  arma::vec c2(K, arma::fill::zeros);
  arma::mat c2_mtx(K,niter+1, arma::fill::zeros);

  if (prior_type=="mixture_normal") {
    c2.fill(0.25);
    c2_mtx.col(0) = c2;
  }
  if (prior_type=="spike_slab") {
    c2.fill(R_NaN);
    c2_mtx.fill(R_NaN);
  }
  List gammaBeta;
  List FW;

  // Gibbs Sampling:
  int iter = 1;

  while (iter <= niter) {
    gammaBeta = sample_gammaBeta_cpp(N, M, K, Z, G, Gamma, beta, sigma_b2, pi_beta);
    Gamma = as<arma::mat>(gammaBeta["Gamma"]);
    beta = as<arma::mat>(gammaBeta["beta"]);
    pi_beta = sample_pi_beta_cpp(M, K, Gamma, prior_pibeta);
    sigma_b2 = sample_sigma_b2_cpp(M, Gamma, beta, prior_sigma2b);

    Z = sample_Z_cpp(N,K,Y,F,W,G,beta,psi);
    psi = sample_psi_cpp(N,P,Y,Z,F,W,prior_psi);

    if (prior_type=="mixture_normal") {
      W = sample_W_cpp(P, K, Y, Z, F, W, psi, sigma_w2, c2);
      F = sample_F_cpp(P, K, W, pi_vec, sigma_w2, c2);
      pi_vec = sample_pi_cpp(P, K, F, prior_pi);
      sigma_w2 = sample_sigma_w2_cpp(K, P, F, W, prior_sigma2w, c2);
      c2 = sample_c2_cpp(K, P, F, W, sigma_w2, prior_c);

      c2_mtx.col(iter) = c2;
    }
    if (prior_type=="spike_slab") {
      FW = sample_FW_spike_slab_cpp(N, P, K, Y, Z, F, W, psi, sigma_w2, pi_vec);
      F = as<arma::mat>(FW["F"]);
      W = as<arma::mat>(FW["W"]);
      pi_vec = sample_pi_cpp(P, K, F, prior_pi);
      sigma_w2 = sample_sigma_w2_spike_slab_cpp(K, F, W, prior_sigma2w);
    }

    // Store samples through iterations:
    Gamma_mtx.slice(iter) = Gamma;
    beta_mtx.slice(iter) = beta;
    pi_beta_mtx.col(iter) = pi_beta;
    Z_mtx.slice(iter) = Z;
    F_mtx.slice(iter) = F;
    W_mtx.slice(iter) = W;
    pi_mtx.col(iter) = pi_vec;
    sigma_w2_mtx.col(iter) = sigma_w2;

    if (verbose) {
      if (iter % 50 == 0) {
        Rprintf("Iteration [%i] finished.\n", iter);
        Rcpp::checkUserInterrupt();
      }
    }
    iter += 1;
  }
  // Save the latest updates in a list:
  List update_list;
  update_list = List::create(Named("Z") = Z, Named("F") = F, Named("W") = W,
                             Named("Gamma") = Gamma, Named("beta") = beta,
                             Named("pi_vec") = pi_vec, Named("pi_beta") = pi_beta,
                             Named("psi") = psi,
                             Named("sigma_w2") = sigma_w2, Named("sigma_b2") = sigma_b2,
                             Named("c2") = c2,
                             Named("niters") = niter, Named("average_niters") = ave_niter);
  // Compute the posterior means:
  List pm_list;
  pm_list = compute_posterior_mean_cpp(Gamma_mtx, beta_mtx, pi_beta_mtx, Z_mtx,
                                       F_mtx, W_mtx, pi_mtx, sigma_w2_mtx, c2_mtx,
                                       niter, ave_niter, prior_type);
  // Compute local false sign rate for each perturbation-gene pair:
  arma::mat lfsr_mat(P,M);
  lfsr_mat = compute_lfsr_cpp(beta_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);
  if (return_samples) {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr") = lfsr_mat,
                        Named("Z_samples") = Z_mtx,
                        Named("F_samples") = F_mtx,
                        Named("W_samples") = W_mtx,
                        Named("Gamma_samples") = Gamma_mtx,
                        Named("beta_samples") = beta_mtx,
                        Named("pi_samples") = pi_mtx,
                        Named("pi_beta_samples") = pi_beta_mtx,
                        Named("sigma_w2_samples") = sigma_w2_mtx,
                        Named("c2_samples") = c2_mtx,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  } else {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr") = lfsr_mat,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  }
}

//' @export
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List restart_gsfa_gibbs_cpp(arma::mat Y, arma::mat G,
                            arma::mat Z, arma::mat F, arma::mat W,
                            arma::mat Gamma, arma::mat beta,
                            arma::vec pi_vec, arma::vec pi_beta,
                            arma::vec psi, arma::vec sigma_w2, arma::vec sigma_b2,
                            arma::vec c2,
                            List prior_params,
                            String prior_type="mixture_normal",
                            int niter=500, int ave_niter=200, int lfsr_niter=200,
                            bool verbose=true, bool return_samples=true){
  // Preparation:
  int N = Y.n_rows;
  int P = Y.n_cols;
  int K = Z.n_cols;
  int M = G.n_cols;

  arma::vec prior_pi(2);
  prior_pi(0) = as<double>(prior_params["prior_s"]);
  prior_pi(1) = as<double>(prior_params["prior_r"]);

  arma::vec prior_pibeta(2);
  prior_pibeta(0) = as<double>(prior_params["prior_sb"]);
  prior_pibeta(1) = as<double>(prior_params["prior_rb"]);

  arma::vec prior_psi(2);
  prior_psi(0) = as<double>(prior_params["prior_gp"]);
  prior_psi(1) = as<double>(prior_params["prior_hp"]);

  arma::vec prior_sigma2w(2);
  prior_sigma2w(0) = as<double>(prior_params["prior_gw"]);
  prior_sigma2w(1) = as<double>(prior_params["prior_hw"]);

  arma::vec prior_sigma2b(2);
  prior_sigma2b(0) = as<double>(prior_params["prior_gb"]);
  prior_sigma2b(1) = as<double>(prior_params["prior_hb"]);

  arma::vec prior_c(2);
  prior_c(0) = as<double>(prior_params["prior_gc"]);
  prior_c(1) = as<double>(prior_params["prior_hc"]);

  arma::cube Z_mtx = arma::zeros<arma::cube>(N,K,niter+1);
  Z_mtx.slice(0) = Z;
  arma::cube F_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  F_mtx.slice(0) = F;
  arma::cube W_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  W_mtx.slice(0) = W;

  arma::cube Gamma_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma_mtx.slice(0) = Gamma;
  arma::cube beta_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta_mtx.slice(0) = beta;

  arma::mat pi_mtx = arma::zeros<arma::mat>(K,niter+1);
  pi_mtx.col(0) = pi_vec;
  arma::mat pi_beta_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta_mtx.col(0) = pi_beta;
  arma::mat sigma_w2_mtx = arma::zeros<arma::mat>(K,niter+1);
  sigma_w2_mtx.col(0) = sigma_w2;
  arma::mat c2_mtx(K,niter+1, arma::fill::zeros);

  if (prior_type=="mixture_normal") {
    c2.fill(0.25);
    c2_mtx.col(0) = c2;
  }
  if (prior_type=="spike_slab") {
    c2.fill(R_NaN);
    c2_mtx.fill(R_NaN);
  }

  // Gibbs Sampling:
  List gammaBeta;
  List FW;

  int iter = 1;

  while (iter <= niter) {
    gammaBeta = sample_gammaBeta_cpp(N, M, K, Z, G, Gamma, beta, sigma_b2, pi_beta);
    Gamma = as<arma::mat>(gammaBeta["Gamma"]);
    beta = as<arma::mat>(gammaBeta["beta"]);
    pi_beta = sample_pi_beta_cpp(M, K, Gamma, prior_pibeta);
    sigma_b2 = sample_sigma_b2_cpp(M, Gamma, beta, prior_sigma2b);

    Z = sample_Z_cpp(N,K,Y,F,W,G,beta,psi);
    psi = sample_psi_cpp(N,P,Y,Z,F,W,prior_psi);

    if (prior_type=="mixture_normal") {
      W = sample_W_cpp(P, K, Y, Z, F, W, psi, sigma_w2, c2);
      F = sample_F_cpp(P, K, W, pi_vec, sigma_w2, c2);
      pi_vec = sample_pi_cpp(P, K, F, prior_pi);
      sigma_w2 = sample_sigma_w2_cpp(K, P, F, W, prior_sigma2w, c2);
      c2 = sample_c2_cpp(K, P, F, W, sigma_w2, prior_c);

      c2_mtx.col(iter) = c2;
    }
    if (prior_type=="spike_slab") {
      FW = sample_FW_spike_slab_cpp(N, P, K, Y, Z, F, W, psi, sigma_w2, pi_vec);
      F = as<arma::mat>(FW["F"]);
      W = as<arma::mat>(FW["W"]);
      pi_vec = sample_pi_cpp(P, K, F, prior_pi);
      sigma_w2 = sample_sigma_w2_spike_slab_cpp(K, F, W, prior_sigma2w);
    }
    // Store samples through iterations:
    Gamma_mtx.slice(iter) = Gamma;
    beta_mtx.slice(iter) = beta;
    pi_beta_mtx.col(iter) = pi_beta;
    Z_mtx.slice(iter) = Z;
    F_mtx.slice(iter) = F;
    W_mtx.slice(iter) = W;
    pi_mtx.col(iter) = pi_vec;
    sigma_w2_mtx.col(iter) = sigma_w2;

    if (verbose) {
      if (iter % 50 == 0) {
        Rprintf("Iteration [%i] finished.\n", iter);
        Rcpp::checkUserInterrupt();
      }
    }
    iter += 1;
  }

  // Save the latest updates in a list:
  List update_list;
  update_list = List::create(Named("Z") = Z, Named("F") = F, Named("W") = W,
                             Named("Gamma") = Gamma, Named("beta") = beta,
                             Named("pi_vec") = pi_vec, Named("pi_beta") = pi_beta,
                             Named("psi") = psi,
                             Named("sigma_w2") = sigma_w2, Named("sigma_b2") = sigma_b2,
                             Named("c2") = c2,
                             Named("niters") = niter, Named("average_niters") = ave_niter);
  // Compute the posterior means:
  List pm_list;
  pm_list = compute_posterior_mean_cpp(Gamma_mtx, beta_mtx, pi_beta_mtx, Z_mtx,
                                       F_mtx, W_mtx, pi_mtx, sigma_w2_mtx, c2_mtx,
                                       niter, ave_niter, prior_type);
  // Compute local false sign rate for each perturbation-gene pair:
  arma::mat lfsr_mat(P,M);
  lfsr_mat = compute_lfsr_cpp(beta_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);
  if (return_samples) {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr") = lfsr_mat,
                        Named("Z_samples") = Z_mtx,
                        Named("F_samples") = F_mtx,
                        Named("W_samples") = W_mtx,
                        Named("Gamma_samples") = Gamma_mtx,
                        Named("beta_samples") = beta_mtx,
                        Named("pi_samples") = pi_mtx,
                        Named("pi_beta_samples") = pi_beta_mtx,
                        Named("sigma_w2_samples") = sigma_w2_mtx,
                        Named("c2_samples") = c2_mtx,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  } else {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr") = lfsr_mat,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  }
}
