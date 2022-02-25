#include <RcppArmadillo.h>
#include "sampling_functions.h"

using namespace Rcpp;

// FUNCTION DEFINITIONS
// ---------------------

//' Perform guided sparse factor analysis (GSFA) given a normalized expression
//' matrix and a perturbation matrix using Gibbs sampling;
//' samples come from two groups, and associations between factors and
//' perturbations are estimated within each group separately.
//'
//' @param Y sample by gene numeric matrix
//' @param G sample by perturbation numeric matrix
//' @param group binary vector indicating the group each sample belongs to
//' @param K number of factors to infer in the model
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
List gsfa_gibbs_2groups_cpp(arma::mat Y, arma::mat G, arma::vec group, int K,
                            String prior_type="mixture_normal",
                            String initialize="svd",
                            Rcpp::Nullable<Rcpp::NumericMatrix> Z_given=R_NilValue,
                            double prior_s=50, double prior_r=0.1,
                            double prior_sb=20, double prior_rb=0.2,
                            double prior_gp=1, double prior_hp=1,
                            double prior_gb=1, double prior_hb=1,
                            double prior_gw=1, double prior_hw=1,
                            double prior_gc=3, double prior_hc=0.5,
                            int niter=500, int ave_niter=200, int lfsr_niter=200,
                            bool verbose=true, bool return_samples=true){
  // Organizing specified hyperparameters:
  List prior_params = List::create(Named("prior_s") = prior_s, Named("prior_r") = prior_r,
                                   Named("prior_sb") = prior_sb, Named("prior_rb") = prior_rb,
                                   Named("prior_gp") = prior_gp, Named("prior_hp") = prior_hp,
                                   Named("prior_gb") = prior_gb, Named("prior_hb") = prior_hb,
                                   Named("prior_gw") = prior_gw, Named("prior_hw") = prior_hw,
                                   Named("prior_gc") = prior_gc, Named("prior_hc") = prior_hc);

  arma::vec prior_pi(2);
  prior_pi(0) = prior_s;
  prior_pi(1) = prior_r;

  arma::vec prior_pibeta(2);
  prior_pibeta(0) = prior_sb;
  prior_pibeta(1) = prior_rb;

  arma::vec prior_psi(2);
  prior_psi(0) = prior_gp;
  prior_psi(1) = prior_hp;

  arma::vec prior_sigma2b(2);
  prior_sigma2b(0) = prior_gb;
  prior_sigma2b(1) = prior_hb;

  arma::vec prior_sigma2w(2);
  prior_sigma2w(0) = prior_gw;
  prior_sigma2w(1) = prior_hw;

  arma::vec prior_c(2);
  prior_c(0) = prior_gc;
  prior_c(1) = prior_hc;

  // Initialization:
  int N = Y.n_rows;
  int P = Y.n_cols;
  int M = G.n_cols;
  arma::uvec index0 = arma::find(group==0);
  arma::uvec index1 = arma::find(group==1);
  int N0 = index0.n_elem;
  int N1 = index1.n_elem;
  Rprintf("There are %i samples in group 0, ", N0);
  Rprintf("and %i samples in group 1.\n", N1);
  arma::mat Y0 = Y.rows(index0);
  arma::mat Y1 = Y.rows(index1);
  arma::mat G0 = G.rows(index0);
  arma::mat G1 = G.rows(index1);

  arma::vec pi_vec(K);
  pi_vec.fill(0.2);

  arma::vec pi_beta0(M);
  pi_beta0.fill(0.2);
  arma::vec pi_beta1(M);
  pi_beta1.fill(0.2);

  arma::vec psi(P);
  psi.fill(1.0);

  arma::vec sigma_b20(M);
  sigma_b20.fill(1.0);
  arma::vec sigma_b21(M);
  sigma_b21.fill(1.0);

  arma::vec sigma_w2(K);
  sigma_w2.fill(1.0);

  arma::vec c2(K, arma::fill::zeros);
  arma::mat c2_mtx(K, niter+1, arma::fill::zeros);
  if (prior_type=="mixture_normal") {
    c2.fill(0.25);
    c2_mtx.col(0) = c2;
  }
  if (prior_type=="spike_slab") {
    c2.fill(R_NaN);
    c2_mtx.fill(R_NaN);
  }

  List intial_ZFW_params;
  if (initialize=="svd"){
    intial_ZFW_params = initialize_svd(K,Y);
  } else {
    intial_ZFW_params = initialize_random(K,Y);
  }
  arma::mat Z = as<arma::mat>(intial_ZFW_params["Z"]);
  arma::mat F = as<arma::mat>(intial_ZFW_params["F"]);
  arma::mat W = as<arma::mat>(intial_ZFW_params["W"]);

  arma::mat Z0 = Z.rows(index0);
  arma::mat Z1 = Z.rows(index1);

  List gammaBeta0;
  gammaBeta0 = initialize_gammaBeta(M, K, G0, Z0);
  arma::mat Gamma0 = as<arma::mat>(gammaBeta0["Gamma"]);
  arma::mat beta0 = as<arma::mat>(gammaBeta0["beta"]);

  List gammaBeta1;
  gammaBeta1 = initialize_gammaBeta(M, K, G1, Z1);
  arma::mat Gamma1 = as<arma::mat>(gammaBeta1["Gamma"]);
  arma::mat beta1 = as<arma::mat>(gammaBeta1["beta"]);

  // Storing the initial values of parameters as samples at iteration 0:
  arma::cube Z_mtx = arma::zeros<arma::cube>(N,K,niter+1);
  Z_mtx.slice(0) = Z;
  arma::cube F_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  F_mtx.slice(0) = F;
  arma::cube W_mtx = arma::zeros<arma::cube>(P,K,niter+1);
  W_mtx.slice(0) = W;

  arma::cube beta0_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta0_mtx.slice(0) = beta0;
  arma::cube Gamma0_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma0_mtx.slice(0) = Gamma0;
  arma::mat pi_beta0_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta0_mtx.col(0) = pi_beta0;

  arma::cube beta1_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta1_mtx.slice(0) = beta1;
  arma::cube Gamma1_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma1_mtx.slice(0) = Gamma1;
  arma::mat pi_beta1_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta1_mtx.col(0) = pi_beta1;

  arma::mat pi_mtx = arma::zeros<arma::mat>(K,niter+1);
  pi_mtx.col(0) = pi_vec;
  arma::mat sigma_w2_mtx = arma::zeros<arma::mat>(K,niter+1);
  sigma_w2_mtx.col(0) = sigma_w2;

  // Gibbs Sampling:
  List FW;
  int iter = 1;

  while (iter <= niter) {
    gammaBeta0 = sample_gammaBeta_cpp(N0, M, K, Z0, G0, Gamma0, beta0, sigma_b20, pi_beta0);
    Gamma0 = as<arma::mat>(gammaBeta0["Gamma"]);
    beta0 = as<arma::mat>(gammaBeta0["beta"]);
    gammaBeta1 = sample_gammaBeta_cpp(N1, M, K, Z1, G1, Gamma1, beta1, sigma_b21, pi_beta1);
    Gamma1 = as<arma::mat>(gammaBeta1["Gamma"]);
    beta1 = as<arma::mat>(gammaBeta1["beta"]);

    pi_beta0 = sample_pi_beta_cpp(M, K, Gamma0, prior_pibeta);
    sigma_b20 = sample_sigma_b2_cpp(M, Gamma0, beta0, prior_sigma2b);
    pi_beta1 = sample_pi_beta_cpp(M, K, Gamma1, prior_pibeta);
    sigma_b21 = sample_sigma_b2_cpp(M, Gamma1, beta1, prior_sigma2b);

    Z0 = sample_Z_cpp(N0, K, Y0, F, W, G0, beta0, psi);
    Z1 = sample_Z_cpp(N1, K, Y1, F, W, G1, beta1, psi);
    Z.rows(index0) = Z0;
    Z.rows(index1) = Z1;

    psi = sample_psi_cpp(N, P, Y, Z, F, W, prior_psi);
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

    // Storing samples throughout iterations:
    Gamma0_mtx.slice(iter) = Gamma0;
    beta0_mtx.slice(iter) = beta0;
    pi_beta0_mtx.col(iter) = pi_beta0;

    Gamma1_mtx.slice(iter) = Gamma1;
    beta1_mtx.slice(iter) = beta1;
    pi_beta1_mtx.col(iter) = pi_beta1;

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
  // Save the latest updates for future Gibbs sampling to pick up from:
  List update_list;
  update_list = List::create(Named("Z") = Z, Named("F") = F, Named("W") = W,
                             Named("Gamma0") = Gamma0, Named("Gamma1") = Gamma1,
                             Named("beta0") = beta0, Named("beta1") = beta1,
                             Named("pi_vec") = pi_vec,
                             Named("pi_beta0") = pi_beta0, Named("pi_beta1") = pi_beta1,
                             Named("psi") = psi,
                             Named("sigma_w2") = sigma_w2,
                             Named("sigma_b20") = sigma_b20, Named("sigma_b21") = sigma_b21,
                             Named("c2") = c2,
                             Named("niters") = niter, Named("average_niters") = ave_niter);
  // Compute the posterior means:
  List pm_list;
  pm_list = compute_posterior_mean_2groups_cpp(Gamma0_mtx, beta0_mtx, pi_beta0_mtx,
                                               Gamma1_mtx, beta1_mtx, pi_beta1_mtx,
                                               Z_mtx, F_mtx, W_mtx, pi_mtx,
                                               sigma_w2_mtx, c2_mtx,
                                               niter, ave_niter, prior_type);
  // Compute local false sign rate for each perturbation-gene pair and each of the sample groups:
  arma::mat lfsr0_mat(P,M);
  lfsr0_mat = compute_lfsr_cpp(beta0_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);
  arma::mat lfsr1_mat(P,M);
  lfsr1_mat = compute_lfsr_cpp(beta1_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);

  if (return_samples) {
    // Save samples at each iteration on top of everything else:
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr0") = lfsr0_mat,
                        Named("lfsr1") = lfsr1_mat,
                        Named("Z_samples") = Z_mtx,
                        Named("F_samples") = F_mtx,
                        Named("W_samples") = W_mtx,
                        Named("Gamma0_samples") = Gamma0_mtx,
                        Named("beta0_samples") = beta0_mtx,
                        Named("pi_beta0_samples") = pi_beta0_mtx,
                        Named("Gamma1_samples") = Gamma1_mtx,
                        Named("beta1_samples") = beta1_mtx,
                        Named("pi_beta1_samples") = pi_beta1_mtx,
                        Named("pi_samples") = pi_mtx,
                        Named("sigma_w2_samples") = sigma_w2_mtx,
                        Named("c2_samples") = c2_mtx,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  } else {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr0") = lfsr0_mat,
                        Named("lfsr1") = lfsr1_mat,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  }
}

//' @export
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List restart_gibbs_2groups_cpp(arma::mat Y, arma::mat G, arma::vec group,
                               arma::mat Z, arma::mat F, arma::mat W,
                               arma::mat Gamma0, arma::mat Gamma1,
                               arma::mat beta0, arma::mat beta1,
                               arma::vec pi_vec,
                               arma::vec pi_beta0, arma::vec pi_beta1,
                               arma::vec psi, arma::vec sigma_w2,
                               arma::vec sigma_b20, arma::vec sigma_b21,
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
  arma::uvec index0 = arma::find(group==0);
  arma::uvec index1 = arma::find(group==1);
  int N0 = index0.n_elem;
  int N1 = index1.n_elem;
  Rprintf("There are %i samples in group 0, ", N0);
  Rprintf("and %i samples in group 1.\n", N1);
  arma::mat Y0 = Y.rows(index0);
  arma::mat Y1 = Y.rows(index1);
  arma::mat G0 = G.rows(index0);
  arma::mat G1 = G.rows(index1);


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

  arma::cube beta0_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta0_mtx.slice(0) = beta0;
  arma::cube Gamma0_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma0_mtx.slice(0) = Gamma0;
  arma::mat pi_beta0_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta0_mtx.col(0) = pi_beta0;

  arma::cube beta1_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  beta1_mtx.slice(0) = beta1;
  arma::cube Gamma1_mtx = arma::zeros<arma::cube>(M,K,niter+1);
  Gamma1_mtx.slice(0) = Gamma1;
  arma::mat pi_beta1_mtx = arma::zeros<arma::mat>(M,niter+1);
  pi_beta1_mtx.col(0) = pi_beta1;

  arma::mat pi_mtx = arma::zeros<arma::mat>(K,niter+1);
  pi_mtx.col(0) = pi_vec;
  arma::mat sigma_w2_mtx = arma::zeros<arma::mat>(K,niter+1);
  sigma_w2_mtx.col(0) = sigma_w2;

  arma::mat c2_mtx(K,niter+1, arma::fill::zeros);

  if (prior_type=="mixture_normal") {
    c2_mtx.col(0) = c2;
  }
  if (prior_type=="spike_slab") {
    c2_mtx.fill(R_NaN);
  }

  // Gibbs Sampling:
  arma::mat Z0 = Z.rows(index0);
  arma::mat Z1 = Z.rows(index1);
  List gammaBeta0;
  List gammaBeta1;
  List FW;
  int iter = 1;

  while (iter <= niter) {
    gammaBeta0 = sample_gammaBeta_cpp(N0, M, K, Z0, G0, Gamma0, beta0, sigma_b20, pi_beta0);
    Gamma0 = as<arma::mat>(gammaBeta0["Gamma"]);
    beta0 = as<arma::mat>(gammaBeta0["beta"]);
    gammaBeta1 = sample_gammaBeta_cpp(N1, M, K, Z1, G1, Gamma1, beta1, sigma_b21, pi_beta1);
    Gamma1 = as<arma::mat>(gammaBeta1["Gamma"]);
    beta1 = as<arma::mat>(gammaBeta1["beta"]);

    pi_beta0 = sample_pi_beta_cpp(M, K, Gamma0, prior_pibeta);
    sigma_b20 = sample_sigma_b2_cpp(M, Gamma0, beta0, prior_sigma2b);
    pi_beta1 = sample_pi_beta_cpp(M, K, Gamma1, prior_pibeta);
    sigma_b21 = sample_sigma_b2_cpp(M, Gamma1, beta1, prior_sigma2b);

    Z0 = sample_Z_cpp(N0, K, Y0, F, W, G0, beta0, psi);
    Z1 = sample_Z_cpp(N1, K, Y1, F, W, G1, beta1, psi);
    Z.rows(index0) = Z0;
    Z.rows(index1) = Z1;

    psi = sample_psi_cpp(N,P,Y,Z,F,W,prior_psi);
    W = sample_W_cpp(P, K, Y, Z, F, W, psi, sigma_w2, c2);
    F = sample_F_cpp(P, K,  W, pi_vec, sigma_w2, c2);
    pi_vec = sample_pi_cpp(P, K, F, prior_pi);
    sigma_w2 = sample_sigma_w2_cpp(K, P, F, W, prior_sigma2w, c2);
    c2 = sample_c2_cpp(K, P, F, W, sigma_w2, prior_c);

    // Store samples throughout iterations:
    Gamma0_mtx.slice(iter) = Gamma0;
    beta0_mtx.slice(iter) = beta0;
    pi_beta0_mtx.col(iter) = pi_beta0;

    Gamma1_mtx.slice(iter) = Gamma1;
    beta1_mtx.slice(iter) = beta1;
    pi_beta1_mtx.col(iter) = pi_beta1;

    Z_mtx.slice(iter) = Z;
    F_mtx.slice(iter) = F;
    W_mtx.slice(iter) = W;
    pi_mtx.col(iter) = pi_vec;
    sigma_w2_mtx.col(iter) = sigma_w2;
    c2_mtx.col(iter) = c2;

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
                             Named("Gamma0") = Gamma0, Named("Gamma1") = Gamma1,
                             Named("beta0") = beta0, Named("beta1") = beta1,
                             Named("pi_vec") = pi_vec,
                             Named("pi_beta0") = pi_beta0, Named("pi_beta1") = pi_beta1,
                             Named("psi") = psi,
                             Named("sigma_w2") = sigma_w2,
                             Named("sigma_b20") = sigma_b20, Named("sigma_b21") = sigma_b21,
                             Named("c2") = c2,
                             Named("niters") = niter, Named("average_niters") = ave_niter);
  // Compute the posterior means:
  List pm_list;
  pm_list = compute_posterior_mean_2groups_cpp(Gamma0_mtx, beta0_mtx, pi_beta0_mtx,
                                               Gamma1_mtx, beta1_mtx, pi_beta1_mtx,
                                               Z_mtx, F_mtx, W_mtx, pi_mtx,
                                               sigma_w2_mtx, c2_mtx,
                                               niter, ave_niter, prior_type);
  // Compute local false sign rate for each perturbation-gene pair and each of the sample groups:
  arma::mat lfsr0_mat(P,M);
  lfsr0_mat = compute_lfsr_cpp(beta0_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);
  arma::mat lfsr1_mat(P,M);
  lfsr1_mat = compute_lfsr_cpp(beta1_mtx, W_mtx, F_mtx, lfsr_niter, prior_type);

  if (return_samples) {
    // Save samples at each iteration on top of everything else
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr0") = lfsr0_mat,
                        Named("lfsr1") = lfsr1_mat,
                        Named("Z_samples") = Z_mtx,
                        Named("F_samples") = F_mtx,
                        Named("W_samples") = W_mtx,
                        Named("Gamma0_samples") = Gamma0_mtx,
                        Named("beta0_samples") = beta0_mtx,
                        Named("pi_beta0_samples") = pi_beta0_mtx,
                        Named("Gamma1_samples") = Gamma1_mtx,
                        Named("beta1_samples") = beta1_mtx,
                        Named("pi_beta1_samples") = pi_beta1_mtx,
                        Named("pi_samples") = pi_mtx,
                        Named("sigma_w2_samples") = sigma_w2_mtx,
                        Named("c2_samples") = c2_mtx,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  } else {
    return List::create(Named("updates") = update_list,
                        Named("posterior_means") = pm_list,
                        Named("lfsr0") = lfsr0_mat,
                        Named("lfsr1") = lfsr1_mat,
                        Named("prior_params") = prior_params,
                        Named("prior_type") = prior_type);
  }
}
