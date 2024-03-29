#ifndef INCLUDE_SAMPLING_FUNS
#define INCLUDE_SAMPLING_FUNS

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec mvrnormArma(arma::vec mu, arma::mat sigma);

void sample_Z(int N,int K, arma::mat Y, arma::mat F, arma::mat W,
              arma::mat G, arma::mat beta, arma::vec psi, arma::mat& Z);

void sample_GammaBeta(int N, int M, int K,
                      arma::mat Z, arma::mat G,
                      arma::mat& Gamma, arma::mat& beta,
                      arma::vec sigma_b2, arma::vec pi_beta);

void sample_W(int P, int K, arma::mat Y, arma::mat Z,
              arma::mat F, arma::mat& W,
              arma::vec psi, arma::vec sigma_w2, arma::vec c2);

void sample_F(int P, int K, arma::mat W, arma::mat& F,
              arma::vec pi_vec, arma::vec sigma_w2, arma::vec c2);

void sample_FW_spike_slab(int N, int P, int K,
                          arma::mat Y, arma::mat Z,
                          arma::mat& F, arma::mat& W,
                          arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec);

void sample_psi(int N, int P,
                arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                arma::vec prior_psi, arma::vec& psi);

void sample_pi(int P, int K, arma::mat F, arma::vec prior_pi, arma::vec& pi_vec);

void sample_pi_beta(int M, int K,
                    arma::mat Gamma, arma::vec prior_pibeta,
                    arma::vec& pi_beta);

void sample_sigma_w2(int K, int P,
                     arma::mat F, arma::mat W,
                     arma::vec prior_sigma2w, arma::vec c2,
                     arma::vec& sigma_w2_vec);

void sample_sigma_w2_spike_slab(int K, arma::mat F, arma::mat W,
                                arma::vec prior_sigma2w, arma::vec& sigma_w2_vec);

void sample_c2(int K, int P,
               arma::mat F, arma::mat W,
               arma::vec sigma2w, arma::vec prior_c,
               arma::vec& c2);

void sample_sigma_b2(int M, arma::mat Gamma, arma::mat beta,
                     arma::vec prior_sigma2b, arma::vec& sigma_b2_vec);

void initialize_random(int K, arma::mat Y,
                       arma::mat& Z, arma::mat& F, arma::mat& W);

void initialize_svd(int K, arma::mat Y,
                    arma::mat& Z, arma::mat& F, arma::mat& W);

void initialize_given_Z(int K, arma::mat Y,
                        arma::mat& Z, arma::mat& F, arma::mat& W);

void initialize_GammaBeta(int M, int K, arma::mat G, arma::mat Z,
                          arma::mat& beta, arma::mat& Gamma);

Rcpp::List compute_posterior_mean_cpp(arma::cube Gamma_mtx, arma::cube beta_mtx,
                                      arma::mat pi_beta_mtx, arma::cube Z_mtx,
                                      arma::cube F_mtx, arma::cube W_mtx,
                                      arma::mat pi_mtx, arma::mat sigma_w2_mtx,
                                      arma::mat c2_mtx,
                                      int niter, int ave_niter,
                                      Rcpp::String prior_type);

Rcpp::List compute_posterior_mean_2groups_cpp(arma::cube Gamma0_mtx, arma::cube beta0_mtx,
                                              arma::mat pi_beta0_mtx, arma::cube Gamma1_mtx,
                                              arma::cube beta1_mtx, arma::mat pi_beta1_mtx,
                                              arma::cube Z_mtx, arma::cube F_mtx,
                                              arma::cube W_mtx, arma::mat pi_mtx,
                                              arma::mat sigma_w2_mtx, arma::mat c2_mtx,
                                              int niter, int ave_niter,
                                              Rcpp::String prior_type);

void compute_lfsr_cpp(arma::cube beta_mtx, arma::cube W_mtx, arma::cube F_mtx,
                      arma::mat& lfsr_mat, arma::mat& total_effect,
                      int use_niter, Rcpp::String prior_type);

arma::cube calibrate_beta_vs_negctrl(arma::cube beta_mtx, int neg_ctrl_index);

#endif
