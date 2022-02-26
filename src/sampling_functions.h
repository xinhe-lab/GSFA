#ifndef INCLUDE_SAMPLING_FUNS
#define INCLUDE_SAMPLING_FUNS

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec mvrnormArma(arma::vec mu, arma::mat sigma);

arma::mat sample_Z_cpp(int N, int K,
                       arma::mat Y, arma::mat F, arma::mat W,
                       arma::mat G, arma::mat beta, arma::vec psi);

void sample_Z(int N,int K, arma::mat Y, arma::mat F, arma::mat W,
              arma::mat G, arma::mat beta, arma::vec psi, arma::mat& Z);

Rcpp::List sample_gammaBeta_cpp(int N, int M, int K,
                          arma::mat Z, arma::mat G,
                          arma::mat Gamma, arma::mat beta,
                          arma::vec sigma_b2, arma::vec pi_beta);

void sample_GammaBeta(int N, int M, int K,
                      arma::mat Z, arma::mat G,
                      arma::mat& Gamma, arma::mat& beta,
                      arma::vec sigma_b2, arma::vec pi_beta);

arma::mat sample_W_cpp(int P, int K,
                       arma::mat Y, arma::mat Z,
                       arma::mat F, arma::mat W,
                       arma::vec psi, arma::vec sigma_w2, arma::vec c2);

void sample_W(int P, int K, arma::mat Y, arma::mat Z,
              arma::mat F, arma::mat& W,
              arma::vec psi, arma::vec sigma_w2, arma::vec c2);

arma::mat sample_F_cpp(int P, int K,
                       arma::mat W, arma::vec pi_vec,
                       arma::vec sigma_w2, arma::vec c2);

void sample_F(int P, int K, arma::mat W, arma::mat& F,
              arma::vec pi_vec, arma::vec sigma_w2, arma::vec c2);

Rcpp::List sample_FW_spike_slab_cpp(int N, int P, int K,
                              arma::mat Y, arma::mat Z,
                              arma::mat F, arma::mat W,
                              arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec);

void sample_FW_spike_slab(int N, int P, int K,
                          arma::mat Y, arma::mat Z,
                          arma::mat& F, arma::mat& W,
                          arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec);

arma::vec sample_psi_cpp(int N, int P,
                         arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                         arma::vec prior_psi);

void sample_psi(int N, int P,
                arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                arma::vec prior_psi, arma::vec& psi);

arma::vec sample_pi_cpp(int P, int K, arma::mat F, arma::vec prior_pi);

void sample_pi(int P, int K, arma::mat F, arma::vec prior_pi, arma::vec& pi_vec);

arma::vec sample_pi_beta_cpp(int M, int K,
                             arma::mat Gamma, arma::vec prior_pibeta);

void sample_pi_beta(int M, int K,
                    arma::mat Gamma, arma::vec prior_pibeta,
                    arma::vec& pi_beta);

arma::vec sample_sigma_w2_cpp(int K, int P,
                              arma::mat F, arma::mat W,
                              arma::vec prior_sigma2w, arma::vec c2);

void sample_sigma_w2(int K, int P,
                     arma::mat F, arma::mat W,
                     arma::vec prior_sigma2w, arma::vec c2,
                     arma::vec& sigma_w2_vec);

arma::vec sample_sigma_w2_spike_slab_cpp(int K, arma::mat F, arma::mat W,
                                         arma::vec prior_sigma2w);

void sample_sigma_w2_spike_slab(int K, arma::mat F, arma::mat W,
                                arma::vec prior_sigma2w, arma::vec& sigma_w2_vec);

arma::vec sample_c2_cpp(int K, int P,
                        arma::mat F, arma::mat W,
                        arma::vec sigma2w, arma::vec prior_c);

void sample_c2(int K, int P,
               arma::mat F, arma::mat W,
               arma::vec sigma2w, arma::vec prior_c,
               arma::vec& c2);

arma::vec sample_sigma_b2_cpp(int M, arma::mat Gamma, arma::mat beta,
                              arma::vec prior_sigma2b);

void sample_sigma_b2(int M, arma::mat Gamma, arma::mat beta,
                     arma::vec prior_sigma2b, arma::vec& sigma_b2_vec);

void initialize_random(int K, arma::mat Y,
                       arma::mat& Z, arma::mat& F, arma::mat& W);

void initialize_svd(int K, arma::mat Y,
                    arma::mat& Z, arma::mat& F, arma::mat& W);

void initialize_given_Z(int K, arma::mat Y,
                        arma::mat& Z, arma::mat& F, arma::mat& W);

Rcpp::List initialize_gammaBeta(int M, int K, arma::mat G, arma::mat Z);

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

arma::mat compute_lfsr_cpp(arma::cube beta_mtx, arma::cube W_mtx, arma::cube F_mtx,
                           int use_niter, Rcpp::String prior_type);

#endif
