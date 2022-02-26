#include "sampling_functions.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void move_seed_rbinom();

double Cquantile(arma::vec x, double q);

// FUNCTION DEFINITIONS
// ---------------------

void move_seed_rbinom(){
  double tmp = R::rbinom(1, 0.5);
}

double Cquantile(arma::vec x, double q) {
  vec y = x;
  y = sort(y);
  return y(floor(x.n_elem * q));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
  // Uses the Cholesky decomposition
  int ncols = sigma.n_cols;
  vec randvec = randn(ncols);
  mat decomp_mat = chol(sigma);
  return mu + decomp_mat.t() * randvec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_Z(int N,int K, arma::mat Y, arma::mat F, arma::mat W,
              arma::mat G, arma::mat beta, arma::vec psi, arma::mat& Z){
  mat Psi_inv = diagmat(1/psi); // 1/(psi_x_vector+0.000001) when necessary
  vec I(K);
  I.ones();
  mat WTP = W.t() * Psi_inv;
  mat Sigma = inv_sympd(WTP * W + diagmat(I));
  vec mu_i(K);
  vec tmp_i(K);
  for (int i=0; i<N; i++) {
    colvec Y_i = conv_to< colvec >::from(Y.row(i));
    tmp_i = (G.row(i) * beta).t();
    mu_i = Sigma * (WTP * Y_i + tmp_i);
    Z.row(i) = conv_to< rowvec >::from(mvnrnd(mu_i, Sigma));
  }
  // mvnrnd(mu,sigma,n) returns n column vectors
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_GammaBeta(int N, int M, int K,
                      arma::mat Z, arma::mat G,
                      arma::mat& Gamma, arma::mat& beta,
                      arma::vec sigma_b2, arma::vec pi_beta){
  double log_qgamma;
  double qgamma;
  mat mu = zeros<mat>(M,K);
  mat L = zeros<mat>(M,K);
  rowvec sum_G2 = sum(square(G),0);

  for (int k=0; k<K; k++) {
    for (int m=0; m<M; m++) {
      L(m,k) = 1.0 / (sum_G2(m) + 1.0/sigma_b2(m));
      double sum_RG = 0;
      for (int i=0; i<N; i++){
        sum_RG = sum_RG + G(i,m) * (Z(i,k) - dot(G.row(i), beta.col(k)) + G(i,m)*beta(m,k));
      }
      mu(m,k) = L(m,k) * sum_RG;

      log_qgamma = std::log(L(m,k)/sigma_b2(m)) / 2 +
        mu(m,k) * mu(m,k) / (L(m,k) * 2) +
        std::log(pi_beta(m)) - std::log(1-pi_beta(m));
      if (log_qgamma > 30){
        Gamma(m,k) = 1;
        // R::rbinom(1, p) does not perform sampling and move the random seed
        // forward when p > 1-E-17 (the exact boundary depends on the machine)
        // Need to manually add a rbinom sampling step to move the random
        // seed along and ensure reproducibility
        move_seed_rbinom();
      } else {
        qgamma = 1.0 / (std::exp(-log_qgamma) + 1);
        Gamma(m,k) = R::rbinom(1, qgamma);
      }

      if (Gamma(m,k) == 1){
        beta(m,k) = R::rnorm(mu(m,k), sqrt(L(m,k)));
      } else {
        beta(m,k) = 0;
      }
    }
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_W(int P, int K, arma::mat Y, arma::mat Z,
              arma::mat F, arma::mat& W,
              arma::vec psi, arma::vec sigma_w2, arma::vec c2){

  mat log_qF = zeros<mat>(P,K);
  mat lambda = zeros<mat>(P,K);
  mat nu = zeros<mat>(P,K);
  mat ZTZ = Z.t()*Z;

  for (int j=0; j<P; j++) {
    vec diag_vec = zeros<vec>(K);
    for (int k=0; k<K; k++) {
      diag_vec(k) = psi(j) / (sigma_w2(k) * (F(j,k) + (1-F(j,k))*c2(k)));
    }
    mat D_j = diagmat(diag_vec);
    mat Lambda_j = inv_sympd(ZTZ + D_j);
    vec nu_j = Lambda_j * Z.t() * Y.col(j);
    W.row(j) = conv_to< rowvec >::from(mvnrnd(nu_j, Lambda_j * psi(j)));
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_F(int P, int K, arma::mat W, arma::mat& F,
              arma::vec pi_vec, arma::vec sigma_w2, arma::vec c2){
  double log_ratio;
  double ratio;
  for (int j=0; j<P; j++) {
    for (int k=0; k<K; k++) {
      log_ratio = std::log(pi_vec(k)) - std::log(1-pi_vec(k)) +
        std::log(c2(k)) / 2. +
        pow(W(j,k), 2)/(2*sigma_w2(k)) * (1.0/c2(k) - 1);
      if (log_ratio > 30){
        F(j,k) = 1;
        // R::rbinom(1, p) does not perform sampling and move the random seed
        // forward when p > 1-E-17 (the exact boundary depends on the machine)
        // Need to manually add a rbinom sampling step to move the random
        // seed along and ensure reproducibility
        move_seed_rbinom();
      } else {
        ratio = 1.0 / (std::exp(-log_ratio) + 1);
        F(j,k) = R::rbinom(1, ratio);
      }
    }
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_FW_spike_slab(int N, int P, int K,
                          arma::mat Y, arma::mat Z,
                          arma::mat& F, arma::mat& W,
                          arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec){

  mat log_qF = zeros<mat>(P,K);
  mat lambda = zeros<mat>(P,K);
  mat nu = zeros<mat>(P,K);
  rowvec sum_Z2 = sum(square(Z),0);

  for (int j=0; j<P; j++) {
    for (int k=0; k<K; k++) {
      lambda(j,k) = 1.0 / (sum_Z2(k)/psi(j) + 1/sigma_w2(k));
      double sum_ZB = 0;
      for (int i=0; i<N; i++){
        sum_ZB = sum_ZB + Z(i,k) * (Y(i,j) - dot(Z.row(i), W.row(j)) + Z(i,k)*W(j,k));
      }
      nu(j,k) = lambda(j,k) * sum_ZB / psi(j);

      log_qF(j,k) = std::log(lambda(j,k)/sigma_w2(k))/2 + nu(j,k)*nu(j,k)/(lambda(j,k)*2) +
        std::log(pi_vec(k)) - std::log(1-pi_vec(k));
      double qF = 1.0 / (std::exp(-log_qF(j,k)) + 1);

      F(j,k) = R::rbinom(1, qF);
      if (F(j,k)==1){
        W(j,k) = R::rnorm(nu(j,k), sqrt(lambda(j,k)));
      } else {
        W(j,k) = 0;
      }
    }
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_psi(int N, int P,
                arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                arma::vec prior_psi, arma::vec& psi){
  mat E = Y-Z*W.t();
  rowvec sum_E2 = sum(E%E,0);
  double a,b;
  a = prior_psi(0) + N/2.0;
  for (int j=0; j<P; j++) {
    b = prior_psi(1) + sum_E2(j)/2.0;
    psi(j) = 1.0/R::rgamma(a, 1.0/b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_pi(int P, int K, arma::mat F, arma::vec prior_pi, arma::vec& pi_vec){
  rowvec sum_F = sum(F,0); // col-wise sum returns a row vector
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_pi(0) * prior_pi(1) + sum_F(k);
    b = prior_pi(0) * (1-prior_pi(1)) + P - sum_F(k);
    pi_vec(k) = R::rbeta(a, b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_pi_beta(int M, int K,
                    arma::mat Gamma, arma::vec prior_pibeta,
                    arma::vec& pi_beta){
  colvec sum_Gamma = sum(Gamma,1); // row-wise sum returns a column vector
  double a,b;
  for (int m=0; m<M; m++) {
    a = prior_pibeta(0)*prior_pibeta(1) + sum_Gamma(m);
    b = prior_pibeta(0)*(1-prior_pibeta(1)) + K - sum_Gamma(m);
    pi_beta(m) = R::rbeta(a, b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_sigma_w2(int K, int P,
                     arma::mat F, arma::mat W,
                     arma::vec prior_sigma2w, arma::vec c2,
                     arma::vec& sigma_w2_vec){
  double a,b;
  for (int k=0; k<K; k++) {
    double sum_W2 = sum(W.col(k)%W.col(k) / (F.col(k)+(1-F.col(k))*c2(k)));
    a = prior_sigma2w(0) + P/2.0;
    b = prior_sigma2w(1) + sum_W2/2.0;
    sigma_w2_vec(k) = 1.0/R::rgamma(a, 1.0/b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_sigma_w2_spike_slab(int K, arma::mat F, arma::mat W,
                                arma::vec prior_sigma2w, arma::vec& sigma_w2_vec){
  rowvec sum_W2 = sum(W%W,0); // column-wise sum
  rowvec sum_F = sum(F,0);
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_sigma2w(0) + sum_F(k)/2.0;
    b = prior_sigma2w(1) + sum_W2(k)/2.0;
    sigma_w2_vec(k) = 1.0/R::rgamma(a, 1.0/b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_c2(int K, int P,
               arma::mat F, arma::mat W,
               arma::vec sigma2w, arma::vec prior_c,
               arma::vec& c2){
  rowvec sum_F = sum(F,0);
  rowvec sum_W2 = sum(W%W%(1-F), 0); // column-wise sum
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_c(0) + (P-sum_F(k))/2.0;
    b = prior_c(1) + sum_W2(k)/(sigma2w(k)*2.0);
    c2(k) = 1.0/R::rgamma(a, 1.0/b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void sample_sigma_b2(int M, arma::mat Gamma, arma::mat beta,
                     arma::vec prior_sigma2b, arma::vec& sigma_b2_vec){
  colvec sum_beta2 = sum(beta%beta,1); // row-wise sum
  colvec sum_Gamma = sum(Gamma,1);
  double a,b;
  for (int m=0; m<M; m++) {
    a = prior_sigma2b(0) + sum_Gamma(m)/2.0;
    b = prior_sigma2b(1) + sum_beta2(m)/2.0;
    sigma_b2_vec(m) = 1.0/R::rgamma(a, 1.0/b);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void initialize_random(int K, arma::mat Y,
                       arma::mat& Z, arma::mat& F, arma::mat& W){
  Rprintf("Initializing Z and W with random values.\n");
  int P = Y.n_cols;
  mat U = randn(P,K);
  U = U * 2.0;

  vec tmp = Rcpp::rbinom(P*K,1,0.1);
  F = mat(tmp);
  F.reshape(P,K);

  W = F % U;
  Z = Y * W * inv(W.t()*W);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void initialize_svd(int K, arma::mat Y,
                    arma::mat& Z, arma::mat& F, arma::mat& W){
  Rprintf("Initializing Z and W with SVD.\n");
  mat svd_U;
  vec d_vec;
  mat svd_V;
  bool SVDsuccess = false;
  while(SVDsuccess == false) {
    SVDsuccess = svd(svd_U, d_vec, svd_V, Y);
    if(SVDsuccess == false){
      Y += 1e-4;
    }
  }
  mat matd = zeros<mat>(K, K);
  matd.diag() = sqrt(d_vec(span(0, K-1)));
  Z = svd_U.cols(0, K-1)*matd;
  Z = (Z - mean(vectorise(Z))) / stddev(vectorise(Z));
  mat U = Y.t() * Z * inv(Z.t()*Z);
  // Impose sparsity on W:
  double threshold = Cquantile(vectorise(abs(U)), 0.2);
  umat init_bool = abs(U) > threshold;
  F = conv_to<mat>::from(init_bool);
  W = U % F;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void initialize_given_Z(int K, arma::mat Y,
                        arma::mat& Z, arma::mat& F, arma::mat& W){
  Rprintf("Initializing Z and W with user provided Z.\n");
  if(K != Z.n_cols){
    throw("The number of factors provided for Z does not match K!");
  }
  mat U = Y.t() * Z * inv(Z.t()*Z);
  // Impose sparsity on W:
  double threshold = Cquantile(vectorise(abs(U)), 0.2);
  umat init_bool = abs(U) > threshold;
  F = conv_to<mat>::from(init_bool);
  W = U % F;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void initialize_GammaBeta(int M, int K, arma::mat G, arma::mat Z,
                          arma::mat& beta, arma::mat& Gamma){
  vec tmp_beta;
  for (int m=0; m<M; m++){
    vec G_col = G.col(m);
    vec tmp_beta = Z.t() * G_col / (sum(G_col%G_col)+1e-4);
    beta.row(m) = conv_to< rowvec >::from(tmp_beta);
  }
  umat init_bool = abs(beta) > Cquantile(vectorise(abs(beta)), 0.5);
  Gamma = conv_to< mat >::from(init_bool);
  beta = beta % Gamma;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List compute_posterior_mean_cpp(arma::cube Gamma_mtx, arma::cube beta_mtx,
                                arma::mat pi_beta_mtx, arma::cube Z_mtx,
                                arma::cube F_mtx, arma::cube W_mtx,
                                arma::mat pi_mtx, arma::mat sigma_w2_mtx,
                                arma::mat c2_mtx,
                                int niter=200, int ave_niter=100,
                                String prior_type="mixture_normal"){
  int N = Z_mtx.n_rows;
  int P = F_mtx.n_rows;
  int M = beta_mtx.n_rows;
  int K = Z_mtx.n_cols;

  mat Gamma_pm = zeros<mat>(M,K);
  mat beta_pm = zeros<mat>(M,K);
  vec pi_beta_pm = zeros<vec>(M);
  mat Z_pm = zeros<mat>(N,K);
  mat F_pm = zeros<mat>(P,K);
  mat W_pm = zeros<mat>(P,K);
  vec pi_pm = zeros<vec>(K);
  vec sigma_w2_pm = zeros<vec>(K);
  vec c2_pm = zeros<vec>(K);

  if (niter >= ave_niter){
    int start_iter = niter-ave_niter+1;
    cube Gamma_slice = Gamma_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    Gamma_pm = sum(Gamma_slice,2) / ave_niter;

    cube beta_slice = beta_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    beta_pm = sum(beta_slice,2) / ave_niter;

    mat pi_beta_slice = pi_beta_mtx(span(0, M-1), span(start_iter, niter));
    pi_beta_pm = sum(pi_beta_slice,1) / ave_niter;

    cube Z_slice = Z_mtx(span(0, N-1), span(0, K-1), span(start_iter, niter));
    Z_pm = sum(Z_slice,2) / ave_niter;

    cube F_slice = F_mtx(span(0, P-1), span(0, K-1), span(start_iter, niter));
    F_pm = sum(F_slice,2) / ave_niter;

    cube W_slice = W_mtx(span(0, P-1), span(0, K-1), span(start_iter, niter));
    W_pm = sum(W_slice,2) / ave_niter;

    mat pi_slice = pi_mtx(span(0, K-1), span(start_iter, niter));
    pi_pm = sum(pi_slice,1) / ave_niter;

    mat sigma_w2_slice = sigma_w2_mtx(span(0, K-1), span(start_iter, niter));
    sigma_w2_pm = sum(sigma_w2_slice,1) / ave_niter;

    if (prior_type=="mixture_normal") {
      mat c2_slices = c2_mtx(span(0, K-1), span(start_iter, niter));
      c2_pm = sum(c2_slices,1) / ave_niter;
    }
  } else {
    warning("Total number of iterations < specified number of iterations to average over, \
            returning the last samples as the estimates instead of posterior means over iterations.");
    Gamma_pm = Gamma_mtx.slice(niter);
    beta_pm = beta_mtx.slice(niter);
    pi_beta_pm = pi_beta_mtx.col(niter);
    Z_pm = Z_mtx.slice(niter);
    F_pm = F_mtx.slice(niter);
    W_pm = W_mtx.slice(niter);
    pi_pm = pi_mtx.col(niter);
    sigma_w2_pm = sigma_w2_mtx.col(niter);
    if (prior_type=="mixture_normal") {
      c2_pm = c2_mtx.col(niter);
    }
  }

  if (prior_type=="mixture_normal") {
    return List::create(Named("Z_pm") = Z_pm,
                        Named("F_pm") = F_pm,
                        Named("W_pm") = W_pm,
                        Named("Gamma_pm") = Gamma_pm,
                        Named("beta_pm") = beta_pm,
                        Named("pi_pm") = pi_pm,
                        Named("pi_beta_pm") = pi_beta_pm,
                        Named("sigma_w2_pm") = sigma_w2_pm,
                        Named("c2_pm") = c2_pm);
  } else {
    return List::create(Named("Z_pm") = Z_pm,
                        Named("F_pm") = F_pm,
                        Named("W_pm") = W_pm,
                        Named("Gamma_pm") = Gamma_pm,
                        Named("beta_pm") = beta_pm,
                        Named("pi_pm") = pi_pm,
                        Named("pi_beta_pm") = pi_beta_pm,
                        Named("sigma_w2_pm") = sigma_w2_pm);
  }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List compute_posterior_mean_2groups_cpp(arma::cube Gamma0_mtx, arma::cube beta0_mtx,
                                        arma::mat pi_beta0_mtx, arma::cube Gamma1_mtx,
                                        arma::cube beta1_mtx, arma::mat pi_beta1_mtx,
                                        arma::cube Z_mtx, arma::cube F_mtx,
                                        arma::cube W_mtx, arma::mat pi_mtx,
                                        arma::mat sigma_w2_mtx, arma::mat c2_mtx,
                                        int niter=200, int ave_niter=100,
                                        String prior_type="mixture_normal"){
  int N = Z_mtx.n_rows;
  int P = F_mtx.n_rows;
  int M = beta0_mtx.n_rows;
  int K = Z_mtx.n_cols;

  mat Gamma0_pm = zeros<mat>(M,K);
  mat beta0_pm = zeros<mat>(M,K);
  vec pi_beta0_pm = zeros<vec>(M);

  mat Gamma1_pm = zeros<mat>(M,K);
  mat beta1_pm = zeros<mat>(M,K);
  vec pi_beta1_pm = zeros<vec>(M);

  mat Z_pm = zeros<mat>(N,K);
  mat F_pm = zeros<mat>(P,K);
  mat W_pm = zeros<mat>(P,K);
  vec pi_pm = zeros<vec>(K);
  vec sigma_w2_pm = zeros<vec>(K);
  vec c2_pm = zeros<vec>(K);

  if (niter >= ave_niter){
    int start_iter = niter-ave_niter+1;
    cube Gamma0_slice = Gamma0_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    Gamma0_pm = sum(Gamma0_slice,2) / ave_niter;

    cube beta0_slice = beta0_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    beta0_pm = sum(beta0_slice,2) / ave_niter;

    mat pi_beta0_slice = pi_beta0_mtx(span(0, M-1), span(start_iter, niter));
    pi_beta0_pm = sum(pi_beta0_slice,1) / ave_niter;

    cube Gamma1_slice = Gamma1_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    Gamma1_pm = sum(Gamma1_slice,2) / ave_niter;

    cube beta1_slice = beta1_mtx(span(0, M-1), span(0, K-1), span(start_iter, niter));
    beta1_pm = sum(beta1_slice,2) / ave_niter;

    mat pi_beta1_slice = pi_beta1_mtx(span(0, M-1), span(start_iter, niter));
    pi_beta1_pm = sum(pi_beta1_slice,1) / ave_niter;

    cube Z_slice = Z_mtx(span(0, N-1), span(0, K-1), span(start_iter, niter));
    Z_pm = sum(Z_slice,2) / ave_niter;

    cube F_slice = F_mtx(span(0, P-1), span(0, K-1), span(start_iter, niter));
    F_pm = sum(F_slice,2) / ave_niter;

    cube W_slice = W_mtx(span(0, P-1), span(0, K-1), span(start_iter, niter));
    W_pm = sum(W_slice,2) / ave_niter;

    mat pi_slice = pi_mtx(span(0, K-1), span(start_iter, niter));
    pi_pm = sum(pi_slice,1) / ave_niter;

    mat sigma_w2_slice = sigma_w2_mtx(span(0, K-1), span(start_iter, niter));
    sigma_w2_pm = sum(sigma_w2_slice,1) / ave_niter;
    if (prior_type=="mixture_normal") {
      mat c2_slices = c2_mtx(span(0, K-1), span(start_iter, niter));
      c2_pm = sum(c2_slices,1) / ave_niter;
    }
  } else {
    warning("Total number of iterations < specified number of iterations to average over, \
            returning the last samples as the estimates instead of posterior means over iterations.");
    Gamma0_pm = Gamma0_mtx.slice(niter);
    beta0_pm = beta0_mtx.slice(niter);
    pi_beta0_pm = pi_beta0_mtx.col(niter);
    Gamma1_pm = Gamma1_mtx.slice(niter);
    beta1_pm = beta1_mtx.slice(niter);
    pi_beta1_pm = pi_beta1_mtx.col(niter);
    Z_pm = Z_mtx.slice(niter);
    F_pm = F_mtx.slice(niter);
    W_pm = W_mtx.slice(niter);
    pi_pm = pi_mtx.col(niter);
    sigma_w2_pm = sigma_w2_mtx.col(niter);
    if (prior_type=="mixture_normal") {
      c2_pm = c2_mtx.col(niter);
    }
  }

  if (prior_type=="mixture_normal") {
    return List::create(Named("Z_pm") = Z_pm,
                        Named("F_pm") = F_pm,
                        Named("W_pm") = W_pm,
                        Named("Gamma0_pm") = Gamma0_pm, Named("Gamma1_pm") = Gamma1_pm,
                        Named("beta0_pm") = beta0_pm, Named("beta1_pm") = beta1_pm,
                        Named("pi_pm") = pi_pm,
                        Named("pi_beta0_pm") = pi_beta0_pm, Named("pi_beta1_pm") = pi_beta1_pm,
                        Named("sigma_w2_pm") = sigma_w2_pm,
                        Named("c2_pm") = c2_pm);
  } else{
    return List::create(Named("Z_pm") = Z_pm,
                        Named("F_pm") = F_pm,
                        Named("W_pm") = W_pm,
                        Named("Gamma0_pm") = Gamma0_pm, Named("Gamma1_pm") = Gamma1_pm,
                        Named("beta0_pm") = beta0_pm, Named("beta1_pm") = beta1_pm,
                        Named("pi_pm") = pi_pm,
                        Named("pi_beta0_pm") = pi_beta0_pm, Named("pi_beta1_pm") = pi_beta1_pm,
                        Named("sigma_w2_pm") = sigma_w2_pm);
  }

}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat compute_lfsr_cpp(arma::cube beta_mtx, arma::cube W_mtx, arma::cube F_mtx,
                           int use_niter=100, String prior_type="mixture_normal"){
  int M = beta_mtx.n_rows;
  int K = beta_mtx.n_cols;
  int P = W_mtx.n_rows;
  int niter = W_mtx.n_slices;

  mat lfsr_mat = zeros<mat>(P,M);
  mat beta_m = zeros<mat>(K,niter);
  mat W_j = zeros<mat>(K,niter);
  mat F_j = zeros<mat>(K,niter);
  for (int m=0; m<M; m++){
    beta_m = beta_mtx.row(m);
    for (int j=0; j<P; j++){
      W_j = W_mtx.row(j);
      if (prior_type=="mixture_normal") {
        F_j = F_mtx.row(j);
        W_j = W_j % F_j; // Set W to zero if F=0 under the normal-mixture prior
      }
      vec bw_prod(use_niter);
      for (int i=0; i<use_niter; i++){
        int slice_indx = niter-use_niter+i;
        bw_prod(i) = sum(W_j.col(slice_indx) % beta_m.col(slice_indx));
      }
      vec sign_count(2);
      sign_count(0) = sum(bw_prod <= 0);
      sign_count(1) = sum(bw_prod >= 0);
      lfsr_mat(j,m) = sign_count.min() / double(use_niter);
    }
  }
  return lfsr_mat;
}
