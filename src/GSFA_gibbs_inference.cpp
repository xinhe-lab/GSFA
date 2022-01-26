#include <RcppArmadillo.h>
using namespace Rcpp;

// You can source this function into an R session using the Rcpp::sourceCpp function
// (or via the Source button on the editor toolbar).

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double Cquantile(arma::vec x, double q) {
  arma::vec y = x;
  y = sort(y);
  return y(floor(x.n_elem * q));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
  // Uses the Cholesky decomposition
  int ncols = sigma.n_cols;
  arma::vec randvec = arma::randn(ncols);
  arma::mat decomp_mat = arma::chol(sigma);
  return mu +  decomp_mat.t() * randvec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat sample_Z_cpp(int N,int K,arma::mat Y,arma::mat F,arma::mat W,
                       arma::mat G,arma::mat beta,arma::vec psi){
  arma::mat Z = arma::zeros<arma::mat>(N,K);
  // arma::mat W = F%U;
  arma::mat Psi_inv = diagmat(1/psi); // 1/(psi_x_vector+0.000001) when necessary
  arma::vec I(K);
  I.ones();
  arma::mat WTP = W.t() * Psi_inv;
  arma::mat Sigma = inv_sympd(WTP * W + diagmat(I));
  arma::vec mu_i(K);
  arma::vec tmp_i(K);
  for (int i=0; i<N; i++) {
    arma::colvec Y_i = arma::conv_to< arma::colvec >::from(Y.row(i));
    tmp_i = (G.row(i) * beta).t();
    mu_i = Sigma * (WTP * Y_i + tmp_i);
    Z.row(i) = arma::conv_to< arma::rowvec >::from(mvnrnd(mu_i, Sigma));
  }
  // mvnrnd(mu,sigma,n) returns n column vectors
  return Z;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List sample_gammaBeta_cpp(int N, int M, int K, arma::mat Z, arma::mat G,
                          arma::mat Gamma, arma::mat beta,
                          arma::vec sigma_b2, arma::vec pi_beta){
  arma::mat log_qgamma = arma::zeros<arma::mat>(M,K);
  arma::mat mu = arma::zeros<arma::mat>(M,K);
  arma::mat L = arma::zeros<arma::mat>(M,K);
  arma::rowvec sum_G2 = sum(square(G),0);

  for (int k=0; k<K; k++) {
    for (int m=0; m<M; m++) {
      L(m,k) = 1.0 / (sum_G2(m) + 1.0/sigma_b2(m));
      double sum_RG = 0;
      for (int i=0; i<N; i++){
        sum_RG = sum_RG + G(i,m) * (Z(i,k) - dot(G.row(i), beta.col(k)) + G(i,m)*beta(m,k));
      }
      mu(m,k) = L(m,k) * sum_RG;

      log_qgamma(m,k) = std::log(L(m,k)/sigma_b2(m)) / 2 +
        mu(m,k) * mu(m,k) / (L(m,k) * 2) +
        std::log(pi_beta(m)) - std::log(1-pi_beta(m));
      double qgamma = 1.0 / (std::exp(-log_qgamma(m,k)) + 1);

      Gamma(m,k) = R::rbinom(1, qgamma);
      if (Gamma(m,k) == 1){
        beta(m,k) = R::rnorm(mu(m,k), sqrt(L(m,k)));
      } else {
        beta(m,k) = 0;
      }
    }
  }
  return List::create(Named("Gamma") = Gamma, Named("beta") = beta);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat sample_W_cpp(int P, int K, arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                       arma::vec psi, arma::vec sigma_w2, arma::vec c2){

  arma::mat log_qF = arma::zeros<arma::mat>(P,K);
  arma::mat lambda = arma::zeros<arma::mat>(P,K);
  arma::mat nu = arma::zeros<arma::mat>(P,K);
  arma::mat ZTZ = Z.t()*Z;

  for (int j=0; j<P; j++) {
    arma::vec diag_vec = arma::zeros<arma::vec>(K);
    for (int k=0; k<K; k++) {
      diag_vec(k) = psi(j) / (sigma_w2(k) * (F(j,k) + (1-F(j,k))*c2(k)));
    }
    arma::mat D_j = diagmat(diag_vec);
    arma::mat Lambda_j = inv_sympd(ZTZ + D_j);
    arma::vec nu_j = Lambda_j * Z.t() * Y.col(j);
    W.row(j) = arma::conv_to< arma::rowvec >::from(mvnrnd(nu_j, Lambda_j * psi(j)));
  }
  return W;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat sample_F_cpp(int P, int K, arma::mat W, arma::vec pi_vec, arma::vec sigma_w2, arma::vec c2){
  arma::mat F = arma::zeros<arma::mat>(P,K);
  for (int j=0; j<P; j++) {
    for (int k=0; k<K; k++) {
      double ratio = pi_vec(k) / (1-pi_vec(k)) * sqrt(c2(k)) *
        arma::trunc_exp(pow(W(j,k), 2)/(2*sigma_w2(k)) * (1.0/c2(k) - 1));
      if (arma::is_finite(ratio)){
        F(j,k) = R::rbinom(1, ratio/(1+ratio));
      } else {
        F(j,k) = 1;
      }
    }
  }
  return F;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List sample_FW_spike_slab_cpp(int N, int P, int K, arma::mat Y, arma::mat Z, arma::mat F, arma::mat W,
                              arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec){

  arma::mat log_qF = arma::zeros<arma::mat>(P,K);
  arma::mat lambda = arma::zeros<arma::mat>(P,K);
  arma::mat nu = arma::zeros<arma::mat>(P,K);
  arma::rowvec sum_Z2 = sum(square(Z),0);

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
        // U(j,k) = R::rnorm(nu(j,k),sqrt(lambda(j,k)));
      } else {
        W(j,k) = 0;
        // U(j,k) = R::rnorm(0,sqrt(sigma_w2(k)));
      }
    }
  }
  return List::create(Named("F") = F, Named("W") = W);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_psi_cpp(int N,int P,arma::mat Y,arma::mat Z,arma::mat F,arma::mat W,arma::vec prior_psi){
  arma::vec psi = arma::zeros<arma::vec>(P);
  // arma::mat W = F%U;
  arma::mat E = Y-Z*W.t();
  arma::rowvec sum_E2 = sum(E%E,0);
  double a,b;
  a = prior_psi(0) + N/2.0;
  for (int j=0; j<P; j++) {
    b = prior_psi(1) + sum_E2(j)/2.0;
    psi(j) = 1.0/R::rgamma(a, 1.0/b);
  }
  return psi;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_pi_cpp(int P, int K, arma::mat F, arma::vec prior_pi){
  arma::vec pi_vec = arma::zeros<arma::vec>(K);
  arma::rowvec sum_F = sum(F,0); // colsum returns a row vector!
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_pi(0)*prior_pi(1) + sum_F(k);
    b = prior_pi(0)*(1-prior_pi(1)) + P - sum_F(k);
    pi_vec(k) = R::rbeta(a, b);
  }
  return pi_vec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_pi_beta_cpp(int M, int K, arma::mat Gamma, arma::vec prior_pibeta){
  arma::vec pi_beta = arma::zeros<arma::vec>(M);
  arma::colvec sum_Gamma = sum(Gamma,1); // rowsum returns a column vector!
  double a,b;
  for (int m=0; m<M; m++) {
    a = prior_pibeta(0)*prior_pibeta(1) + sum_Gamma(m);
    b = prior_pibeta(0)*(1-prior_pibeta(1)) + K - sum_Gamma(m);
    pi_beta(m) = R::rbeta(a, b);
  }
  return pi_beta;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_sigma_w2_cpp(int K, int P, arma::mat F, arma::mat W, arma::vec prior_sigma2w, arma::vec c2){
  arma::vec sigma_w2_vec = arma::zeros<arma::vec>(K);
  double a,b;
  for (int k=0; k<K; k++) {
    double sum_W2 = sum(W.col(k)%W.col(k) / (F.col(k)+(1-F.col(k))*c2(k)));
    a = prior_sigma2w(0) + P/2.0;
    b = prior_sigma2w(1) + sum_W2/2.0;
    sigma_w2_vec(k) = 1.0/R::rgamma(a, 1.0/b);
  }
  return sigma_w2_vec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_sigma_w2_spike_slab_cpp(int K, arma::mat F, arma::mat W, arma::vec prior_sigma2w){
  arma::vec sigma_w2_vec = arma::zeros<arma::vec>(K);
  // arma::mat W = F%U;
  arma::rowvec sum_W2 = sum(W%W,0); // colsum
  arma::rowvec sum_F = sum(F,0);
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_sigma2w(0) + sum_F(k)/2.0;
    b = prior_sigma2w(1) + sum_W2(k)/2.0;
    sigma_w2_vec(k) = 1.0/R::rgamma(a, 1.0/b);
  }
  return sigma_w2_vec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_c2_cpp(int K, int P, arma::mat F, arma::mat W, arma::vec sigma2w, arma::vec prior_c){
  arma::vec c2 = arma::zeros<arma::vec>(K);
  arma::rowvec sum_F = sum(F,0);
  arma::rowvec sum_W2 = sum(W%W%(1-F), 0); // colsum
  double a,b;
  for (int k=0; k<K; k++) {
    a = prior_c(0) + (P-sum_F(k))/2.0;
    b = prior_c(1) + sum_W2(k)/(sigma2w(k)*2.0);
    c2(k) = 1.0/R::rgamma(a, 1.0/b);
  }
  return c2;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec sample_sigma_b2_cpp(int M, arma::mat Gamma, arma::mat beta, arma::vec prior_sigma2b){
  arma::vec sigma_b2_vec = arma::zeros<arma::vec>(M);
  arma::colvec sum_beta2 = sum(beta%beta,1); // rowsum
  arma::colvec sum_Gamma = sum(Gamma,1);
  double a,b;
  for (int m=0; m<M; m++) {
    a = prior_sigma2b(0) + sum_Gamma(m)/2.0;
    b = prior_sigma2b(1) + sum_beta2(m)/2.0;
    sigma_b2_vec(m) = 1.0/R::rgamma(a, 1.0/b);
  }
  return sigma_b2_vec;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List initialize_random(int K, arma::mat Y){
  int P = Y.n_cols;
  arma::mat init_U = arma::randn(P,K);
  init_U = init_U * 2.0;

  arma::vec tmp = Rcpp::rbinom(P*K,1,0.1);
  arma::mat init_F = arma::mat(tmp);
  init_F.reshape(P,K);

  arma::mat init_W = init_F % init_U;
  arma::mat init_Z = Y * init_W * inv(init_W.t()*init_W);
  return List::create(Named("W") = init_W, Named("F") = init_F, Named("Z") = init_Z);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List initialize_svd(int K, arma::mat Y){
  arma::mat svd_U;
  arma::vec d_vec;
  arma::mat svd_V;
  bool SVDsuccess = false;
  while(SVDsuccess == false) {
    SVDsuccess = svd(svd_U, d_vec, svd_V, Y);
    if(SVDsuccess == false){
      Y += 1e-4;
    }
  }
  arma::mat matd = arma::zeros<arma::mat>(K, K);
  matd.diag() = sqrt(d_vec(arma::span(0, K-1)));
  arma::mat init_Z = svd_U.cols(0, K-1)*matd;
  init_Z = (init_Z - mean(vectorise(init_Z))) / stddev(vectorise(init_Z));
  arma::mat init_U = Y.t() * init_Z * inv(init_Z.t()*init_Z);
  // Impose sparsity on W:
  double threshold = Cquantile(vectorise(abs(init_U)), 0.2);
  arma::umat init_bool = abs(init_U) > threshold;
  arma::mat init_F = arma::conv_to<arma::mat>::from(init_bool);
  arma::mat init_W = init_U % init_F;
  return List::create(Named("W") = init_W, Named("F") = init_F, Named("Z") = init_Z);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List initialize_given_Z(int K, arma::mat Y, arma::mat init_Z){
  if(K != init_Z.n_cols){
    throw("Number of factors provided in init_Z does not match K!");
  }
  arma::mat init_U = Y.t() * init_Z * inv(init_Z.t()*init_Z);
  // Impose sparsity on W:
  double threshold = Cquantile(vectorise(abs(init_U)), 0.2);
  arma::umat init_bool = abs(init_U) > threshold;
  arma::mat init_F = arma::conv_to<arma::mat>::from(init_bool);
  arma::mat init_W = init_U % init_F;
  return List::create(Named("W") = init_W, Named("F") = init_F, Named("Z") = init_Z);
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

  arma::mat Gamma_pm = arma::zeros<arma::mat>(M,K);
  arma::mat beta_pm = arma::zeros<arma::mat>(M,K);
  arma::vec pi_beta_pm = arma::zeros<arma::vec>(M);
  arma::mat Z_pm = arma::zeros<arma::mat>(N,K);
  arma::mat F_pm = arma::zeros<arma::mat>(P,K);
  arma::mat W_pm = arma::zeros<arma::mat>(P,K);
  arma::vec pi_pm = arma::zeros<arma::vec>(K);
  arma::vec sigma_w2_pm = arma::zeros<arma::vec>(K);
  arma::vec c2_pm = arma::zeros<arma::vec>(K);

  if (niter >= ave_niter){
    int start_iter = niter-ave_niter+1;
    arma::cube Gamma_slice = Gamma_mtx(arma::span(0, M-1), arma::span(0, K-1), arma::span(start_iter, niter));
    Gamma_pm = sum(Gamma_slice,2) / ave_niter;

    arma::cube beta_slice = beta_mtx(arma::span(0, M-1), arma::span(0, K-1), arma::span(start_iter, niter));
    beta_pm = sum(beta_slice,2) / ave_niter;

    arma::mat pi_beta_slice = pi_beta_mtx(arma::span(0, M-1), arma::span(start_iter, niter));
    pi_beta_pm = sum(pi_beta_slice,1) / ave_niter;

    arma::cube Z_slice = Z_mtx(arma::span(0, N-1), arma::span(0, K-1), arma::span(start_iter, niter));
    Z_pm = sum(Z_slice,2) / ave_niter;

    arma::cube F_slice = F_mtx(arma::span(0, P-1), arma::span(0, K-1), arma::span(start_iter, niter));
    F_pm = sum(F_slice,2) / ave_niter;

    arma::cube W_slice = W_mtx(arma::span(0, P-1), arma::span(0, K-1), arma::span(start_iter, niter));
    W_pm = sum(W_slice,2) / ave_niter;

    arma::mat pi_slice = pi_mtx(arma::span(0, K-1), arma::span(start_iter, niter));
    pi_pm = sum(pi_slice,1) / ave_niter;

    arma::mat sigma_w2_slice = sigma_w2_mtx(arma::span(0, K-1), arma::span(start_iter, niter));
    sigma_w2_pm = sum(sigma_w2_slice,1) / ave_niter;

    if (prior_type=="mixture_normal") {
      arma::mat c2_slices = c2_mtx(arma::span(0, K-1), arma::span(start_iter, niter));
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
arma::mat compute_lfsr_cpp(arma::cube beta_mtx, arma::cube W_mtx, arma::cube F_mtx,
                           int use_niter=100, String prior_type="mixture_normal"){
  int M = beta_mtx.n_rows;
  int K = beta_mtx.n_cols;
  int P = W_mtx.n_rows;
  int niter = W_mtx.n_slices;

  arma::mat lfsr_mat = arma::zeros<arma::mat>(P,M);
  arma::mat beta_m = arma::zeros<arma::mat>(K,niter);
  arma::mat W_j = arma::zeros<arma::mat>(K,niter);
  arma::mat F_j = arma::zeros<arma::mat>(K,niter);
  for (int m=0; m<M; m++){
    beta_m = beta_mtx.row(m);
    for (int j=0; j<P; j++){
      W_j = W_mtx.row(j);
      if (prior_type=="mixture_normal") {
        F_j = F_mtx.row(j);
        W_j = W_j % F_j; // Set W to zero if F=0, since sparsity is not enforced
      }
      arma::vec bw_prod(use_niter);
      for (int i=0; i<use_niter; i++){
        int slice_indx = niter-use_niter+i;
        bw_prod(i) = sum(W_j.col(slice_indx) % beta_m.col(slice_indx));
      }
      arma::vec sign_count(2);
      sign_count(0) = sum(bw_prod <= 0);
      sign_count(1) = sum(bw_prod >= 0);
      lfsr_mat(j,m) = sign_count.min() / double(use_niter);
    }
  }
  return lfsr_mat;
}

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
  arma::mat c2_mtx(K,niter+1, arma::fill::zeros); // = arma::zeros<arma::mat>(K,niter+1);

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
                                       niter, ave_niter, prior_type=prior_type);
  // Compute local false sign rate for each perturbation-gene pair:
  arma::mat lfsr_mat(P,M);
  lfsr_mat = compute_lfsr_cpp(beta_mtx, W_mtx, F_mtx, lfsr_niter, prior_type=prior_type);
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
                                       niter, ave_niter, prior_type=prior_type);
  // Compute local false sign rate for each perturbation-gene pair:
  arma::mat lfsr_mat(P,M);
  lfsr_mat = compute_lfsr_cpp(beta_mtx, W_mtx, F_mtx, lfsr_niter, prior_type=prior_type);
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
