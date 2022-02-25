// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gsfa_gibbs_cpp
List gsfa_gibbs_cpp(arma::mat Y, arma::mat G, int K, String prior_type, String initialize, double prior_s, double prior_r, double prior_sb, double prior_rb, double prior_gp, double prior_hp, double prior_gb, double prior_hb, double prior_gw, double prior_hw, double prior_gc, double prior_hc, int niter, int ave_niter, int lfsr_niter, bool verbose, bool return_samples);
RcppExport SEXP _GSFA_gsfa_gibbs_cpp(SEXP YSEXP, SEXP GSEXP, SEXP KSEXP, SEXP prior_typeSEXP, SEXP initializeSEXP, SEXP prior_sSEXP, SEXP prior_rSEXP, SEXP prior_sbSEXP, SEXP prior_rbSEXP, SEXP prior_gpSEXP, SEXP prior_hpSEXP, SEXP prior_gbSEXP, SEXP prior_hbSEXP, SEXP prior_gwSEXP, SEXP prior_hwSEXP, SEXP prior_gcSEXP, SEXP prior_hcSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP lfsr_niterSEXP, SEXP verboseSEXP, SEXP return_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< String >::type initialize(initializeSEXP);
    Rcpp::traits::input_parameter< double >::type prior_s(prior_sSEXP);
    Rcpp::traits::input_parameter< double >::type prior_r(prior_rSEXP);
    Rcpp::traits::input_parameter< double >::type prior_sb(prior_sbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_rb(prior_rbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gp(prior_gpSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hp(prior_hpSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gb(prior_gbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hb(prior_hbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gw(prior_gwSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hw(prior_hwSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gc(prior_gcSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hc(prior_hcSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< int >::type lfsr_niter(lfsr_niterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type return_samples(return_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(gsfa_gibbs_cpp(Y, G, K, prior_type, initialize, prior_s, prior_r, prior_sb, prior_rb, prior_gp, prior_hp, prior_gb, prior_hb, prior_gw, prior_hw, prior_gc, prior_hc, niter, ave_niter, lfsr_niter, verbose, return_samples));
    return rcpp_result_gen;
END_RCPP
}
// restart_gsfa_gibbs_cpp
List restart_gsfa_gibbs_cpp(arma::mat Y, arma::mat G, arma::mat Z, arma::mat F, arma::mat W, arma::mat Gamma, arma::mat beta, arma::vec pi_vec, arma::vec pi_beta, arma::vec psi, arma::vec sigma_w2, arma::vec sigma_b2, arma::vec c2, List prior_params, String prior_type, int niter, int ave_niter, int lfsr_niter, bool verbose, bool return_samples);
RcppExport SEXP _GSFA_restart_gsfa_gibbs_cpp(SEXP YSEXP, SEXP GSEXP, SEXP ZSEXP, SEXP FSEXP, SEXP WSEXP, SEXP GammaSEXP, SEXP betaSEXP, SEXP pi_vecSEXP, SEXP pi_betaSEXP, SEXP psiSEXP, SEXP sigma_w2SEXP, SEXP sigma_b2SEXP, SEXP c2SEXP, SEXP prior_paramsSEXP, SEXP prior_typeSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP lfsr_niterSEXP, SEXP verboseSEXP, SEXP return_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_beta(pi_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_w2(sigma_w2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_b2(sigma_b2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< List >::type prior_params(prior_paramsSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< int >::type lfsr_niter(lfsr_niterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type return_samples(return_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(restart_gsfa_gibbs_cpp(Y, G, Z, F, W, Gamma, beta, pi_vec, pi_beta, psi, sigma_w2, sigma_b2, c2, prior_params, prior_type, niter, ave_niter, lfsr_niter, verbose, return_samples));
    return rcpp_result_gen;
END_RCPP
}
// gsfa_gibbs_2groups_cpp
List gsfa_gibbs_2groups_cpp(arma::mat Y, arma::mat G, arma::vec group, int K, String prior_type, String initialize, Rcpp::Nullable<Rcpp::NumericMatrix> Z_given, double prior_s, double prior_r, double prior_sb, double prior_rb, double prior_gp, double prior_hp, double prior_gb, double prior_hb, double prior_gw, double prior_hw, double prior_gc, double prior_hc, int niter, int ave_niter, int lfsr_niter, bool verbose, bool return_samples);
RcppExport SEXP _GSFA_gsfa_gibbs_2groups_cpp(SEXP YSEXP, SEXP GSEXP, SEXP groupSEXP, SEXP KSEXP, SEXP prior_typeSEXP, SEXP initializeSEXP, SEXP Z_givenSEXP, SEXP prior_sSEXP, SEXP prior_rSEXP, SEXP prior_sbSEXP, SEXP prior_rbSEXP, SEXP prior_gpSEXP, SEXP prior_hpSEXP, SEXP prior_gbSEXP, SEXP prior_hbSEXP, SEXP prior_gwSEXP, SEXP prior_hwSEXP, SEXP prior_gcSEXP, SEXP prior_hcSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP lfsr_niterSEXP, SEXP verboseSEXP, SEXP return_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< String >::type initialize(initializeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Z_given(Z_givenSEXP);
    Rcpp::traits::input_parameter< double >::type prior_s(prior_sSEXP);
    Rcpp::traits::input_parameter< double >::type prior_r(prior_rSEXP);
    Rcpp::traits::input_parameter< double >::type prior_sb(prior_sbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_rb(prior_rbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gp(prior_gpSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hp(prior_hpSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gb(prior_gbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hb(prior_hbSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gw(prior_gwSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hw(prior_hwSEXP);
    Rcpp::traits::input_parameter< double >::type prior_gc(prior_gcSEXP);
    Rcpp::traits::input_parameter< double >::type prior_hc(prior_hcSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< int >::type lfsr_niter(lfsr_niterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type return_samples(return_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(gsfa_gibbs_2groups_cpp(Y, G, group, K, prior_type, initialize, Z_given, prior_s, prior_r, prior_sb, prior_rb, prior_gp, prior_hp, prior_gb, prior_hb, prior_gw, prior_hw, prior_gc, prior_hc, niter, ave_niter, lfsr_niter, verbose, return_samples));
    return rcpp_result_gen;
END_RCPP
}
// restart_gibbs_2groups_cpp
List restart_gibbs_2groups_cpp(arma::mat Y, arma::mat G, arma::vec group, arma::mat Z, arma::mat F, arma::mat W, arma::mat Gamma0, arma::mat Gamma1, arma::mat beta0, arma::mat beta1, arma::vec pi_vec, arma::vec pi_beta0, arma::vec pi_beta1, arma::vec psi, arma::vec sigma_w2, arma::vec sigma_b20, arma::vec sigma_b21, arma::vec c2, List prior_params, String prior_type, int niter, int ave_niter, int lfsr_niter, bool verbose, bool return_samples);
RcppExport SEXP _GSFA_restart_gibbs_2groups_cpp(SEXP YSEXP, SEXP GSEXP, SEXP groupSEXP, SEXP ZSEXP, SEXP FSEXP, SEXP WSEXP, SEXP Gamma0SEXP, SEXP Gamma1SEXP, SEXP beta0SEXP, SEXP beta1SEXP, SEXP pi_vecSEXP, SEXP pi_beta0SEXP, SEXP pi_beta1SEXP, SEXP psiSEXP, SEXP sigma_w2SEXP, SEXP sigma_b20SEXP, SEXP sigma_b21SEXP, SEXP c2SEXP, SEXP prior_paramsSEXP, SEXP prior_typeSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP lfsr_niterSEXP, SEXP verboseSEXP, SEXP return_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma0(Gamma0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma1(Gamma1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_beta0(pi_beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_beta1(pi_beta1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_w2(sigma_w2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_b20(sigma_b20SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_b21(sigma_b21SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< List >::type prior_params(prior_paramsSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< int >::type lfsr_niter(lfsr_niterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type return_samples(return_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(restart_gibbs_2groups_cpp(Y, G, group, Z, F, W, Gamma0, Gamma1, beta0, beta1, pi_vec, pi_beta0, pi_beta1, psi, sigma_w2, sigma_b20, sigma_b21, c2, prior_params, prior_type, niter, ave_niter, lfsr_niter, verbose, return_samples));
    return rcpp_result_gen;
END_RCPP
}
// Cquantile
double Cquantile(arma::vec x, double q);
RcppExport SEXP _GSFA_Cquantile(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(Cquantile(x, q));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::vec mvrnormArma(arma::vec mu, arma::mat sigma);
RcppExport SEXP _GSFA_mvrnormArma(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// sample_Z_cpp
arma::mat sample_Z_cpp(int N, int K, arma::mat Y, arma::mat F, arma::mat W, arma::mat G, arma::mat beta, arma::vec psi);
RcppExport SEXP _GSFA_sample_Z_cpp(SEXP NSEXP, SEXP KSEXP, SEXP YSEXP, SEXP FSEXP, SEXP WSEXP, SEXP GSEXP, SEXP betaSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_Z_cpp(N, K, Y, F, W, G, beta, psi));
    return rcpp_result_gen;
END_RCPP
}
// sample_gammaBeta_cpp
List sample_gammaBeta_cpp(int N, int M, int K, arma::mat Z, arma::mat G, arma::mat Gamma, arma::mat beta, arma::vec sigma_b2, arma::vec pi_beta);
RcppExport SEXP _GSFA_sample_gammaBeta_cpp(SEXP NSEXP, SEXP MSEXP, SEXP KSEXP, SEXP ZSEXP, SEXP GSEXP, SEXP GammaSEXP, SEXP betaSEXP, SEXP sigma_b2SEXP, SEXP pi_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_b2(sigma_b2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_beta(pi_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_gammaBeta_cpp(N, M, K, Z, G, Gamma, beta, sigma_b2, pi_beta));
    return rcpp_result_gen;
END_RCPP
}
// sample_W_cpp
arma::mat sample_W_cpp(int P, int K, arma::mat Y, arma::mat Z, arma::mat F, arma::mat W, arma::vec psi, arma::vec sigma_w2, arma::vec c2);
RcppExport SEXP _GSFA_sample_W_cpp(SEXP PSEXP, SEXP KSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP FSEXP, SEXP WSEXP, SEXP psiSEXP, SEXP sigma_w2SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_w2(sigma_w2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(sample_W_cpp(P, K, Y, Z, F, W, psi, sigma_w2, c2));
    return rcpp_result_gen;
END_RCPP
}
// sample_F_cpp
arma::mat sample_F_cpp(int P, int K, arma::mat W, arma::vec pi_vec, arma::vec sigma_w2, arma::vec c2);
RcppExport SEXP _GSFA_sample_F_cpp(SEXP PSEXP, SEXP KSEXP, SEXP WSEXP, SEXP pi_vecSEXP, SEXP sigma_w2SEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_w2(sigma_w2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(sample_F_cpp(P, K, W, pi_vec, sigma_w2, c2));
    return rcpp_result_gen;
END_RCPP
}
// sample_FW_spike_slab_cpp
List sample_FW_spike_slab_cpp(int N, int P, int K, arma::mat Y, arma::mat Z, arma::mat F, arma::mat W, arma::vec psi, arma::vec sigma_w2, arma::vec pi_vec);
RcppExport SEXP _GSFA_sample_FW_spike_slab_cpp(SEXP NSEXP, SEXP PSEXP, SEXP KSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP FSEXP, SEXP WSEXP, SEXP psiSEXP, SEXP sigma_w2SEXP, SEXP pi_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_w2(sigma_w2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pi_vec(pi_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_FW_spike_slab_cpp(N, P, K, Y, Z, F, W, psi, sigma_w2, pi_vec));
    return rcpp_result_gen;
END_RCPP
}
// sample_psi_cpp
arma::vec sample_psi_cpp(int N, int P, arma::mat Y, arma::mat Z, arma::mat F, arma::mat W, arma::vec prior_psi);
RcppExport SEXP _GSFA_sample_psi_cpp(SEXP NSEXP, SEXP PSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP FSEXP, SEXP WSEXP, SEXP prior_psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_psi(prior_psiSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_psi_cpp(N, P, Y, Z, F, W, prior_psi));
    return rcpp_result_gen;
END_RCPP
}
// sample_pi_cpp
arma::vec sample_pi_cpp(int P, int K, arma::mat F, arma::vec prior_pi);
RcppExport SEXP _GSFA_sample_pi_cpp(SEXP PSEXP, SEXP KSEXP, SEXP FSEXP, SEXP prior_piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_pi(prior_piSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_pi_cpp(P, K, F, prior_pi));
    return rcpp_result_gen;
END_RCPP
}
// sample_pi_beta_cpp
arma::vec sample_pi_beta_cpp(int M, int K, arma::mat Gamma, arma::vec prior_pibeta);
RcppExport SEXP _GSFA_sample_pi_beta_cpp(SEXP MSEXP, SEXP KSEXP, SEXP GammaSEXP, SEXP prior_pibetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_pibeta(prior_pibetaSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_pi_beta_cpp(M, K, Gamma, prior_pibeta));
    return rcpp_result_gen;
END_RCPP
}
// sample_sigma_w2_cpp
arma::vec sample_sigma_w2_cpp(int K, int P, arma::mat F, arma::mat W, arma::vec prior_sigma2w, arma::vec c2);
RcppExport SEXP _GSFA_sample_sigma_w2_cpp(SEXP KSEXP, SEXP PSEXP, SEXP FSEXP, SEXP WSEXP, SEXP prior_sigma2wSEXP, SEXP c2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_sigma2w(prior_sigma2wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c2(c2SEXP);
    rcpp_result_gen = Rcpp::wrap(sample_sigma_w2_cpp(K, P, F, W, prior_sigma2w, c2));
    return rcpp_result_gen;
END_RCPP
}
// sample_sigma_w2_spike_slab_cpp
arma::vec sample_sigma_w2_spike_slab_cpp(int K, arma::mat F, arma::mat W, arma::vec prior_sigma2w);
RcppExport SEXP _GSFA_sample_sigma_w2_spike_slab_cpp(SEXP KSEXP, SEXP FSEXP, SEXP WSEXP, SEXP prior_sigma2wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_sigma2w(prior_sigma2wSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_sigma_w2_spike_slab_cpp(K, F, W, prior_sigma2w));
    return rcpp_result_gen;
END_RCPP
}
// sample_c2_cpp
arma::vec sample_c2_cpp(int K, int P, arma::mat F, arma::mat W, arma::vec sigma2w, arma::vec prior_c);
RcppExport SEXP _GSFA_sample_c2_cpp(SEXP KSEXP, SEXP PSEXP, SEXP FSEXP, SEXP WSEXP, SEXP sigma2wSEXP, SEXP prior_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma2w(sigma2wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_c(prior_cSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_c2_cpp(K, P, F, W, sigma2w, prior_c));
    return rcpp_result_gen;
END_RCPP
}
// sample_sigma_b2_cpp
arma::vec sample_sigma_b2_cpp(int M, arma::mat Gamma, arma::mat beta, arma::vec prior_sigma2b);
RcppExport SEXP _GSFA_sample_sigma_b2_cpp(SEXP MSEXP, SEXP GammaSEXP, SEXP betaSEXP, SEXP prior_sigma2bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_sigma2b(prior_sigma2bSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_sigma_b2_cpp(M, Gamma, beta, prior_sigma2b));
    return rcpp_result_gen;
END_RCPP
}
// initialize_random
List initialize_random(int K, arma::mat Y);
RcppExport SEXP _GSFA_initialize_random(SEXP KSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_random(K, Y));
    return rcpp_result_gen;
END_RCPP
}
// initialize_svd
List initialize_svd(int K, arma::mat Y);
RcppExport SEXP _GSFA_initialize_svd(SEXP KSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_svd(K, Y));
    return rcpp_result_gen;
END_RCPP
}
// initialize_given_Z
List initialize_given_Z(int K, arma::mat Y, arma::mat init_Z);
RcppExport SEXP _GSFA_initialize_given_Z(SEXP KSEXP, SEXP YSEXP, SEXP init_ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init_Z(init_ZSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_given_Z(K, Y, init_Z));
    return rcpp_result_gen;
END_RCPP
}
// initialize_gammaBeta
List initialize_gammaBeta(int M, int K, arma::mat G, arma::mat Z);
RcppExport SEXP _GSFA_initialize_gammaBeta(SEXP MSEXP, SEXP KSEXP, SEXP GSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(initialize_gammaBeta(M, K, G, Z));
    return rcpp_result_gen;
END_RCPP
}
// compute_posterior_mean_cpp
List compute_posterior_mean_cpp(arma::cube Gamma_mtx, arma::cube beta_mtx, arma::mat pi_beta_mtx, arma::cube Z_mtx, arma::cube F_mtx, arma::cube W_mtx, arma::mat pi_mtx, arma::mat sigma_w2_mtx, arma::mat c2_mtx, int niter, int ave_niter, String prior_type);
RcppExport SEXP _GSFA_compute_posterior_mean_cpp(SEXP Gamma_mtxSEXP, SEXP beta_mtxSEXP, SEXP pi_beta_mtxSEXP, SEXP Z_mtxSEXP, SEXP F_mtxSEXP, SEXP W_mtxSEXP, SEXP pi_mtxSEXP, SEXP sigma_w2_mtxSEXP, SEXP c2_mtxSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP prior_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Gamma_mtx(Gamma_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type beta_mtx(beta_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pi_beta_mtx(pi_beta_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Z_mtx(Z_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type F_mtx(F_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W_mtx(W_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pi_mtx(pi_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_w2_mtx(sigma_w2_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type c2_mtx(c2_mtxSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_posterior_mean_cpp(Gamma_mtx, beta_mtx, pi_beta_mtx, Z_mtx, F_mtx, W_mtx, pi_mtx, sigma_w2_mtx, c2_mtx, niter, ave_niter, prior_type));
    return rcpp_result_gen;
END_RCPP
}
// compute_posterior_mean_2groups_cpp
List compute_posterior_mean_2groups_cpp(arma::cube Gamma0_mtx, arma::cube beta0_mtx, arma::mat pi_beta0_mtx, arma::cube Gamma1_mtx, arma::cube beta1_mtx, arma::mat pi_beta1_mtx, arma::cube Z_mtx, arma::cube F_mtx, arma::cube W_mtx, arma::mat pi_mtx, arma::mat sigma_w2_mtx, arma::mat c2_mtx, int niter, int ave_niter, String prior_type);
RcppExport SEXP _GSFA_compute_posterior_mean_2groups_cpp(SEXP Gamma0_mtxSEXP, SEXP beta0_mtxSEXP, SEXP pi_beta0_mtxSEXP, SEXP Gamma1_mtxSEXP, SEXP beta1_mtxSEXP, SEXP pi_beta1_mtxSEXP, SEXP Z_mtxSEXP, SEXP F_mtxSEXP, SEXP W_mtxSEXP, SEXP pi_mtxSEXP, SEXP sigma_w2_mtxSEXP, SEXP c2_mtxSEXP, SEXP niterSEXP, SEXP ave_niterSEXP, SEXP prior_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Gamma0_mtx(Gamma0_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type beta0_mtx(beta0_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pi_beta0_mtx(pi_beta0_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Gamma1_mtx(Gamma1_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type beta1_mtx(beta1_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pi_beta1_mtx(pi_beta1_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Z_mtx(Z_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type F_mtx(F_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W_mtx(W_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pi_mtx(pi_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_w2_mtx(sigma_w2_mtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type c2_mtx(c2_mtxSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type ave_niter(ave_niterSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_posterior_mean_2groups_cpp(Gamma0_mtx, beta0_mtx, pi_beta0_mtx, Gamma1_mtx, beta1_mtx, pi_beta1_mtx, Z_mtx, F_mtx, W_mtx, pi_mtx, sigma_w2_mtx, c2_mtx, niter, ave_niter, prior_type));
    return rcpp_result_gen;
END_RCPP
}
// compute_lfsr_cpp
arma::mat compute_lfsr_cpp(arma::cube beta_mtx, arma::cube W_mtx, arma::cube F_mtx, int use_niter, String prior_type);
RcppExport SEXP _GSFA_compute_lfsr_cpp(SEXP beta_mtxSEXP, SEXP W_mtxSEXP, SEXP F_mtxSEXP, SEXP use_niterSEXP, SEXP prior_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type beta_mtx(beta_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W_mtx(W_mtxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type F_mtx(F_mtxSEXP);
    Rcpp::traits::input_parameter< int >::type use_niter(use_niterSEXP);
    Rcpp::traits::input_parameter< String >::type prior_type(prior_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_lfsr_cpp(beta_mtx, W_mtx, F_mtx, use_niter, prior_type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GSFA_gsfa_gibbs_cpp", (DL_FUNC) &_GSFA_gsfa_gibbs_cpp, 22},
    {"_GSFA_restart_gsfa_gibbs_cpp", (DL_FUNC) &_GSFA_restart_gsfa_gibbs_cpp, 20},
    {"_GSFA_gsfa_gibbs_2groups_cpp", (DL_FUNC) &_GSFA_gsfa_gibbs_2groups_cpp, 24},
    {"_GSFA_restart_gibbs_2groups_cpp", (DL_FUNC) &_GSFA_restart_gibbs_2groups_cpp, 25},
    {"_GSFA_Cquantile", (DL_FUNC) &_GSFA_Cquantile, 2},
    {"_GSFA_mvrnormArma", (DL_FUNC) &_GSFA_mvrnormArma, 2},
    {"_GSFA_sample_Z_cpp", (DL_FUNC) &_GSFA_sample_Z_cpp, 8},
    {"_GSFA_sample_gammaBeta_cpp", (DL_FUNC) &_GSFA_sample_gammaBeta_cpp, 9},
    {"_GSFA_sample_W_cpp", (DL_FUNC) &_GSFA_sample_W_cpp, 9},
    {"_GSFA_sample_F_cpp", (DL_FUNC) &_GSFA_sample_F_cpp, 6},
    {"_GSFA_sample_FW_spike_slab_cpp", (DL_FUNC) &_GSFA_sample_FW_spike_slab_cpp, 10},
    {"_GSFA_sample_psi_cpp", (DL_FUNC) &_GSFA_sample_psi_cpp, 7},
    {"_GSFA_sample_pi_cpp", (DL_FUNC) &_GSFA_sample_pi_cpp, 4},
    {"_GSFA_sample_pi_beta_cpp", (DL_FUNC) &_GSFA_sample_pi_beta_cpp, 4},
    {"_GSFA_sample_sigma_w2_cpp", (DL_FUNC) &_GSFA_sample_sigma_w2_cpp, 6},
    {"_GSFA_sample_sigma_w2_spike_slab_cpp", (DL_FUNC) &_GSFA_sample_sigma_w2_spike_slab_cpp, 4},
    {"_GSFA_sample_c2_cpp", (DL_FUNC) &_GSFA_sample_c2_cpp, 6},
    {"_GSFA_sample_sigma_b2_cpp", (DL_FUNC) &_GSFA_sample_sigma_b2_cpp, 4},
    {"_GSFA_initialize_random", (DL_FUNC) &_GSFA_initialize_random, 2},
    {"_GSFA_initialize_svd", (DL_FUNC) &_GSFA_initialize_svd, 2},
    {"_GSFA_initialize_given_Z", (DL_FUNC) &_GSFA_initialize_given_Z, 3},
    {"_GSFA_initialize_gammaBeta", (DL_FUNC) &_GSFA_initialize_gammaBeta, 4},
    {"_GSFA_compute_posterior_mean_cpp", (DL_FUNC) &_GSFA_compute_posterior_mean_cpp, 12},
    {"_GSFA_compute_posterior_mean_2groups_cpp", (DL_FUNC) &_GSFA_compute_posterior_mean_2groups_cpp, 15},
    {"_GSFA_compute_lfsr_cpp", (DL_FUNC) &_GSFA_compute_lfsr_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GSFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
