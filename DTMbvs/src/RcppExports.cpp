// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// DTMbvs
List DTMbvs(int iterations, int thin, String prior, arma::mat x, arma::vec branch_location, List node_children_pointer, arma::mat branch_counts, arma::mat subtree_counts, arma::mat alpha, arma::cube phi, arma::cube zeta, double sigma2_alpha, double sigma2_phi, double aa, double bb, arma::cube Omega, arma::cube G, arma::cube Var, arma::mat S, double v0, double v1, double a_G, double b_G, double pie, double lambda);
RcppExport SEXP _DTMbvs_DTMbvs(SEXP iterationsSEXP, SEXP thinSEXP, SEXP priorSEXP, SEXP xSEXP, SEXP branch_locationSEXP, SEXP node_children_pointerSEXP, SEXP branch_countsSEXP, SEXP subtree_countsSEXP, SEXP alphaSEXP, SEXP phiSEXP, SEXP zetaSEXP, SEXP sigma2_alphaSEXP, SEXP sigma2_phiSEXP, SEXP aaSEXP, SEXP bbSEXP, SEXP OmegaSEXP, SEXP GSEXP, SEXP VarSEXP, SEXP SSEXP, SEXP v0SEXP, SEXP v1SEXP, SEXP a_GSEXP, SEXP b_GSEXP, SEXP pieSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< String >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type branch_location(branch_locationSEXP);
    Rcpp::traits::input_parameter< List >::type node_children_pointer(node_children_pointerSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type branch_counts(branch_countsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type subtree_counts(subtree_countsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_alpha(sigma2_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_phi(sigma2_phiSEXP);
    Rcpp::traits::input_parameter< double >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< double >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Var(VarSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< double >::type a_G(a_GSEXP);
    Rcpp::traits::input_parameter< double >::type b_G(b_GSEXP);
    Rcpp::traits::input_parameter< double >::type pie(pieSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(DTMbvs(iterations, thin, prior, x, branch_location, node_children_pointer, branch_counts, subtree_counts, alpha, phi, zeta, sigma2_alpha, sigma2_phi, aa, bb, Omega, G, Var, S, v0, v1, a_G, b_G, pie, lambda));
    return rcpp_result_gen;
END_RCPP
}
// dmbvs_ss
List dmbvs_ss(arma::mat XX, arma::mat YY, arma::vec alpha, arma::vec beta, arma::vec mu_al, arma::vec sig_al, arma::mat mu_be, arma::mat sig_be, double aa_hp, double bb_hp, arma::vec prop_per_alpha, arma::vec prop_per_beta, int GG, int thin, int Log);
RcppExport SEXP _DTMbvs_dmbvs_ss(SEXP XXSEXP, SEXP YYSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP mu_alSEXP, SEXP sig_alSEXP, SEXP mu_beSEXP, SEXP sig_beSEXP, SEXP aa_hpSEXP, SEXP bb_hpSEXP, SEXP prop_per_alphaSEXP, SEXP prop_per_betaSEXP, SEXP GGSEXP, SEXP thinSEXP, SEXP LogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_al(mu_alSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig_al(sig_alSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_be(mu_beSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_be(sig_beSEXP);
    Rcpp::traits::input_parameter< double >::type aa_hp(aa_hpSEXP);
    Rcpp::traits::input_parameter< double >::type bb_hp(bb_hpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prop_per_alpha(prop_per_alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prop_per_beta(prop_per_betaSEXP);
    Rcpp::traits::input_parameter< int >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type Log(LogSEXP);
    rcpp_result_gen = Rcpp::wrap(dmbvs_ss(XX, YY, alpha, beta, mu_al, sig_al, mu_be, sig_be, aa_hp, bb_hp, prop_per_alpha, prop_per_beta, GG, thin, Log));
    return rcpp_result_gen;
END_RCPP
}
// dmbvs_gibbs
List dmbvs_gibbs(arma::mat XX, arma::mat YY, arma::vec alpha, arma::vec beta, arma::vec mu_al, arma::vec sig_al, arma::mat mu_be, arma::mat sig_be, double aa_hp, double bb_hp, arma::vec prop_per_alpha, arma::vec prop_per_beta, int GG, int thin, int Log);
RcppExport SEXP _DTMbvs_dmbvs_gibbs(SEXP XXSEXP, SEXP YYSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP mu_alSEXP, SEXP sig_alSEXP, SEXP mu_beSEXP, SEXP sig_beSEXP, SEXP aa_hpSEXP, SEXP bb_hpSEXP, SEXP prop_per_alphaSEXP, SEXP prop_per_betaSEXP, SEXP GGSEXP, SEXP thinSEXP, SEXP LogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_al(mu_alSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig_al(sig_alSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_be(mu_beSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_be(sig_beSEXP);
    Rcpp::traits::input_parameter< double >::type aa_hp(aa_hpSEXP);
    Rcpp::traits::input_parameter< double >::type bb_hp(bb_hpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prop_per_alpha(prop_per_alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prop_per_beta(prop_per_betaSEXP);
    Rcpp::traits::input_parameter< int >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type Log(LogSEXP);
    rcpp_result_gen = Rcpp::wrap(dmbvs_gibbs(XX, YY, alpha, beta, mu_al, sig_al, mu_be, sig_be, aa_hp, bb_hp, prop_per_alpha, prop_per_beta, GG, thin, Log));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DTMbvs_DTMbvs", (DL_FUNC) &_DTMbvs_DTMbvs, 25},
    {"_DTMbvs_dmbvs_ss", (DL_FUNC) &_DTMbvs_dmbvs_ss, 15},
    {"_DTMbvs_dmbvs_gibbs", (DL_FUNC) &_DTMbvs_dmbvs_gibbs, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_DTMbvs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
