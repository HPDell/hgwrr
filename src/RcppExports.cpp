// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "hlmgwr.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


// hlmgwr_backfitting_maximum_likelihood
RcppExport SEXP _hgwrr_backfitting_maximum_likelihood(
    SEXP gSEXP, 
    SEXP xSEXP, 
    SEXP zSEXP, 
    SEXP ySEXP, 
    SEXP uSEXP, 
    SEXP groupSEXP, 
    SEXP bwSEXP,
    SEXP alphaSEXP, 
    SEXP eps_iterSEXP, 
    SEXP eps_gradientSEXP, 
    SEXP max_itersSEXP, 
    SEXP max_retriesSEXP,
    SEXP ml_typeSEXP, 
    SEXP verboseSEXP
) {
BEGIN_RCPP
    Rcpp::List rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type gParam(gSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xParam(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type zParam(zSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yParam(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type uParam(uSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type groupParam(groupSEXP);
    Rcpp::traits::input_parameter< const double& >::type bwParams(bwSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps_iter(eps_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps_gradient(eps_gradientSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type max_iters(max_itersSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type max_retries(max_retriesSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type ml_type(ml_typeSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type verbose(verboseSEXP);
    arma::mat g = arma::mat(gParam);
    arma::mat x = arma::mat(xParam);
    arma::mat z = arma::mat(zParam);
    arma::vec y = arma::vec(yParam);
    arma::mat u = arma::mat(uParam);
    arma::uvec group = arma::conv_to<arma::uvec>::from(arma::vec(groupParam)) - 1;
    double bw = double(bwParams);
    HLMGWROptions options(alpha, eps_iter, eps_gradient, max_iters, max_retries, verbose, ml_type);
    HLMGWRParams hgwr_result = backfitting_maximum_likelihood(g, x, z, y, u, group, bw, options);
    rcpp_result_gen = List::create(
        Named("gamma") = hgwr_result.gamma,
        Named("beta") = hgwr_result.beta,
        Named("mu") = hgwr_result.mu,
        Named("D") = hgwr_result.D
    );
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hgwrr_backfitting_maximum_likelihood", (DL_FUNC) &_hgwrr_backfitting_maximum_likelihood, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_hgwrr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
