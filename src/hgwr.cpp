// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "hlmgwr.h"
#include "utils.h"

using namespace std;
using namespace Rcpp;
using namespace hgwr;

// [[Rcpp::export]]
List hgwr_bfml(
    const arma::mat& g,
    const arma::mat& x,
    const arma::mat& z,
    const arma::vec& y,
    const arma::mat& u,
    const arma::vec& group,
    double bw,
    int bw_optim,
    size_t kernel,
    double alpha,
    double eps_iter,
    double eps_gradient,
    size_t max_iters,
    size_t max_retries,
    size_t ml_type,
    size_t verbose
) {
    arma::uvec mgroup = arma::conv_to<arma::uvec>::from(group) - 1;
    auto mkernel = HGWR::KernelType(size_t(kernel));
    HGWR::Options options { alpha, eps_iter, eps_gradient, max_iters, max_retries, verbose, ml_type };
    HGWR algorithm(g, x, z, y, u, mgroup, mkernel, options);
    if (bw_optim < 0) {
        algorithm.set_bw(bw);
    } else {
        algorithm.set_bw_optim(true);
        algorithm.set_bw_criterion_type(HGWR::BwOptimCriterionType(bw_optim));
    }
    algorithm.set_printer(&prcout);
    auto hgwr_result = algorithm.fit();

    return List::create(
        Named("gamma") = hgwr_result.gamma,
        Named("beta") = hgwr_result.beta,
        Named("mu") = hgwr_result.mu,
        Named("D") = hgwr_result.D,
        Named("sigma") = hgwr_result.sigma,
        Named("bw") = hgwr_result.bw,
        Named("logLik") = algorithm.get_loglik(),
        Named("trS") = algorithm.get_trS(),
        Named("var_beta") = algorithm.get_var_beta()
    );
}
