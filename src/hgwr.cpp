#include <Rcpp.h>
#include <armadillo>
#include "hlmgwr.h"
#include "utils.h"

using namespace std;
using namespace Rcpp;
using namespace hgwr;

// [[Rcpp::export]]
List hgwr_bfml(
    const NumericMatrix& g,
    const NumericMatrix& x,
    const NumericMatrix& z,
    const NumericVector& y,
    const NumericMatrix& u,
    const NumericVector& group,
    double bw,
    size_t kernel,
    double alpha,
    double eps_iter,
    double eps_gradient,
    size_t max_iters,
    size_t max_retries,
    size_t ml_type,
    size_t verbose
) {
    arma::mat mg = myas(g);
    arma::mat mx = myas(x);
    arma::mat mz = myas(z);
    arma::vec my = myas(y);
    arma::mat mu = myas(u);
    arma::uvec mgroup = arma::conv_to<arma::uvec>::from(myas(group)) - 1;
    auto mkernel = HGWR::KernelType(size_t(kernel));
    HGWR::Options options { alpha, eps_iter, eps_gradient, max_iters, max_retries, verbose, ml_type };
    HGWR algorithm(mg, mx, mz, my, mu, mgroup, mkernel, bw, options);
    algorithm.set_printer(&prcout);
    auto hgwr_result = algorithm.fit();

    return List::create(
        Named("gamma") = mywrap(hgwr_result.gamma),
        Named("beta") = mywrap(hgwr_result.beta),
        Named("mu") = mywrap(hgwr_result.mu),
        Named("D") = mywrap(hgwr_result.D),
        Named("sigma") = hgwr_result.sigma,
        Named("bw") = hgwr_result.bw
    );
}
