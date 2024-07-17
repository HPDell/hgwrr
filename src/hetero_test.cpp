// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <optional>
#include "hlmgwr.h"
#include "utils.h"
#include "progress.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

vec kernel_bisquare_ada(const vec& d, double bw, uword focus) {
    double b = hgwr::HGWR::actual_bw(d, bw);
    vec wi = (1 - (d % d) / (b * b)) % (1 - (d % d) / (b * b)) % (d <= b);
    vec w = wi / sum(wi);
    return w;
}

mat denreg_poly(
    const mat& x,
    const mat& uv,
    size_t poly = 2,
    double bw = 10,
    int kernel = 0
) {
    uword ndp = uv.n_rows, dim = uv.n_cols;
    mat g(arma::size(x));
    for (uword i = 0; i < x.n_rows; i++) {
        mat duv = (uv.each_row() - uv.row(i));
        vec d = sqrt(sum(duv % duv, 1));
        mat wi = kernel_bisquare_ada(d, bw, i);
        mat U(ndp, dim * poly + 1, fill::ones);
        for (size_t p = 0; p < poly; p++) {
            U.cols(p * dim + 1, p * dim + dim) = pow(duv, p + 1);
        }
        mat Utw = trans(U.each_col() % wi);
        mat UtwU = Utw * U;
        for (uword k = 0; k < x.n_cols; k++) {
            vec r = solve(UtwU, Utw * x.col(k));
            g(i, k) = r(0);
        }
    }
    return g;
}

mat denreg_poly(
    const mat& x,
    const mat& uv,
    const cube& L,
    size_t poly = 2,
    double bw = 10,
    int kernel = 0
) {
    mat g(arma::size(x));
    for (uword i = 0; i < x.n_rows; i++) {
        for (uword k = 0; k < x.n_cols; k++) {
            vec r = L.slice(i) * x.col(k);
            g(i, k) = r(0);
        }
    }
    return g;
}

// [[Rcpp::export]]
List spatial_hetero_perm(
    const arma::mat& x,
    const arma::mat& uv,
    int poly = 2,
    int resample = 5000,
    double bw = 10,
    int kernel = 0,
    int verbose = 0
) {
    bool precalc_dw = false;
    uword ndp = uv.n_rows, dim = uv.n_cols, ncL = dim * poly + 1;
    cube L(ncL, ndp, ndp);
    if (ndp < 4096) {
        precalc_dw = true;
        if (verbose > 0) Rcout << "* Calculating spatial weights in advance" << "\n";
        for (uword i = 0; i < ndp; i++)
        {
            mat duv = uv.each_row() - uv.row(i);
            vec d = sqrt(sum(duv % duv, 1));
            vec wi = kernel_bisquare_ada(d, bw, i);
            mat U(ndp, ncL, arma::fill::ones);
            for (size_t p = 0; p < poly; p++) {
                U.cols(p * dim + 1, p * dim + dim) = pow(duv, p + 1);
            }
            mat Utw = trans(U.each_col() % wi);
            L.slice(i) = inv(Utw * U) * Utw;
        }
        if (verbose > 0) Rcout << "* Testing with pre-calculated spatial weights" << "\n";
    } else {
        if (verbose > 0) Rcout << "* Testing without pre-calculated spatial weights" << "\n";
    }
    mat r0 = precalc_dw ? denreg_poly(x, uv, L, poly, bw, kernel) : denreg_poly(x, uv, poly, bw, kernel);
    rowvec stat0 = var(r0, 0, 0);
    mat stats(resample, x.n_cols);
    ProgressBar p(resample, verbose > 0);
    p.display();
    for (size_t i = 0; i < resample; i++) {
        mat xi(size(x));
        for (size_t c = 0; c < x.n_cols; c++)
        {
            xi.col(c) = shuffle(x.col(c));
        }
        mat ri = precalc_dw ? denreg_poly(xi, uv, L, poly, bw, kernel) : denreg_poly(xi, uv, poly, bw, kernel);
        stats.row(i) = var(ri, 0, 0);
        p.tic();
    }
    return Rcpp::List::create(
        Rcpp::Named("t") = stats,
        Rcpp::Named("t0") = stat0
    );
}
