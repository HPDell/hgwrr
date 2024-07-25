// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <optional>
#include <algorithm>
#include "hlmgwr.h"
#include "utils.h"
#include "progress.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace hgwr;

HGWR::GWRKernelFunctionSquared KERNEL[2] = {
    HGWR::gwr_kernel_bisquare2,
    HGWR::gwr_kernel_gaussian2
};

uword factor(uword x) {
    uword f = 1;
    for (uword i = 1; i <= x; i++)
    {
        f *= i;
    }
    return f;
}

umat cart_prod(const umat& a, const uvec& b) {
    return join_rows(repelem(a, b.n_rows, 1), repmat(b, a.n_rows, 1));
}

umat poly_alpha(uword k, uword p) {
    umat ialpha = linspace<uvec>(0, p, p + 1);
    umat all_alpha = ialpha;
    for (size_t i = 1; i < k; i++)
    {
        all_alpha = cart_prod(all_alpha, ialpha);
    }
    uvec sum_alpha = sum(all_alpha, 1);
    return all_alpha.rows(find((sum_alpha > 0) && (sum_alpha <= k)));
}

mat poly_items(const mat& x, const mat& alpha, uword p) {
    umat ialpha = linspace<uvec>(0, p, p + 1);
    umat all_alpha = ialpha;
    uword k = x.n_cols;
    for (size_t i = 1; i < k; i++)
    {
        all_alpha = cart_prod(all_alpha, ialpha);
    }
    uvec sum_alpha = sum(all_alpha, 1);
    mat fact_alpha = alpha;
    fact_alpha.transform([](uword i){ return factor(i); });
    vec div = prod(fact_alpha, 1);
    mat px(x.n_rows, alpha.n_rows + 1, arma::fill::ones);
    for (size_t i = 0; i < alpha.n_rows; i++)
    {
        px.col(i + 1) = prod(arma::pow(x.each_row(), alpha.row(i)), 1) / div(i);
    }
    return px;
}

mat denreg_poly(
    const mat& x,
    const mat& uv,
    const mat& alpha,
    size_t poly = 2,
    double bw = 10,
    int kernel = 0
) {
    mat g(arma::size(x));
    for (uword i = 0; i < x.n_rows; i++) {
        mat duv = (uv.each_row() - uv.row(i));
        mat U = poly_items(duv, alpha, poly);
        vec d = sqrt(sum(duv % duv, 1));
        double b = HGWR::actual_bw(d, bw);
        mat wi = (*(KERNEL + kernel))(d % d, b * b);
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
    const cube& L
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
    uword ndp = uv.n_rows;
    mat alpha = arma::conv_to<mat>::from(poly_alpha(uv.n_cols, poly));
    cube L;
    if (ndp < 4096) {
        precalc_dw = true;
        L = cube(alpha.n_rows + 1, ndp, ndp);
        if (verbose > 0) Rcout << "* Calculating spatial weights in advance" << "\n";
        for (uword i = 0; i < ndp; i++)
        {
            mat duv = uv.each_row() - uv.row(i);
            mat U = poly_items(duv, alpha, poly);
            vec d = sqrt(sum(duv % duv, 1));
            double b = HGWR::actual_bw(d, bw);
            mat wi = (*(KERNEL + kernel))(d % d, b * b);
            mat Utw = trans(U.each_col() % wi);
            L.slice(i) = inv(Utw * U) * Utw;
        }
        if (verbose > 0) Rcout << "* Testing with pre-calculated spatial weights" << "\n";
    } else {
        if (verbose > 0) Rcout << "* Testing without pre-calculated spatial weights" << "\n";
    }
    mat r0 = precalc_dw ? denreg_poly(x, uv, L) : denreg_poly(x, uv, alpha, poly, bw, kernel);
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
        mat ri = precalc_dw ? denreg_poly(xi, uv, L) : denreg_poly(xi, uv, alpha, poly, bw, kernel);
        stats.row(i) = var(ri, 0, 0);
        p.tic();
    }
    return Rcpp::List::create(
        Rcpp::Named("t") = stats,
        Rcpp::Named("t0") = stat0
    );
}
