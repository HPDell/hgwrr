#include <Rcpp.h>
#include <armadillo>
#include "utils.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

vec kernel_bisquare_ada(const mat& uv, double bw, uword focus) {
    mat duv = uv.each_row() - uv.row(focus);
    vec d2 = sum(duv % duv, 1);
    vec ds = sort(d2);
    // vec dsd = join_cols(vec({ 0 }), diff(ds));
    // vec dsc = ds(find(dsd > 0));
    // double b = uword(bw) <= dsc.n_elem ? dsc(uword(bw) - 1) : dsc(dsc.n_elem - 1);
    double b = ds(uword(bw) - 1);
    double b2 = b * b;
    vec wi = (1 - d2 / b2) % (1 - d2 / b2) % (d2 <= b2);
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
        mat wi = kernel_bisquare_ada(uv, bw, i);
        mat duv = (uv.each_row() - uv.row(i));
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

// [[Rcpp::export]]
List spatial_hetero_perm(
    const NumericMatrix& rx,
    const NumericMatrix& ruv,
    int poly = 2,
    int resample = 5000,
    double bw = 10,
    int kernel = 0
) {
    mat x = myas(rx);
    mat uv = myas(ruv);
    mat r0 = denreg_poly(x, uv, poly, bw, kernel);
    rowvec stat0 = var(r0, 0, 0);
    mat stats(resample, x.n_cols);
    for (size_t i = 0; i < resample; i++) {
        mat xi = shuffle(x, 0);
        mat ri = denreg_poly(xi, uv, poly, bw, kernel);
        stats.row(i) = var(ri, 0, 0);
    }
    return Rcpp::List::create(
        Rcpp::Named("t") = mywrap(stats),
        Rcpp::Named("t0") = mywrap(stat0)
    );
}
