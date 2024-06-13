#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <armadillo>

void prcout(std::string message)
{
    Rcpp::Rcout << message;
}

arma::mat myas(const Rcpp::NumericMatrix& rmat)
{
    return arma::mat(rmat.begin(), rmat.nrow(), rmat.ncol());
}

arma::vec myas(const Rcpp::NumericVector& rvec)
{
    return arma::vec(rvec.begin(), rvec.size());
}

SEXP mywrap(const arma::mat& amat)
{
    Rcpp::RObject x = Rcpp::wrap(amat.begin(), amat.end());
    x.attr("dim") = Rcpp::Dimension(amat.n_rows, amat.n_cols);
    return x;
}

SEXP mywrap(const arma::vec& avec)
{
    Rcpp::RObject x = Rcpp::wrap(avec.begin(), avec.end());
    return x;
}

#endif // UTILS_H