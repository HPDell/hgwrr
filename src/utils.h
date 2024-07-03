#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <armadillo>

inline void prcout(const std::string& message)
{
    Rcpp::Rcout << message;
}

inline arma::mat myas(const Rcpp::NumericMatrix& rmat)
{
    return arma::mat(rmat.begin(), rmat.nrow(), rmat.ncol());
}

inline arma::vec myas(const Rcpp::NumericVector& rvec)
{
    return arma::vec(rvec.begin(), rvec.size());
}

inline SEXP mywrap(const arma::mat& amat)
{
    Rcpp::RObject x = Rcpp::wrap(amat.begin(), amat.end());
    x.attr("dim") = Rcpp::Dimension(amat.n_rows, amat.n_cols);
    return x;
}

inline SEXP mywrap(const arma::vec& avec)
{
    Rcpp::RObject x = Rcpp::wrap(avec.begin(), avec.end());
    return x;
}

#endif // UTILS_H