#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <armadillo>

inline void prcout(const std::string& message)
{
    Rcpp::Rcout << message;
}

inline void prcancel()
{
    Rcpp::checkUserInterrupt();
}

#endif // UTILS_H