#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <armadillo>

inline void prcout(const std::string& message)
{
    Rcpp::Rcout << message;
}

#endif // UTILS_H