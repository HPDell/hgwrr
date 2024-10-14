#' HGWR: Hierarchical and Geographically Weighted Regression
#' 
#' An R and C++ implementation of Hierarchical and Geographically Weighted
#' Regression (HGWR) model is provided in this package. This model divides
#' coefficients into three types: local fixed effects, global fixed effects,
#' and random effects. If data have spatial hierarchical structures
#' (especially are overlapping on some locations),
#' it is worth trying this model to reach better fitness.
#' 
#' @docType package
#' @name hgwrr-package
#' @details \packageDESCRIPTION{hgwrr}
#' @author Yigong Hu, Richard Harris, Richard Timmerman
#' @note Acknowledgement:
#' We gratefully acknowledge support from China Scholarship Council.
#' @references Hu, Y., Lu, B., Ge, Y., Dong, G., 2022.
#' Uncovering spatial heterogeneity in real estate prices via
#' combined hierarchical linear model and geographically weighted regression.
#' Environment and Planning B: Urban Analytics and City Science.
#' \doi{10.1177/23998083211063885}
#'
#' @useDynLib hgwrr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
NULL