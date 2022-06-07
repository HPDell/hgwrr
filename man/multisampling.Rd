\name{multisampling}
\alias{multisampling}
\alias{Multisampling}
\docType{data}
\title{Simulated Spatial Multisampling Data (DataFrame)}
\description{
A simulation data of spatial hierarchical structure and
samples overlapping on certain locations.
}
\usage{data(multisampling)}
\format{
A list of two items called "data" and "coord".
Item "data" is a data frame with 484 observations at 16 locations
on the following 6 variables.
\describe{
    \item{\code{y}}{a numeric vector, dependent variable \eqn{y}}
    \item{\code{g1}}{a numeric vector, group level independent variable \eqn{g_1}}
    \item{\code{g2}}{a numeric vector, group level independent variable \eqn{g_2}}
    \item{\code{z1}}{a numeric vector, sample level independent variable \eqn{z_1}}
    \item{\code{x1}}{a numeric vector, sample level independent variable \eqn{x_1}}
    \item{\code{group}}{a numeric vector, group id of each sample}
}
where \code{g1} and \code{g2} are used to estimate local fixed effects;
\code{x1} is used to estimate global fixed effects
and \code{z1} is used to estimate random effects.
}
\author{Yigong Hu \email{yigong.hu@bristol.ac.uk}}
\examples{
data(multisampling)
hgwr(formula = y ~ g1 + g2 + x1 + (z1 | group),
     data = multisampling$data,
     local.fixed = c("g1", "g2"),
     coords = multisampling$coords,
     bw = 10)
}
