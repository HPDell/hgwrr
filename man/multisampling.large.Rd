\name{multisampling.large}
\alias{multisampling.large}
\docType{data}
\title{Large Scale Simulated Spatial Multisampling Data (DataFrame)}
\description{
A large scale simulation data of spatial hierarchical structure and
samples overlapping on certain locations.
}
\usage{data(multisampling)}
\format{
A list of three items called "data", "coords" and "beta".
Item "data" is a data frame with 13862 observations at 200 locations
and the following 6 variables.
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
\dontrun{
data(multisampling.large)
hgwr(formula = y ~ L(g1 + g2) + x1 + (z1 | group),
     data = multisampling.large$data,
     coords = multisampling.large$coords,
     bw = 32, kernel = "bisquared")
}
}
