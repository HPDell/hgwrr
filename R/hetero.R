#' Test the spatial heterogeneity in data based on permutation.
#' 
#' @param x A matrix of data to be tested. Each column is a variable.
#' @param coords A matrix of coordinates.
#' @param \dots Additional arguments.
#' @param resample The total times of resampling with replacement. Default to 5000.
#' @param poly The number of polynomial terms used by the polynomial estimator. Default to 2.
#' @param bw The adaptive bandwidth used by the polynomial estimator. Default to 10.
#' @param kernel The kernel function used by the polynomial estimator.
#' @param verbose The verbosity level. Default to 0.
#' 
#' @return A `shgt` object of permutation-test results with the following items:
#' \describe{
#'  \item{\code{vars}}{The names of variables.}
#'  \item{\code{t0}}{The value of the statistics (variance of density estimation) on original values.}
#'  \item{\code{t}}{The value of the same statistics on permuted values.}
#'  \item{\code{p}}{The p-value for each variable.}
#' }
#' 
#' @examples
#' data(multisampling.large)
#' spatial_hetero_test(multisampling.large$beta, multisampling.large$coords)
#' 
#' @export
spatial_hetero_test <- function(
    x,
    coords,
    ...,
    resample = 5000,
    poly = 2,
    bw = 10,
    kernel = c("bisquared", "gaussian"),
    verbose = 0
) {
    kernel <- match.arg(kernel)
    x <- as.matrix(x)
    coords <- as.matrix(coords)
    if (nrow(x) != nrow(coords)) {
        stop('The rows of x and coords must match.')
    }
    var_names <- colnames(x)
    kernel_id <- switch(kernel, "gaussian" = 0, "bisquared" = 1)
    tv <- spatial_hetero_perm(x, coords, poly, resample, bw, kernel_id, verbose)
    pv <- sapply(seq_along(tv$t0), function(i) {
        with(tv, mean(t[,i] > t0[i]))
    })
    res <- tv
    res$vars <- var_names
    res$p <- pv
    class(res) <- "shgt"
    res
}

#' Print the result of spatial heterogeneity test
#' 
#' @param x A `shgt` object.
#' 
#' @method print shgt
#' @export 
print.shgt <- function(x, ...) {
    if (!inherits(x, 'shgt')) {
        stop("The `x` must be an shgt object.")
    }
    cat("Spatial Heterogeneity Test\n", fill = T)
    show_tbl <- data.frame(
        t0 = as.vector(x$t0),
        t = colMeans(x$t),
        P = x$p,
        stars = sapply(x$p, pv2stars)
    )
    colnames(show_tbl) <- c("t0", "Avg. t", "Pr(t>t0)", "")
    rownames(show_tbl) <- x$vars
    print(show_tbl)
    cat("\n", fill = T)
    cat("Significance levels: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' '\n", fill = T)
}
