#' Test the spatial heterogeneity in data based on permutation.
#'
#' @param x A matrix of data to be tested. Each column is a variable.
#' @param coords A matrix of coordinates.
#' @param \dots Additional arguments.
#' @param resample The total times of resampling with replacement.
#'    Default to 5000.
#' @param poly The number of polynomial terms used by the polynomial estimator.
#'    Default to 2.
#' @param bw The adaptive bandwidth used by the polynomial estimator.
#'    Default to 10.
#' @param kernel The kernel function used by the polynomial estimator.
#' @param verbose The verbosity level. Default to 0.
#'
#' @return A `shgt` object of permutation-test results with the following items:
#' \describe{
#'  \item{\code{vars}}{The names of variables.}
#'  \item{\code{t0}}{The value of the statistics on original values.}
#'  \item{\code{t}}{The value of the same statistics on permuted values.}
#'  \item{\code{p}}{The p-value for each variable.}
#' }
#' Currently, variance is used as the statistics.
#'
#' @examples
#' data(multisampling.large)
#' spatial_hetero_test(multisampling.large$beta, multisampling.large$coords)
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
    stop("The rows of x and coords must match.")
  }
  var_names <- colnames(x)
  kernel_id <- switch(kernel,
    "gaussian" = 0,
    "bisquared" = 1
  )
  tv <- spatial_hetero_bootstrap(x, coords, poly, resample,
                                 bw, kernel_id, verbose)
  pv <- sapply(seq_along(tv$t0), function(i) {
    with(tv, mean(abs(t[, i]) > abs(t0[i])))
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
#' @param \dots Other unused arguments.
#'
#' @method print shgt
#' @export
print.shgt <- function(x, ...) {
  if (!inherits(x, "shgt")) {
    stop("The `x` must be an shgt object.")
  }
  cat("Spatial Heterogeneity Test\n", fill = TRUE)
  show_tbl <- data.frame(
    t0 = matrix2char(as.vector(x$t0)),
    t = apply(
      x$t, 2,
      function(r) sprintf("[%s,%s]", matrix2char(min(r)), matrix2char(max(r)))
    ),
    P = matrix2char(x$p),
    stars = sapply(x$p, pv2stars)
  )
  colnames(show_tbl) <- c("t0", "t", "Pr(t>t0)", "")
  rownames(show_tbl) <- x$vars
  print(show_tbl)
  cat("\n", fill = TRUE)
  cat("Significance levels: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' '\n",
      fill = TRUE)
}

#' Generic method to test spatial heterogeneity
#'
#' @param x The data to be tested.
#' @inheritParams spatial_hetero_test
#' @export
spatial_hetero <- function(x, ...) UseMethod("spatial_hetero")

#' @describeIn spatial_hetero
#' For the matrix, `coords` is necessary.
#'
#' @param coords The coordinates used for testing.
#' @method spatial_hetero matrix
spatial_hetero.matrix <- function(x, coords, ...) {
  if (!inherits(x, "matrix")) {
    stop("Argument x is not a matrix")
  }
  if (!(is.numeric(x) || is.integer(x))) {
    stop("Only support numeric or integer matrix")
  }
  call <- match.call(spatial_hetero_test, expand.dots = TRUE)
  eval.parent(call)
}

#' @describeIn spatial_hetero
#' For the matrix, `coords` is necessary.
#'
#' @param coords The coordinates used for testing.
#' @importFrom sf st_centroid st_coordinates st_drop_geometry
#' @method spatial_hetero sf
spatial_hetero.sf <- function(x, ...) {
  if (!inherits(x, "matrix")) {
    stop("Argument x is not a matrix")
  }
  if (!(is.numeric(x) || is.integer(x))) {
    stop("Only support numeric or integer matrix")
  }
  coords <- sf::st_coordinates(sf::st_centroid(data))
  x_nogeo <- sf::st_drop_geometry(x)
  x_numerical <- x_nogeo[vapply(x_nogeo, function(x) is.numeric(x), FALSE)]
  call <- match.call(spatial_hetero_test, expand.dots = TRUE)
  call[["x"]] <- x_numerical
  call[["coords"]] <- coords
  eval.parent(call)
}

#' @describeIn spatial_hetero
#' Default behavior.
#' @method spatial_hetero default
spatial_hetero.default <- function(x, ...) stop("Method not implemented")
