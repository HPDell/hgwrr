#' Generic method to test spatial heterogeneity
#'
#' @param x The data to be tested.
#' @export
spatial_hetero_test <- function(x, ...) UseMethod("spatial_hetero_test")

#' @describeIn spatial_hetero_test
#' Default behavior.
#' @method spatial_hetero_test default
spatial_hetero_test.default <- function(x, ...) stop("Method not implemented")

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
#' @return A `spahetbootres` object of permutation-test results with the following items:
#' \describe{
#'  \item{\code{vars}}{The names of variables.}
#'  \item{\code{t0}}{The value of the statistics on original values.}
#'  \item{\code{t}}{The value of the same statistics on permuted values.}
#'  \item{\code{p}}{The p-value for each variable.}
#' }
#' Currently, variance is used as the statistics.
#'
#' @export
spatial_hetero_test_data <- function(
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
  class(res) <- "spahetbootres"
  res
}

#' Print the result of spatial heterogeneity test
#'
#' @param x A `spahetbootres` object.
#' @param \dots Other unused arguments.
#'
#' @method print spahetbootres
#' @export
print.spahetbootres <- function(x, ...) {
  if (!inherits(x, "spahetbootres")) {
    stop("The `x` must be an spahetbootres object.")
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

#' @describeIn spatial_hetero_test
#' For the matrix, `coords` is necessary.
#'
#' @inheritDotParams spatial_hetero_test_data resample:verbose
#' @param coords The coordinates used for testing.
#' Accepts a matrix or vector.
#' For matrix, it needs to have the same number of rows as `x`.
#' For vector, it indicates the columns in `x`
#' and the actual coordinates will be taken from `x`.
#'
#' @method spatial_hetero_test matrix
#' @export
spatial_hetero_test.matrix <- function(x, coords, ...) {
  if (!inherits(x, "matrix")) {
    stop("Argument x is not a matrix")
  }
  if (!(is.numeric(x) || is.integer(x))) {
    stop("Only support numeric or integer matrix")
  }
  if (is.matrix(coords)) {
    if (nrow(coords) != nrow(x))
      stop("`coords` needs to have the same number of rows as `x`")
  }
  call <- match.call(spatial_hetero_test_data, expand.dots = TRUE)
  call[[1]] <- spatial_hetero_test_data
  if (is.vector(coords)) {
    call[["x"]] <- x[, -coords]
    call[["coords"]] <- x[, coords]
  }
  eval.parent(call)
}

#' @describeIn spatial_hetero_test
#'
#' @inheritParams spatial_hetero_test.matrix
#'
#' @method spatial_hetero_test numeric
#' @export
spatial_hetero_test.numeric <- function(x, coords, ...) {
  if (!is.numeric(x) && !is.vector(x)) {
    stop("Argument x needs to be a numeric vector")
  }
  x_mat <- matrix(x, ncol = 1)
  spatial_hetero_test.matrix(x_mat, coords, ...)
}

#' @describeIn spatial_hetero_test
#'
#' @inheritParams spatial_hetero_test.matrix
#'
#' @method spatial_hetero_test vector
#' @export
spatial_hetero_test.vector <- function(x, coords, ...) {
  if (!inherits(x, "vector")) {
    stop("Argument x is not a vector")
  }
  x_mat <- matrix(x, ncol = 1)
  spatial_hetero_test.matrix(x_mat, coords, ...)
}

#' @describeIn spatial_hetero_test
#' For data frame, `coords` is necessary.
#'
#' @inheritParams spatial_hetero_test.matrix
#'
#' @method spatial_hetero_test data.frame
#' @export
spatial_hetero_test.data.frame <- function(x, coords, ...) {
  if (!inherits(x, "data.frame")) {
    stop("Argument x is not a data.frame")
  }
  x_numerical <- as.matrix(x[vapply(x, function(x) is.numeric(x), FALSE)])
  if (length(x_numerical) < 1) {
    stop("There needs at least one numeric column in `x`")
  }
  if (is.vector(coords)) {
    coords <- x[coords]
  }
  spatial_hetero_test.matrix(x_numerical, coords, ...)
}

#' @describeIn spatial_hetero_test
#' For the `sf` object, coordinates of centroids are used.
#' Only the numerical columns are tested.
#'
#' @inheritDotParams spatial_hetero_test_data resample:verbose
#'
#' @importFrom sf st_centroid st_coordinates st_drop_geometry
#'
#' @method spatial_hetero_test sf
#' @export
spatial_hetero_test.sf <- function(x, ...) {
  if (!inherits(x, "sf")) {
    stop("Argument x is not an sf object")
  }
  coords <- sf::st_coordinates(sf::st_centroid(x))
  x_nogeo <- sf::st_drop_geometry(x)
  spatial_hetero_test.data.frame(x_nogeo, coords, ...)
}

#' @importFrom stats var
stat_glsw <- function(g, sd) diag(var(g / sd))

#' @describeIn hgwr
#' Test the spatial heterogeneity with bootstrapping.
#'
#' @param x An `hgwrm` object
#' @param round The number of times to sampling from model.
#' @param statistic A function used to calculate the statistics
#' on the original data and bootstrapped data.
#' Default to the variance of standardlised GLSW estimates.
#' @param parallel If TRUE, use `furrr` package to parallel.
#' @param verbose Override the verbose setting
#' in the `hgwrm` model with this value.
#' @param \dots Arguments to be used in re-fitting the model
#' with function [hgwrr::hgwr].
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats fitted
#' @method spatial_hetero_test hgwrm
#' @export
spatial_hetero_test.hgwrm <- function(
  x,
  round = 99,
  statistic = stat_glsw,
  parallel = FALSE,
  verbose = 0,
  ...
) {
  t0 <- with(x, statistic(gamma, gamma_se))
  args <- as.list(x$call[-1])
  args$coords <- x$coords
  args$f_test <- FALSE
  args$max_iters <- 100
  args$bw <- x$bw
  args$verbose <- verbose
  if (...length() > 0) {
    settings <- list(...)
    for (key in names(settings)) {
      args[[key]] <- settings[[key]]
    }
  }
  yhat <- fitted(x)
  covar0 <- x$D
  sigma0 <- x$sigma
  data0 <- x$frame
  resp <- x$effects$response
  ge <- x$effects$group
  groups <- unique(data0[[ge]])
  z0 <- x$frame.parsed$z
  z_group <- lapply(groups, function(gr) z0[data0[[ge]] == gr, ])
  pgb <- NULL
  if (requireNamespace("progressr", quietly = TRUE)) {
    pgb <- progressr::progressor(round)
  }
  worker <- function(i, pgb = NULL) {
    ei <- do.call(c, lapply(z_group, function(zi) {
      MASS::mvrnorm(mu = rep(0, nrow(zi)),
                    Sigma = zi %*% covar0 %*% t(zi) + diag(sigma0, nrow(zi)))
    }))
    ysi <- ei + yhat
    datai <- data0
    datai[[resp]] <- ysi
    argsi <- args
    argsi$data <- datai
    modeli <- do.call(hgwrr::hgwr_fit, argsi)
    stati <- with(modeli, statistic(gamma, gamma_se))
    if (!is.null(pgb)) pgb()
    stati
  }
  if (is.logical(parallel) && parallel) {
    if (requireNamespace("furrr", quietly = TRUE)) {
      t <- do.call(rbind, furrr::future_map(
        seq_len(round), worker,
        pgb = pgb,
        .options = furrr::furrr_options(seed = TRUE)
      ))
    } else {
      stop("Package \"furrr\" must be installed to parallel",
           call. = FALSE)
    }
  } else {
    t <- do.call(rbind, lapply(seq_len(round), worker, pgb = pgb))
  }
  pv <- sapply(seq_along(t0), function(i) {
    mean(abs(t[, i]) > abs(t0[i]))
  })
  res <- list(
    t0 = t0,
    t = t,
    vars = x$effects$glsw,
    pv = pv
  )
  class(res) <- "spahetbootres"
  res
}
