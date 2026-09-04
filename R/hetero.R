#' Generic method to test spatial heterogeneity
#'
#' @param x The data to be tested.
#' @export
spatial_hetero_test <- function(x, ...) UseMethod("spatial_hetero_test")

#' @describeIn spatial_hetero_test
#' Default behavior.
#' @method spatial_hetero_test default
#' @export
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
#' Takes `x` as values of a series variables stored by column,
#' and `coords` as coordinates for each row in `x`.
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
#' Takes `x` as values of the variable,
#' and `coords` as coordinates for each element in `x`.
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
#' Takes `x` as variable values (each column is a variable),
#' and `coords` as coordinates for each row in `x`.
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

.effect_terms <- function(intercept, effects) {
  c(if (isTRUE(intercept != 0)) "1" else "0", effects)
}

.effectwise_null_formula <- function(x, effect) {
  glsw <- setdiff(x$effects$glsw, c("Intercept", effect))
  fixed <- setdiff(x$effects$fixed, "Intercept")
  slr <- setdiff(x$effects$slr, "Intercept")
  if (!identical(effect, "Intercept")) {
    fixed <- unique(c(fixed, effect))
  }

  glsw_intercept <- x$intercept$glsw
  if (identical(effect, "Intercept")) {
    glsw_intercept <- 0
  }
  if (length(glsw) == 0 && glsw_intercept == 0) {
    stop(
      "The effect-wise null model must retain at least one GLSW effect.",
      call. = FALSE
    )
  }

  stats::as.formula(sprintf(
    "%s ~ %s + L(%s) + (%s | %s)",
    x$effects$response,
    paste(.effect_terms(x$intercept$fixed, fixed), collapse = " + "),
    paste(.effect_terms(glsw_intercept, glsw), collapse = " + "),
    paste(.effect_terms(x$intercept$slr, slr), collapse = " + "),
    x$effects$group
  ))
}

.simulate_hgwr_error <- function(model, legacy_sigma = FALSE) {
  data0 <- model$frame
  ge <- model$effects$group
  groups <- unique(data0[[ge]])
  z0 <- model$frame.parsed$z
  sigma2 <- if (legacy_sigma) model$sigma else model$sigma^2
  simulated <- numeric(nrow(data0))
  for (gr in groups) {
    index <- which(data0[[ge]] == gr)
    zi <- as.matrix(z0[index, , drop = FALSE])
    simulated[index] <- MASS::mvrnorm(
      mu = rep(0, nrow(zi)),
      Sigma = zi %*% model$D %*% t(zi) + diag(sigma2, nrow(zi))
    )
  }
  simulated
}

#' @describeIn hgwr
#' Test the spatial heterogeneity with bootstrapping.
#'
#' @param x An `hgwrm` object
#' @param round The number of times to sample from the model.
#' @param statistic A function used to calculate the statistics on the
#' original data and bootstrapped data. Defaults to the variance of
#' standardised GLSW estimates.
#' @param null Bootstrap null-generating mechanism. `"fitted"` preserves the
#' historical behaviour centred on the fitted spatially varying model.
#' `"effectwise"` fits a spatially constant null model for each tested effect.
#' @param effects Character vector naming GLSW effects to test. `NULL` tests
#' all GLSW effects.
#' @param refit_bw Whether to re-select the bandwidth in every bootstrap fit.
#' @param p_correction Bootstrap p-value convention. `"none"` returns the
#' empirical proportion; `"plus_one"` adds one to numerator and denominator.
#' @param fail_on_nonconvergence Whether a non-converged null or bootstrap
#' HGWR fit should stop the test. Defaults to `FALSE` for compatibility.
#' @param parallel If TRUE, use `furrr` package to parallelise resamples.
#' @param verbose Override the verbose setting in the `hgwrm` model.
#' @param \dots Arguments to be used in re-fitting the model with
#' [hgwrr::hgwr()].
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats as.formula fitted
#' @method spatial_hetero_test hgwrm
#' @export
spatial_hetero_test.hgwrm <- function(
  x,
  round = 99,
  statistic = stat_glsw,
  null = c("fitted", "effectwise"),
  effects = NULL,
  refit_bw = FALSE,
  p_correction = c("none", "plus_one"),
  fail_on_nonconvergence = FALSE,
  parallel = FALSE,
  verbose = 0,
  ...
) {
  null <- match.arg(null)
  p_correction <- match.arg(p_correction)
  if (!inherits(x, "hgwrm")) {
    stop("Argument x must be an hgwrm object.", call. = FALSE)
  }
  if (length(round) != 1 || !is.numeric(round) || round < 1) {
    stop("round must be a positive integer.", call. = FALSE)
  }
  round <- as.integer(round)

  all_t0 <- with(x, statistic(gamma, gamma_se))
  names(all_t0) <- x$effects$glsw
  if (is.null(effects)) {
    effects <- x$effects$glsw
  }
  if (!is.character(effects) || length(effects) < 1 ||
      any(!effects %in% x$effects$glsw)) {
    stop("effects must name one or more GLSW effects in x.", call. = FALSE)
  }
  effects <- unique(effects)
  effect_index <- match(effects, x$effects$glsw)
  t0 <- all_t0[effect_index]

  args <- as.list(x$call[-1])
  args$coords <- x$coords
  args$data <- x$frame
  args$f_test <- FALSE
  args$max_iters <- 100
  args$max_retries <- 10
  args$bw <- if (refit_bw) "CV" else x$bw
  args$verbose <- verbose
  if (...length() > 0) {
    settings <- list(...)
    for (key in names(settings)) {
      args[[key]] <- settings[[key]]
    }
  }

  pgb <- NULL
  if (requireNamespace("progressr", quietly = TRUE)) {
    nticks <- round * ifelse(null == "effectwise", length(effects), 1)
    pgb <- progressr::progressor(nticks)
  }

  map_bootstrap <- function(null_model, selected_index, legacy_sigma) {
    worker <- function(i, pgb = NULL) {
      datai <- x$frame
      datai[[x$effects$response]] <-
        fitted(null_model) + .simulate_hgwr_error(null_model, legacy_sigma)
      argsi <- args
      argsi$data <- datai
      modeli <- do.call(hgwrr::hgwr_fit, argsi)
      if (isTRUE(fail_on_nonconvergence) && identical(modeli$converged, FALSE)) {
        stop(
          "A bootstrap HGWR fit reached its iteration or retry limit.",
          call. = FALSE
        )
      }
      stati <- with(modeli, statistic(gamma, gamma_se))[selected_index]
      if (!is.null(pgb)) pgb()
      stati
    }
    if (is.logical(parallel) && parallel) {
      if (!requireNamespace("furrr", quietly = TRUE)) {
        stop("Package \"furrr\" must be installed to parallel", call. = FALSE)
      }
      draws <- furrr::future_map(
        seq_len(round), worker,
        pgb = pgb,
        .options = furrr::furrr_options(seed = TRUE)
      )
    } else {
      draws <- lapply(seq_len(round), worker, pgb = pgb)
    }
    matrix(unlist(draws), nrow = round, byrow = TRUE)
  }

  if (null == "fitted") {
    # The historical fitted-model bootstrap uses one simulated data set per
    # round and evaluates every requested GLSW effect on that same refit.
    t <- map_bootstrap(x, effect_index, legacy_sigma = TRUE)
  } else {
    t <- vapply(seq_along(effect_index), function(j) {
      null_args <- args
      null_args$formula <- .effectwise_null_formula(x, effects[j])
      null_model <- do.call(hgwrr::hgwr_fit, null_args)
      if (isTRUE(fail_on_nonconvergence) && identical(null_model$converged, FALSE)) {
        stop(
          "An effectwise null HGWR fit reached its iteration or retry limit.",
          call. = FALSE
        )
      }
      as.numeric(map_bootstrap(null_model, effect_index[j], legacy_sigma = FALSE))
    }, numeric(round))
  }
  if (length(effects) == 1) {
    t <- matrix(t, ncol = 1)
  }
  colnames(t) <- effects

  exceedances <- vapply(seq_along(t0), function(i) {
    if (p_correction == "plus_one") {
      sum(abs(t[, i]) >= abs(t0[i]))
    } else {
      sum(abs(t[, i]) > abs(t0[i]))
    }
  }, numeric(1))
  if (p_correction == "plus_one") {
    pv <- (1 + exceedances) / (round + 1)
  } else {
    pv <- exceedances / round
  }
  names(pv) <- effects

  res <- list(
    t0 = t0,
    t = t,
    vars = effects,
    pv = pv,
    null = null,
    round = round,
    refit_bw = refit_bw,
    p_correction = p_correction,
    fail_on_nonconvergence = fail_on_nonconvergence
  )
  class(res) <- "spahetbootres"
  res
}
