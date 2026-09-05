data(mulsam.test)

fit_optimizer_fixture <- function(alpha, ml_type = "D_Only") {
  hgwr(
    y ~ L(g1 + g2) + x1 + (z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = alpha,
    eps_iter = 1e-5,
    eps_gradient = 1e-5,
    max_iters = 200,
    max_retries = 10,
    ml_type = ml_type
  )
}

test_that("covariance optimizer returns a positive-definite matrix and diagnostics", {
  model <- fit_optimizer_fixture(0.01)
  expect_true(all(is.finite(model$D)))
  expect_true(isTRUE(model$min_eigen_D > 0))
  expect_equal(model$min_eigen_D, min(eigen(model$D, symmetric = TRUE)$values), tolerance = 1e-8)
  expect_type(model$ml_converged, "logical")
  expect_type(model$ml_status, "integer")
  expect_true(is.numeric(model$ml_iterations))
  expect_true(is.numeric(model$ml_failures))
  expect_false(isTRUE(all.equal(model$D, diag(nrow(model$D)), tolerance = 1e-6)))
})

test_that("D-only optimizer is insensitive to reasonable initial step sizes", {
  models <- lapply(c(0.001, 0.01, 0.1), fit_optimizer_fixture)
  expect_true(all(vapply(models, function(x) all(is.finite(x$D)), logical(1))))
  expect_true(all(vapply(models, function(x) x$min_eigen_D > 0, logical(1))))
  objectives <- vapply(models, function(x) -x$logLik / nrow(x$frame), numeric(1))
  expect_lt(max(objectives) - min(objectives), 1e-4)
  reference <- models[[2]]$D
  relative_difference <- vapply(
    models,
    function(x) norm(x$D - reference, type = "F") / max(norm(reference, type = "F"), 1e-12),
    numeric(1)
  )
  expect_true(all(relative_difference < 0.05))
})

test_that("joint D-beta optimizer stays finite and positive definite", {
  model <- fit_optimizer_fixture(0.01, ml_type = "D_Beta")
  expect_true(all(is.finite(model$D)))
  expect_true(all(is.finite(model$beta)))
  expect_true(model$min_eigen_D > 0)
  expect_true(is.finite(model$logLik))
})
