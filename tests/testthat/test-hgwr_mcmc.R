data(multisampling)
m <- NULL

test_that("hgwr fit", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = 10,
      eps_iter = 1e-3,
      max_iters = 20,
      global_options = mcmc_options(
        mcmc_iters = 120,
        mcmc_burnin = 40,
        save_draws = TRUE
      )
    )
  })
})

test_that("local intercept is constrained without projecting random effects", {
  group_weights <- as.numeric(table(m$frame.parsed$group))
  expect_equal(weighted.mean(m$gamma[, "Intercept"], group_weights), 0,
               tolerance = 1e-8)
  expect_gt(abs(weighted.mean(m$mu[, "Intercept"], group_weights)), 1e-6)
})

test_that("fit diagnostics and final conditional draws are returned", {
  expect_named(
    m$fit_diagnostics,
    c(
      "converged", "stop_reason", "outer_iterations",
      "total_inner_updates", "final_inner_iterations", "inner_converged",
      "all_inner_converged", "inner_failures",
      "max_inner_iterations_used",
      "max_parameter_change", "optimizer_status",
      "optimizer_iterations", "optimizer_measure", "optimizer_failures",
      "optimizer_restarts", "bw_optimizer_status",
      "bw_optimizer_failures", "bw_lower", "bw_upper",
      "bw_objective", "bw_evaluations",
      "effective_neighbors", "spd_corrections", "inverse_jitter_uses",
      "pseudoinverse_uses", "cholesky_jitter_uses"
    )
  )
  expect_gte(nrow(m$conditional_draws$beta), 80)
  expect_equal(nrow(m$conditional_draws$D), nrow(m$conditional_draws$beta))
  expect_equal(length(m$conditional_draws$sigma2), nrow(m$conditional_draws$beta))
  expect_true(all(is.finite(m$conditional_draws$beta)))
  expect_true(m$fit_diagnostics$converged)
  expect_identical(m$fit_diagnostics$stop_reason, "tolerance_reached")
  expect_lte(
    m$fit_diagnostics$max_parameter_change,
    sqrt(1e-3)
  )
})

test_that("conditional-mode comparator is deterministic and converges internally", {
  fit_cm <- function() hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = 10,
    eps_iter = 1e-3,
    max_iters = 20,
    global_options = cm_options(max_inner_iters = 100, inner_tolerance = 1e-6)
  )
  cm_a <- expect_no_error(fit_cm())
  cm_b <- expect_no_error(fit_cm())
  expect_equal(cm_a$beta, cm_b$beta)
  expect_equal(cm_a$D, cm_b$D)
  expect_equal(cm_a$gamma, cm_b$gamma)
  expect_true(cm_a$fit_diagnostics$inner_converged)
  expect_true(cm_a$fit_diagnostics$all_inner_converged)
  expect_identical(cm_a$fit_diagnostics$inner_failures, 0)
  expect_lte(cm_a$fit_diagnostics$final_inner_iterations, 100)
  expect_lte(cm_a$fit_diagnostics$max_inner_iterations_used, 100)
  expect_gte(
    cm_a$fit_diagnostics$max_inner_iterations_used,
    cm_a$fit_diagnostics$final_inner_iterations
  )
})

test_that("invalid prior settings are rejected", {
  expect_error(
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = 10,
      global_options = mcmc_options(
        mcmc_iters = 20,
        mcmc_burnin = 10,
        d_prior_df_offset = 1
      )
    ),
    "d_prior_df_offset"
  )
})

test_that("hgwr bandwidth optimisation", {
  fit_bw <- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = "CV",
      eps_iter = 1e-3,
      max_iters = 10,
      global_options = mcmc_options(
        mcmc_iters = 80,
        mcmc_burnin = 30
      )
    )
  })
  expect_true(
    fit_bw$fit_diagnostics$stop_reason %in% c("tolerance_reached", "max_iters")
  )
  expect_equal(fit_bw$fit_diagnostics$bw_optimizer_failures, 0)
  expect_equal(fit_bw$bw, as.integer(fit_bw$bw))
  expect_gt(fit_bw$fit_diagnostics$bw_evaluations, 0)
})
