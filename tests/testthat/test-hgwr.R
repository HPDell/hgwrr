data(mulsam.test)
m <- NULL

test_that("HGWR reports bounded-iteration diagnostics", {
  model <- hgwr(
    y ~ L(g1 + g2) + x1 + (z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = 1e-8,
    max_iters = 20,
    max_retries = 10
  )
  expect_true(is.numeric(model$iterations))
  expect_true(is.numeric(model$retries))
  expect_type(model$converged, "logical")
  expect_lte(model$iterations, 20L)
  expect_lte(model$retries, 10L)
})

test_that("hgwr fit", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = 10,
      alpha = 1e-8
    )
  })
})

test_that("hgwr fit no intercept", {
  data(mulsam.test)
  expect_no_error(hgwr(
    formula = y ~ L(0 + g1 + g2) + x1 + (z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = 1e-8
  ))
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + 0 + x1 + (z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = 1e-8
  ))
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (0 + z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr fit sf", {
  data(mulsam.test)
  ms_sf <- with(mulsam.test, cbind(coords[data$group, ], data))
  ms_sf <- sf::st_as_sf(ms_sf, coords = 1:2)
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = ms_sf,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr parse formula", {
  expect_setequal(m$effects$glsw, c("Intercept", "g1", "g2"))
  expect_setequal(m$effects$fixed, c("Intercept", "x1"))
  expect_setequal(m$effects$slr, c("Intercept", "z1"))
  expect_equal(m$effects$group, c("group"))
  expect_equal(m$effects$response, c("y"))
})

test_that("hgwr s3 methods", {
  expect_no_error(print(m))
  expect_no_error(print(m, table.style = "md"))
  expect_no_error(coef(m))
  expect_no_error(coefficients(m))
  expect_no_error(fitted(m))
  expect_no_error(residuals(m))
  expect_no_error(summary(m))
  expect_no_error(print(summary(m), table.style = "md"))
})

test_that("hgwr data.frame coords check", {
  expect_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = mulsam.test$data,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr bandwidth optimisation", {
  expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
})

test_that("hgwr summary with spatial heterogeneity test", {
  expect_no_error({
    summary(m, test_hetero = TRUE)
  })
})

test_that("single-model F3 experiment returns aligned diagnostics", {
  model <- hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = mulsam.test$data,
    coords = mulsam.test$coords,
    bw = 10,
    alpha = 1e-8,
    f_test = TRUE
  )
  variants <- c(
    "f_test", "f_test_scale", "f_test_nuisance",
    "f_test_combined", "f_test_orthogonal"
  )
  for (variant in variants) {
    expect_identical(names(model[[variant]]), model$effects$glsw)
    expect_true(all(vapply(model[[variant]], function(x) {
      length(x) == 4L && all(is.finite(x)) && x[[4L]] >= 0 && x[[4L]] <= 1
    }, logical(1))))
  }
  expect_identical(
    rownames(model$f_test_diagnostics),
    model$effects$glsw
  )
  expect_true(all(is.finite(model$f_test_diagnostics)))
  expect_lt(max(model$f_test_diagnostics[, "map_identity_error"]), 1e-8)
})

test_that("hgwr s3 methods with no random slop or intercept", {
  m_rn <- NULL
  expect_no_error({
    m_rn <<- hgwr(
      formula = y ~ L(g1 + g2) + x1 + (1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_rn))
  expect_no_error(print(m_rn, table.style = "md"))
  expect_no_error(coef(m_rn))
  expect_no_error(coefficients(m_rn))
  expect_no_error(fitted(m_rn))
  expect_no_error(residuals(m_rn))
  expect_no_error(summary(m_rn))
  expect_no_error(print(summary(m_rn), table.style = "md"))

  expect_no_error({
    m_rn <<- hgwr(
      formula = y ~ L(g1 + g2) + x1 + (0 + z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_rn))
  expect_no_error(print(m_rn, table.style = "md"))
  expect_no_error(coef(m_rn))
  expect_no_error(coefficients(m_rn))
  expect_no_error(fitted(m_rn))
  expect_no_error(residuals(m_rn))
  expect_no_error(summary(m_rn))
  expect_no_error(print(summary(m_rn), table.style = "md"))
})

test_that("hgwr s3 methods with no fixed slop or intercept", {
  m_fn <- NULL
  expect_no_error({
    m_fn <<- hgwr(
      formula = y ~ L(g1 + g2) + 1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_fn))
  expect_no_error(print(m_fn, table.style = "md"))
  expect_no_error(coef(m_fn))
  expect_no_error(coefficients(m_fn))
  expect_no_error(fitted(m_fn))
  expect_no_error(residuals(m_fn))
  expect_no_error(summary(m_fn))
  expect_no_error(print(summary(m_fn), table.style = "md"))

  expect_no_error({
    m_fn <<- hgwr(
      formula = y ~ L(g1 + g2) + 0 + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_fn))
  expect_no_error(print(m_fn, table.style = "md"))
  expect_no_error(coef(m_fn))
  expect_no_error(coefficients(m_fn))
  expect_no_error(fitted(m_fn))
  expect_no_error(residuals(m_fn))
  expect_no_error(summary(m_fn))
  expect_no_error(print(summary(m_fn), table.style = "md"))
})

test_that("hgwr s3 methods with no GLSW slop", {
  m_gn <- NULL
  expect_no_error({
    m_gn <<- hgwr(
      formula = y ~ L(1) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_gn))
  expect_no_error(print(m_gn, table.style = "md"))
  expect_no_error(coef(m_gn))
  expect_no_error(coefficients(m_gn))
  expect_no_error(fitted(m_gn))
  expect_no_error(residuals(m_gn))
  expect_no_error(summary(m_gn))
  expect_no_error(print(summary(m_gn), table.style = "md"))

  expect_no_error({
    m_gn <<- hgwr(
      formula = y ~ L(0 + g1 + g2) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = "CV",
      alpha = 1e-8
    )
  })
  expect_no_error(print(m_gn))
  expect_no_error(print(m_gn, table.style = "md"))
  expect_no_error(coef(m_gn))
  expect_no_error(coefficients(m_gn))
  expect_no_error(fitted(m_gn))
  expect_no_error(residuals(m_gn))
  expect_no_error(summary(m_gn))
  expect_no_error(print(summary(m_gn), table.style = "md"))
})

test_that("hgwr order data", {
  m0 <- expect_no_error({
    hgwr(
      formula = y ~ L(0 + g1 + g2) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = 10,
      alpha = 0.01
    )
  })
  set.seed(1)
  data_perm <- with(
    mulsam.test,
    data[order(runif(nrow(data))), ]
  )
  m_perm <- expect_no_error({
    hgwr(
      formula = y ~ L(0 + g1 + g2) + x1 + (z1 | group),
      data = data_perm,
      coords = mulsam.test$coords,
      bw = 10,
      alpha = 0.01
    )
  })
  expect_equal(m_perm$gamma, m0$gamma)
  expect_equal(m_perm$beta, m0$beta)
  expect_equal(m_perm$mu, m0$mu)
  expect_equal(m_perm$D, m0$D)
  expect_equal(m_perm$sigma, m0$sigma)
  expect_equal(m_perm$bw, m0$bw)
  expect_equal(m_perm$gamma_se, m0$gamma_se)
  expect_equal(m_perm$logLik, m0$logLik)
  expect_equal(m_perm$trS, m0$trS)
  expect_equal(m_perm$var_beta, m0$var_beta)
})
