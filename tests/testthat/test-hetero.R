library(testthat)

test_that("spatial heterogeneity: vector matrix data.frame", {
  data(mulsam.test)
  g <- with(mulsam.test, {
    aggregate(data[c("g1", "g2")], by = list(data$group), mean)
  })[,-1]
  expect_no_error({
    spatial_hetero_test_data(g, as.matrix(mulsam.test$coords))
  })
  expect_no_error({
    spatial_hetero_test(as.matrix(g), mulsam.test$coords)
  })
  expect_no_error({
    spatial_hetero_test(g[["g1"]], mulsam.test$coords)
  })
  expect_no_error({
    spatial_hetero_test(g, mulsam.test$coords)
  })
})

test_that("spatial heterogeneity: sf", {
  data(wuhan.hp)
  g <- aggregate(wuhan.hp, list(wuhan.hp$group), mean)[1:16, -1]
  g <- g[c("d.Commercial", "d.GreenLand")]
  expect_no_error({
    spatial_hetero_test(g)
  })
})

test_that("spatial heterogeneity: HGWR", {
  data(mulsam.test)
  m <- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = mulsam.test$data,
      coords = mulsam.test$coords,
      bw = 10,
      alpha = 1e-8
    )
  })
  expect_no_error({
    spatial_hetero_test(m)
  })

  fullmodel_formula <- hgwrr:::.fullmodel_null_formula(m)
  fullmodel_terms <- hgwrr:::parse_formula(fullmodel_formula)
  expect_setequal(fullmodel_terms$fixed.effects, c("g1", "g2", "x1"))
  expect_length(fullmodel_terms$glsw.effects, 0)

  set.seed(123)
  fullmodel <- spatial_hetero_test(
    m,
    round = 2,
    effects = c("g1", "g2"),
    null = "fullmodel",
    p_correction = "plus_one"
  )
  expect_s3_class(fullmodel, "spahetbootres")
  expect_identical(fullmodel$null, "fullmodel")
  expect_identical(fullmodel$vars, c("g1", "g2"))
  expect_identical(dim(fullmodel$t), c(2L, 2L))
  expect_true(all(fullmodel$pv >= 1 / 3))

  parsed <- m$frame.parsed
  random_fit <- rowSums(parsed$z * m$mu[parsed$group, ])
  expect_equal(
    fitted(m) - hgwrr:::.hgwr_null_mean(m),
    random_fit
  )

  frame_before <- m$frame
  coords_before <- m$coords
  model_with_other_mu <- m
  model_with_other_mu$mu <- model_with_other_mu$mu + 100
  set.seed(321)
  bootstrap_data <- hgwrr:::.bootstrap_hgwr_data(
    m$frame, m, m$effects$response, legacy_sigma = FALSE
  )
  set.seed(321)
  bootstrap_data_other_mu <- hgwrr:::.bootstrap_hgwr_data(
    m$frame, model_with_other_mu, m$effects$response, legacy_sigma = FALSE
  )
  expect_identical(
    bootstrap_data[[m$effects$response]],
    bootstrap_data_other_mu[[m$effects$response]]
  )
  unchanged <- setdiff(names(m$frame), m$effects$response)
  expect_identical(bootstrap_data[unchanged], m$frame[unchanged])
  expect_identical(m$frame, frame_before)
  expect_identical(m$coords, coords_before)

  set.seed(123)
  legacy <- spatial_hetero_test(
    m,
    round = 3,
    effects = "g2",
    null = "fitted",
    p_correction = "none"
  )
  expect_s3_class(legacy, "spahetbootres")
  expect_identical(legacy$vars, "g2")
  expect_identical(dim(legacy$t), c(3L, 1L))
  expect_no_error(print(legacy))

  set.seed(123)
  legacy_both <- spatial_hetero_test(
    m,
    round = 2,
    effects = c("g1", "g2"),
    null = "fitted"
  )
  expect_identical(dim(legacy_both$t), c(2L, 2L))
  expect_identical(colnames(legacy_both$t), c("g1", "g2"))

  set.seed(123)
  corrected <- spatial_hetero_test(
    m,
    round = 3,
    effects = "g2",
    null = "effectwise",
    p_correction = "plus_one"
  )
  expect_s3_class(corrected, "spahetbootres")
  expect_identical(corrected$null, "effectwise")
  expect_identical(corrected$vars, "g2")
  expect_gte(corrected$pv[[1]], 1 / 4)
  expect_lte(corrected$pv[[1]], 1)

  set.seed(123)
  refitted_bw <- spatial_hetero_test(
    m,
    round = 1,
    effects = "g2",
    null = "effectwise",
    refit_bw = TRUE,
    max_iters = 100,
    max_retries = 10
  )
  expect_true(refitted_bw$refit_bw)
  expect_identical(dim(refitted_bw$t), c(1L, 1L))

  expect_error(
    spatial_hetero_test(m, effects = "missing"),
    "effects must name"
  )
})
