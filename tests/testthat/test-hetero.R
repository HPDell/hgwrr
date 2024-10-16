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
})
