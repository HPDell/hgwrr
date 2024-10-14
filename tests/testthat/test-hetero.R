library(testthat)

test_that("spatial heterogeneity: vector matrix data.frame", {
  data(multisampling)
  g <- with(multisampling, {
    aggregate(data[c("g1", "g2")], by = list(data$group), mean)
  })[,-1]
  expect_no_error({
    spatial_hetero_test_data(g, as.matrix(multisampling$coords))
  })
  expect_no_error({
    spatial_hetero_test(as.matrix(g), multisampling$coords)
  })
  expect_no_error({
    spatial_hetero_test(g[["g1"]], multisampling$coords)
  })
  expect_no_error({
    spatial_hetero_test(g, multisampling$coords)
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
  data(multisampling)
  m <- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = 10,
      alpha = 1e-8
    )
  })
  expect_no_error({
    spatial_hetero_test(m)
  })
})
