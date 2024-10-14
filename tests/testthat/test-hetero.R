library(testthat)

test_that("spatial heterogeneity: data.frame", {
  data(multisampling)
  expect_no_error({
    with(multisampling, spatial_hetero_test(data, coords))
  })
})

test_that("spatial heterogeneity: sf", {
  data(wuhan.hp)
  expect_no_error({
    spatial_hetero_test(wuhan.hp)
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
