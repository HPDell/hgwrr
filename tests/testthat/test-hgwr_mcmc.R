data(multisampling)
m <- NULL

test_that("hgwr fit", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = 10,
      global_options = mcmc_options(),
      verbose = 1
    )
  })
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
