data(multisampling)
m <- NULL

test_that("hgwr fit", {
  expect_no_error({
    summary(hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = "AIC",
      alpha = 1e-4,
      verbose = 2
    ))
  })
})

test_that("hgwr multiscale fit (CV)", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = "CV",
      multiscale = TRUE,
      alpha = 1e-2,
      verbose = 1
    )
  })
  summary(m)
})

test_that("hgwr multiscale fit (AIC)", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = "AIC",
      multiscale = TRUE,
      alpha = 1e-2,
      verbose = 1
    )
  })
  summary(m)
})