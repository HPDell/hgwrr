data(multisampling)
m <- NULL

test_that("hgwr fit", {
  m <<- expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
      bw = 10,
      alpha = 1e-8
    )
  })
})

test_that("hgwr fit no intercept", {
  data(multisampling)
  expect_no_error(hgwr(
    formula = y ~ L(0 + g1 + g2) + x1 + (z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = 10,
    alpha = 1e-8
  ))
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + 0 + x1 + (z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = 10,
    alpha = 1e-8
  ))
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (0 + z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr fit sf", {
  data(multisampling)
  ms_sf <- with(multisampling, cbind(coords[data$group, ], data))
  ms_sf <- sf::st_as_sf(ms_sf, coords = 1:2)
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = ms_sf,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr parse formula", {
  expect_setequal(m$effects$local.fixed, c("Intercept", "g1", "g2"))
  expect_setequal(m$effects$global.fixed, c("Intercept", "x1"))
  expect_setequal(m$effects$random, c("Intercept", "z1"))
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
    data = multisampling$data,
    bw = 10,
    alpha = 1e-8
  ))
})

test_that("hgwr bandwidth optimisation", {
  expect_no_error({
    hgwr(
      formula = y ~ L(g1 + g2) + x1 + (z1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
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

test_that("hgwr s3 methods with no random slop or intercept", {
  m_rn <- NULL
  expect_no_error({
    m_rn <<- hgwr(
      formula = y ~ L(g1 + g2) + x1 + (1 | group),
      data = multisampling$data,
      coords = multisampling$coords,
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
      data = multisampling$data,
      coords = multisampling$coords,
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
      data = multisampling$data,
      coords = multisampling$coords,
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
      data = multisampling$data,
      coords = multisampling$coords,
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
      data = multisampling$data,
      coords = multisampling$coords,
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
      data = multisampling$data,
      coords = multisampling$coords,
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
