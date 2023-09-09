m <- NULL

test_that("hgwr fit", {
  data(multisampling)
  m <<- expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = 10
  ))
})

test_that("hgwr fit sf", {
  data(multisampling)
  ms_sf <- with(multisampling, cbind(coords[data$group, ], data))
  ms_sf <- sf::st_as_sf(ms_sf, coords = 1:2)
  expect_no_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = ms_sf,
    bw = 10
  ))
})

test_that("hgwr parse formula", {
  expect_setequal(m$effects$local.fixed, c("g1", "g2"))
  expect_setequal(m$effects$global.fixed, c("x1"))
  expect_setequal(m$effects$random, c("z1"))
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
  data(multisampling)
  expect_error(hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = multisampling$data,
    bw = 10
  ))
})
