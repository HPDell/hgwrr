library(testthat)

test_that("Multisampling", {
    data(multisampling)
    m <- expect_no_error({hgwr(
        formula = y ~ L(g1 + g2) + x1 + (z1 | group),
        data = multisampling$data,
        coords = multisampling$coords,
        bw = 10,
        alpha = 1e-8
    )})
    expect_no_error({
        summary(m, test_hetero = T)
    })
    expect_no_error({
        summary(m, test_hetero = list(kernel = "gaussian"))
    })
})