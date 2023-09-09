d <- data.frame(
  years = 1:5,
  number = rnorm(5),
  level = c("top", "mid", "low", "mid", "top"),
  distance = factor(c("far", "near", "near", "far", "far")),
  sold = c(TRUE, TRUE, FALSE, FALSE, TRUE)
)

test_that("make dummy", {
    dummy_d <- expect_no_error(make.dummy(d))
    expect_setequal(
        names(dummy_d),
        c("years", "number", "level.top", "level.mid",
          "distance.far", "sold")
    )
})

