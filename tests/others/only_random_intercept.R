data(mulsam.test)
m <- hgwr(
  formula = y ~ L(g1 + g2) + x1 + (1 | group),
  data = mulsam.test$data,
  coords = mulsam.test$coords,
  alpha = 1e-8,
  verbose = 1
)
