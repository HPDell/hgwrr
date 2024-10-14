data(multisampling)
m <- hgwr(
  formula = y ~ L(g1 + g2) + x1 + (1 | group),
  data = multisampling$data,
  coords = multisampling$coords,
  alpha = 1e-8,
  verbose = 1
)
