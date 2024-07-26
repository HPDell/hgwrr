data(multisampling.large)
m <- hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = multisampling.large$data,
    coords = multisampling.large$coords,
    alpha = 1e-8,
    verbose = 1
)
