devtools::load_all()
data(multisampling)
hgwr(
    formula = y ~ L(g1 + g2) + x1 + (z1 | group),
    data = multisampling$data,
    coords = multisampling$coords,
    bw = "CV",
    multiscale = TRUE,
    alpha = 1e-2,
    verbose = 1
)