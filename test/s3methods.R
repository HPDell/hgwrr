library(hgwrr)
load("R/data/multisampling.rda")
hgwr_formula <- y ~ g1 + g2 + x1 + (z1 | group)
result <- hgwr(hgwr_formula, multisampling$data, c("g1", "g2"), multisampling$coord, 10)
result

coef(result)
coefficients(result)
