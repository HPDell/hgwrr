library(hgwrr)

data_dir <- Sys.getenv("TEST_DATA_DIR")
data_group <- read.csv(file.path(data_dir, "hlmgwr_group.csv"), header = F)[[1]] + 1
data_X <- read.csv(file.path(data_dir, "hlmgwr_x.csv"), header = F)[, -1]
data_G <- read.csv(file.path(data_dir, "hlmgwr_g.csv"), header = F)[data_group, -1]
data_Z <- read.csv(file.path(data_dir, "hlmgwr_z.csv"), header = F)[, -1]
data_u <- read.csv(file.path(data_dir, "hlmgwr_u.csv"), header = F)
data_y <- read.csv(file.path(data_dir, "hlmgwr_y.csv"), header = F)

hgwr_data <- cbind(data_y, data_G, data_X, data_Z, data_group)
colnames(hgwr_data) <- c("y", "g1", "g2", "x1", "z1", "group")
rownames(hgwr_data) <- NULL
hgwr_coords <- as.matrix(data_u)
hgwr_formula <- y ~ g1 + g2 + x1 + (z1 | group)

result <- hgwr(hgwr_formula, hgwr_data, c("g1", "g2"), hgwr_coords, "CV", kernel = "bisquared", verbose = 1, max_retries = 10000)
fitted <- result$fitted
y <- hgwr_data$y
residual <- y - fitted
rsquared <- 1 - sum(residual^2) / sum((y - mean(y))^2)
cat("Rsquared:", rsquared, "\n")
