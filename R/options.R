ml_options <- function(
  alpha = 0.01,
  eps_gradient = 1e-6,
  ml_type = c("D_Only", "D_Beta")
) list(
  method = "ML",
  alpha = 0.01,
  eps_gradient = 1e-6,
  ml_type = match.arg(ml_type)
)

mcmc_options <- function(
  mcmc_iters = 10000,
  mcmc_burnin = 2000
) list(
  method = "MCMC",
  mcmc_iters = 10000,
  mcmc_burnin = 2000
)