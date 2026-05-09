ml_options <- function(
  alpha = 0.01,
  eps_gradient = 1e-6,
  ml_type = c("D_Only", "D_Beta")
) list(
  method = "ML",
  alpha = alpha,
  eps_gradient = eps_gradient,
  ml_type = match.arg(ml_type)
)

mcmc_options <- function(
  mcmc_iters = 5000,
  mcmc_burnin = 400
) list(
  method = "MCMC",
  mcmc_iters = mcmc_iters,
  mcmc_burnin = mcmc_burnin
)