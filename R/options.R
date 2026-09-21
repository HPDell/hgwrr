#' Maximum likelihood options
#'
#' @param alpha Initial step size used by the maximum likelihood optimizer.
#' @param eps_gradient Gradient convergence threshold for the optimizer.
#' @param ml_type Maximum likelihood variant. `"D_Only"` estimates the random
#'   effect covariance with `beta` updated by GLS; `"D_Beta"` estimates `D` and
#'   `beta` jointly.
#'
#' @return A list passed to [hgwr()] through `global_options`.
#'
#' @export
ml_options <- function(
  alpha = 0.01,
  eps_gradient = 1e-6,
  ml_type = c("D_Only", "D_Beta")
) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("`alpha` must be a finite positive scalar.")
  }
  if (length(eps_gradient) != 1L || !is.finite(eps_gradient) ||
      eps_gradient <= 0) {
    stop("`eps_gradient` must be a finite positive scalar.")
  }
  list(
    method = "ML",
    alpha = alpha,
    eps_gradient = eps_gradient,
    ml_type = match.arg(ml_type)
  )
}

#' Markov Chain Monte Carlo estimator options
#'
#' @param mcmc_iters Total MCMC iterations within each backfitting iteration.
#' @param mcmc_burnin Number of initial MCMC iterations discarded as burn-in.
#' @param beta_prior_variance Scalar variance multiplier in
#'   `beta | sigma^2 ~ N(0, sigma^2 beta_prior_variance I)`.
#' @param d_prior_df_offset Degrees-of-freedom offset in
#'   `D ~ IW(q + d_prior_df_offset, S0)`. It must exceed 1 so the prior mean
#'   exists.
#' @param d_prior_mean Diagonal value of the inverse-Wishart prior mean.
#' @param sigma_prior_shape,sigma_prior_scale Shape and scale in
#'   `sigma^2 ~ Inv-Gamma(shape, scale)`.
#' @param save_draws Save retained draws from the last conditional chain. These
#'   draws condition on the final spatial surface and are not draws from a
#'   joint HGWR posterior.
#' @param inner_tolerance Parameter-change tolerance used only by conditional
#'   mode backfitting. It is ignored for stochastic MCMC updates.
#'
#' @return A list passed to [hgwr()] through `global_options`.
#'
#' @export
mcmc_options <- function(
  mcmc_iters = 5000,
  mcmc_burnin = 400,
  beta_prior_variance = 100,
  d_prior_df_offset = 5,
  d_prior_mean = 1,
  sigma_prior_shape = 0.001,
  sigma_prior_scale = 0.001,
  inner_tolerance = 1e-6,
  save_draws = FALSE
) list(
  method = "MCMC",
  mcmc_iters = mcmc_iters,
  mcmc_burnin = mcmc_burnin,
  beta_prior_variance = beta_prior_variance,
  d_prior_df_offset = d_prior_df_offset,
  d_prior_mean = d_prior_mean,
  sigma_prior_shape = sigma_prior_shape,
  sigma_prior_scale = sigma_prior_scale,
  inner_tolerance = inner_tolerance,
  save_draws = save_draws
)

#' Conditional-mode backfitting options
#'
#' Uses the same conditional hierarchical model and priors as [mcmc_options()],
#' but replaces random draws with coordinate-wise conditional modes. This
#' deterministic comparator separates prior regularization from Monte Carlo
#' averaging; it is not a joint MAP estimator.
#'
#' @param max_inner_iters Maximum number of conditional-mode coordinate updates
#'   within each backfitting iteration.
#' @inheritParams mcmc_options
#'
#' @return A list passed to [hgwr()] through `global_options`.
#'
#' @export
cm_options <- function(
  max_inner_iters = 100,
  inner_tolerance = 1e-6,
  beta_prior_variance = 100,
  d_prior_df_offset = 5,
  d_prior_mean = 1,
  sigma_prior_shape = 0.001,
  sigma_prior_scale = 0.001
) list(
  method = "CM",
  mcmc_iters = max_inner_iters,
  mcmc_burnin = 0L,
  beta_prior_variance = beta_prior_variance,
  d_prior_df_offset = d_prior_df_offset,
  d_prior_mean = d_prior_mean,
  sigma_prior_shape = sigma_prior_shape,
  sigma_prior_scale = sigma_prior_scale,
  inner_tolerance = inner_tolerance,
  save_draws = FALSE
)

#' @rdname cm_options
#' @param inner_iters Deprecated alias for `max_inner_iters`.
#' @export
map_options <- function(
  inner_iters = 100,
  beta_prior_variance = 100,
  d_prior_df_offset = 5,
  d_prior_mean = 1,
  sigma_prior_shape = 0.001,
  sigma_prior_scale = 0.001,
  inner_tolerance = 1e-6
) {
  warning("`map_options()` is deprecated; use `cm_options()`.", call. = FALSE)
  cm_options(
    max_inner_iters = inner_iters,
    inner_tolerance = inner_tolerance,
    beta_prior_variance = beta_prior_variance,
    d_prior_df_offset = d_prior_df_offset,
    d_prior_mean = d_prior_mean,
    sigma_prior_shape = sigma_prior_shape,
    sigma_prior_scale = sigma_prior_scale
  )
}
