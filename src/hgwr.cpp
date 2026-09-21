// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "hlmgwr.h"
#include "utils.h"

using namespace std;
using namespace Rcpp;
using namespace hgwr;

List fit_diagnostics(HGWR& algorithm)
{
    return List::create(
        Named("converged") = algorithm.get_converged(),
        Named("stop_reason") = algorithm.get_stop_reason(),
        Named("outer_iterations") = algorithm.get_outer_iterations(),
        Named("total_inner_updates") = algorithm.get_total_inner_updates(),
        Named("final_inner_iterations") = algorithm.get_final_inner_iterations(),
        Named("inner_converged") = algorithm.get_inner_converged(),
        Named("all_inner_converged") = algorithm.get_all_inner_converged(),
        Named("inner_failures") = algorithm.get_inner_failures(),
        Named("max_inner_iterations_used") =
            algorithm.get_max_inner_iterations_used(),
        Named("max_parameter_change") = algorithm.get_max_parameter_change(),
        Named("optimizer_status") = algorithm.get_optimizer_status(),
        Named("optimizer_iterations") = algorithm.get_optimizer_iterations(),
        Named("optimizer_measure") = algorithm.get_optimizer_measure(),
        Named("optimizer_failures") = algorithm.get_optimizer_failures(),
        Named("optimizer_restarts") = algorithm.get_optimizer_restarts(),
        Named("bw_optimizer_status") = algorithm.get_bw_optimizer_status(),
        Named("bw_optimizer_failures") = algorithm.get_bw_optimizer_failures(),
        Named("bw_lower") = algorithm.get_bw_lower(),
        Named("bw_upper") = algorithm.get_bw_upper(),
        Named("bw_objective") = algorithm.get_bw_objective(),
        Named("bw_evaluations") = algorithm.get_bw_evaluations(),
        Named("effective_neighbors") = std::ceil(algorithm.get_bw()),
        Named("spd_corrections") = algorithm.get_spd_corrections(),
        Named("inverse_jitter_uses") = algorithm.get_inverse_jitter_uses(),
        Named("pseudoinverse_uses") = algorithm.get_pseudoinverse_uses(),
        Named("cholesky_jitter_uses") = algorithm.get_cholesky_jitter_uses()
    );
}

// [[Rcpp::export]]
List hgwr_bfml(
    const arma::mat& g,
    const arma::mat& x,
    const arma::mat& z,
    const arma::vec& y,
    const arma::mat& u,
    const arma::vec& group,
    double bw,
    int bw_optim,
    size_t kernel,
    double alpha,
    double eps_iter,
    double eps_gradient,
    size_t max_iters,
    size_t max_retries,
    size_t ml_type,
    bool f_test,
    size_t verbose
) {
    arma::uvec mgroup = arma::conv_to<arma::uvec>::from(group) - 1;
    auto mkernel = HGWR::KernelType(size_t(kernel));
    HGWR::Options options { alpha, eps_iter, eps_gradient, max_iters, max_retries, verbose, ml_type };
    HGWR algorithm(g, x, z, y, u, mgroup, mkernel, options);
    if (bw_optim < 0) {
        algorithm.set_bw(bw);
    } else {
        algorithm.set_bw_optim(true);
        algorithm.set_bw_criterion_type(HGWR::BwOptimCriterionType(bw_optim));
    }
    algorithm.set_printer(&prcout);
    algorithm.set_canceler(&prcancel);
    auto hgwr_result = algorithm.fit(f_test);

    auto results = List::create(
        Named("gamma") = hgwr_result.gamma,
        Named("beta") = hgwr_result.beta,
        Named("mu") = hgwr_result.mu,
        Named("D") = hgwr_result.D,
        Named("sigma") = hgwr_result.sigma,
        Named("bw") = hgwr_result.bw,
        Named("gamma_se") = algorithm.get_gamma_se(),
        Named("logLik") = algorithm.get_loglik(),
        Named("trS") = algorithm.get_trS(),
        Named("var_beta") = algorithm.get_var_beta(),
        Named("edf") = algorithm.edf(),
        Named("enp") = algorithm.enp(),
        Named("fit_diagnostics") = fit_diagnostics(algorithm)
    );

    if (f_test)
    {
        auto ftest_result = algorithm.test_glsw();
        results["f_test"] = ftest_result;
    }

    return results;
}

// [[Rcpp::export]]
List hgwr_mcmc(
    const arma::mat& g,
    const arma::mat& x,
    const arma::mat& z,
    const arma::vec& y,
    const arma::mat& u,
    const arma::vec& group,
    double bw,
    int bw_optim,
    size_t kernel,
    double eps_iter,
    size_t max_iters,
    size_t max_retries,
    size_t mcmc_iters,
    size_t mcmc_burnin,
    double beta_prior_variance,
    double d_prior_df_offset,
    double d_prior_mean,
    double sigma_prior_shape,
    double sigma_prior_scale,
    double inner_tolerance,
    bool save_draws,
    bool conditional_mode,
    bool f_test,
    size_t verbose
) {
    arma::uvec mgroup = arma::conv_to<arma::uvec>::from(group) - 1;
    auto mkernel = HGWR::KernelType(size_t(kernel));
    HGWR::Options options;
    options.eps_iter = eps_iter;
    options.max_iters = max_iters;
    options.max_retries = max_retries;
    options.verbose = verbose;
    HGWR algorithm(g, x, z, y, u, mgroup, mkernel, options);
    if (bw_optim < 0) {
        algorithm.set_bw(bw);
    } else {
        algorithm.set_bw_optim(true);
        algorithm.set_bw_criterion_type(HGWR::BwOptimCriterionType(bw_optim));
    }
    algorithm.set_printer(&prcout);
    algorithm.set_canceler(&prcancel);
    HGWR::MonteCarloOptions mc_options;
    mc_options.iters = mcmc_iters;
    mc_options.burnin = mcmc_burnin;
    mc_options.beta_prior_variance = beta_prior_variance;
    mc_options.d_prior_df_offset = d_prior_df_offset;
    mc_options.d_prior_mean = d_prior_mean;
    mc_options.sigma_prior_shape = sigma_prior_shape;
    mc_options.sigma_prior_scale = sigma_prior_scale;
    mc_options.inner_tolerance = inner_tolerance;
    mc_options.save_draws = save_draws;
    mc_options.conditional_mode = conditional_mode;
    auto hgwr_result = algorithm.fit_mcmc_backfitting(f_test, mc_options);

    auto results = List::create(
        Named("gamma") = hgwr_result.gamma,
        Named("beta") = hgwr_result.beta,
        Named("mu") = hgwr_result.mu,
        Named("D") = hgwr_result.D,
        Named("sigma") = hgwr_result.sigma,
        Named("bw") = hgwr_result.bw,
        Named("gamma_se") = algorithm.get_gamma_se(),
        Named("logLik") = algorithm.get_loglik(),
        Named("trS") = algorithm.get_trS(),
        Named("var_beta") = algorithm.get_var_beta(),
        Named("edf") = algorithm.edf(),
        Named("enp") = algorithm.enp(),
        Named("fit_diagnostics") = fit_diagnostics(algorithm)
    );

    if (save_draws && !conditional_mode)
    {
        results["conditional_draws"] = List::create(
            Named("beta") = algorithm.get_last_beta_draws(),
            Named("D") = algorithm.get_last_D_draws(),
            Named("sigma2") = algorithm.get_last_sigma2_draws()
        );
    }

    if (f_test)
    {
        auto ftest_result = algorithm.test_glsw();
        results["f_test"] = ftest_result;
    }

    return results;
}
