#' Hierarchical and Geographically Weighted Regression
#'
#' A Hierarchical Linear Model (HLM) with local fixed effects.
#'
#' @param formula A formula.
#' Its structure is similar to \code{\link[lme4]{lmer}} function
#' in **lme4** package.
#' @param data A DataFrame.
#' @param local.fixed A character vector.
#' It contains names of local fixed effects.
#' @param coords A 2-column matrix.
#' It consists of coordinates for each group.
#' @param bw A numeric value. It is the value of bandwidth.
#' In this stage this function only support adaptive bandwidth.
#' And its unit must be the number of nearest neighbours.
#' @param alpha A numeric value. It is the size of the first trial step in
#' maximum likelihood algorithm.
#' @param eps_iter A numeric value. Terminate threshold of back-fitting.
#' @param eps_gradient A numeric value. Terminate threshold of
#' maximum likelihood algorithm.
#' @param max_iters An integer value. The maximum of iteration.
#' @param max_retries An integer value. If the algorithm tends to be diverge,
#' it stops automatically after trying *max_retires* times.
#' @param ml_type An integer value. Represent which maximum likelihood
#' algorithm is used. Possible values are:
#' \describe{
#'  \item{\code{HGWR_ML_TYPE_D_ONLY}}{Only \eqn{D} is specified by maximum likelihood.}
#'  \item{\code{HGWR_ML_TYPE_D_BETA}}{Both \eqn{D} and \eqn{beta} is specified by maximum likelihood.}
#' }
#' @param verbose An integer value. Determine the log level.
#' Possible values are:
#' \describe{
#'  \item{0}{no log is printed.}
#'  \item{1}{only logs in back-fitting are printed.}
#'  \item{2}{all logs are printed.}
#' }
#'
#' @return A list consists of \eqn{\gamma}, \eqn{\beta}, \eqn{D}, \eqn{\mu}
#' of a HGWR model.
#'
#' @examples
#' data(multisampling)
#' hgwr(formula = y ~ g1 + g2 + x1 + (z1 | group),
#'      data = multisampling$data,
#'      local.fixed = c("g1", "g2"),
#'      coords = multisampling$coord,
#'      bw = 10)
#'
hgwr <- function(formula, data, local.fixed, coords, bw,
                 alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
                 max_iters = 1e6, max_retries = 10,
                 ml_type = HGWR_ML_TYPE_D_ONLY, verbose = 0) {
    ### Extract variables
    model_desc <- parse.formula(formula)
    y <- as.vector(data[[model_desc$response]])
    group <- as.vector(as.integer(data[[model_desc$group]]))
    group.unique <- unique(group)
    z <- as.matrix(cbind(1, data[model_desc$random.effects]))
    fe <- model_desc$fixed.effects
    lfe <- fe[fe %in% local.fixed]
    gfe <- fe[!(fe %in% local.fixed)]
    x <- as.matrix(cbind(1, data[gfe]))
    g <- as.matrix(cbind(1, aggregate(data[lfe], list(group), mean)[,-1]))
    hgwr_result <- .hgwr_bml(g, x, z, y, as.matrix(coords), group, bw,
                             alpha, eps_iter, eps_gradient,
                             as.integer(max_iters), as.integer(max_retries),
                             as.integer(ml_type), as.integer(verbose))
    gamma <- hgwr_result$gamma
    beta <- t(hgwr_result$beta)
    mu <- hgwr_result$mu
    D <- hgwr_result$D
    intercept <- gamma[,1] + mu[,1] + beta[,1]
    coefficients <- as.data.frame(cbind(intercept, gamma[,-1], beta[,-1], mu[,-1], group.unique))
    colnames(coefficients) <- c("Intercept", lfe, gfe, model_desc$random.effects, model_desc$group)
    fitted <- rowSums(g * gamma)[group] + x %*% t(beta) + rowSums(z * mu[group,])
    list(
       coefficients = coefficients,
       random.effects.var = D,
       fitted = fitted
    )
}

#' HGWR maximum likelihood algorithm type: D only.
#'
#' @family HGWR ML types
HGWR_ML_TYPE_D_ONLY <- as.integer(0)

#' HGWR maximum likelihood algorithm type: D and beta.
#'
#' @family HGWR ML types
HGWR_ML_TYPE_D_BETA <- as.integer(1)
