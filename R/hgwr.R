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
    sigma <- hgwr_result$sigma
    fitted <- rowSums(g * gamma)[group] + x %*% t(beta) + rowSums(z * mu[group,])
    result <- list(
        gamma = gamma,
        beta = beta,
        mu = mu,
        D = D,
        fitted = fitted,
        sigma = sigma,
        model.effects = list(
            global.fixed = gfe,
            local.fixed = lfe,
            random = model_desc$random.effects,
            group = model_desc$group
        ),
        call = match.call(),
        frame = data,
        groups = group.unique
    )
    class(result) <- "hgwrm"
    result
}

#' HGWR maximum likelihood algorithm type: D only.
#'
#' @family HGWR ML types
HGWR_ML_TYPE_D_ONLY <- as.integer(0)

#' HGWR maximum likelihood algorithm type: D and beta.
#'
#' @family HGWR ML types
HGWR_ML_TYPE_D_BETA <- as.integer(1)


print.table.md <- function(x) {
    x.length <- apply(x, c(1, 2), nchar)
    x.length.max <- apply(x.length, 2, max)
    x.fmt <- sprintf("%%%ds", x.length.max)
    for(c in 1:ncol(x)) {
        if(x.length.max[c] > 0)
            cat("|", sprintf(x.fmt[c], x[1, c]), " ")
    }
    cat("|\n")
    for(c in 1:ncol(x)) {
        if(x.length.max[c] > 0)
            cat("|", sprintf(paste(rep("-", x.length.max[c]), collapse = "")), " ")
    }
    cat("|\n")
    for (r in 2:nrow(x)) {
        for (c in 1:ncol(x)) {
            if(x.length.max[c] > 0)
                cat("|", sprintf(x.fmt[c], x[r, c]), " ")
        }
        cat("|\n")
    }
}

print.hgwrm <- function(x, decimal.fmt = "%.6f", ...) {
    matrix2char <- function(m, fmt = decimal.fmt) {
        apply(m, c(1, 2), function(x) { sprintf(fmt, x) })
    }
    if (class(x) != "hgwrm") {
        stop("It's not a hgwm object.")
    }
    cat("Hierarchical and geographically weighted regression model", "\n")
    cat("=========================================================", "\n")
    cat("Formula:", deparse(x$call[[2]]), "\n")
    cat(" Method:", "Back-fitting and Maximum likelihood", "\n")
    cat("   Data:", deparse(x$call[[3]]), "\n")
    cat("\n")
    effects <- x$model.effects
    cat("Global Fixed Effects", "\n")
    cat("-------------------", "\n")
    beta_str <- rbind(
        c("Intercept", effects$global.fixed),
        matrix2char(x$beta)
    )
    print.table.md(beta_str)
    cat("\n")
    cat("Local Fixed Effects", "\n")
    cat("-------------------", "\n")
    gamma_fivenum <- t(apply(x$gamma, 2, fivenum))
    gamma_str <- rbind(
        c("Coefficient", "Min", "1st Quartile", "Median", "3rd Quartile", "Max"),
        cbind(c("Intercept", effects$local.fixed), matrix2char(gamma_fivenum))
    )
    print.table.md(gamma_str)
    cat("\n")
    cat("Random Effects", "\n")
    cat("--------------", "\n")
    random_effects <- effects$random
    random_corr_cov <- x$sigma * x$sigma * x$D
    random_stddev <- sqrt(diag(random_corr_cov))
    random_corr <- t(random_corr_cov / random_stddev) / random_stddev
    diag(random_corr) <- 1
    random_corr_str <- matrix2char(random_corr)
    random_corr_str[!lower.tri(random_corr)] <- ""
    random_corr_str <- rbind("", random_corr_str)
    random_corr_str[1, 1] <- "Corr"
    random_dev_str <- cbind(
        "", c("Intercept", random_effects), matrix2char(matrix(random_stddev, ncol = 1))
    )
    random_dev_str[1, 1] <- effects$group
    random_dev_str <- rbind(
        c("Groups", "Name", "Std.Dev."),
        random_dev_str
    )
    random_residual_str <- cbind(
        matrix(c("Residual", "", sprintf(decimal.fmt, x$sigma)), nrow = 1),
        matrix("", nrow = 1, ncol = ncol(random_corr))
    )
    random_str <- rbind(
        cbind(random_dev_str, random_corr_str),
        random_residual_str
    )
    print.table.md(random_str)
    cat("\n")
    cat("Other Information", "\n")
    cat("-----------------", "\n")
    cat("Number of Obs:", nrow(x$frame), "\n")
    cat("       Groups:", effects$group, ",", nrow(x$mu), "\n")
}

coef.hgwrm <- function(x, ...) {
    if (class(x) != "hgwrm") {
        stop("It's not a hgwm object.")
    }
    gamma <- x$gamma
    beta <- matrix(x$beta, nrow = length(x$groups), ncol = ncol(x$beta), byrow = T)
    mu <- x$mu
    intercept <- gamma[,1] + mu[,1] + beta[,1]
    effects <- x$model.effects
    coef <- as.data.frame(cbind(intercept, gamma[,-1], beta[,-1], mu[,-1], x$groups))
    colnames(coef) <- c("Intercept", effects$global.fixed, effects$local.fixed, effects$random, effects$group)
    coef
}
