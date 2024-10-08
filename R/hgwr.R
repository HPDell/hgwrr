#' Hierarchical and Geographically Weighted Regression
#'
#' A Hierarchical Linear Model (HLM) with group-level geographically weighted effects.
#'
#' @param formula A formula.
#' Its structure is similar to \code{\link[lme4]{lmer}} function
#' in **lme4** package.
#' Models can be specified with the following form:
#' ```r
#' response ~ L(local.fixed) + global.fixed + (random | group)
#' ```
#' For more information, please see the `formula` subsection in details.
#' @param data The data.
#' @param \dots Further arguments for the specified type of `data`.
#' @param bw A numeric value. It is the value of bandwidth or `"CV"`.
#' In this stage this function only support adaptive bandwidth.
#' And its unit must be the number of nearest neighbours.
#' If `"CV"` is specified, the algorithm will automatically select an
#' optimized bandwidth value.
#' @param kernel A character value. It specify which kernel function is used
#' in GWR part. Possible values are
#' \describe{
#'  \item{\code{gaussian}}{Gaussian kernel function \eqn{k(d)=\exp\left(-\frac{d^2}{b^2}\right)}}
#'  \item{\code{bisquared}}{Bi-squared kernel function. If \eqn{d<b} then \eqn{k(d)=\left(1-\frac{d^2}{b^2}\right)^2} else \eqn{k(d)=0}}
#' }
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
#'  \item{\code{D_Only}}{Only \eqn{D} is specified by maximum likelihood.}
#'  \item{\code{D_Beta}}{Both \eqn{D} and \eqn{beta} is specified by maximum likelihood.}
#' }
#' @param f_test A logical value. Determine whether to do F test on GLSW effects.
#' If `f_test=TURE`, there will be a `f_test` item in the returned object
#' showing the F test for each GLSW effect.
#' @param verbose An integer value. Determine the log level.
#' Possible values are:
#' \describe{
#'  \item{0}{no log is printed.}
#'  \item{1}{only logs in back-fitting are printed.}
#'  \item{2}{all logs are printed.}
#' }
#'
#' @return A list describing the model with following fields.
#' \describe{
#'  \item{\code{gamma}}{Coefficients of group-level spatially weighted effects.}
#'  \item{\code{beta}}{Coefficients of fixed effects.}
#'  \item{\code{mu}}{Coefficients of sample-level random effects.}
#'  \item{\code{D}}{Variance-covariance matrix of sample-level random effects.}
#'  \item{\code{sigma}}{Variance of errors.}
#'  \item{\code{effects}}{A list including names of all effects.}
#'  \item{\code{call}}{Calling of this function.}
#'  \item{\code{frame}}{The DataFrame object sent to this call.}
#'  \item{\code{frame.parsed}}{Variables extracted from the data.}
#'  \item{\code{groups}}{Unique group labels extracted from the data.}
#'  \item{\code{f_test}}{A list of F test for GLSW effects. Only exists when `f_test=TRUE`. Each item contains the F value, degrees of freedom in the numerator, degrees of freedom in the denominator, and \eqn{p} value of \eqn{F>F_\alpha}.}
#' }
#' 
#' @details  
#' ## Effect Specification in Formula
#' In the HGWR model, there are three types of effects specified by the
#' `formula` argument:
#' \describe{
#'  \item{Group-level spatially weighted (GLSW, aka. local fixed) effects}{Effects wrapped by functional symbol `L`.}
#'  \item{Sample-level random (SLR) effects}{Effects specified outside the functional symbol `L` but to the left of symbol `|`.}
#'  \item{Fixed effects}{Other effects}
#' }
#' For example, the following formula in the example of this function below is written as
#' ```r
#' y ~ L(g1 + g2) + x1 + (z1 | group)
#' ```
#' where `g1` and `g2` are GLSW effects,
#' `x1` is the fixed effects,
#' and `z1` is the SLR effects grouped by the group indicator `group`.
#' Note that SLR effects can only be specified once!
#'
#' @examples
#' data(multisampling)
#' hgwr(formula = y ~ L(g1 + g2) + x1 + (z1 | group),
#'      data = multisampling$data,
#'      coords = multisampling$coords,
#'      bw = 10)
#' 
#' mod_Ftest <-hgwr(
#'  formula = y ~ L(g1 + g2) + x1 + (z1 | group),
#'  data = multisampling$data,
#'  coords = multisampling$coords,
#'  bw = 10
#' )
#' summary(mod_Ftest)
#' 
#' 
#' @importFrom stats aggregate
#' 
#' @export 
hgwr <- function(
    formula, data, ..., bw = "CV",
    kernel = c("gaussian", "bisquared"),
    alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
    max_iters = 1e6, max_retries = 1e6,
    ml_type = c("D_Only", "D_Beta"), f_test = FALSE, verbose = 0
) {
    UseMethod("hgwr", data)
}

#' @rdname hgwr
#' @method hgwr sf
#' 
#' @importFrom sf st_centroid st_coordinates
#' @importFrom stats aggregate
#' 
#' @export 
hgwr.sf <- function(
    formula, data, ..., bw = "CV",
    kernel = c("gaussian", "bisquared"),
    alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
    max_iters = 1e6, max_retries = 1e6,
    ml_type = c("D_Only", "D_Beta"), f_test = FALSE, verbose = 0
) {
    ### Generate group-level coordinates by taking means
    data_coords <- sf::st_coordinates(sf::st_centroid(data))
    data <- sf::st_drop_geometry(data)
    group <- data[[parse.formula(formula)$group]]
    group_unique <- unique(group)
    group_index <- as.vector(match(group, group_unique))
    group_coords <- aggregate(data_coords, by = list(group_index), FUN = mean)[,-1]
    mc0 <- mc <- match.call(expand.dots = TRUE)
    mc[[1]] <- as.name("hgwr_fit")
    mc[["data"]] <- data
    mc[["coords"]] <- group_coords
    mev <- eval.parent(mc)
    mev$call <- mc0
    mev
}

#' @rdname hgwr
#' @method hgwr data.frame
#' 
#' @param coords A 2-column matrix.
#' It consists of coordinates for each group.
#' 
#' @export 
hgwr.data.frame <- function(
    formula, data, ..., coords, bw = "CV",
    kernel = c("gaussian", "bisquared"),
    alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
    max_iters = 1e6, max_retries = 1e6,
    ml_type = c("D_Only", "D_Beta"), f_test = FALSE, verbose = 0
) {
    mc0 <- mc <- match.call(expand.dots = TRUE)
    mc[[1]] <- as.name("hgwr_fit")
    mev <- eval.parent(mc)
    mev$call <- mc0
    mev
}

#' @describeIn hgwr Fit a HGWR model
#' @export 
hgwr_fit <- function(
    formula, data, coords, bw = c("CV", "AIC"),
    kernel = c("gaussian", "bisquared"),
    alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
    max_iters = 1e6, max_retries = 1e6,
    ml_type = c("D_Only", "D_Beta"), f_test = FALSE, verbose = 0
) {
    ### Extract variables
    kernel <- match.arg(kernel)
    kernel_index <- switch(
        kernel,
        "gaussian" = 0L,
        "bisquared" = 1L
    )
    ml_type <- switch(
        match.arg(ml_type),
        "D_Only" = 0L,
        "D_Beta" = 1L
    )
    model_desc <- parse.formula(formula)
    y <- as.vector(data[[model_desc$response]])
    group <- data[[model_desc$group]]
    group_unique <- unique(group)
    group_index <- as.vector(match(group, group_unique))
    re <- model_desc$random.effects
    if (length(model_desc$random.effects) > 0) {
        z <- as.matrix(make.dummy(data[model_desc$random.effects]))
        if (model_desc$intercept$random > 0) {
            re <- c("Intercept", re)
            z <- cbind(model_desc$intercept$random, z)
        }
    } else {
        if (model_desc$intercept$random > 0) {
            re <- c("Intercept", re)
            z <- matrix(rep(model_desc$intercept$random, times = nrow(data)), ncol = 1)
        } else stop("Please provide a SLR effect (including intercept) or use other models.")
    }
    gfe <- model_desc$fixed.effects
    if (length(model_desc$fixed.effects) > 0) {
        x <- as.matrix(make.dummy(data[model_desc$fixed.effects]))
        if (model_desc$intercept$fixed > 0) {
            gfe <- c("Intercept", gfe)
            x <- cbind(model_desc$intercept$fixed, x)
        }
    } else {
        if (model_desc$intercept$fixed > 0) {
            gfe <- c("Intercept", gfe)
            x <- matrix(rep(model_desc$intercept$fixed, times = nrow(data)), ncol = 1)
        }
        else stop("Please provide a fixed effect (including intercept) or use other models.")
    }
    lfe <- model_desc$local.fixed.effects
    if (length(model_desc$local.fixed.effects) > 0) {
        g <- as.matrix(aggregate(make.dummy(data[model_desc$local.fixed.effects]), list(group), mean)[,-1])
        if (model_desc$intercept$local > 0) {
            lfe <- c("Intercept", lfe)
            g <- cbind(model_desc$intercept$local, g)
        }
    } else {
        if (model_desc$intercept$local > 0) {
            lfe <- c("Intercept", lfe)
            g <- matrix(rep(model_desc$intercept$local, times = length(group_unique)), ncol = 1)
        }
        else stop("Please provide a GLSW effect (including intercept) or use other models.")
    }
    ### Get bandwidth value
    if (is.character(bw)) {
        bw <- match.arg(bw)
        bw_value <- NA_real_
        optim_bw <- TRUE
        bw_criterion <- switch(bw, "CV" = 0L, "AIC" = 1L)
    } else if (is.numeric(bw) || is.integer(bw)) {
        bw_value <- bw
        optim_bw <- FALSE
        bw_criterion <- -1L
    } else {
        bw_value <- Inf
        optim_bw <- TRUE
        bw_criterion <- 0L
    }

    ### Call C
    hgwr_result <- tryCatch({
        hgwr_bfml(
            g, x, z, y, as.matrix(coords), group_index, bw_value, bw_criterion, kernel_index,
            alpha, eps_iter, eps_gradient,
            as.integer(max_iters), as.integer(max_retries),
            as.integer(ml_type), f_test, as.integer(verbose)
        )
    }, error = function(e) {
        stop("Error occurred when estimating HGWR parameters.")
    })
    rownames(hgwr_result$beta) <- c(gfe)
    colnames(hgwr_result$gamma) <- c(lfe)
    colnames(hgwr_result$mu) <- c(re)
    if (optim_bw)
        bw_value <- hgwr_result$bw

    ### Prepare Return Result
    result <- c(hgwr_result, list(
        effects = list(
            global.fixed = gfe,
            local.fixed = lfe,
            random = re,
            group = model_desc$group,
            response = model_desc$response
        ),
        intercept = model_desc$intercept,
        frame = data,
        frame.parsed = list(
            y = y,
            x = x,
            g = g,
            z = z,
            group = group_index
        ),
        groups = group_unique,
        coords = coords
    ))
    result$bw <- bw_value
    class(result) <- "hgwrm"
    result
}

#' Get estimated coefficients.
#'
#' @param object An `hgwrm` object returned by [hgwr()].
#' @param \dots Parameter received from other functions.
#'
#' @return A \code{DataFrame} object consists of all estimated coefficients.
#'
#' @seealso [hgwr()], [summary.hgwrm()], [fitted.hgwrm()] and [residuals.hgwrm()].
#' 
#' @export
coef.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }
    gamma <- object$gamma
    beta <- matrix(object$beta, nrow = length(object$groups), ncol = length(object$beta), byrow = T)
    mu <- object$mu
    intercept <- matrix(0, length(object$groups), 1)
    if (object$intercept$local > 0) {
        intercept <- intercept + gamma[,1]
        gamma <- gamma[,-1]
    }
    if (object$intercept$fixed > 0) {
        intercept <- intercept + beta[,1]
        beta <- beta[,-1]
    }
    if (object$intercept$random > 0) {
        intercept <- intercept + mu[,1]
        mu <- mu[,-1]
    }
    effects <- object$effects
    coef <- as.data.frame(cbind(intercept, gamma, beta, mu))
    coef_names <- c(
        effects$local.fixed[effects$local.fixed != "Intercept"],
        effects$global.fixed[effects$global.fixed != "Intercept"],
        effects$random[effects$random != "Intercept"]
    )
    if (any(unlist(object$intercept) > 0)) {
        colnames(coef) <- c("Intercept", coef_names)
    } else {
        colnames(coef) <- coef_names
    }
    coef
}

#' Get fitted response.
#'
#' @inheritParams coef.hgwrm
#'
#' @return A vector consists of fitted response values.
#'
#' @seealso [hgwr()], [summary.hgwrm()], [coef.hgwrm()] and [residuals.hgwrm()].
#' 
#' @export 
fitted.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }
    xf <- object$frame.parsed
    rowSums(xf$g * object$gamma)[xf$group] +
        as.vector(xf$x %*% object$beta) +
        rowSums(xf$z * object$mu[xf$group,])
}

#' Get residuals.
#'
#' @inheritParams coef.hgwrm
#'
#' @return A vector consists of residuals.
#'
#' @seealso [hgwr()], [summary.hgwrm()], [coef.hgwrm()] and [fitted.hgwrm()].
#' 
#' @export 
#' 
residuals.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }
    object$frame.parsed$y - fitted.hgwrm(object)
}

#' Summary an `hgwrm` object.
#'
#' @param object An `hgwrm` object returned from [hgwr()].
#' @param \dots Other arguments passed from other functions.
#' @param test_hetero Logical/list value.
#' Whether to test the spatial heterogeneity of GLSW effects.
#' If it is set to `FALSE`, the test will not be executed.
#' If it is set to `TRUE`, the test will be executed with default parameters (see details below).
#' It accepts a list to enable the test with specified parameters.
#' @param verbose An Integer value to control whether additional messages
#' during testing spatial heterogeneity should be reported.
#'
#' @return A list containing summary informations of this `hgwrm` object
#' with the following fields.
#' \describe{
#'  \item{\code{diagnostic}}{A list of diagnostic information.}
#'  \item{\code{random.stddev}}{The standard deviation of random effects.}
#'  \item{\code{random.corr}}{The correlation matrix of random effects.}
#'  \item{\code{residuals}}{The residual vector.}
#' }
#' 
#' @details The parameters used to perform test of spatial heterogeneity are
#' \describe{
#'  \item{\code{bw}}{Bandwidth (unit: number of nearest neighbours) used to make spatial kernel density estimation. Default: `10`.}
#'  \item{\code{poly}}{The number of polynomial terms used in the local polynomial estimation. Default: `2`.}
#'  \item{\code{resample}}{Total resampling times. Default: `5000`.}
#'  \item{\code{kernel}}{The kernel function used in the local polynomial estimation. Options are `"gaussian"` and `"bisquared"`. Default: `"bisquared"`.}
#' }
#' 
#' @examples 
#' data(multisampling)
#' m <- hgwr(
#'  formula = y ~ L(g1 + g2) + x1 + (z1 | group),
#'  data = multisampling$data,
#'  coords = multisampling$coords,
#'  bw = 10
#' )
#' summary(m)
#' summary(m, test_hetero = TRUE)
#' summary(m, test_hetero = list(kernel = "gaussian"))
#' 
#' @importFrom stats AIC logLik pt sd
#' @seealso [hgwr()].
#' @export 
summary.hgwrm <- function(object, ..., test_hetero = FALSE, verbose = 0) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }

    res <- as.list(object)

    ### Diagnostics
    #### R-squared
    y <- object$frame[[object$effects$response]]
    tss <- sum((y - mean(y))^2)
    x_residuals <- residuals.hgwrm(object)
    rss <- sum(x_residuals^2)
    rsquared <- 1 - rss / tss
    #### AIC
    logLik_object <- stats::logLik(object)
    aic <- stats::AIC(logLik_object)
    #### Record results
    res$diagnostic <- list(
        rsquared = rsquared,
        logLik = object$logLik,
        AIC = aic
    )

    ### Significance test
    significance <- list()
    #### Beta
    edf <- object$edf
    se_beta <- sqrt(object$var_beta)
    t_beta <- abs(object$beta) / se_beta
    p_beta <- (1 - stats::pt(t_beta, edf)) * 2
    significance$beta <- data.frame(
        est = object$beta,
        sd = se_beta,
        tv = t_beta,
        pv = p_beta
    )
    #### Gamma different from 0
    t_gamma <- object$gamma / object$gamma_se
    p_gamma <- stats::pt(-abs(object$gamma / object$gamma_se), object$edf) * 2
    significance$gamma <- list(
        tv = t_gamma,
        pv = p_gamma
    )
    #### Gamma spatial heterogeneity
    is_test_hetero <- FALSE
    if (length(test_hetero) == 1 && is.logical(test_hetero) && test_hetero == TRUE) {
        is_test_hetero <- TRUE
    } else if (is.list(test_hetero)) {
        is_test_hetero <- TRUE
    }
    if (is_test_hetero) {
        bw <- 10L
        resample <- 5000L
        poly <- 2L
        kernel <- "bisquared"
        if (is.list(test_hetero)) {
            bw <- ifelse("bw" %in% names(test_hetero), test_hetero$bw, bw)
            resample <- ifelse("resample" %in% names(test_hetero), test_hetero$resample, resample)
            poly <- ifelse("poly" %in% names(test_hetero), test_hetero$poly, poly)
            kernel <- ifelse("kernel" %in% names(test_hetero), test_hetero$kernel, kernel)
        }
        kernel <- match.arg(kernel, c("gaussian", "bisquared"))
        kernel_id <- switch(kernel, "gaussian" = 0, "bisquared" = 1)
        sd_gamma <- apply(object$gamma, 2, stats::sd)
        t_gamma <- spatial_hetero_perm(
            object$gamma,
            as.matrix(object$coords),
            poly = poly,
            resample = resample,
            bw = bw,
            kernel = kernel_id,
            verbose = as.integer(verbose)
        )
        t_gamma$t <- sqrt(t_gamma$t)
        t_gamma$t0 <- sqrt(t_gamma$t0)
        pv <- sapply(seq_along(t_gamma$t0), function(i) {
            with(t_gamma, mean(t[,i] > t0[i]))
        })
        significance$gamma$hetero <- list(
            sd = sd_gamma,
            t = t_gamma$t,
            t0 = t_gamma$t0,
            pv = pv
        )
    }
    #### Save results
    res$significance <- significance

    ### Random effects
    random_corr_cov <- object$sigma * object$sigma * object$D
    random_stddev <- sqrt(diag(random_corr_cov))
    random_corr <- t(random_corr_cov / random_stddev) / random_stddev
    diag(random_corr) <- 1
    res$random.stddev <- random_stddev
    res$random.corr <- random_corr

    ### Residuals
    res$residuals <- x_residuals

    ### return
    class(res) <- "summary.hgwrm"
    res
}

#' Print description of a `hgwrm` object.
#'
#' @param x An `hgwrm` object returned by [hgwr()].
#' @param decimal.fmt The format string passing to [base::sprintf()].
#' @inheritDotParams print.table.md
#'
#' @return No return.
#'
#' @examples
#' data(multisampling)
#' model <- hgwr(formula = y ~ L(g1 + g2) + x1 + (z1 | group),
#'               data = multisampling$data,
#'               coords = multisampling$coords,
#'               bw = 10)
#' print(model)
#' print(model, table.style = "md")
#'
#' @seealso [summary.hgwrm()], [print.table.md()].
#' 
#' @importFrom stats fivenum
#' 
#' @export 
#'
print.hgwrm <- function(x, decimal.fmt = "%.6f", ...) {
    if (!inherits(x, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }

    ### Basic Information
    cat("Hierarchical and geographically weighted regression model", fill = T)
    cat("=========================================================", fill = T)
    cat("Formula:", deparse(x$call[[2]]), fill = T)
    cat(" Method:", "Back-fitting and Maximum likelihood", fill = T)
    cat("   Data:", deparse(x$call[[3]]), fill = T)
    cat("\n")
    effects <- x$effects
    cat("Fixed Effects", fill = T)
    cat("-------------", fill = T)
    beta_str <- rbind(
        effects$global.fixed,
        matrix2char(t(x$beta), decimal.fmt)
    )
    print.table.md(beta_str, ...)
    cat("\n")
    cat("Group-level Spatially Weighted Effects", fill = T)
    cat("--------------------------------------", fill = T)
    cat("Bandwidth:", x$bw, "(nearest neighbours)", fill = T)
    cat("\n")
    cat("Coefficient estimates:", fill = T)
    gamma_fivenum <- t(apply(x$gamma, 2, fivenum))
    gamma_str <- rbind(
        c("Coefficient", "Min", "1st Quartile", "Median", "3rd Quartile", "Max"),
        cbind(effects$local.fixed, matrix2char(gamma_fivenum, decimal.fmt))
    )
    print.table.md(gamma_str, ...)
    cat("\n")
    cat("Sample-level Random Effects", fill = T)
    cat("---------------------------", fill = T)
    re <- x$effects$random
    x <- summary.hgwrm(x)
    random_stddev <- x$random.stddev
    random_dev_str <- cbind(
        "", re, matrix2char(matrix(random_stddev, ncol = 1), decimal.fmt)
    )
    random_dev_str[1, 1] <- effects$group
    random_dev_str <- rbind(
        c("Groups", "Name", "Std.Dev."),
        random_dev_str
    )
    random_residual_str <- matrix(c("Residual", "", sprintf(decimal.fmt, x$sigma)), nrow = 1)
    random_str <- rbind(
        random_dev_str,
        random_residual_str
    )
    if (length(re) > 1) {
        random_corr <- x$random.corr
        random_corr_str <- matrix2char(random_corr, decimal.fmt)
        random_corr_str[!lower.tri(random_corr)] <- ""
        random_corr_str <- rbind("", random_corr_str)
        random_corr_str[1, 1] <- "Corr"
        random_str <- cbind(
            random_str, rbind(
                random_corr_str,
                matrix("", nrow = nrow(random_residual_str), ncol = ncol(random_corr_str))
            )
        )
        print.table.md(random_str, ...)
    }
    print.table.md(random_str, ...)
    cat("\n")
    cat("Other Information", fill = T)
    cat("-----------------", fill = T)
    cat("Number of Obs:", nrow(x$frame), fill = T)
    cat("       Groups:", effects$group, ",", nrow(x$mu), fill = T)
}

#' Print summary of an `hgwrm` object.
#'
#' @param x An object returned from [summary.hgwrm()].
#' @inherit print.hgwrm
#' @inheritDotParams print.table.md
#'
#' @return No return.
#'
#' @examples
#' data(multisampling)
#' model <- hgwr(formula = y ~ L(g1 + g2) + x1 + (z1 | group),
#'               data = multisampling$data,
#'               coords = multisampling$coords,
#'               bw = 10)
#' summary(model)
#' 
#' @export 
#'
print.summary.hgwrm <- function(x, decimal.fmt = "%.6f", ...) {
    if (!inherits(x, "summary.hgwrm")) {
        stop("It's not a summary.hgwrm object.")
    }
    args = list(...)
    is_md = ifelse(!is.null(args$table.style), args$table.style == "md", FALSE)

    ### Call information
    cat("Hierarchical and geographically weighted regression model", fill = T)
    cat("=========================================================", fill = T)
    cat("Formula:", deparse(x$call[[2]]), fill = T)
    cat(" Method:", "Back-fitting and Maximum likelihood", fill = T)
    cat("   Data:", deparse(x$call[[3]]), fill = T)
    cat("\n")

    ### Parameter Estimates
    cat("Parameter Estimates", fill = T)
    cat("-------------------", fill = T)
    cat("Fixed effects:", fill = T)
    gfe <- x$effects$global.fixed
    pv_gfe <- x$significance$beta$pv
    stars <- vapply(pv_gfe, pv2stars, rep(" ", n = length(pv_gfe)))
    print.table.md(rbind(
        c("", "Estimated", "Sd. Err", "t.val", "Pr(>|t|)", ""),
        as.matrix(cbind(variable = gfe, matrix2char(as.matrix(x$significance$beta), fmt = decimal.fmt), stars = stars))
    ), ...)
    cat("\n")
    cat("Bandwidth:", x$bw, "(nearest neighbours)", fill = T)
    cat("\n")
    cat("GLSW effects:", fill = T)
    lfe <- x$effects$local.fixed
    s_gamma <- t(apply(x$significance$gamma$pv, 2, function(p) { c(
        "***" = mean(p < 0.001),
        "**" = mean(p >= 0.001 & p < 0.01),
        "*" = mean(p >= 0.01 & p < 0.05),
        "." = mean(p >= 0.05 & p < 0.1)
    ) })) * 100
    rownames(s_gamma) <- lfe
    gamma_stats <- cbind(
        gamma_mean <- matrix2char(as.vector(colMeans(x$gamma)), decimal.fmt),
        gamma_se_mean <- matrix2char(as.vector(colMeans(x$gamma_se)), decimal.fmt),
        matrix2char(s_gamma, fmt = "%.1f%%")
    )
    gamma_stats_name <- c("Mean Est.", "Mean Sd.", colnames(s_gamma))
    gamma_str <- rbind(
        c("", gamma_stats_name),
        cbind(lfe, gamma_stats)
    )
    print.table.md(gamma_str, ...)
    cat("\n")
    if (!is.null(x$f_test)) {
        cat("GLSW effect F test:", fill = T)
        f_test_res <- do.call(rbind, x$f_test)
        f_test_stars <-  vapply(f_test_res[,4], FUN = pv2stars, rep(" ", n = nrow(f_test_res)))
        f_test_str <- rbind(
            c("", "F value", "Num. D.F.", "Den. D.F.", "Pr(>F)", ""),
            cbind(lfe, matrix2char(f_test_res, decimal.fmt), f_test_stars)
        )
        print.table.md(f_test_str, ...)
        cat("\n")
    }
    if (!is.null(x$significance$gamma$hetero)) {
        cat("GLSW effect spatial heterogeneity:", fill = T)
        h_gamma <- x$significance$gamma$hetero
        pv_hetero <- h_gamma$pv
        hetero_stats <- matrix2char(cbind(
            gamma_sd = as.vector(h_gamma$sd),
            t0 = as.vector(h_gamma$t0),
            t.min = apply(h_gamma$t, 2, min),
            t.max = apply(h_gamma$t, 2, max),
            pv = pv_hetero
        ), decimal.fmt)
        hetero_stats <- cbind(
            hetero_stats,
            stars = vapply(pv_hetero, pv2stars, rep(" ", n = length(pv_hetero)))
        )
        hetero_stats_name <- c("Sd. Est.", "t0", "Min. t", "Max. t", "Pr(t>t0)", "")
        hetero_str <- rbind(
            c("", hetero_stats_name),
            cbind(lfe, hetero_stats)
        )
        print.table.md(hetero_str, ...)
        cat("\n")
    }
    cat("SLR effects:", fill = T)
    re <- x$effects$random
    random_stddev <- x$random.stddev
    random_dev_str <- cbind(
        "", re, matrix2char(cbind(colMeans(x$mu), matrix(random_stddev, ncol = 1)), decimal.fmt)
    )
    random_dev_str[1, 1] <- x$effects$group
    random_dev_str <- rbind(
        c("Groups", "Name", "Mean", "Std.Dev."),
        random_dev_str
    )
    random_residual_str <- matrix(c("Residual", "", sprintf(decimal.fmt, c(mean(x$residuals), x$sigma))), nrow = 1)
        random_str <- rbind(
            random_dev_str,
            random_residual_str
        )
    if (length(re) > 1) {
        random_corr <- x$random.corr
        random_corr_str <- matrix2char(random_corr, decimal.fmt)
        random_corr_str[!lower.tri(random_corr)] <- ""
        random_corr_str <- rbind("", random_corr_str)
        random_corr_str[1, 1] <- "Corr"
        random_str <- cbind(
            random_str, rbind(
                random_corr_str,
                matrix("", nrow = nrow(random_residual_str), ncol = ncol(random_corr_str))
            )
        )
        print.table.md(random_str, ...)
    }
    print.table.md(random_str, ...)
    cat("\n", fill = T)

    ### Diagnostics
    cat("Diagnostics", fill = T)
    cat("-----------", fill = T)
    diagnostic_chr <- cbind(
        names(x$diagnostic),
        matrix2char(unlist(x$diagnostic), decimal.fmt)
    )
    print.table.md(diagnostic_chr, ...)
    cat("\n")

    ### Residuals
    cat("Scaled Residuals", fill = T)
    cat("----------------", fill = T)
    resiudal_fivenum <- fivenum(x$residuals)
    residual_fivenum_mat <- matrix(resiudal_fivenum, nrow = 1)
    residual_fivenum_chr <- rbind(
        c("Min", "1Q", "Median", "3Q", "Max"),
        matrix2char(residual_fivenum_mat, decimal.fmt)
    )
    print.table.md(residual_fivenum_chr, ...)
    cat("\n")
    cat("Other Information", fill = T)
    cat("-----------------", fill = T)
    cat("Number of Obs:", nrow(x$frame), fill = T)
    cat("       Groups:", x$effects$group, ",", nrow(x$mu), fill = T)
}

#' Log likelihood function
#' 
#' @param object An `hgwrm` object.
#' @param \dots Additional arguments.
#' 
#' @return An `logLik` instance used for S3 method `logLik()`.
#' 
#' @method logLik hgwrm
#' @export 
logLik.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not an hgwrm object.")
    }
    n <- object$frame.parsed$y
    p <- length(object$beta)
    q <- ncol(object$mu)
    enp_gwr <- object$trS[1]
    enp_hlm <- p + q * (q + 1) / 2 + 1
    enp <- enp_gwr + enp_hlm
    val <- object$logLik
    attr(val, "df") <- object$enp
    attr(val, "nall") <- n
    attr(val, "nobs") <- n
    class(val) <- "logLik"
    val
}
