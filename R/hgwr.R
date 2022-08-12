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
#'  \item{\code{gamma}}{Coefficients of local fixed effects.}
#'  \item{\code{beta}}{Coefficients of global fixed effects.}
#'  \item{\code{mu}}{Coefficients of random effects.}
#'  \item{\code{D}}{Variance-covariance matrix of random effects.}
#'  \item{\code{sigma}}{Variance of errors.}
#'  \item{\code{effects}}{A list including names of all effects.}
#'  \item{\code{call}}{Calling of this function.}
#'  \item{\code{frame}}{The DataFrame object sent to this call.}
#'  \item{\code{frame.parsed}}{Variables extracted from the data.}
#'  \item{\code{groups}}{Unique group labels extracted from the data.}
#' }
#'
#' @examples
#' data(multisampling)
#' hgwr(formula = y ~ g1 + g2 + x1 + (z1 | group),
#'      data = multisampling$data,
#'      local.fixed = c("g1", "g2"),
#'      coords = multisampling$coords,
#'      bw = 10)
#'
hgwr <- function(
    formula, data, local.fixed, coords, bw,
    kernel = c("gaussian", "bisquared"),
    alpha = 0.01, eps_iter = 1e-6, eps_gradient = 1e-6,
    max_iters = 1e6, max_retries = 10,
    ml_type = c("D_Only", "D_Beta"), verbose = 0
) {
    ### Extract variables
    kernel <- match.arg(kernel)
    ml_type <- switch(
        match.arg(ml_type),
        "D_Only" = 0L,
        "D_Beta" = 1L
    )
    model_desc <- parse.formula(formula)
    y <- as.vector(data[[model_desc$response]])
    group <- as.vector(as.integer(data[[model_desc$group]]))
    group_unique <- unique(group)
    group_index <- match(group, group_unique)
    z <- as.matrix(cbind(1, data[model_desc$random.effects]))
    fe <- model_desc$fixed.effects
    lfe <- fe[fe %in% local.fixed]
    gfe <- fe[!(fe %in% local.fixed)]
    x <- as.matrix(cbind(1, data[gfe]))
    g <- as.matrix(cbind(1, aggregate(data[lfe], list(group), mean)[,-1]))
    hgwr_result <- .hgwr_bml(
        g, x, z, y, as.matrix(coords), group_index, bw, kernel,
        alpha, eps_iter, eps_gradient,
        as.integer(max_iters), as.integer(max_retries),
        as.integer(ml_type), as.integer(verbose)
    )
    gamma <- hgwr_result$gamma
    beta <- t(hgwr_result$beta)
    mu <- hgwr_result$mu
    D <- hgwr_result$D
    sigma <- hgwr_result$sigma
    result <- list(
        gamma = gamma,
        beta = beta,
        mu = mu,
        D = D,
        sigma = sigma,
        effects = list(
            global.fixed = gfe,
            local.fixed = lfe,
            random = model_desc$random.effects,
            group = model_desc$group,
            response = model_desc$response
        ),
        call = match.call(),
        frame = data,
        frame.parsed = list(
            y = y,
            x = x,
            g = g,
            z = z,
            group = group
        ),
        groups = group_unique
    )
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
coef.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }
    gamma <- object$gamma
    beta <- matrix(object$beta, nrow = length(object$groups), ncol = ncol(object$beta), byrow = T)
    mu <- object$mu
    intercept <- gamma[,1] + mu[,1] + beta[,1]
    effects <- object$effects
    coef <- as.data.frame(cbind(intercept, gamma[,-1], beta[,-1], mu[,-1], object$groups))
    colnames(coef) <- c("Intercept", effects$local.fixed, effects$global.fixed, effects$random, effects$group)
    coef
}

#' Get fitted reponse.
#'
#' @inheritParams coef.hgwrm
#'
#' @return A vector consists of fitted response values.
#'
#' @seealso [hgwr()], [summary.hgwrm()], [coef.hgwrm()] and [residuals.hgwrm()].
fitted.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }
    xf <- object$frame.parsed
    rowSums(xf$g * object$gamma)[xf$group] +
        as.vector(xf$x %*% t(object$beta)) +
        rowSums(xf$z * object$mu[xf$group,])
}

#' Get residuals.
#'
#' @inheritParams coef.hgwrm
#'
#' @return A vector consists of residuals.
#'
#' @seealso [hgwr()], [summary.hgwrm()], [coef.hgwrm()] and [fitted.hgwrm()].
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
#' @seealso [hgwr()].
#'
summary.hgwrm <- function(object, ...) {
    if (!inherits(object, "hgwrm")) {
        stop("It's not a hgwrm object.")
    }

    res <- as.list(object)

    ### Diagnostics
    y <- object$frame[[object$effects$response]]
    tss <- sum((y - mean(y))^2)
    x_residuals <- residuals.hgwrm(object)
    rss <- sum(x_residuals^2)
    rsquared <- 1 - rss / tss
    res$diagnostic <- list(
        rsquared = rsquared
    )

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

#' Print a character matrix as a table.
#'
#' @param x A character matrix.
#' @param col.sep Column seperator. Default to `""`.
#' @param header.sep Header seperator. Default to `"-"`.
#' @param row.begin Character at the beginning of each row.
#' Default to `col.sep`.
#' @param row.end Character at the ending of each row.
#' Default to `col.sep`.
#' @param table.style Name of pre-defined style.
#' Possible values are `"plain"`, `"md"` or `"latex"`. Default to `"plain"`.
#' @param \dots Additional style control arguments.
#'
#' @return No return.
#'
#' @details
#' When `table.style` is specified, `col.sep`, `header.sep`, `row.begin`
#' and `row.end` would not take effects.
#' Because this function will automatically set their values.
#' For each possible value of `table.style`, its corresponding style settings
#' are shown in the following table.
#' \tabular{llll}{
#'                   \tab \strong{\code{plain}} \tab \strong{\code{md}} \tab \strong{\code{latex}} \cr
#' \code{col.sep}    \tab \code{""}             \tab \code{"|"}         \tab \code{"&"}            \cr
#' \code{header.sep} \tab \code{""}             \tab \code{"-"}         \tab \code{""}             \cr
#' \code{row.begin}  \tab \code{""}             \tab \code{"|"}         \tab \code{""}             \cr
#' \code{row.end}    \tab \code{""}             \tab \code{"|"}         \tab \code{"\\\\"}
#' }
#'
#' In this function, characters are right padded by spaces.
#'
#' @seealso [print.hgwrm()], [summary.hgwrm()].
print.table.md <- function(x, col.sep = "", header.sep = "",
                           row.begin = "", row.end = "",
                           table.style = c("plain", "md", "latex"), ...) {
    if (!missing(table.style)) {
        table.style <- match.arg(table.style)
        if (table.style == "md") {
            col.sep <- "|"
            header.sep <- "-"
            row.begin <- "|"
            row.end <- "|"
        } else if (table.style == "latex") {
            col.sep <- "&"
            header.sep <- ""
            row.begin <- ""
            row.end <- "\\\\"
        } else if (table.style == "plain") {
            col.sep <- ""
            header.sep <- ""
            row.begin <- ""
            row.end <- ""
        } else {
           stop("Unknown table.style.")
        }
    }
    if (nchar(header.sep) > 1) {
       stop("Currently only 1 character header.sep is supported.")
    }
    ### Print table
    x.length <- apply(x, c(1, 2), nchar)
    x.length.max <- apply(x.length, 2, max)
    x.fmt <- sprintf("%%%ds", x.length.max)
    for(c in 1:ncol(x)) {
        if(x.length.max[c] > 0)
            cat(ifelse(c == 1, row.begin, col.sep),
                sprintf(x.fmt[c], x[1, c]), "")
    }
    cat(paste0(row.end, "\n"))
    if (nchar(header.sep) > 0) {
        for(c in 1:ncol(x)) {
            if(x.length.max[c] > 0) {
                header.sep.full <- paste(rep("-", x.length.max[c]),
                                         collapse = "")
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(header.sep.full), "")
            }
        }
        cat(paste0(row.end, "\n"))
    }
    for (r in 2:nrow(x)) {
        for (c in 1:ncol(x)) {
            if(x.length.max[c] > 0)
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(x.fmt[c], x[r, c]), "")
        }
        cat(paste0(row.end, "\n"))
    }
}

#' Convert a numeric matrix to character matrix according to a format string.
#'
#' @param m A numeric matrix.
#' @param fmt Format string. Passing to [base::sprintf()].
#'
#' @seealso [base::sprintf()], [print.hgwrm()], [print.summary.hgwrm()].
matrix2char <- function(m, fmt = "%.6f") {
    mc <- NULL
    if ("array" %in% class(m)) {
        mc <- apply(m, seq(length(dim(m))), function(x) { sprintf(fmt, x) })
    } else {
        mc <- sprintf(fmt, m)
    }
    mc
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
#' model <- hgwr(formula = y ~ g1 + g2 + x1 + (z1 | group),
#'               data = multisampling$data,
#'               local.fixed = c("g1", "g2"),
#'               coords = multisampling$coords,
#'               bw = 10)
#' print(model)
#' print(model, table.style = "md")
#'
#' @seealso [summary.hgwrm()], [print.table.md()].
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
    cat("Global Fixed Effects", fill = T)
    cat("-------------------", fill = T)
    beta_str <- rbind(
        c("Intercept", effects$global.fixed),
        matrix2char(x$beta)
    )
    print.table.md(beta_str, ...)
    cat("\n")
    cat("Local Fixed Effects", fill = T)
    cat("-------------------", fill = T)
    gamma_fivenum <- t(apply(x$gamma, 2, fivenum))
    gamma_str <- rbind(
        c("Coefficient", "Min", "1st Quartile", "Median", "3rd Quartile", "Max"),
        cbind(c("Intercept", effects$local.fixed), matrix2char(gamma_fivenum))
    )
    print.table.md(gamma_str, ...)
    cat("\n")
    cat("Random Effects", fill = T)
    cat("--------------", fill = T)
    x_summary <- summary.hgwrm(x)
    random_stddev <- x_summary$random.stddev
    random_corr <- x_summary$random.corr
    random_corr_str <- matrix2char(random_corr)
    random_corr_str[!lower.tri(random_corr)] <- ""
    random_corr_str <- rbind("", random_corr_str)
    random_corr_str[1, 1] <- "Corr"
    random_dev_str <- cbind(
        "", c("Intercept", x$effects$random), matrix2char(matrix(random_stddev, ncol = 1))
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
#' model <- hgwr(formula = y ~ g1 + g2 + x1 + (z1 | group),
#'               data = multisampling$data,
#'               local.fixed = c("g1", "g2"),
#'               coords = multisampling$coords,
#'               bw = 10)
#' summary(model)
#'
print.summary.hgwrm <- function(x, decimal.fmt = "%.6f", ...) {
    if (!inherits(x, "summary.hgwrm")) {
        stop("It's not a summary.hgwrm object.")
    }

    ### Call information
    cat("Hierarchical and geographically weighted regression model", fill = T)
    cat("=========================================================", fill = T)
    cat("Formula:", deparse(x$call[[2]]), fill = T)
    cat(" Method:", "Back-fitting and Maximum likelihood", fill = T)
    cat("   Data:", deparse(x$call[[3]]), fill = T)
    cat("\n")

    ### Diagnostics
    cat("Diagnostics", fill = T)
    cat("-----------", fill = T)
    rsquared <- x$diagnostic$rsquared
    diagnostic_mat <- matrix(c(rsquared), nrow = 1, ncol = 1)
    diagnostic_chr <- rbind(
        c("Rsquared"),
        matrix2char(diagnostic_mat, decimal.fmt)
    )
    print.table.md(diagnostic_chr, ...)
    cat("\n")

    ### Residuals
    cat("Scaled residuals", fill = T)
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
