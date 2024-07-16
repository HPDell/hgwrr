#' Spatial heterogeneity test
#' 
#' @export
spatial_hetero_test <- function(
    x,
    coords,
    ...,
    resample = 5000,
    poly = 2,
    bw = 10,
    verbose = 0
) {
    var_names <- colnames(x)
    tv <- spatial_hetero_perm(x, coords, poly, resample, bw, 0, verbose)
    pv <- sapply(seq_along(tv$t0), function(i) {
        with(tv, mean(t[,i] > t0[i]))
    })
    res <- tv
    res$vars <- var_names
    res$p <- pv
    class(res) <- "shgt"
    res
}

#' Print the result of spatial heterogeneity test
#' 
#' @method print shgt
#' @export 
print.shgt <- function(x, ...) {
    cat("Spatial Heterogeneity Test\n", fill = T)
    cat("\n", fill = T)
    show_tbl <- data.frame(
        t0 = as.vector(x$t0),
        t = colMeans(x$t),
        P = x$p,
        stars = sapply(x$p, pv2stars)
    )
    colnames(show_tbl) <- c("t0", "Avg. t", "Pr(t>t0)", "")
    rownames(show_tbl) <- x$vars
    print(show_tbl)
}
