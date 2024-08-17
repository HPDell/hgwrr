

#' Print a character matrix as a table.
#'
#' @param x A character matrix.
#' @param col.sep Column separator. Default to `""`.
#' @param header.sep Header separator. Default to `"-"`.
#' If `header.sep` only contains one character, it will be repeated for each column.
#' If it contains more than one character, it will be printed below the first row.
#' @param row.begin Character at the beginning of each row.
#' Default to `col.sep`.
#' @param row.end Character at the ending of each row.
#' Default to `col.sep`.
#' @param table.before Characters to be printed before the table.
#' @param table.after Characters to be printed after the table.
#' @param table.style Name of pre-defined style.
#' Possible values are `"plain"`, `"md"`, `"latex"`, or `"booktabs"`. Default to `"plain"`.
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
print.table.md <- function(
    x, col.sep = "", header.sep = "",
    row.begin = "", row.end = "",
    table.before = NA_character_, table.after = NA_character_,
    table.style = c("plain", "md", "latex", "booktabs"),
    ...
) {
    ### Special characters
    if (!missing(table.style)) {
        if (table.style == "md") {
            x <- apply(x, c(1, 2), function(t) {
                gsub("([\\|*])", "\\\\\\1", t)
            })
        } else if (table.style == "latex" || table.style == "booktabs") {
            x <- apply(x, c(1, 2), function(t) {
                gsub("([&%])", "\\\\\\1", t)
            })
        }
    }
    ### Process formats
    x.length <- apply(x, c(1, 2), nchar)
    x.length.max <- apply(x.length, 2, max)
    x.fmt <- sprintf("%%%ds", x.length.max)
    if (!missing(table.style)) {
        table.style <- match.arg(table.style)
        if (table.style == "md") {
            col.sep <- ifelse(missing(col.sep), "|", col.sep)
            header.sep <- ifelse(missing(header.sep), "-", header.sep)
            row.begin <- ifelse(missing(row.begin), "|", row.begin)
            row.end <- ifelse(missing(row.end), "|", row.end)
            table.before <- ifelse(missing(table.before), "", table.before)
        } else if (table.style == "latex") {
            col.sep <- ifelse(missing(col.sep), "&", col.sep)
            header.sep <- ifelse(missing(header.sep), "", header.sep)
            row.begin <- ifelse(missing(row.begin), "", row.begin)
            row.end <- ifelse(missing(row.end), "\\\\", row.end)
            table.before <- ifelse(missing(table.before), sprintf("\\begin{tabular}{%s}", paste(rep("l", times = ncol(x)), collapse = "")), table.before)
            table.after <- ifelse(missing(table.after), "\\end{tabular}", table.after)
        } else if (table.style == "plain") {
            col.sep <- ifelse(missing(col.sep), "", col.sep)
            header.sep <- ifelse(missing(header.sep), "", header.sep)
            row.begin <- ifelse(missing(row.begin), "", row.begin)
            row.end <- ifelse(missing(row.end), "", row.end)
        } else if (table.style == "booktabs") {
            col.sep <- ifelse(missing(col.sep), "&", col.sep)
            header.sep <- ifelse(missing(header.sep), "\\midrule", header.sep)
            row.begin <- ifelse(missing(row.begin), "", row.begin)
            row.end <- ifelse(missing(row.end), "\\\\", row.end)
            table.before <- ifelse(missing(table.before), sprintf("\\begin{tabular}{%s}\n\\toprule", paste(rep("l", times = ncol(x)), collapse = "")), table.before)
            table.after <- ifelse(missing(table.after), "\\bottomrule\n\\end{tabular}", table.after)
        } else {
           stop("Unknown table.style.")
        }
    }
    ### Print table
    if (!is.na(table.before)) cat(table.before, fill = T)
    for(c in 1:ncol(x)) {
        if(x.length.max[c] > 0)
            cat(ifelse(c == 1, row.begin, col.sep),
                sprintf(x.fmt[c], x[1, c]), "")
    }
    cat(paste0(row.end, "\n"))
    if (nchar(header.sep) == 1) {
        for(c in 1:ncol(x)) {
            if(x.length.max[c] > 0) {
                header.sep.full <- paste(rep("-", x.length.max[c]),
                                         collapse = "")
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(header.sep.full), "")
            }
        }
        cat(paste0(row.end, "\n"))
    } else if (nchar(header.sep) > 1) {
        cat(paste0(header.sep, "\n"))
    }
    for (r in 2:nrow(x)) {
        for (c in 1:ncol(x)) {
            if(x.length.max[c] > 0)
                cat(ifelse(c == 1, row.begin, col.sep),
                    sprintf(x.fmt[c], x[r, c]), "")
        }
        cat(paste0(row.end, "\n"))
    }
    if (!is.na(table.after)) cat(table.after, fill = T)
}

#' Convert a numeric matrix to character matrix according to a format string.
#'
#' @param m A numeric matrix.
#' @param fmt Format string. Passing to [base::sprintf()].
#'
#' @seealso [base::sprintf()], [print.hgwrm()], [print.summary.hgwrm()].
#' @noRd 
matrix2char <- function(m, fmt = "%.6f") {
    mc <- NULL
    if ("array" %in% class(m)) {
        mc <- apply(m, seq(length(dim(m))), function(x) { sprintf(fmt, x) })
    } else {
        mc <- sprintf(fmt, m)
    }
    mc
}

#' Convert p values to stars
#' 
#' @param p The p value
#' @return Stars representing significance level. See details below.
#' 
#' @details 
#' - `***` for $p<=0.001$;
#' - `**` for $0.001<p<=0.01$;
#' - `*` for $0.01<p<=0.05$;
#' - `.` for $0.05<p<=0.1$;
#' - ` ` for $0.1<p$.
#' 
#' @noRd
pv2stars <- function(p) {
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1) "."
    else " "
}