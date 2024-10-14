#' Print a character matrix as a table.
#'
#' @param x A character matrix.
#' @param col_sep Column separator. Default to `""`.
#' @param header_sep Header separator. Default to `"-"`.
#' If `header_sep` only contains one character,
#' it will be repeated for each column.
#' If it contains more than one character,
#' it will be printed below the first row.
#' @param row_begin Character at the beginning of each row.
#' Default to `col_sep`.
#' @param row_end Character at the ending of each row.
#' Default to `col_sep`.
#' @param table_before Characters to be printed before the table.
#' @param table_after Characters to be printed after the table.
#' @param table_style Name of pre-defined style.
#' Possible values are `"plain"`, `"md"`, `"latex"`, or `"booktabs"`.
#' Default to `"plain"`.
#' @param \dots Additional style control arguments.
#'
#' @return No return.
#'
#' @details
#' When `table_style` is specified, `col_sep`, `header_sep`, `row_begin`
#' and `row_end` would not take effects.
#' Because this function will automatically set their values.
#' For each possible value of `table_style`, its corresponding style settings
#' are shown in the following table.
#' \tabular{llll}{
#'                   \tab \strong{\code{plain}} \tab \strong{\code{md}} \tab \strong{\code{latex}} \cr
#' \code{col_sep}    \tab \code{""}             \tab \code{"|"}         \tab \code{"&"}            \cr
#' \code{header_sep} \tab \code{""}             \tab \code{"-"}         \tab \code{""}             \cr
#' \code{row_begin}  \tab \code{""}             \tab \code{"|"}         \tab \code{""}             \cr
#' \code{row_end}    \tab \code{""}             \tab \code{"|"}         \tab \code{"\\\\"}
#' }
#'
#' In this function, characters are right padded by spaces.
#'
#' @seealso [print.hgwrm()], [summary.hgwrm()].
print_table_md <- function(
  x, col_sep = "", header_sep = "",
  row_begin = "", row_end = "",
  table_before = NA_character_, table_after = NA_character_,
  table_style = c("plain", "md", "latex", "booktabs"),
  ...
) {
  ### Special characters
  if (!missing(table_style)) {
    if (table_style == "md") {
      x <- apply(x, c(1, 2), function(t) {
        gsub("([\\|*])", "\\\\\\1", t)
      })
    } else if (table_style == "latex" || table_style == "booktabs") {
      x <- apply(x, c(1, 2), function(t) {
        gsub("([&%])", "\\\\\\1", t)
      })
    }
  }
  ### Process formats
  x_length <- apply(x, c(1, 2), nchar)
  x_length_max <- apply(x_length, 2, max)
  x_fmt <- sprintf("%%%ds", x_length_max)
  if (!missing(table_style)) {
    table_style <- match.arg(table_style)
    if (table_style == "md") {
      col_sep <- ifelse(missing(col_sep), "|", col_sep)
      header_sep <- ifelse(missing(header_sep), "-", header_sep)
      row_begin <- ifelse(missing(row_begin), "|", row_begin)
      row_end <- ifelse(missing(row_end), "|", row_end)
      table_before <- ifelse(missing(table_before), "", table_before)
    } else if (table_style == "latex") {
      col_sep <- ifelse(missing(col_sep), "&", col_sep)
      header_sep <- ifelse(missing(header_sep), "", header_sep)
      row_begin <- ifelse(missing(row_begin), "", row_begin)
      row_end <- ifelse(missing(row_end), "\\\\", row_end)
      table_before <- ifelse(
        missing(table_before),
        sprintf(
          "\\begin{tabular}{%s}",
          paste(rep("l", times = ncol(x)), collapse = "")
        ),
        table_before
      )
      table_after <- ifelse(missing(table_after), "\\end{tabular}", table_after)
    } else if (table_style == "plain") {
      col_sep <- ifelse(missing(col_sep), "", col_sep)
      header_sep <- ifelse(missing(header_sep), "", header_sep)
      row_begin <- ifelse(missing(row_begin), "", row_begin)
      row_end <- ifelse(missing(row_end), "", row_end)
    } else if (table_style == "booktabs") {
      col_sep <- ifelse(missing(col_sep), "&", col_sep)
      header_sep <- ifelse(missing(header_sep), "\\midrule", header_sep)
      row_begin <- ifelse(missing(row_begin), "", row_begin)
      row_end <- ifelse(missing(row_end), "\\\\", row_end)
      table_before <- ifelse(
        missing(table_before),
        sprintf(
          "\\begin{tabular}{%s}\n\\toprule",
          paste(rep("l", times = ncol(x)), collapse = "")
        ),
        table_before
      )
      table_after <- ifelse(
        missing(table_after),
        "\\bottomrule\n\\end{tabular}",
        table_after
      )
    } else {
      stop("Unknown table.style.")
    }
  }
  ### Print table
  if (!is.na(table_before)) cat(table_before, fill = TRUE)
  for (c in seq_len(ncol(x))) {
    if (x_length_max[c] > 0) {
      cat(
        ifelse(c == 1, row_begin, col_sep),
        sprintf(x_fmt[c], x[1, c]), ""
      )
    }
  }
  cat(paste0(row_end, "\n"))
  if (nchar(header_sep) == 1) {
    for (c in seq_len(ncol(x))) {
      if (x_length_max[c] > 0) {
        header_sep_full <- paste(rep("-", x_length_max[c]),
          collapse = ""
        )
        cat(
          ifelse(c == 1, row_begin, col_sep),
          sprintf(header_sep_full), ""
        )
      }
    }
    cat(paste0(row_end, "\n"))
  } else if (nchar(header_sep) > 1) {
    cat(paste0(header_sep, "\n"))
  }
  for (r in 2:nrow(x)) {
    for (c in seq_len(ncol(x))) {
      if (x_length_max[c] > 0) {
        cat(
          ifelse(c == 1, row_begin, col_sep),
          sprintf(x_fmt[c], x[r, c]), ""
        )
      }
    }
    cat(paste0(row_end, "\n"))
  }
  if (!is.na(table_after)) cat(table_after, fill = TRUE)
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
    mc <- apply(m, seq_along(dim(m)), function(x) {
      sprintf(fmt, x)
    })
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
  if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else if (p < 0.1) {
    "."
  } else {
    " "
  }
}
