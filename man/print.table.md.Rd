% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print.table.md}
\alias{print.table.md}
\title{Print a character matrix as a table.}
\usage{
\method{print}{table.md}(
  x,
  col.sep = "",
  header.sep = "",
  row.begin = "",
  row.end = "",
  table.style = c("plain", "md", "latex"),
  ...
)
}
\arguments{
\item{x}{A character matrix.}

\item{col.sep}{Column separator. Default to \code{""}.}

\item{header.sep}{Header separator. Default to \code{"-"}.}

\item{row.begin}{Character at the beginning of each row.
Default to \code{col.sep}.}

\item{row.end}{Character at the ending of each row.
Default to \code{col.sep}.}

\item{table.style}{Name of pre-defined style.
Possible values are \code{"plain"}, \code{"md"} or \code{"latex"}. Default to \code{"plain"}.}

\item{\dots}{Additional style control arguments.}
}
\value{
No return.
}
\description{
Print a character matrix as a table.
}
\details{
When \code{table.style} is specified, \code{col.sep}, \code{header.sep}, \code{row.begin}
and \code{row.end} would not take effects.
Because this function will automatically set their values.
For each possible value of \code{table.style}, its corresponding style settings
are shown in the following table.
\tabular{llll}{
\tab \strong{\code{plain}} \tab \strong{\code{md}} \tab \strong{\code{latex}} \cr
\code{col.sep}    \tab \code{""}             \tab \code{"|"}         \tab \code{"&"}            \cr
\code{header.sep} \tab \code{""}             \tab \code{"-"}         \tab \code{""}             \cr
\code{row.begin}  \tab \code{""}             \tab \code{"|"}         \tab \code{""}             \cr
\code{row.end}    \tab \code{""}             \tab \code{"|"}         \tab \code{"\\\\"}
}

In this function, characters are right padded by spaces.
}
\seealso{
\code{\link[=print.hgwrm]{print.hgwrm()}}, \code{\link[=summary.hgwrm]{summary.hgwrm()}}.
}
