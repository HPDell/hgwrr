#' Make Dummy Variables
#'
#' Function \code{make_dummy}
#' converts categorical variables in a data frame to dummy variables.
#'
#' @param data The data frame from which dummy variables need to be extracted.
#'
#' @return The data frame with extracted dummy variables.
#'
#' @examples
#' make_dummy(iris["Species"])
#'
#' @export
#'
make_dummy <- function(data) {
  as.data.frame(Reduce(cbind, Map(make_dummy_extract, data, names(data))))
}

#' @description
#' Function \code{make_dummy_extract}
#' converts a column to dummy variables if necessary
#' and assign appropriate names.
#' See the "detail" section for further information.
#' Users can define their own functions to allow the model
#' deal with some types of variables properly.
#'
#' @param col A vector to extract dummy variables.
#' @param name The vector's name.
#'
#' @examples
#' make_dummy_extract(iris$Species, "Species")
#'
#' @rdname make_dummy
#'
#' @export
#'
make_dummy_extract <- function(col, name) {
  UseMethod("make_dummy_extract")
}

#' @details
#' If \code{col} is a character vector,
#' the function will get unique values of its elements
#' and leave out the last one.
#' Then, all the unique values are combined with the \code{name} argument
#' as names of new columns.
#'
#' @examples
#' make_dummy_extract(c("top", "mid", "low", "mid", "top"), "level")
#'
#' @rdname make_dummy
#' @method make_dummy_extract character
#' @export
#'
make_dummy_extract.character <- function(col, name) {
  lev <- unique(col)
  lev <- lev[-length(lev)]
  dummy <- Reduce(c, lapply(lev, function(x) {
    list(as.numeric(col == x))
  }))
  lev_names <- paste(name, lev, sep = ".")
  names(dummy) <- lev_names
  as.data.frame(dummy)
}

#' @details
#' If \code{col} is a factor vector,
#' the function will get its levels and leave out the last one.
#' Then, all level labels are combined with the \code{name} argument
#' as names of new columns.
#'
#' @examples
#' make_dummy_extract(factor(c("far", "near", "near")), "distance")
#'
#' @rdname make_dummy
#' @method make_dummy_extract factor
#' @export
#'
make_dummy_extract.factor <- function(col, name) {
  lev <- levels(col)
  lev <- lev[-length(lev)]
  dummy <- Reduce(c, lapply(lev, function(x) {
    list(as.numeric(col == x))
  }))
  lev_names <- paste(name, lev, sep = ".")
  names(dummy) <- lev_names
  as.data.frame(dummy)
}

#' @details
#' If \code{col} is a logical vector,
#' the function will convert it to a numeric vector
#' with value `TRUE` mapped to `1` and `FALSE` to `0`.
#'
#' @examples
#' make_dummy_extract(c(TRUE, TRUE, FALSE), "sold")
#'
#' @rdname make_dummy
#' @method make_dummy_extract logical
#' @export
#'
make_dummy_extract.logical <- function(col, name) {
  dummy <- list(as.numeric(col == TRUE))
  names(dummy) <- c(name)
  as.data.frame(dummy)
}

#' @details
#' If \code{col} is of other types,
#' the default behaviour for extracting dummy variables is
#' just to copy the original value and try to convert it to numeric values.
#'
#' @rdname make_dummy
#' @method make_dummy_extract default
#' @export
#'
make_dummy_extract.default <- function(col, name) {
  dummy <- list(as.numeric(col))
  names(dummy) <- c(name)
  as.data.frame(dummy)
}
