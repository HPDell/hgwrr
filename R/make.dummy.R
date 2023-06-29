#' Make Dummy Variables
#' 
#' Function \code{make.dummy}
#' converts categorical variables in a data frame to dummy variables.
#' 
#' @param data The data frame from which dummy variables need to be extracted.
#' 
#' @return The data frame with extracted dummy variables.
#' 
#' @examples
#' make.dummy(iris["Species"])
#' 
#' @export 
#' 
make.dummy <- function(data) {
    as.data.frame(Reduce(cbind, Map(make.dummy.extract, data, names(data))))
}

#' @description 
#' Function \code{make.dummy.extract}
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
#' make.dummy.extract(iris$Species, "Species")
#' 
#' @rdname make.dummy
#' 
#' @export 
#' 
make.dummy.extract <- function(col, name) {
    UseMethod("make.dummy.extract")
}

#' @details 
#' If \code{col} is a character vector,
#' the function will get unique values of its elements
#' and leave out the last one.
#' Then, all the unique values are combined with the \code{name} argument
#' as names of new columns.
#' 
#' @examples
#' make.dummy.extract(c("top", "mid", "low", "mid", "top"), "level")
#' 
#' @rdname make.dummy
#' @method make.dummy.extract character
#' @export 
#' 
make.dummy.extract.character <- function(col, name) {
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
#' make.dummy.extract(factor(c("far", "near", "near")), "distance")
#' 
#' @rdname make.dummy
#' @method make.dummy.extract factor
#' @export 
#' 
make.dummy.extract.factor <- function(col, name) {
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
#' make.dummy.extract(c(T, T, F), "sold")
#' 
#' @rdname make.dummy
#' @method make.dummy.extract logical
#' @export 
#' 
make.dummy.extract.logical <- function(col, name) {
    dummy <- list(as.numeric(col == TRUE))
    names(dummy) <- c(name)
    as.data.frame(dummy)
}

#' @details
#' If \code{col} is of other types,
#' the default behaviour for extracting dummy variables is
#' just to copy the original value and try to convert it to numeric values.
#' 
#' @rdname make.dummy
#' @method make.dummy.extract default
#' @export 
#' 
make.dummy.extract.default <- function(col, name) {
    dummy <- list(as.numeric(col))
    names(dummy) <- c(name)
    as.data.frame(dummy)
}
