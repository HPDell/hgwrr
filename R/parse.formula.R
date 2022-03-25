#' Parse a HGWR formula.
#'
#' This function accept an R formula object and extract names of
#' the group variable, local fixed effects, global fixed effects
#' and random effects.
#'
#' @param formula A `formula` object.
#'      Its format is much like the formula used in \code{\link[lme4]{lmer}}
#'      from package ''**lme4**''.
#'
#' @return A `list` consists of:
#'      - `response`: name of dependent (response) variable.
#'      - `group`: name of group variable.
#'      - `random.effects`: a vector of names of random effects.
#'      - `fixed.effects`: a vector of names of fixed effects.
#'
parse.formula <- function(formula) {
    model <- list()
    root <- as.formula(formula)
    stack <- list(root)
    re <- list()
    fe <- list()
    random_mode <- FALSE
    random_start_length <- 0
    while(length(stack) > 0) {
        cur <- tail(stack, 1)[[1]]
        stack <- stack.pop(stack)
        if (length(cur) > 1) {
            cur_symbol <- cur[[1]]
            if (cur_symbol == '~') {
                ### Response and Causes
                if (length(cur) == 3) {
                    ### Extract Response
                    model$response <- as.character(cur[[2]])
                    ### Push right part into stack
                    stack <- stack.push(stack, cur[[3]])
                }
                else stop("Error in formula: cannot extract response.")
            } else if (cur_symbol == '+') {
                ### Left hand
                stack <- stack.push(stack, cur[[2]])
                if (length(cur) > 2)
                    stack <- stack.push(stack, cur[[3]])
            } else if (cur_symbol == '(') {
                random_mode <- TRUE
                random_start_length <- length(stack)
                stack <- stack.push(stack, cur[[2]])
            } else if (cur_symbol == '|') {
                model$group <- as.character(cur[[3]])
                stack <- stack.push(stack, cur[[2]])
            } else stop("Error in formula: unrecognized symbol.")
        } else {
            if (random_mode) {
                re <- c(re, cur)
                if (length(stack) == random_start_length) {
                    random_mode <- FALSE
                }
            } else fe <- c(fe, cur)
        }
    }
    model$random.effects <- as.character(re)
    model$fixed.effects <- as.character(fe)
    model
}

#' Push stack
#'
#' @param s A `list`, `vector` or any other object which works with
#'      function \code{\link[base]{c}}
#' @param x An object which can be appended to `s`.
#'
stack.push <- function(s, x) {
    c(s, x)
}

#' Pop stack
#'
#' @param s A `list`, `vector` or any other object which works with
#'      function \code{\link[base]{c}}
#'
stack.pop <- function(s) {
    s[-length(s)]
}
