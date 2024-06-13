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
#' @importFrom stats as.formula
#' @importFrom utils tail
#' 
#' @noRd
#'
parse.formula <- function(formula) {
    model <- list(
        intercept = list(
            random = TRUE,
            fixed = TRUE,
            local = TRUE
        )
    )
    root <- as.formula(formula)
    stack <- list(root)
    re <- list()
    fe <- list()
    le <- list()
    local_mode <- FALSE
    random_mode <- FALSE
    random_start_length <- 0
    local_start_length <- 0
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
                stack <- stack.push(stack, cur[[2]])
            } else if (cur_symbol == '|') {
                if (local_mode) stop("Error in formula: cannot set random effects for local fixed effects.")
                if (length(re) > 0) stop("Error in formula: only can set random effects once.")
                model$group <- as.character(cur[[3]])
                random_mode <- TRUE
                random_start_length <- length(stack)
                stack <- stack.push(stack, cur[[2]])
            } else if (cur_symbol == 'L') {
                if (random_mode) stop("Error in formula: cannot set local fixed effects for random effects.")
                local_mode <- TRUE
                local_start_length <- length(stack)
                stack <- stack.push(stack, cur[[2]])
            } else stop("Error in formula: unrecognized symbol.")
        } else {
            if (random_mode) {
                if (inherits(cur, "numeric") && cur == 0) model$intercept$random = FALSE
                else re <- c(re, cur)
                if (length(stack) == random_start_length) {
                    random_mode <- FALSE
                }
            } else if (local_mode) {
                if (inherits(cur, "numeric") && cur == 0) model$intercept$local = FALSE
                else le <- c(le, cur)
                if (length(stack) == local_start_length) {
                    local_mode <- FALSE
                }
            } else {
                if (inherits(cur, "numeric") && cur == 0) model$intercept$fixed = FALSE
                else fe <- c(fe, cur)
            }
        }
    }
    model$random.effects <- rev(as.character(re))
    model$fixed.effects <- rev(as.character(fe))
    model$local.fixed.effects <- rev(as.character(le))
    model
}

#' Push stack
#'
#' @param s A `list`, `vector` or any other object which works with
#'      function \code{\link[base]{c}}
#' @param x An object which can be appended to `s`.
#' 
#' @rdname parse.formula
#' 
#' @noRd 
#'
stack.push <- function(s, x) c(s, x)

#' Pop stack
#'
#' @param s A `list`, `vector` or any other object which works with
#'      function \code{\link[base]{c}}
#' 
#' @rdname parse.formula
#' 
#' @noRd 
#'
stack.pop <- function(s) s[-length(s)]
