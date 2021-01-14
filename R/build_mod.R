#' Build a usable model matrix from a formula and data
#'
#' @param formula formula to use
#' @param data data to use
#' @param subset subset of the data
#' @param na.action what to do with NAs
#' @param contrasts optional linear contrasts
#'
#' @return a remoteGLS object
build_mod <- function(formula, data, subset,
                      # weights, offset,  ## can't handle these yet
                      na.action, contrasts = NULL
){
  cl <- match.call() # function call
  mf <- match.call(expand.dots = FALSE) # don't expand ...
  m <- match(c("formula", "data", "subset",
               "weights", "na.action", "offset"),
             names(mf), 0L) # match arguments provided by call
  mf <- mf[c(1L, m)] #function name, plus arguments matched
  mf$drop.unused.levels <- TRUE # show that we dropped levels
  mf[[1L]] <- quote(stats::model.frame) # rename the function call
  mf <- eval(mf, parent.frame()) # evaluate the model frame with the data
  mt <- attr(mf, "terms") # model terms
  y <- model.response(mf, "numeric") # response (vector)
  # w <- as.vector(model.weights(mf)) # model.weights
  # offset <- model.offset(mf) # model offset
  ny <- length(y)
  x <- model.matrix(mt, mf, contrasts) # create model matrix
  return(list(call = cl,
              formula = formula,
              response = y,
              response.name = names(mf)[1],
              var.names = names(mf),
              mod.mat = x))
}
