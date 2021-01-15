## fitGLS2

#' fit a remotePARTS GLS model
#'
#' @param formula formula to build the model with
#' @param data object containing the data
#' @param V Variance matrix
#' @param nugget nugget
#' @param form.0 null model formula (default: "y ~ 1")
#' @param save_xx should xx be saved?
#' @param threads how many threads (functionality not yet available)
#' @param contrasts optional linear contrasts to use
#' @param ... additional arguments passed to \code{optimize_nugget}
#'
#' @return a remoteGLS object
#' @export
#'
#' @examples TBA
fitGLS2 <- function(formula, data, V, nugget = 0,
                    form.0 = NULL,
                    save_xx = FALSE,
                    threads = 1,
                    contrasts = NULL,
                    ...){

  ## Parse formula arguments to make model matrix ----
  call <- match.call() # function call
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
  if (is.matrix(y)){stop("response is a matrix: must be a vector")}
  ny <- length(y)
  X <- model.matrix(mt, mf, contrasts) # create model matrix
  rm(mf) # delete the large model frame from memory

  ## Handle missing nugget (NULL or NA) ----
  if (is.null(nugget) || is.na(nugget)){
    nugget = optimize_nugget(X, V, y, ...)
  }

  ## Build null model----
  if (is.null(form.0)){
    form.0 = formula(y ~ 1)
  } else {
    form.0 = formula(form.0)
  }
  # conditionally assign X0
  X0 <- if (missing(data) || is.null(data))
    model.matrix(form.0)
  else
    model.matrix(form.0, data)

  ## coerce to matrices ----
  # X = as(X, "matrix")
  # V = as(V, "matrix")
  # # y = as(y, "matrix")
  # X0 = as(X0, "matrix")

  ## error handling ----
  stopifnot(all(is.double(X), is.double(V), is.double(y), is.double(X0)))
  stopifnot(all.equal(nrow(X), nrow(V), nrow(X0), ny))
  stopifnot(all(check_posdef(V)))

  GLS <- remoteGLS(form = formula)
  GLS$model.info$call <- call
  GLS$nugget = nugget

  ## Run GLS ----
  .Call(`_remoteSTAR_fitGLS2_cpp`, GLS, X, V, y, X0, nugget, save_xx, threads)

  # add in p values
  GLS$pval.t <- 2 * pt(abs(GLS$tstat), df = GLS$dft, lower.tail = F)
  GLS$pval.F <- pf(GLS$Fstat, df1 = GLS$df.F[1], df2 = GLS$df.F[2], lower.tail = F)

  ## Return ----
  return(GLS)


  # return(list(test.out = list(X = X, y = y, X0 = X0),
  #             out = out))

}

# fitGLS2(formula = test.y ~ 0 + test.land, V = test.V)
#
#
# ## arguments
# formula = formula(test.y ~ 0 + test.land)
# data = NULL
# subset = NULL
# contrast = NULL
# test.model =
# ## Build empty remote GLS (right now the C++ function doesn't use this)
# test.GLS = remoteGLS(test.model$formula)
# # test.GLS$model.info$call <- call
#
#
# out <- .Call(`_remoteSTAR_fitGLS_cpp`,
#              y = test.model$response,
#              X = test.model$mod.mat,
#              V = V,
#              X0 = X0,
#              nugget = nugget,
#              save_xx = save_xx,
#              threads = threads)
# out
