# Optimize GLS

## GLS function to optimize
#' function to optimize GLS parameters r, a, and nugget
#' @param par parameters to optimize
#' @param y response vector
#' @param modmat predictor model matrix
#' @rdname optimize_GLS
optim_GLS_func <- function(par, y, modmat, D, verbose = FALSE,
                           V.meth = "exponential-power"){

  if (V.meth == "exponential-power"){
    spcor = c(r = par["r"], a = par["a"])
  } else {
    spcor = c(r = par["r"])
  }
  r = unname(par["r"])
  a = unname(par["a"])
  nug = unname(par["nug"])

  ## fit a variance matrix
  V = fitV(Dist = D/max(D), spatialcor = spcor, method = V.meth)

  ## Don't use a non-positive definitive matrix
  if(any(is.nan(V))){return(NA)}
  if(!all(check_posdef(V))){return(NA)}

  ## calculate the log-liklihood
  LL = fitGLS2(y ~ 0 + modmat, V = V, nugget = nug, LL_only = TRUE)

  ## Print the values, if asked
  if(verbose){ print(c(r = r, a = a, nugget = nug, LL = LL)) }

  ## return the negative log-liklihood
  return(-LL)
}


#' Fit a GLS by estimating r, a, and nugget
#'
#' @rdname optimize_GLS
#'
#' @param formula model formula
#' @param D Distance matrix
#' @param V.meth method passed to \code{fitV()} default: "exponential-power"
#' @param nugget NA: find maximum liklihood nugget
#' @param spcor NA: find maximum liklihood spatial correlation
#' @param verbose should the optimizer steps be printed to the console?
#' @param data optional data to search for objects in formula
#' @param contrasts possible contrasts object
#' @param ... additional arguments passed to \code{fitGLS2()}
#'
#' @seealso [fitGLS2()]
#'
#' @return a remoteGLS object
#' @export
#'
#' @examples
#' set.seed(916)
#'
#' ## load Alaska 3000 data
#' data("ndvi_AK3000")
#'
#' ## take a random subset of 100 pixels (to make example fast)
#' subsamp = sample.int(n = nrow(ndvi_AK3000), size = 100)
#'
#' ## subset the data: we now have 100 pixels, latitude, longitude, and land class
#' df = ndvi_AK3000[subsamp, c("lng", "lat", "land")] # subset the data
#'
#' ## simulate a response variable kappa
#' df$kappa = .5*df$lat + .2*df$lng + rnorm(100) #simulate response variable
#'
#' ## calculate distance matrix
#' D = geosphere::distm(df[, c("lng", "lat")])/1000 #distance in km
#'
#' ## fit the GLS, including ML spatial correlation and nugget
#' optimize_GLS(kappa ~ 0 + land, data = df, D = D, verbose = FALSE)
optimize_GLS <- function(formula, D, V.meth = "exponential-power",
                         nugget = NA, spcor = NA, verbose = FALSE,
                         contrasts = NULL, data, ...){

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
  modmat <- model.matrix(mt, mf, contrasts) # create model matrix
  rm(mf) # delete the large model frame from memory

  pars.start = c(r = .5, a = 1, nug = 0)

  # values to be fit are NA
  pars.fixed = pars.start
  pars.fixed["r"] = spcor[1]
  if(V.meth == "exponential-power"){ pars.fixed["a"] = spcor[2] }
  pars.fixed["nug"] = nugget

  pars.full = pars.start
  pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]

  ## parameters to be passed to optim()
  pars = pars.full[is.na(pars.fixed)]

  if(verbose){
    cat("optimizing parameters:\n")
  }
  # if(missing(data)){
  #   opt.df = NULL
  # } else {opt.df = data}
  opt.out = optim(pars, fn = optim_GLS_func, y = y, modmat = modmat,
                  D = D, verbose = verbose)
  if(verbose){cat("\n")}

  r.ml = opt.out$par["r"]
  a.ml = opt.out$par["a"]
  nug.ml = opt.out$par["nug"]
  if(V.meth == "exponential-power"){
    spcor.ml = c(r.ml, a.ml)
  } else {
    spcor.ml = r.ml
  }

  V =  fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)

  GLS.out = fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml, ...)
  GLS.out$model.info$call = match.call()

  return(GLS.out)
}

# GLS.test <- optimize_GLS(y = y, modmat = modmat, D = D, V.meth = "exponential-power",
#                          nugget = NA, spcor = NA, verbose = TRUE, form0 = y ~ 1,
#                          no_F = TRUE)
# GLS.test2 <- optimize_GLS(y = y, modmat = modmat, D = D, V.meth = "exponential-power",
#                          nugget = NA, spcor = NA, verbose = FALSE, form0 = y ~ 1,
#                          no_F = FALSE)
