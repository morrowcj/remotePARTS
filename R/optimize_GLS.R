# Optimize GLS

## GLS function to optimize
#' function to optimize GLS parameters r, a, and nugget
#' @param in.pars parameters to optimize
#' @param y response vector
#' @param modmat predictor model matrix
#' @param D distance matrix
#' @param constrain do parameters need to be constrained internally?
#' @rdname optimize_GLS
#' @details When \code{constrain = FALSE}, it is assumed that parameter values
#'  given by \code{in.pars} are already constrained to their desired range. When
#'  When \code{constrain = TRUE}, parameter values are internally transformed.
#'
#'  r and nugget are transformed with an inverse logit
#'  (i.e., \code{function(x){1 / (1 + exp(-x))}}) and \code{a} is transformed
#'  with an exponential function (i.e., exp(a)).
#'
optim_GLS_func <- function(in.pars = c(r = 0.01, a = 1, nug = 0.1),
                               y, modmat, D, verbose = FALSE,
                               V.meth = "exponential-power",
                               constrain = FALSE # should pars be constrained?
){
  # define inverse logit function
  inv_logit = function(l) {1 / (1 + exp(-l))}

  ## constrain parameters, if needed
  if(constrain){
    in.pars["r"] = inv_logit(in.pars["r"])
    if(V.meth == "exponential-power"){in.pars["a"] = exp(in.pars["a"])}
    in.pars["nug"] = inv_logit(in.pars["nug"])
  }

  ## Extract spatial parameters
  if (V.meth == "exponential-power"){
    # force values between 0,1
    spcor = c(r = in.pars["r"], a = in.pars["a"])
    nug = in.pars["nug"]
  } else if (V.meth == "exponential") {
    spcor = c(r = in.pars["r"])
    nug = in.pars["nug"]
  }

  ## Fit V matrix
  V = fitV(Dist = D/max(D), spatialcor = spcor, method = V.meth)

  ## Check that V is valid covariance matrix
  if(any(is.nan(V)) || !all(check_posdef(V))){
    LL = NA # return NA if not valid
    # LL = -9e100 ## return lowest finite LL
  } else {
    ## Calculate log-liklihood
    LL = fitGLS2(y ~ 0 + modmat, V = V, nugget = nug, LL_only = TRUE)
  }
  ## Print the values, if asked
  if(verbose){ print(c(r = unname(spcor[1]),
                       a = unname(spcor[2]),
                       nugget = unname(nug),
                       LL = LL)) }

  # if(constrain){warning("constrain = TRUE: Parameter values need to be back-transformed to be valid.")}

  return(-LL)
}

#' Fit a GLS by estimating r, a, and nugget
#'
#' @param formula model formula
#' @param coords 2-column matrix of x and y coordinates
#' @param dist_f distance function to use
#' @param V.meth covariance method for fitting V
#' @param nugget spatial nugget. If NA, a nugget will be estimated
#' @param spcor spatial correlation paramters either NA, a single value (r) or
#' a named vector: c(r = ..., a = ...). if spcor = NA, these parameters will be
#' estimated
#' @param verbose should additional info be printed?
#' @param data data object for which formula tries to match
#' @param pars.start default starting values for the spcor
#' @param save_xx save cross correlation stats?
#' @param ret.GLS return GLS?
#' @param algorithm method passed to optimizer (default = 'Nelder-Mead')
#' see \code{?optim()} for available methods
#' @param debug boolean: debug mode?
#' @param ... additional arguments passed to \code{fitGLS2()}
#'
#' @rdname optimize_GLS
#'
#' @seealso [fitGLS2()]
#'
#' @return a remoteGLS object
#'
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
#' coords = df[, c("lng", "lat")]
#' D = geosphere::distm(coords)/1000 #distance in km
#'
#' ## fit the GLS, including ML spatial correlation and nugget
#' optimize_GLS(kappa ~ 0 + land, data = df, coords = coords, verbose = FALSE)
#'
#' ## other examples
#'
#' # CG algorithm
#' optimize_GLS(kappa ~ 0 + land, data = df, coords = coords, algorithm = "CG")
#'
#' # L-BFGS-B algorithm
#' optimize_GLS(kappa ~ 0 + land, data = df, coords = coords, algorithm = "L-BFGS-B")
#'
#' # exponential-power covariance function
#' optimize_GLS(kappa ~ 0 + land, data = df, coords = coords, V.meth = "exponential-power")
#'
#' # Exponential-power and L-BFGS-B (NOT RUN)
#' ## Produces Error: can't use exponential-power and L-BFGS-B together
#' if (FALSE)
#'   optimize_GLS(kappa ~ 0 + land, data = df, coords = coords, debug = TRUE,
#'                V.meth = "exponential-power", algorithm = "L-BFGS-B")
#' }
optimize_GLS <- function(formula, coords, dist_f = "dist_km",
                         V.meth = "exponential", nugget = NA,
                         spcor = NA, verbose = FALSE, data,
                         pars.start = c(r = 0.1, a = 1, nug = 0.2),
                         save_xx = FALSE, ret.GLS = FALSE,
                         algorithm = "Nelder-Mead",
                         debug = FALSE,
                         ...){
  # Model Setup ----
  ## formula handling
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  if (is.matrix(y)) {
    stop("response is a matrix: must be a vector")
  }
  ny <- length(y)
  modmat <- model.matrix(mt, mf)
  rm(mf)

  # initial error handling ----
  ## unsupported method
  if (!V.meth %in% c("exponential", "exponential-power")){
    stop("Unsupported V.meth: must be 'exponential' or 'exponential-power'")
  }
  ## incompatible methods
  if (V.meth == "exponential-power" & algorithm == "L-BFGS-B"){
    stop("Cannot use algorithm 'L-BFGS-B' with V.meth 'exponential-power'")
  }

  # paramter setup ----
  pars.fixed = pars.start
  pars.fixed["r"] = spcor[1]
  pars.fixed["nug"] = nugget
  if (V.meth == "exponential-power") {
    pars.fixed["a"] = spcor[2]
  }
  pars.full = pars.start
  pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
  pars = pars.full[is.na(pars.fixed)]

  # Calculate distance ----
  D_func <- match.fun(dist_f)
  D = D_func(coords)

  # Algorithm-specific optimization ----
  if (algorithm == "L-BFGS-B"){
    constrain = FALSE
    opt.out <- tryCatch(
      ## try running optim()
      expr = {optim(par = pars, fn = optim_GLS_func, y = y,
                    V.meth = V.meth,
                    modmat = modmat, D = D, verbose = verbose,
                    constrain = constrain, method = "L-BFGS-B",
                    lower = c(r = 1e-100, nug = 1e-100),
                    upper = c(r = 1, nug = 1))},
      ## catch any errors
      error = function(e){return(list(failed = TRUE, error = e))}
    )

  } else {
    logit = function(p) {log(p / (1 - p))}
    inv_logit = function(l) {1 / (1 + exp(-l))}

    ## prepare paramters for constraint
    pars["r"] = logit(pars["r"])
    if (V.meth == "exponential-power") {
      pars["a"] = log(pars["a"])
    }
    pars["nug"] = logit(pars["nug"])

    opt.out <- tryCatch(
      ## try running optim()
      expr = {optim(par = pars, fn = optim_GLS_func, y = y,
                    V.meth = V.meth,
                    modmat = modmat, D = D, verbose = verbose,
                    constrain = TRUE, method = algorithm)},
      ## catch any errors
      error = function(e){return(list(failed = TRUE, error = e))})

    ## back-transform results
    opt.out$par["r"] = inv_logit(opt.out$par["r"])
    opt.out$par["nug"] = inv_logit(opt.out$par["nug"])
    if (V.meth == "exponential-power"){
      opt.out$par["a"] = exp(opt.out$par["a"])
    }
  }

  # Error handling (optim) ----
  if (!is.null(opt.out$failed)){
    message("Optimization was not able to complete successfully.")
    message("Here's the original error message from optim():")
    message(opt.out$error)
    if (V.meth == "exponential-power"){
      message("If problems persist, try using V.meth = 'exponential' instead.")
      message("'exponential-power' often fails to create a valid covariance matrix.")
    }
    stop("Failed to optimize parameters")
  }

  # collect results ----
  return.list <- list()

  spcor.ml = c(r = unname(opt.out$par["r"]))
  if (V.meth == "exponential-power"){
    spcor.ml["a"] = unname(opt.out$par["a"])
  }
  nug.ml = unname(opt.out$par["nug"])


  if (debug) {
    return.list$in.pars = pars
    return.list$optim.out = opt.out
    return.list$spcor.ml = spcor.ml
    return.list$nug.ml = nug.ml
  }

  if (ret.GLS) {
    V = fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)
    GLS.out = fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml,
                      save_xx = save_xx, ...)
    GLS.out$model.info$call = call
    return.list$GLS = GLS.out
  }

  return.list$pars = c(spcor.ml, nug = nug.ml)
  return.list$scaled.range = unname(spcor.ml["r"] * max(D))

  return(return.list)
}
