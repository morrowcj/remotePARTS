# Optimize GLS

## GLS function to optimize
#' function to optimize GLS parameters r, a, and nugget
#' @param par parameters to optimize
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
optim_GLS_func_NEW <- function(in.pars = c(r = 0.01, a = 1, nug = 0.1),
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
#' @param contrasts linear contrasts to apply
#' @param data data object for which formula tries to match
#' @param pars.start default starting values for the spcor
#' @param save_xx save cross correlation stats?
#' @param ret.GLS return GLS?
#' @param ... additional arguments passed to \code{fitGLS2()}
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
#' D = geosphere::distm(df[, c("lng", "lat")])/1000 #distance in km
#'
#' ## fit the GLS, including ML spatial correlation and nugget
#' optimize_GLS(kappa ~ 0 + land, data = df, coords = df[, c("lng", "lat")], verbose = FALSE)
optimize_GLS <- function (formula, coords, dist_f = "dist_km",
                          V.meth = "exponential-power", nugget = NA,
                          spcor = NA, verbose = FALSE, contrasts = NULL, data,
                          pars.start = c(r = 0.1, a = 1, nug = 1e-9), save_xx = F,
                          ret.GLS = TRUE, ...)
{
  # Setup ----
  ## logit function
  logit = function(p) {log(p / (1 - p))}
  inv_logit = function(l) {1 / (1 + exp(-l))}
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
  modmat <- model.matrix(mt, mf, contrasts)
  rm(mf)
  ## setup optimizer parameters
  pars.fixed = pars.start
  pars.fixed["r"] = spcor[1]
  pars.fixed["nug"] = nugget
  if (V.meth == "exponential-power") {
    pars.fixed["a"] = spcor[2]
  }
  pars.full = pars.start
  pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
  pars = pars.full[is.na(pars.fixed)]
  ## setup upper and lower bounds for spatial parameters
  lowr.all = c(r = 1e-9, a = 1e-9, nug = 1e-9)
  uppr.all = c(r = 1, a = 3, nug = 1)
  lowr = lowr.all[is.na(pars.fixed)]
  uppr = uppr.all[is.na(pars.fixed)]

  # Run Optimizer ----
  if (verbose) {
    cat("optimizing parameters:\n")
  }
  ## calculate D once
  D_func <- match.fun(dist_f)
  D = D_func(coords)
  ## optimize over the custom function
  # opt.out = optim(par = pars,
  #                 fn = remotePARTS:::optim_GLS_func,
  #                 y = y, V.meth = V.meth,
  #                 modmat = modmat, D = D,
  #                 verbose = verbose,
  #                 control = list(pgtol = 1e-6), # tolerance
  #                 lower = lowr, upper = uppr, # bounds
  #                 method = "L-BFGS-B") # fast optimizer that allows bounds

  ## transform the parameters
  pars["r"] = log(pars["r"])
  if(V.meth == "exponential-power"){pars["a"] = log(pars["a"])}
  pars["nug"] = logit(pars["nug"])

  ## use alternate optimized function
  if(V.meth == "exponential-power"){
      opt.out = optim(par = pars,
                  fn = optim_GLS_func,
                  y = y, V.meth = V.meth,
                  modmat = modmat, D = D,
                  verbose = verbose,
                  control = list(pgtol = 1e-6),
                  method = "Nelder-Mead") # fast optimizer that allows bounds

  } else {
    opt.out = optim(par = pars,
                    fn = optim_GLS_func,
                    y = y, V.meth = V.meth,
                    modmat = modmat, D = D,
                    verbose = verbose,
                    control = list(pgtol = 1e-6),
                    lower = lowr, upper = uppr,
                    method = "L-BFGS-B") # fast optimizer that allows bounds
  }


  ## add blank line after trace
  if (verbose) {
    cat("\n")
  }

  # Collect parameters ----
  r.ml = opt.out$par["r"] # range param
  a.ml = opt.out$par["a"] # shape param
  nug.ml = opt.out$par["nug"] # nugget
  # r.ml = exp(opt.out$par["r"]) # range param
  # a.ml = exp(opt.out$par["a"]) # shape param
  # nug.ml = inv_logit(opt.out$par["nug"]) # nugget
  ## r and/or a
  if (V.meth == "exponential-power") {
    spcor.ml = c(r.ml, a.ml)
  }
  else {
    spcor.ml = r.ml
  }

  # Fit GLS one last time with estimated parameters ----
  if(ret.GLS){
    V = fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)
    GLS.out = fitGLS2(formula = y ~ 0 + modmat, V = V,
                                    nugget = nug.ml, save_xx=save_xx, ...)
    GLS.out$model.info$call = call
    ret.list = list(GLS = GLS.out,
                    spatial.pars = c(r.ml, a.ml, nug.ml),
                    scaled.range = unname(r.ml * max(D)))
  } else {
    ret.list = list(spatial.pars = c(r.ml, a.ml, nug.ml),
                    scaled.range = unname(r.ml * max(D)))
  }

  # Return statement ----
  return(ret.list)
}


#' OLD Fit a GLS by estimating r, a, and nugget
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
#'
optimize_GLS_OLD <- function(formula, D, V.meth = "exponential-power",
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
