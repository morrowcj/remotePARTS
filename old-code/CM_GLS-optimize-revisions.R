optim_GLS_func <- function(par, y, modmat, D, verbose = FALSE,
                           V.meth = "exponential-power"){

  ## logit function
  logit = function(p) {log(p / (1 - p))}
  inv_logit = function(l) {1 / (1 + exp(-l))}

  if (V.meth == "exponential-power"){
    spcor = c(r = exp(par["r"]), a = exp(par["a"]))
  } else {
    spcor = c(r = exp(par["r"]))
  }
  r = exp(unname(par["r"]))
  a = exp(unname(par["a"]))
  nug = inv_logit(unname(par["nug"]))

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


exp_fun = function(d, r){
  return(exp(-d/r))
}

exppow_fun = function(d, r, a){
  return(exp(-(d/r)^a))
}

logit = function(p) {log(p / (1 - p))}
inv_logit = function(l) {1 / (1 + exp(-l))}

## Testing Setup ----
set.seed(510)
# load Alaska 3000 data
data("ndvi_AK3000")
# take a random subset of 100 pixels (to make example fast)
subsamp = sample.int(n = nrow(ndvi_AK3000), size = 100)
# subset the data: we now have 100 pixels, latitude, longitude, and land class
df = ndvi_AK3000[subsamp, c("lng", "lat", "land")] # subset the data
# simulate a response variable kappa
df$kappa = .5*df$lat + .2*df$lng + rnorm(100) #simulate response variable
# calculate distance matrix
D = geosphere::distm(df[, c("lng", "lat")])/1000 #distance in km
mf = model.frame(kappa ~ 0 + land, data = df)
modmat = model.matrix(kappa ~ 0 + land, data = df)
y = model.response(mf, "numeric")

## Testing
r = seq(from = 1e-100, to = 1, length.out = 4)
a = c(seq(1e-100, 2, length.out = 4), 1.99999, 9e100)
nug = seq(from = 0, to = 1, length.out = 4)

in.df = expand.grid(r.in = r, a.in = a, nug.in = nug)
in.df$LL.exp_pwr = NA
in.df$LL.exp = NA

for(i in 1:nrow(in.df)){
  in.pars = c(r = in.df[i, "r.in"], a = in.df[i, "a.in"], nug = in.df[i, "nug.in"])
  in.df[i, "LL.exp_pwr"] = optim_GLS_func_NEW(in.pars, y = y, D = D, modmat = modmat, V.meth = "exponential-power")
  in.df[i, "LL.exp"] = optim_GLS_func_NEW(in.pars, y = y, D = D, modmat = modmat, V.meth = "exponential")
}

in.df

library(ggplot2)
ggplot(in.df, aes(y = is.na(LL.exp_pwr), x = nug.in)) +
  geom_point() +
  facet_grid(round(r.in,2) ~ round(a.in, 2))

ggplot(in.df, aes(y = LL.exp_pwr, x = nug.in)) +
  geom_point() +
  facet_grid(round(r.in,2) ~ round(a.in, 2))

# geom_contour_filled()

# curve(exp_fun(d = x, r = 0) ,from = 0, to = 1) # bad
# curve(exp_fun(d = x, r = 1e-9) ,from = 0, to = 1) # OK
# curve(exp_fun(d = x, r = 1) ,from = 0, to = 1) # OK
#
# curve(exppow_fun(d = x, r = .1, a = 0), from = 0, to = 1) # OK?
# curve(exppow_fun(d = x, r = .1, a = Inf), from = 0, to = 1) # OK?
# curve(exppow_fun(d = x, r = 0, a = Inf), from = 0, to = 1) # bad
# curve(exppow_fun(d = x, r = .1, a = 2), from = 0, to = 1) # bad

## Exponential-power (Needs constraining)
opt.exppow <- optim(par = c(r = logit(.1), a = log(1), nug = logit(.2)),
      fn = optim_GLS_func_NEW, y = y, V.meth = "exponential-power",
      modmat = modmat, D = D, verbose = TRUE, constrain = TRUE,
      method = "Nelder-Mead")
opt.exppow$par
opt.exppow$par.correct = c(r = inv_logit(opt.exppow$par["r"]),
                           a = exp(opt.exppow$par["a"]),
                           nug = inv_logit(opt.exppow$par["nug"]))
round(opt.exppow$par.correct, 4)

## Exponential pre-constrained
opt.exp <- optim(par = c(r = .1, nug = .2),
                 fn = optim_GLS_func_NEW, y = y, V.meth = "exponential",
                 modmat = modmat, D = D, verbose = TRUE, method = "L-BFGS-B",
                 lower = c(r = 1e-100, nug = 1e-100), upper = c(r = 1, nug = 1))
opt.exp$par

## Exponential (Needs constraining)
opt.exp2 <- optim(par = c(r = logit(.1), a = log(1), nug = logit(.2)),
                    fn = optim_GLS_func_NEW, y = y, V.meth = "exponential",
                    modmat = modmat, D = D, verbose = TRUE, constrain = TRUE,
                    method = "Nelder-Mead")
opt.exp2$par
opt.exp2$par.correct = c(r = inv_logit(opt.exp2$par["r"]),
                           a = NA,
                           nug = inv_logit(opt.exp2$par["nug"]))
round(opt.exp2$par.correct, 4)




## New optimize_GLS function
optimize_GLS_NEW = function(formula, coords, dist_f = "dist_km",
                            V.meth = "exponential", nugget = NA,
                            spcor = NA, verbose = FALSE, data,
                            pars.start = c(r = 0.1, a = 1, nug = 0.2),
                            save_xx = FALSE, ret.GLS = TRUE,
                            algorithm = "Nelder-Mead",
                            debug = TRUE,
                            ...){
  inv_logit = function(l) {1 / (1 + exp(-l))}
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

  # Optimizer setup: dependent upon algorithm and V.meth ----
  ## prepare parameters to pass to optim
  pars.fixed = pars.start
  pars.fixed["r"] = spcor[1]
  pars.fixed["nug"] = nugget
  if (V.meth == "exponential-power") {
    pars.fixed["a"] = spcor[2]
  }
  pars.full = pars.start
  pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
  pars = pars.full[is.na(pars.fixed)]

  ## exponential-power method
  if (V.meth == "exponential-power"){

    ### check that an incompatible algorithm isn't used
    if (algorithm == "L-BFGS-B"){
      stop("cannot use algorithm = 'L-BFGS-B' with V.meth = 'exponential-power'")
    }
    ### transform paramters for constraining
    constrain = TRUE
    pars["r"] = logit(pars["r"])
    pars["a"] = log(pars["a"])
    pars["nug"] = logit(pars["nug"])

  ## exponential method
  } else if (V.meth == "exponential"){

    ### L-BFGS-B algorithm
    if (algorithm == "L-BFGS-B"){
      constrain = FALSE
    ### All other algorithms
    } else {
      constrain = TRUE
      pars["r"] = logit(pars["r"])
      pars["nug"] = logit(pars["nug"])
    }
  ## incompatible V.meth
  } else {
    stop("V.meth must be either 'exponential' or 'exponential-power'")
  }

  # Calculate Distance matrix ----
  D_func <- match.fun(dist_f)
  D = D_func(coords)

  # Run Optimizer ----
  if (algorithm == "L-BFGS-B"){
    opt.out <- tryCatch(
      ## try running optim()
      expr = {optim(par = pars, fn = optim_GLS_func_NEW, y = y,
                    V.meth = V.meth,
                    modmat = modmat, D = D, verbose = verbose,
                    constrain = constrain, method = algorithm,
                    lower = c(r = 1e-100, nug = 1e-100),
                    upper = c(r = 1, nug = 1))},
      ## catch any errors
      error = function(e){
        message("Optimizer was not able to complete successfully.")
        message("Here's the original error message from optim():")
        message(e, "\n")
        return(NA)
      }
    )

  } else {

    opt.out <- tryCatch(
      ## try running optim()
      expr = {optim(par = pars, fn = optim_GLS_func_NEW, y = y,
                    V.meth = V.meth,
                    modmat = modmat, D = D, verbose = verbose,
                    constrain = constrain, method = algorithm)},
      ## catch any errors
      error = function(e){
        message("Optimizer was not able to complete successfully.")
        message("Here's the original error message from optim():")
        message(e, "\n")
        return(NA)
      }
    )
  }

  # Output ----
  ## additional error handling
  if (any(is.na(opt.out))) {
    if(V.meth == "exponential-power"){
      message("try using V.meth = 'exponential' instead.")
      message("'exponential-power' often fails to create a valid covariance matrix.")
    }
    stop("Failed to optimize parameters.")
  ## Return results, if no error exists
  } else {

    if(constrain){
      opt.out$par["r"] = inv_logit(opt.out$par["r"])
      opt.out$par["nug"] = inv_logit(opt.out$par["nug"])
      opt.out$par["a"] = exp(opt.out$par["a"])
      # opt.out$par = c(r = inv_logit(opt.out$par["r"]),
      #                 a = exp(opt.out$par["a"]),
      #                 nug = inv_logit(opt.out$par["nug"]))
    }

    if(V.meth == "exponential-power"){
      spcor.ml = c(r = unname(opt.out$par["r"]),
                   a = unname(opt.out$par["a"]))
    } else if (V.meth == "exponential"){
      spcor.ml = c(r = unname(opt.out$par["r"]))
    }

    nug.ml = unlist(opt.out$par["nug"])
    scaled.range = unname(spcor.ml["r"] * max(D))

    if (debug){
      return(list(opt.out = opt.out,
                  pars.in = pars,
                  spcor.ml = spcor.ml,
                  nug.ml = nug.ml,
                  pars = c(spcor.ml, nug.ml),
                  scaled.range = scaled.range
                  ))
    } else if (ret.GLS){
      V = fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)
      GLS.out = fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml,
                        save_xx = save_xx, ...)
      GLS.out$model.info$call = call
      return(list(GLS = GLS.out,
                  pars = c(spcor.ml, nug.ml),
                  scaled.range = scaled.range))
    } else {
      return(list(pars = c(spcor.ml, nug.ml),
                  scaled.range = scaled.range))
    }

  }

}

# ## Testing (Part 2) - Contour plot
# r = seq(from = 0, to = 1, length.out = 100)
# # a = c(seq(0, 2, length.out = 20))
# nug = seq(from = 0, to = 1, length.out = 100)
#
# tmp.df = expand.grid(r.in = r, nug.in = nug)
# tmp.df$LL = NA
#
# pb <- txtProgressBar(min = 0, max = nrow(tmp.df), style = 3)
# for(i in 1:nrow(tmp.df)){
#   tmp.pars = c(r = tmp.df[i, "r.in"], nug = tmp.df[i, "nug.in"])
#   tmp.df[i, "LL"] = optim_GLS_func_NEW(tmp.pars, y = y, D = D, modmat = modmat, V.meth = "exponential")
#   setTxtProgressBar(pb, i)
# }
#
# tmp.df %>% ggplot(aes(x = r.in, y = nug.in, z = LL)) +
#   geom_contour_filled(label = "-LL")
