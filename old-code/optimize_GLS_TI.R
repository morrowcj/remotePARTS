optimize_GLS_TI <- function (formula, D, V.meth = "exponential-power", nugget = NA,
                             spcor = NA, verbose = FALSE, contrasts = NULL, data, pars.start = c(r = 0.5, a = 1, nug = 0), save_xx = F)
{
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

    pars.fixed = pars.start
    lowr.all = c(r = 1e-9, a = 1e-9, nug = 0)
    uppr.all = c(r = 1, a = Inf, nug = 1)
    pars.fixed["r"] = spcor[1]
    pars.fixed["nug"] = nugget
    if (V.meth == "exponential-power") {
        pars.fixed["a"] = spcor[2]
    }
    pars.full = pars.start
    pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
    pars = pars.full[is.na(pars.fixed)]
    lowr = lowr.all[is.na(pars.fixed)]
    uppr = uppr.all[is.na(pars.fixed)]

    if (verbose) {
        cat("optimizing parameters:\n")
    }
    # opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat,
    #     D = D, verbose = verbose, control=list(reltol=1e-6), method = "Nelder-Mead")
    opt.out = optim(par = pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth,
                    modmat = modmat, D = D, verbose = TRUE,
                    control = list(pgtol = 1e-6),
                    lower = lowr, upper = uppr,
                    method = "L-BFGS-B")
    # opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat,
    # D = D, verbose = verbose, control=list(reltol=1e-6), method = "L-BFGS-B", lower=c(1e-3,1e-3), upper=c(1,1))
    if (verbose) {
        cat("\n")
    }
    r.ml = opt.out$par["r"]
    a.ml = opt.out$par["a"]
    nug.ml = opt.out$par["nug"]
    if (V.meth == "exponential-power") {
        spcor.ml = c(r.ml, a.ml)
    }
    else {
        spcor.ml = r.ml
    }
    V = remotePARTS:::fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)

    GLS.out = remotePARTS:::fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml, save_xx=save_xx)
    GLS.out$model.info$call = call
    GLS.out$spcor = spcor.ml
    return(GLS.out)
}

logit = function(p) {log(p / (1 - p))}
inv_logit = function(l) {1 / (1 + exp(-l))}


optim_GLS_func2 <- function(par, y, modmat, D, verbose = FALSE,
                           V.meth = "exponential-power"){

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


optimize_GLS_CM <- function (formula, coords, dist_f = "dist_km",
                             V.meth = "exponential-power", nugget = NA,
                             spcor = NA, verbose = FALSE, contrasts = NULL, data,
                             pars.start = c(r = 0.1, a = 1, nug = 1e-9), save_xx = F,
                             ret.GLS = TRUE)
{
    # Setup ----
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
    lowr.all = c(r = 1e-9, a = 1e-9, nug = 0)
    uppr.all = c(r = 1, a = Inf, nug = 1)
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
    pars["a"] = log(pars["a"])
    pars["nug"] = logit(pars["nug"])
    ## use alternate optimized function
    opt.out = optim(par = pars,
                    fn = optim_GLS_func2,
                    y = y, V.meth = V.meth,
                    modmat = modmat, D = D,
                    verbose = verbose,
                    # control = list(pgtol = 1e-6), # tolerance
                    # lower = lowr, upper = uppr, # bounds
                    method = "Nelder-Mead") # fast optimizer that allows bounds

    ## add blank line after trace
    if (verbose) {
        cat("\n")
    }

    # Collect parameters ----
    # r.ml = opt.out$par["r"] # range param
    # a.ml = opt.out$par["a"] # shape param
    # nug.ml = opt.out$par["nug"] # nugget
    r.ml = exp(opt.out$par["r"]) # range param
    a.ml = exp(opt.out$par["a"]) # shape param
    nug.ml = inv_logit(opt.out$par["nug"]) # nugget
    ## r and/or a
    if (V.meth == "exponential-power") {
        spcor.ml = c(r.ml, a.ml)
    }
    else {
        spcor.ml = r.ml
    }

    # Fit GLS one last time with estimated parameters ----
    if(ret.GLS){
        V = remotePARTS::fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)
        GLS.out = remotePARTS:::fitGLS2(formula = y ~ 0 + modmat, V = V,
                                        nugget = nug.ml, save_xx=save_xx)
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

# #################################
# # example
# library(geosphere)
#
# b0 <- 0
# b1 <- 10
# s <- 1
# r <- .2
# a <- 1
# nug <- .05
#
# nSpace <- 40
# xdim <- nSpace
# ydim <- nSpace
# n <- nSpace^2
#
# location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
# colnames(location) <- c("lng", "lat")
# Dist <- distm(location)
# Dist <- Dist/(max(Dist)/2^.5)
#
# # Note that r is in the same units and Dist. It makes sense to standardize Dist to have maximum distance = 1 and then use r < 1.
# V <- (1 - nug) * exp(-(Dist/r)^a) + nug * diag(n)
# Dr <- t(chol(V))
#
# dat <- data.frame(X=location[,"lng"]/max(location[,"lng"]))
#
# dat$Y <- b1*dat$X + Dr %*% rnorm(n=n)
#
# # you can run this with either exponential-power or exponential; the latter is going to be a lot faster.
# mod <- optimize_GLS_TI(Y ~ X, data=dat, V.meth="exponential", D=Dist, verbose=T, pars.start = c(r = r, a = a, nug = nug))
# mod
