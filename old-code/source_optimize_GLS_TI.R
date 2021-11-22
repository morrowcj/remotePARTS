optimize_GLS_TI <- function (formula, D, V.meth = "exponential-power", nugget = NA, spcor = NA, verbose = FALSE, contrasts = NULL, data, pars.start = c(r = 0.5, a = 1, nug = 0), save_xx = T) 
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
    pars.fixed["r"] = log(spcor[1])
    if (V.meth == "exponential-power") {
        pars.fixed["a"] = log(spcor[2])
    }
    pars.fixed["nug"] = log(nugget/(1-nugget))
    pars.full = pars.start
    pars.full["r"] <- log(pars.full["r"])
    pars.full["a"] <- log(pars.full["a"])
    pars.full["nug"] <- log(pars.full["nug"]/(1-pars.full["nug"]))
    pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
    pars = pars.full[is.na(pars.fixed)]
    if (verbose) {
        cat("optimizing parameters:\n")
    }

    opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat, 
        D = D, verbose = verbose, control=list(reltol=1e-6), method = "Nelder-Mead")
    # opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat, 
        # D = D, verbose = verbose, control=list(reltol=1e-6), method = "L-BFGS-B", lower=c(1e-3,1e-3), upper=c(1,1))
    if (verbose) {
        cat("\n")
    }
    r.ml = exp(opt.out$par["r"])
    a.ml = exp(opt.out$par["a"])
    nug.ml = 1/(1+exp(-opt.out$par["nug"]))
    if (V.meth == "exponential-power") {
        spcor.ml = c(r.ml, a.ml)
    }
    else {
        spcor.ml = r.ml
    }
    V = remotePARTS:::fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)

    GLS.out = remotePARTS:::fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml, save_xx=save_xx)
    GLS.out$model.info$call = match.call()
    return(list(GLS.out=GLS.out, spatial.pars=list(r=r.ml, a=a.ml, nug=nug.ml)))
}

#################################
# example
# library(geosphere)

# b0 <- 0
# b1 <- 10
# s <- 1
# r <- .2
# a <- 1
# nug <- .05

# nSpace <- 40
# xdim <- nSpace
# ydim <- nSpace
# n <- nSpace^2

# location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
# colnames(location) <- c("lng", "lat")
# Dist <- distm(location)
# Dist <- Dist/(max(Dist)/2^.5)

# # Note that r is in the same units and Dist. It makes sense to standardize Dist to have maximum distance = 1 and then use r < 1.
# V <- (1 - nug) * exp(-(Dist/r)^a) + nug * diag(n)
# Dr <- t(chol(V))

# dat <- data.frame(X=location[,"lng"]/max(location[,"lng"]))

# dat$Y <- b1*dat$X + Dr %*% rnorm(n=n)

# # you can run this with either exponential-power or exponential; the latter is going to be a lot faster.
# mod <- optimize_GLS_TI(Y ~ X, data=dat, V.meth="exponential", D=Dist, verbose=T, pars.start = c(r = r, a = a, nug = nug))
# mod