## invert_choldec ----

#' inverted cholesky decomposition of a matrix M
#'
#' @description
#'
#' Note: this should not be confused with the inverse of M derived from the
#' cholesky decomposition (i.e. `chol2inv(M)`)
#'
#' @param M numeric matrix for which the inverse cholesky matrix is to be
#' obtained.
#' @param nugget a length 1 numeric vector containing the nugget that should
#' be added to the var-cov matrix. Default is 0.
#' @param debug logical: debug mode (print additional info)?
#'
#' @return the inverse (lower triangle only) of the cholesky decomposition of M
#' @export
#'
#' @examples #TBA
invert_cholR <- function(M, nugget = 0, debug = FALSE){
  stopifnot(nrow(M) == ncol(M))

  n = nrow(M)

  # handle nugget
  if(nugget != 0){
    if (debug) {print("using nugget")}
    M <- (1 - nugget) * M + nugget * diag(n)
  }

  # return the result
  return(t(backsolve(chol(M), diag(n))))
}


## fitGLS ----
## This function calls invert_cholR()

#' fit GLS to remote sensing data
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param V n x n numeric covariance matrix
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#' @param nugget nugget to be added to variance matrix. see `?invert_cholR()`
#'
#' @return list containing relevant output paramters from model fit
#'
#' Note: this output should be converted to an S3 class and values should be
#' extracted by defining `print.method()` and other methods.
#'
#' @export
#'
#' @examples #TBA
fitGLS <- function(X, V, y, X0 = NULL, nugget = 0){
  stopifnot(all.equal(nrow(X), ncol(V), nrow(V), length(y)))
  n <- nrow(X)
  invcholV <- invert_cholR(V, nugget)
  xx <- invcholV %*% X
  yy <- invcholV %*% y
  varX <- crossprod(xx)
  beta <- as.numeric(solve(varX, crossprod(xx, yy)))
  SSE <- as.numeric(crossprod(yy - xx %*% beta))
  MSE <- SSE/(n - ncol(xx))
  varcov <- MSE * solve(varX)
  se <- diag(varcov)^0.5
  t <- beta/se
  df.t <- n - ncol(xx)
  p.t <- 2 * pt(abs(t), df = df.t, lower.tail = F) # compute in R not C++

  logdetV <- -2 * sum(log(diag(invcholV)))
  logLik <- -0.5 * (n * log(2 * pi) + n * log((n-ncol(xx))*MSE/n) + logdetV + n)

  if(is.null(X0) | missing(X0)){
    X0 <- matrix(1, nrow = n, ncol = 1)
  }
  xx0 <- invcholV %*% X0
  varX0 <- crossprod(xx0)
  betahat0 <- as.numeric(solve(varX0, crossprod(xx0, yy)))
  SSE0 <- as.numeric(crossprod(yy - xx0 %*% betahat0))
  df0 <- length(betahat0)
  MSE0 <- SSE0/(n - ncol(xx))
  MSR <- (SSE0 - SSE)/(ncol(xx) - ncol(xx0))
  logLik0 <- -0.5 * (n * log(2 * pi) + n * log((n-df0)*MSE0/n) + logdetV + n)

  varX0 <- t(xx0) %*% xx0
  varcov0 <- MSE0 * solve(varX0)
  se0 <- diag(varcov0)^0.5

  if (ncol(xx) > 1) { # when/why would ncol(xx) == 1 ever ?
    FF <- (n - ncol(xx))/(ncol(xx) - ncol(xx0)) * (SSE0 - SSE)/SSE
    df1.F <- ncol(xx) - ncol(xx0)
    df2.F <- n - ncol(xx)
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  } else {
    FF <- (n - 1) * (SSE0 - SSE)/SSE
    df1.F <- 1
    df2.F <- n - 1
    p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
    df.F <- c(df1.F, df2.F)
  }

  return(list(betahat = beta, VarX = varX, SSE = SSE, MSE = MSE, varcov = varcov,
              SE = se, tstat = t, pval.t = p.t, dft = df.t, logDetV = logdetV,
              logLik = logLik, betahat0 = betahat0, SE0 = se0, SSE0 = SSE0,
              SSR = SSE0 - SSE, MSE0 = MSE0, MSR = MSR, df0 = df0,
              logLik0 = logLik0, Fstat = FF, pval.F = p.F, df.F = df.F
              ))
}

# Rcpp::sourceCpp("Cpp/GLS_Chol.cpp")
## This function can be depracated if I implement t-test and F-test in C++

#' wrapper for fitGLS_cpp
#'
#' @details at present, the C++ implementations cannot calculate pvalues from
#' t tests or F tests. This function is meant to add them after the fact
#'
#' Note: this function should be deprecated and instead an S3 class and a
#' `coef.method()` should be created to obtain the p values from any `starmod.gls`
#' objects
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param V n x n numeric covariance matrix
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#'
#' @return
#' @export
#'
#' @examples #TBA
wrap_glscpp <- function(X, V, y, X0 = NULL){
  if(is.null(X0) | missing(X0)){
    X0 <- matrix(as.double(1), nrow = nrow(X), ncol = 1)
  }
  out <- fitGLS_cpp(X = X, V = V, y = y, X0 = X0)

  out$pval.t <- 2 * pt(abs(out$tstat), df = out$dft, lower.tail = F)
  out$pval.F <- pf(out$Fstat, df1 = out$df.F[1], df2 = out$df.F[2], lower.tail = F)
}

## fitV ----


#' Fit varcov matrix from distance matrix
#'
#' @details
#'
#' Note: currently the "exponential" version yields a non-positive definitive
#' (and therefore invalid) V matrix.
#'
#' @param Dist n x n numeric distance matrix
#' @param spatialcor spatial correlation parameter(s)
#' @param fun function with which to transform the distance matrix into the
#' varcov matrix
#'
#' @return n x n varcov matrix
#' @export
#'
#' @examples #TBA
fitV <- function(Dist, spatialcor, fun = "exponential"){
  switch(fun,
         # exponential (with alias)
         "exponential" = exp(-Dist/spatialcor), ## This version yeilds non positive-definitive matrix!!!
         "exp" = exp(-Dist/spatialcor),
         # exponential power (with alias)
         "exponential-power" = exp(-(Dist/spatialcor[1])^spatialcor[2]),
         "exp-pwr" = exp(-(Dist/spatialcor[1])^spatialcor[2]),
         # taper-spherical (aliases)
         "taper-spherical" = taper_sphere(Dist, spatialcor),
         "taper" = taper_sphere(Dist, spatialcor),
         "sphr" = taper_sphere(Dist, spatialcor))
}

#' Fit varcov matrix from distance matrix
#'
#' @details
#'
#' Note: currently the "exponential" version yields a non-positive definitive
#' (and therefore invalid) V matrix.
#'
#' Additionally, V.fit() is functionally identical to fitV() however, fitV()
#' should be used instead. This function should be phased out
#'
#' @param Dist n x n numeric distance matrix
#' @param spatialcor spatial correlation parameter(s)
#' @param fun function with which to transform the distance matrix into the
#' varcov matrix
#'
#' @return
# @export ## Don't export this function - it is the same as fitV
#'
#' @examples #TBA
V.fit <- function(Dist, spatialcor, FUN = "exponential") {

  if (FUN == "exponential")
    return(exp(-Dist/spatialcor))

  if (FUN == "exponential-power")
    return(exp(-(Dist/spatialcor[1])^spatialcor[2]))

  if (FUN == "taper-spherical")
    return(taper.spherical(Dist, spatialcor))

}

## fitNugget ----
## This function calls fitGLS() and invert_choldec()

#' Obtain the maximum likelihood nugget estimate from a `starmod.gls` model
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param V n x n numeric covariance matrix
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#' @param int length 2 numeric vector specifying the inverval over which to
#' search of the MLE estimator. Default = `c(0, 1)`.
#' @param tol the desired accuracy of the mathematical optimization.
#' see `?optimize()`.
#'
#' @return maximum likelihood estimate of the nugget
#' @export
#'
#' @examples #TBA
fitNugget <-  function(X, V, y, int = c(0,1), tol = .00001){
  N.opt <- optimize(f = function(nug){return(fitGLS(X, V, y, nugget = nug)$logLik)},
                    interval = int, tol = tol, maximum = TRUE)
  if(N.opt$maximum < tol){
    N0.LL <- fitGLS(X, V, y, nugget = 0)$logLik
    if(N0.LL > N.opt$objective){
      N.opt$maximum <- 0
    }
  }
  return(N.opt$maximum)
}

#' Obtain the maximum likelihood nugget estimate from a `starmod.gls.cpp` model
#'
#' @details this function is a wrapper that calls the C++ version of fitGLS()
#' and can be faster than fitNugget(). However, this functionallity is also
#' available in the C++ version of fitGLS()
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param V n x n numeric covariance matrix
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#' @param int length 2 numeric vector specifying the inverval over which to
#' search of the MLE estimator. Default = `c(0, 1)`.
#' @param tol the desired accuracy of the mathematical optimization.
#' see `?optimize()`.
#'
#' @return maximum likelihood estimate of the nugget
#' @export
#'
#' @examples #TBA
fitNugget_Rcpp <-  function(X, V, y, int = c(0,1), tol = .00001){
  N.opt <- optimize(f = function(nug){return(LogLikGLS_cpp(nugget = nug, X, V, y))},
                    interval = int, tol = tol, maximum = TRUE)
  if(N.opt$maximum < tol){
    N0.LL <- LogLikGLS_cpp(nugget = 0, X, V, y)
    if(N0.LL > N.opt$objective){
      N.opt$maximum <- 0
    }
  }
  return(N.opt$maximum)
}

## Partitioned GLS ----
#' fit GLS model by partitioning remote sensing data
#'
#' @details
#'
#' Note: This function is not complete yet. Use the C++ version instead
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param V n x n numeric covariance matrix
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#' @param nugget nugget to be added to variance matrix. see `?invert_cholR()`
#' @param npart integer: number of of partitions to divide the data into
#' @param mincross intiger: minimum number of partition pairs from which to
#' calculate statistics (i.e. )
#' @param nug.int interval of nugget passed to fitGLS()
#' @param nug.tol accuracy of nugget calculation passed to fitGLS()
#'
#' @return list of GLS statistics
#' @export
#'
#' @examples #TBA
fitGLS.partition <- function(X, V, y, X0, nugget = 0, npart = 10, mincross = 5,
                             nug.int = c(0, 1), nug.tol = 0.00001){
  ## Select random subsets according to the number of partitions
  n <- nrow(data) # full data n
  nn <- n - (n%%npart) # n divisible by npart
  n.p <- nn/npart # size of each partition
  shuff <- sample(n)[1:nn] # shuffled rows
  # shuff.mat <- matrix(shuff, nrow = npart)
    ## TBA: handle user-defined partitions?

  ## calculate degrees of freedom
  df2 <- n.p - (ncol(X) - 1)
  df0 <- n.p - (ncol(X0) - 1)
  df1 <- df0 - df2

  ## adjust the minimum number of crossed partitions
  if(mincross > npart | is.na(mincross)|is.null(mincross) | missing(mincross)){
    mincross <- npart
  }

  ## loop through each partition and gather results
  # for(partition in seq_len(npart)){ ## lapply is better for now
  results <- lapply(seq_len(npart), function(partition){
    ## subset the full data according to the partion
    subset <- (partition - 1)*n.p + (seq_len(n.p))
    tmp <- fitGLS(X = X[subset, ], V = V[subset, subset], y = y[subset],
                  X0 = X0[subset, ], nugget = nugget)

      ## TBA: fit V matrix to individual partitions
      ## TBA: allow for non-fixed nugget

    out <- tmp[c("SSR", "SSE", "SSE0","betahat", "betahat0", "SE", "SE0",
                      "Fstat", "pval.F", "logLik", "logLik0")]

    ## include incvhol, xx, and xx0 for the first few subsets
    if(!is.na(mincross) && partition <= mincross){
      out$invcholV <- invert_cholR(V[subset, subset], nugget = nugget)
      out$xx <- tmp$xx
      out$xx0 <- tmp$xx0
    } else{
      out$invcholV <- NULL
      out$xx <- NULL
      out$xx0 <- NULL
    }
    return(out)})

  ## Calculate pairwise cross-partition statistics
return(results)

}

#' fit GLS model by partitioning remote sensing data via Rcpp
#'
#' @param X n x p numeric design matrix for predictor variables
#' @param y length n numeric resposne vector
#' @param X0 n x p0 null numeric design matrix
#' @param Dist distance matrix
#' @param spatcor spatial correlation parameter(s)
#' @param Vfit.fun function to use for Vfit() calculation
#' @param npart number of partitions
#' @param mincross minimum number of parition pairs from which to calculate
#' statistics
#' @param nug.int interval of nugget search
#' @param nug.tol accuracy of nugget estimate
#'
#' @return list of GLS statistics
#' @export
#'
#' @examples #TBA
fitGLS.partition_rcpp <- function(X, y, X0, Dist, spatcor,
                                  Vfit.fun = "exponential-power",
                                  npart = 5, mincross = 4,
                                  nug.int = c(0, 1), nug.tol = .00001){

  ## Select random subsets according to the number of partitions
  n <- nrow(X) # full data n
  nn <- n - (n%%npart) # n divisible by npart
  n.p <- nn/npart # size of each partition
  shuff <- sample(n)[1:nn] # shuffled rows
  partition <- matrix(shuff, nrow = npart)
  # shuff.mat <- matrix(shuff, nrow = npart)
  ## TBA: handle user-defined partitions?

  ## calculate degrees of freedom
  df2 <- n.p - (ncol(X) - 1)
  df0 <- n.p - (ncol(X0) - 1)
  df1 <- df0 - df2

  ## adjust the minimum number of crossed partitions
  if(mincross > npart | is.na(mincross)|is.null(mincross) | missing(mincross)){
    mincross <- npart
  }

  out = lapply(seq_len(npart), function(i){
    yi <- y[partition[, i]]
    Xi <- as.matrix(X[partition[,i], ])
    Xi0 <- as.matrix(X0[partition[, i]])
    # loci <- loc[partition[, i], ]
    Vi <- fitV(Dist[partition[, i], partition[, i]],
                spatialcor = spatcor, fun = Vfit.fun)
    save_xx = ifelse(i <= mincross, TRUE, FALSE)
    return(GLS_worker_cpp(yi, Xi, Vi, Xi0, save_xx = save_xx))
  })

  out.cross = lapply(seq_len(mincross - 1), function(x){
    i = x; j = x+1
    Xij = as.matrix(X[partition[, c(i,j)], ])
    # locij = loc[X[partition[, c(i,j)], ]
    Vij <- fitV(Dist[partition[, c(i,j)], partition[, c(i,j)]], spatialcor = spatcor,
                 fun = Vfit.fun)
    Li <- out[[i]]; Lj = out[[j]]
    res = crosspart_worker_cpp(Li, Lj, Vij, df1, df2)
  })
  return(list("part_results" = out, "crosspart_results" = out.cross))
}
