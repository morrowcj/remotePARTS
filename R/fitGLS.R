## invert_choldec ----

#' Title
#'
#' @param M
#' @param nugget
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
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
## This function calls invert_choldec()

#' Title
#'
#' @param X
#' @param V
#' @param y
#' @param X0
#' @param nugget
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param X
#' @param V
#' @param y
#' @param X0
#'
#' @return
#' @export
#'
#' @examples
wrap_glscpp <- function(X, V, y, X0 = NULL){
  if(is.null(X0) | missing(X0)){
    X0 <- matrix(as.double(1), nrow = nrow(X), ncol = 1)
  }
  out <- fitGLS_cpp(X = X, V = V, y = y, X0 = X0)

  out$pval.t <- 2 * pt(abs(out$tstat), df = out$dft, lower.tail = F)
  out$pval.F <- pf(out$Fstat, df1 = out$df.F[1], df2 = out$df.F[2], lower.tail = F)
}

## fitDistVar ----


#' Title
#'
#' @param Dist
#' @param spatialcor
#' @param fun
#'
#' @return
#' @export
#'
#' @examples
fitDistVar_R <- function(Dist, spatialcor, fun = "exponential"){
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

#' Title
#'
#' @param Dist
#' @param spatialcor
#' @param FUN
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param V
#' @param y
#' @param int
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param V
#' @param y
#' @param int
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param X
#' @param V
#' @param y
#' @param X0
#' @param nugget
#' @param npart
#' @param mincross
#' @param nug.int
#' @param nug.tol
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param y
#' @param X0
#' @param Dist
#' @param spatcor
#' @param Vfit.fun
#' @param npart
#' @param mincross
#' @param nug.int
#' @param nug.tol
#'
#' @return
#' @export
#'
#' @examples
fitGLS.partition_rcpp <- function(X, y, X0, Dist, spatcor,
                                  Vfit.fun = "exponential-power",
                                  npart = 5, mincross = 4,
                                  nug.int = c(0, 1), nug.tol = .00001){

  ## Select random subsets according to the number of partitions
  n <- nrow(data) # full data n
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
    Vi <- V.fit(Dist[partition[, i], partition[, i]],
                spatialcor = spatcor, FUN = Vfit.fun)
    save_xx = ifelse(i <= mincross, TRUE, FALSE)
    return(GLS_worker_cpp(yi, Xi, Vi, Xi0, save_xx = save_xx))
  })

  out.cross = lapply(seq_len(mincross - 1), function(x){
    i = x; j = x+1
    Xij = as.matrix(X[partition[, c(i,j)], ])
    # locij = loc[X[partition[, c(i,j)], ]
    Vij <- V.fit(Dist[partition, c(i,j), partition, c(i,j)], spatialcor = spatcor,
                 FUN = Vfit.fun)
    Li <- out[[i]]; Lj = out[[j]]
    res = crosspart_worker_cpp(Li, Lj, Vij, df1, df2)
  })
  return(list("part_results" = out, "crosspart_results" = out.cross))
}
