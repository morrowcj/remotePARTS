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
#' @param nug.int interval of nugget passed to fitGLS_R()
#' @param nug.tol accuracy of nugget calculation passed to fitGLS_R()
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
    tmp <- fitGLS_R(X = X[subset, ], V = V[subset, subset], y = y[subset],
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


#' Worker function 2 for partitioned GLS
#'
#' @details this is the second worker function for the partitioned GLS analysis.
#'
#' NOTE: currently, there is no native parallel functionality and the partitioned
#' form of the GLS is not implemented entirely in C++. Instead, the R function
#' fitGLS.partition_rcpp() weaves between R and C++ on a single core. While
#' this method is still much faster than the purely R implementation, migration
#' to entirely C++ will greatly improve speed further. This migration requires
#' calculating geographic distances with C++ which I've not yet written.
#'
#' Additionally, there seems to be a memory-related issue with the cpp version
#' of this code. I've
#' successfully used the function when partitions have 100 or fewer rows (too
#' small). However, larger partitions cause a fatal error that causes a crash.
#'
#' @param xxi numeric matrix xx from  partition i
#' @param xxj numeric matrix xx from  partition j
#' @param xxi0 numeric matrix xx0 from  partition i
#' @param xxj0 numeric matrix xx0 from  partition j
#' @param tUinv_i numeric matrix tInvCholV from  partition i
#' @param tUinv_j numeric matrix tInvCholV from  partition j
#' @param Vsub numeric variance matrix for Xij (upper block)
#' @param df1 first degree of freedom
#' @param df2 second degree of freedom
#'
#' @export
#' @examples #TBA
crosspart_worker_R <- function(xxi, xxj, xxi0, xxj0, tUinv_i, tUinv_j,
                             Vsub,
                             # nug_i, nug_j,
                             df1, df2){
  np = nrow(xxi)
  # rescale nuggets
  # nug_i = ifelse(nug_i == 0, 0, (1 - nug_i)/nug_i)
  # nug_j = ifelse(nug_i == 0, 0, (1 - nug_i)/nug_i)

  # variance matrix with nuggets (ARE NEVER USED)
  # Vn <- diag(rep(c(nug_i, nug_j), each = np)) + Vij
  # Vn <- Vn[1:np, (np+1):(2*np)] # upper right block

  # calculate stats
  Rij <- crossprod(t(tUinv_i), tcrossprod(Vsub, tUinv_j))

  Hi <- xxi %*% solve(crossprod(xxi)) %*% t(xxi)
  Hj <- xxj %*% solve(crossprod(xxj)) %*% t(xxj)

  Hi0 <- xxi0 %*% solve(crossprod(xxi0)) %*% t(xxi0)
  Hj0 <- xxj0 %*% solve(crossprod(xxj0)) %*% t(xxj0)

  SiR <- Hi - Hi0
  SjR <- Hj - Hj0

  SiE <- diag(np) - Hi
  SjE <- diag(np) - Hj

  # rSSRij <- (SiR %*% (Rij %*% SjR %*% t(Rij)))/df1
  # rSSEij <- (SiE %*% (Rij %*% SjE %*% t(Rij)))/df2

  rSSRij <- matrix(SiR, nrow=1) %*%
    matrix(Rij %*% SjR %*% t(Rij), ncol=1)/df1

  rSSEij <- matrix(SiE, nrow=1) %*%
    matrix(Rij %*% SjE %*% t(Rij), ncol=1)/df2


  # output
  out_lst <- list("Rij" = Rij,
                  "Hi" = Hi,
                  "Hj" = Hj,
                  "Hi0" = Hi0,
                  "Hj0" = Hj0,
                  "SiR" = SiR,
                  "SjR" = SjR,
                  "rSSRij" = rSSRij,
                  "rSSEij" = rSSEij)
  return(out_lst)
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
#' @param workerB_cpp logical: should the cpp version of worker function be
#' used? this argument is deprecated and was just used to test
#'
#' @return list of GLS statistics
#'
#' @examples #TBA
fitGLS.partition_rcpp <- function(X, y, X0, Dist, spatcor,
                                  Vfit.fun = "exponential-power",
                                  npart = 5, mincross = 4,
                                  nug.int = c(0, 1), nug.tol = .00001,
                                  workerB_cpp = TRUE){

  ## Select random subsets according to the number of partitions
  n <- nrow(X) # full data n
  nn <- n - (n%%npart) # n divisible by npart
  n.p <- nn/npart # size of each partition
  shuff <- sample(n)[1:nn] # shuffled rows
  partition <- matrix(shuff, ncol = npart)
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
    gls.out <- GLS_worker(yi, Xi, Vi, Xi0, save_xx = save_xx)

    #add pvalues
    gls.out$pval.t <- sapply(gls.out$tstat, function(x){
      2 * pt(abs(x), df = gls.out$dft, lower.tail = FALSE)
    })
    gls.out$pval.F <- pf(gls.out$Fstat, gls.out$df.F[1], gls.out$df.F[2],
                         lower.tail = FALSE)

    return(gls.out)
  })

  out.cross = lapply(seq_len(mincross - 1), function(x){
    i = x; j = x + 1
    Xij = as.matrix(X[partition[, c(i,j)], ])
    # locij = loc[X[partition[, c(i,j)], ]
    Vij <- fitV(Dist[partition[, c(i,j)], partition[, c(i,j)]], spatialcor = spatcor,
                 fun = Vfit.fun)
    Vsub <- Vij[1:n.p, (n.p+1):(2*n.p)] # off diaganal element
    Li <- out[[i]]; Lj = out[[j]]

    ## use crosspart worker function
    if(workerB_cpp){
    res = crosspart_worker(xxi = Li$xx, xxj = Lj$xx,
                           xxi0 = Li$xx0, xxj0 = Lj$xx0,
                           nug_i = Li$nugget,
                           nug_j = Lj $nugget,
                           invChol_i = Li$invcholV,
                           invChol_j = Lj$invcholV,
                           Vsub = Vsub,
                           df1 = df1, df2 = df2)
    } else {
      res = crosspart_worker_R(xxi = Li$xx, xxj = Lj$xx,
                                 xxi0 = Li$xx0, xxj0 = Lj$xx0,
                                 tUinv_i = Li$tInvCholV,
                                 tUinv_j = Lj$tInvCholV,
                                 Vsub = Vsub,
                                 df1 = df1, df2 = df2)
    }

    return(res)
  })

  # average statistics
  Fmean <-  numeric(length(out))
  coef <- matrix(NA, ncol = ncol(X), nrow = length(out))
  coef0 <- matrix(NA, ncol = ncol(X0), nrow = length(out))
  for (i in 1:length(out)){
    Fmean[i] = out[[i]]$Fstat
    coef[i, ] = out[[i]]$betahat
    coef0[i, ] = out[[i]]$betahat0
  }
  rSSE <- numeric(length(out.cross))
  rSSR <- numeric(length(out.cross))
  for(i in 1:length(out.cross)){
    rSSR[i] = out.cross[[i]]$rSSRij
    rSSE[i] = out.cross[[i]]$rSSEij
  }
  Fmean = mean(Fmean, na.rm = TRUE)
  rSSR = mean(rSSR, na.rm = TRUE)
  rSSE = mean(rSSE, na.rm = TRUE)
  coef = colMeans(coef, na.rm = TRUE)
  coef0 = colMeans(coef0, na.rm = TRUE)

  return(list("part_results" = out, "crosspart_results" = out.cross,
              "Fmean" = Fmean, "rSSR" = rSSR,
                               "rSSE" = rSSE, "coef" = coef,
                               "coef0" = coef0, "npart" = npart,
                               "np" = n.p, "df1" = df1, "df2" = df2))
}

#' bootsrap test
#'
#' @param Fmean.obs Fmean
#' @param rSSR rSSR
#' @param rSSE rSSE
#' @param df1 df1
#' @param df2 df2
#' @param npart partitions
#' @param nboot bootstraps
#'
#' @return
#' @export
#'
#' @examples #TBA
boot_corF <- function(Fmean.obs, rSSR, rSSE, df1, df2,
                                   npart, nboot = 2000){
  part <- rep(1:npart, each=df1)

  rZ <- rSSR^.5/df1
  v.MSR <- diag(df1) - rZ
  v.MSR <- kronecker(diag(npart),v.MSR) + rZ
  D.MSR <- t(chol(v.MSR, pivot=T))
  if(attr(D.MSR, "rank") < npart){
    rank.MSR <- attr(D.MSR, "rank")
    v.MSR <- diag(df1) - .99/df1
    v.MSR <- kronecker(diag(npart),v.MSR) + rZ
    D.MSR <- t(chol(v.MSR, pivot=T))
  }else{
    rank.MSR <- NA
  }

  v.MSE <- (1-rSSE) * diag(npart) + rSSE
  D.MSE <- t(chol(v.MSE))

  count <- 0
  for(boot in 1:nboot){
    Z1 <- D.MSR %*% rnorm(npart*df1)
    MSR.boot <- aggregate(Z1^2, by=list(part), FUN=sum)[,2]/df1
    MSE.boot <- 1 + D.MSE %*% rnorm(npart, mean=0, sd=(2*df2)^.5/df2)
    if(Fmean.obs < mean(MSR.boot/MSE.boot)) count <- count + 1
  }
  pval = ifelse(count/nboot == 0, 1/nboot, count/nboot)
  return(list(pvalue = pval, nboot = nboot, rank.MSR = rank.MSR))
}

#' Correlated chi-squared test
#'
#' @details performs a correlated chi-squared test on cross-partition GLS
#' results.
#'
#' @param Fmean average partition F value
#' @param rSSR cross-partition regression sum of squares
#' @param df1 first degrees of freedom
#' @param npart number of partitions used to split the data
#'
#' @export
cor_chisq <- function(Fmean, rSSR, df1, npart){
  rZ <- rSSR^.5/df1
  v.MSR <- diag(df1) - rZ
  V.MSR <- kronecker(diag(npart),v.MSR) + rZ
  lambda <- eigen(V.MSR)$values
  pvalue <- suppressWarnings(CompQuadForm::imhof(q = npart * df1 * Fmean, lambda = lambda)$Qq)
  pvalue = ifelse(pvalue <= 1e-06, 1e-06, pvalue) # prevent from being negative/too low
  return(pvalue)
}

#' analytical t-test for partitioned GLS
#'
#' @param coefs the coefficients to test - averaged across partitions
#' @param part.SEs matrix of partition standard errors (columns are partitions)
#' @param rcoef cross-coefficient averaged across partition pairs
#' @param df2 second degrees of freedom from partitioned GLS
#' @param npart number of partitions the data was split into
#'
#' @export
cor_t <- function(coefs, part.SEs, rcoef, df2, npart){

  secoef <- matrix(NA, length(coefs), 1)
  for(i in seq_len(length(coefs))){
    R <- (1 - rcoef[i]) * diag(npart) + rcoef[i] * matrix(1,npart,npart)
    secoef[i,] <- (part.SEs[i,] %*% R %*% part.SEs[i,])^.5/npart
  }

  tscore <- coefs/secoef
  pvalue <- 2 * pt(abs(tscore), df=df2 * npart, lower.tail = FALSE)

  ttest <- cbind(coefs, secoef, tscore, pvalue)
  colnames(ttest) <- c("coefs", "se", "tscore", "P")

  return(p.t = ttest)
}

#' Wrapper for bootsrap test
#'
#' @param starpart output of partitioned GLS model...
#' @param nboot number of bootstraps to run for F-test (NA skips bootstrap)
#'
#' @return
#' @export
#'
#' @examples #TBA
GLS.partition.pvalue <- function(starpart, nboot = 2000){
  if(is.finite(starpart$rSSR) & !is.na(nboot)) {
    p.Fmean <- boot_corF(Fmean.obs = starpart$Fmean,
                                      rSSR = starpart$rSSR,
                                      rSSE = starpart$rSSE,
                                      df1 = starpart$df1,
                                      df2 = starpart$df2,
                                      npart = starpart$npart,
                                      nboot = nboot)
  }else{
    p.Fmean <- list(NA,NA,NA)
  }
  pF.part = sapply(starpart$part_results, function(x)x$pval.F)

  pchisqr = cor_chisq(starpart$Fmean, starpart$rSSR,
                             starpart$df1, starpart$npart)

  pF.hoch <- min(p.adjust(pF.part, "hochberg"))
  pF.homm <- min(p.adjust(pF.part, "hommel"))
  pF.fdr <- min(p.adjust(pF.part, "fdr"))

  pF.adjust = c("hochberg" = pF.hoch, "hommel" = pF.homm, "fdr" = pF.fdr,
                "boot" = p.Fmean$pvalue)

  # p.t <- cor_t(coef = z$coef, se.part = z$se.part,
  #                     rcoef = z$rcoef, z$df2, npart = z$npart)




  return(list("pF.boot" = p.Fmean, "p.chisqr" = pchisqr, "pF.adj" = pF.adjust))
}


