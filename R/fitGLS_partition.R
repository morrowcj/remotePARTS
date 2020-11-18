## Re-working of fitGLS.partition

#' function to calculate partition size or number of partitions
#'
#' @param npix number of pixels in full dataset
#'
#' @examples
#' dat.M <- matrix(rnorm(3000*20), ncol = 20)
#' dat.rast <- raster::raster(dat.M)
sample_partitions <- function(npix, npart = 10, partsize = NA, verbose = TRUE){


  ## check which npart of partsize was given
  no.partsize <- (missing(partsize) || is.na(partsize) | is.null(partsize))
  no.npart <- (missing(npart) || is.na(npart) | is.null(npart))

  ## caclulate partition size
  if(no.partsize){
    if(verbose){print("calculating partsize\n")}
    partsize = (npix - (npix%%npart)) / npart
  }

  ## OR calculate number of partitions
  if(no.npart){
    if(verbose){print("calculating npart\n")}
    npart = floor(npix/part.size)
  }

  if(npart * partsize > npix){
    stop("npart * partsize may not be greater than npix")
  }

  part.mat <- matrix(sample.int(npix, size = npix, replace = FALSE),
                     ncol = npart, nrow = partsize)
  colnames(part.mat) <- paste("part",1:npart, sep = ".")

  return(part.mat)
}


##### NOTE:: This function is the same name as the one in GLS_functions.R and
##### should eventually be changed to just fitGLS.partition()

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
#' @export
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
                             invChol_i = Li$invcholV,
                             invChol_j = Lj$invcholV,
                             nug_i = Li$nugget,
                             nug_j = Li$nugget,
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
