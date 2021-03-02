## invert_chol() helpers ----

#' check if matrix is positive definite
#'
#' @details check if a matrix is 1) square, 2) symmetric, and 3) positive
#' definite
#'
#' @param M numeric matrix
#'
#' @export
check_posdef <- function(M){

  res <- c(sqr = FALSE, sym = FALSE, posdef = FALSE)

  res["sqr"] = nrow(M) == ncol(M)
  # if(!res["sqr"]){
  #   stop("M not square")
  # } # not square

  if(res["sqr"]){

    res["sym"] = isSymmetric(M)
    # if(!res["sym"]){
    #   stop("M not symmetric")
    # } # not symetric

    if(res["sym"]){

      res["posdef"] = all(eigen(M, only.values = TRUE)$values > 1e-8)
      # if(!res["posdef"]){
      #   stop("M not positive definite")
      # }
    }
  }
  return(res)
}

## GLS Helpers ----

## GLS Partition Helpers ----

#' function to calculate partition size or number of partitions
#'
#' @param npix number of pixels in full dataset
#' @param npart number of partitions to create
#' @param partsize size of each partition
#' @param pixels vector of pixel indexes to sample from
#' @param verbose logical: TRUE prints additional info
#' @export
#'
#' @examples
#' # simulate data with 100 pixels and 20 time points
#' dat.M <- matrix(rnorm(100*20), ncol = 20)
#' # 4 partitions (exhaustive)
#' sample_partitions(npix = nrow(dat.M), npart = 4)
#' # partitions with 10 pixels each (exhaustive)
#' sample_partitions(npix = nrow(dat.M), partsize = 10)
#' # 4 partitions each with 10 pixels (non-exhaustive)
#' sample_partitions(npix = nrow(dat.M), npart = 4, partsize = 10)
#'
#' # index of 50 pixels to subset
#' sub.indx <- c(1:10, 21:25, 30:62, 70:71)
#' # 5 partitions (exhaustive) using only the specified pixels
#' sample_partitions(npix = nrow(dat.M), npart = 5, pixels = sub.indx)
sample_partitions <- function(npix, npart = 10, partsize = NA,
                              pixels = NA, verbose = TRUE){

  if(all(!is.na(pixels)) & (length(pixels) > 1)){
    npix = length(pixels)
    from = pixels
  } else {
    from = 1:npix
  }

  ## check which npart of partsize was given
  no.partsize <- (missing(partsize) || is.na(partsize) | is.null(partsize))
  no.npart <- (missing(npart) || is.na(npart) | is.null(npart))

  ## caclulate partition size
  if(no.partsize){
    if(verbose){print("calculating partsize")}
    partsize = (npix - (npix%%npart)) / npart
  }

  ## OR calculate number of partitions
  if(no.npart){
    if(verbose){print("calculating npart")}
    npart = floor(npix/partsize)
  }

  if(npart * partsize > npix){
    stop("npart * partsize may not be greater than npix")
  }

  remainder = npix %% partsize
  samp <- sample(from, size = npix - remainder, replace = FALSE)

  part.mat <- matrix(samp, ncol = npart, nrow = partsize)
  colnames(part.mat) <- paste("part",1:npart, sep = ".")

  return(part.mat)
}

#' calculate degrees of freedom for partitioned GLS
#'
#' @param partsize number of pixels in each partition
#' @param p number of predictors in alternate model
#' @param p0 number of parameters in null model
#'
#' @export
#'
#' @examples
#' calc_dfpart(partsize = 2000, p = 4, p0 = 1)
calc_dfpart <- function(partsize, p, p0){
  stopifnot(length(partsize) == 1)
  df2 = partsize - (p - 1)
  df0 = partsize - (p0 - 1)
  df1 = df0 - df2
  return(c("df1" = df1, "df2" = df2))
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
#' @return list of pvalues, number of bootstrap iterations, and MSR rank
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
  return(c("pval.chisqr" = pvalue))
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
    R <- (1 - rcoef[i]) * diag(npart) + rcoef[i] * matrix(1, npart, npart)
    secoef[i, ] <- (part.SEs[,i] %*% R %*% part.SEs[, i])^.5/npart
  }

  tscore <- coefs/secoef
  pvalue <- 2 * pt(abs(tscore), df=df2 * npart, lower.tail = FALSE)

  ttest <- cbind(coefs, secoef, tscore, pvalue)
  colnames(ttest) <- c("Est", "SE", "t.stat", "pval.t")

  return(p.t = ttest)
}


#' Wrapper for bootsrap test
#'
#' @param part.out output of partitioned GLS model...
#' @param nboot number of bootstraps to run for F-test (NA skips bootstrap)
#'
#' @return list of p values calculated with different methods
#' @export
#'
#' @examples #TBA
GLS.partition.pvalue <- function(part.out, nboot = 2000){
  if(is.finite(part.out$rSSR) & !is.na(nboot)) {
    p.Fmean <- boot_corF(Fmean.obs = part.out$Fmean,
                         rSSR = part.out$rSSR,
                         rSSE = part.out$rSSE,
                         df1 = part.out$df1,
                         df2 = part.out$df2,
                         npart = part.out$npart,
                         nboot = nboot)
  }else{
    p.Fmean <- list(NA,NA,NA)
  }
  pF.part = sapply(part.out$part_results, function(x)x$pval.F)

  pchisqr = cor_chisq(part.out$Fmean, part.out$rSSR,
                      part.out$df1, part.out$npart)

  pF.hoch <- min(p.adjust(pF.part, "hochberg"))
  pF.homm <- min(p.adjust(pF.part, "hommel"))
  pF.fdr <- min(p.adjust(pF.part, "fdr"))

  pF.adjust = c("hochberg" = pF.hoch, "hommel" = pF.homm, "fdr" = pF.fdr,
                "boot" = p.Fmean$pvalue)

  # p.t <- cor_t(coef = z$coef, se.part = z$se.part,
  #                     rcoef = z$rcoef, z$df2, npart = z$npart)




  return(list("pF.boot" = p.Fmean, "p.chisqr" = pchisqr, "pF.adj" = pF.adjust))
}
