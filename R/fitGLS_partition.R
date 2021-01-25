
#' fit GLS model by partitioning remote sensing data via Rcpp
#' @rdname fitGLS_partition
#'
#' @param X n x p numeric design matrix of predictor variables
#' @param y length n numeric response vector
#' @param X0 n x p0 null numeric design matrix
#' @param Dist distance matrix
#' @param spatcor spatial correlation parameter(s)
#' @param Vfit.fun function to use for Vfit() calculation
#' @param npart number of partitions
#' @param mincross minimum number of partition pairs from which to calculate
#' statistics
#' @param nug.int interval of nugget search
#' @param nug.tol accuracy of nugget estimate
#' @param workerB_cpp logical: should the cpp version of worker function be
#' used? this argument is deprecated and was just used to test
#'
#' @return list of GLS statistics
#'
#' @details \code{fitGLS.partition_rcpp()} uses \code{GLS_worker()} to
#' obtain partition-specific statistics and then \code{crosspart_worker()}
#' to obtain cross-partition statistics.
#'
#' Currently, \code{fitGLS.parition_rcpp()} needs the full distance
#' matrix \code{Dist}. \code{fitGLS.partition()} is in development and will
#' not require the full matrix be loaded into memory.
#'
#' @examples
#'
#' @export
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
               spatialcor = spatcor, method = Vfit.fun)
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
                method = Vfit.fun)
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

# # # Test how to call a function, passed as an argument: ----
# # The best way I can think of to make a versatile and memory-limited
# # fitGLS.partition() function is by taking a function as an argument.
# # This function could be user-defined so that any input structure could
# # be used. this user function func() would need the following characteristics
# #
# # 1) the first argument of the function should be a single integer that
# # indicates which partition the function will be creating/revealing.
# #
#
# # 2) the output should be a list that contains at least the following NAMED
# # elements: "X" a design matrix of predictors with exactly n.p rows; "y"
# # the GLS response vector (e.g. residuals from CLS) that is exactly length n.p;
# # Either "D" an n.p x n.p distance matrix OR "V" a n.p x n.p Covariance matrix.
# # if "D" is output instead of "V", an additional object "Vfit.meth" containing
# # method by which "V" should be obtained from "D" (see Vfit()).
# #
# test.data <- as.data.frame(matrix(rnorm(1000*20), ncol = 20))
#
# test.partition = sample_partitions(1000, npart = 4)
#
# part_funct.A <- function(part.i, data, partition){
#   sub.index = as.vector(partition[,part.i])
#   return(data[sub.index, ])
# }
#
# call_part_funct.A <- function(func, part.i, ...){
#   f <- match.fun(func)
#   f(part.i, ...)
# }
