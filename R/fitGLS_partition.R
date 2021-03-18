
#' fit GLS model by partitioning remote sensing data via Rcpp
#' @rdname fitGLS.partition_rcpp
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
#' @examples #TBA
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
                               invChol_i = Li$tInvCholV,
                               invChol_j = Lj$tInvCholV,
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



#' @title Extract a partition from a csv file
#'
#' @param part.i which partition to extract (integer)
#' @param csv.path path to csv file containing data
#' @param part.mat a matrix whose columns contain indices for a partition.
#'
#' @return
#' a list with the following elements:
#'
#' \describe{
#'     \item{\code{$y}}{the response vector of length \code{partsize}
#'     corresponding to parition i}
#'     \item{\code{$X}}{a model matrix for parition i}
#'     \item{\code{$coords}}{a matrix or data frame of spatial coordinates for
#'     partition i. The first column are x coordinates and the second are y
#'     coordinates}
#' }
#'
#' @details \code{part_csv()} is meant to be used with \code{fitGLS.partition()}
#' @seealso [fitGLS.partition()]
#'
#' @export
#' @examples
#' n.pix = 30865 # pixels in AK_ndvi_common-land.csv
#' parts = sample_partitions(npix = n.pix, npart = 4, partsize = 1000)
#'
#' data.file = system.file("extdata", "AK_ndvi_common-land.csv",
#'                         package = "remotePARTS")
#'
#' A.out = part_csv(part.i = 1, csv.path = data.file, part.mat = parts)
#'
#' # dimensions of A.out elements:
#' dim(A.out$V)
#' dim(A.out$X)
#' dim(A.out$X0)
#' length(A.out$y)
part_csv <- function(part.i, csv.path, part.mat){
  partDF = data.frame(part = part.mat[, part.i]) #df with 1 column: current partition

  # use SQL syntax with read.csv.sql to only read specific rows (very fast)
  data = sqldf::read.csv.sql(file = csv.path,
                             sql = "select file.* from file join partDF on file.rowid = partDF.part",
                             dbname = tempfile(),
                             header = TRUE)

    return(list(y = data$cls.coef,
                X = model.matrix(cls.coef ~ 0 + land, data = data),
                coords = data[, c("lng", "lat")]))
}

#' @title Compute distance matrix in kilometers
#'
#' @param coords coordinate matrix
#' @param coords2 optional coordinate matrix
#'
#' @return distance matrix
#'
#' @details this function is simply a wrapper for \code{geosphere::distm()}
#' @seealso [geosphere::distm()]
#'
#' @export
dist_km <- function(coords, coords2 = NULL){
  if(is.null(coords2)){
    return(geosphere::distm(coords)/1000)
  } else {
    return(geosphere::distm(coords, coords2)/1000)
  }
}

#' @title Partitioned GLS
#' @description Fit a partitioned GLS using the PARTS and calculate
#' cross-partition statitics.
#' @rdname fitGLS.partition
#'
#'
#' @param part_f function to partition data. See details for more info
#' @param dist_f function to calculate distance. See details for more info
#' @param V.meth method passed to \code{fitV()}
#' @param spatcor spatial correlation used by \code{fitV()}
#' @param partsize number of pixels in each partition
#' @param npart number of partitions
#' @param mincross number of partition pairs to calculate cross-partition
#' statistics from.
#' @param X0 null model matrix with \code{partsize} rows
#' @param ... additional arguments passed to \code{part_f}
#'
#' @details
#'
#' \code{fitGLS.partition()} calls the function specified by \code{part_f} to
#' get the partitions. \code{part_f} is called \code{npart} times and uses an
#' iterator as the first argument (i.e. \code{part_f(1)},
#' \code{part_f(2)}, ... \code{part_f(npart)}). A GLS is fit to each partition
#' and cross-partition statistics are calculated for each pair of partitions up
#' to the smaller of \code{mincross} and \eqn{\code{npart} - 1}.
#'
#' \code{part_f} can be any function that takes an integer i (from 1 to
#' \code{npart}) as its first argument and whose output is a list with at
#' least the following elements:
#'
#' \describe{
#'     \item{\code{$y}}{the response vector of length \code{partsize}
#'     corresponding to parition i}
#'     \item{\code{$X}}{a model matrix for parition i}
#'     \item{\code{$coords}}{a matrix or data frame of spatial coordinates for
#'     partition i. The first column are x coordinates and the second are y
#'     coordinates}
#' }
#'
#' \code{dist_f} can be any function that returns a distance matrix. This
#' function should be able to calculate pairwise distances from a single
#' coordinate matrix and should also be able to calculate distances from
#' a pair of coordinate matrices (as in \code{distgeo::distm()}).
#' \code{dist_km()} is an example of an appropriate distance function and
#' is the default for \code{GLS.partition}.
#'
#' @seealso [fitGLS()], [sample_partition()], [calc_df()], [dist_km()],
#' [part_csv()], [geosphere::distgeo()],
#'
#' @export
#'
#' @examples
#' n.pix = 30865 # pixels in AK_ndvi_common-land.csv
#' parts = sample_partitions(npix = n.pix, npart = 4, partsize = 1000)
#' data.file = system.file("extdata", "AK_ndvi_common-land.csv",
#'                         package = "remotePARTS")
#'
#' GLS.part = fitGLS.partition(part_f = "part_csv", dist_f = "dist_km",
#'                          partsize = nrow(parts), npart = ncol(parts),
#'                          V.meth = "exponential", spatcor = .5,
#'                          # additional arguments passed to part_csv():
#'                          csv.path = data.file, part.mat = parts)
fitGLS.partition <- function(part_f = "part_csv", dist_f = "dist_km",
                          V.meth = "exponential", spatcor,
                          partsize, npart, mincross = 6, X0, ...){
  # Setup ----
  ## match the input functions
  func <- match.fun(part_f)
  D_func <- match.fun(dist_f)

  ## Default X0
  if (missing(X0)) {
    X0 = model.matrix(rep(0, partsize) ~ 1)
  }

  ## adjust mincross
  if (mincross > (npart - 1)){
    mincross <- npart - 1
  }

  ## Main Loop
  for (i in 1:(npart - 1)) {
    j = i + 1
    if(i <= mincross){xx.print = TRUE}

    # partition i ----
    if (i == 1) {
      ## function output
      out.i <- func(i, ...) # out.i <- func(i, data.file, parts)
      out.j = NULL

      stopifnot("part_f(i)$X must be a matrix" = is.matrix(out.i$X))
      stopifnot("part_f(i)$coords must be a matrix" = is.matrix(out.i$coords))

      ## calculate df
      dfs <- calc_dfpart(partsize, p = ncol(out.i$X), p0 = ncol(X0))
      ## setup output
      betas = matrix(NA, ncol = ncol(out.i$X), nrow = npart,
                     dimnames = list(NULL, colnames(out.i$X)))
      SEs = betas
      t.stats = betas
      pvals.t = betas
      mod.stats = matrix(NA, nrow = npart, ncol = 7,
                         dimnames = list(NULL, c("nugget","logLik", "SSE", "MSE",
                                                 "MSR", "F.stat", "pval.F")))
      cross.SS = matrix(NA, nrow = mincross, ncol = 2,
                           dimnames = list(NULL, c("rSSR", "rSSE")))
      rcoefs = matrix(NA, nrow = mincross, ncol = ncol(out.i$X),
                      dimnames = list(NULL, colnames(out.i$X)))
      ## calculate V
      D.i <- D_func(out.i$coords)
      V.i <- fitV(D.i, spatcor, V.meth)
      ## fit GLS
      gls.i <- GLS_worker(y = out.i$y, X = out.i$X, V = V.i, X0 = X0, save_xx = xx.print)
    } else {
      ## copy out, V, and gls from previous j
      out.i <-  out.j
      V.i <- V.j
      gls.i <- gls.j
    }
    ## fill in the output
    betas[i, ] <- gls.i$betahat
    SEs[i, ] <- gls.i$SE
    t.stats[i, ] <- gls.i$tstat
    pvals.t[i, ] <- gls.i$pval.t
    mod.stats[i, ] <- c(gls.i$nugget, gls.i$logLik, gls.i$SSE, gls.i$MSE, gls.i$MSR,
                        gls.i$Fstat, gls.i$pval.F)

    # partition j ----
    ## function output
    out.j <- func(j, ...) # out.j <- func(j, data.file, parts)
    ## check i and j are comparable
    stopifnot("part_f(i) and part_f(j) must have equal dimensions" = {
      all.equal(dim(out.i$X), dim(out.j$X))
      all.equal(length(out.i$y), length(out.j$y))
    })
    ## Calculate V
    D.j <- D_func(out.j$coords)
    V.j <- fitV(D.j, spatcor, V.meth)
    ## fit GLS
    gls.j <- GLS_worker(y = out.j$y, X = out.j$X, V = V.j, X0 = X0, save_xx = xx.print)
    ## fill in data for the last partition
    if(j == npart){
      betas[j, ] <- gls.j$betahat
      SEs[j, ] <- gls.j$SE
      t.stats[j, ] <- gls.j$tstat
      pvals.t[j, ] <- gls.j$pval.t
      mod.stats[j, ] <- c(gls.j$nugget, gls.j$logLik, gls.j$SSE, gls.j$MSE, gls.j$MSR,
                          gls.j$Fstat, gls.j$pval.F)
    }


    # cross-partition ----
    ## cross-partition varcovar matrix
    D.ij <- D_func(out.i$coords, out.j$coords)
    V.ij <- fitV(D.ij, spatcor, V.meth)

    ## cross-partition statistics
    if(i <= mincross){
      cross.ij = crosspart_worker(xxi = gls.i$xx, xxj = gls.j$xx,
                                  xxi0 = gls.i$xx0, xxj0 = gls.j$xx0,
                                  invChol_i = gls.i$invcholV,
                                  invChol_j = gls.j$invcholV,
                                  Vsub = V.ij,
                                  nug_i = gls.i$nugget, nug_j = gls.j$nugget,
                                  df1 = dfs[1], df2 = dfs[2])
      cross.SS[i, ] = c(cross.ij$rSSRij, cross.ij$rSSEij)
      rcoefs[i, ] = as.vector(cross.ij$rcoefij)
    }

  }
    # Overall Statistics ----
    fmean = mean(mod.stats[, "F.stat"])
    rSSRmean = mean(cross.SS[, "rSSR"])
    rSSEmean = mean(cross.SS[, "rSSE"])
    rcoefmean = colMeans(rcoefs)
    coefmean = colMeans(betas)

    # Output ----
    out.list <- list(call = match.call(),
                     part.stats = list("coefficients" = betas,
                                       "SEs" = SEs,
                                       "t.stats" = t.stats,
                                       "pvals.t" = pvals.t,
                                       "mod.stats" = mod.stats),
                     cross.stats = list(cross.SS = cross.SS, rcoefs = rcoefs),
                     overall.stats = list("dfs" = dfs,
                                          "coefmean" = coefmean,
                                          "rcoefmean" = rcoefmean,
                                          "meanstats" = c("fmean" = fmean,
                                                          "rSSRmean" = rSSRmean,
                                                          "rSSEmean" = rSSEmean))
    )
    class(out.list) <- "remoteGLS.parts"
    return(out.list)
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

## OLD first attempt at fitGLS.partition() ----
# #' fit GLS model by partitioning remote sensing data
# #'
# #' @details
# #'
# #' Note: This function is not complete yet. Use the C++ version instead
# #'
# #' @param X n x p numeric design matrix for predictor variables
# #' @param V n x n numeric covariance matrix
# #' @param y length n numeric resposne vector
# #' @param X0 n x p0 null numeric design matrix
# #' @param nugget nugget to be added to variance matrix. see `?invert_cholR()`
# #' @param npart integer: number of of partitions to divide the data into
# #' @param mincross intiger: minimum number of partition pairs from which to
# #' calculate statistics (i.e. )
# #' @param nug.int interval of nugget passed to fitGLS_R()
# #' @param nug.tol accuracy of nugget calculation passed to fitGLS_R()
# #'
# #' @return list of GLS statistics
# #' @export
# #'
# #' @examples #TBA
# fitGLS.partition <- function(X, V, y, X0, nugget = 0, npart = 10, mincross = 5,
#                              nug.int = c(0, 1), nug.tol = 0.00001){
#   ## Select random subsets according to the number of partitions
#   n <- nrow(data) # full data n
#   nn <- n - (n%%npart) # n divisible by npart
#   n.p <- nn/npart # size of each partition
#   shuff <- sample(n)[1:nn] # shuffled rows
#   # shuff.mat <- matrix(shuff, nrow = npart)
#   ## TBA: handle user-defined partitions?
#
#   ## calculate degrees of freedom
#   df2 <- n.p - (ncol(X) - 1)
#   df0 <- n.p - (ncol(X0) - 1)
#   df1 <- df0 - df2
#
#   ## adjust the minimum number of crossed partitions
#   if(mincross > npart | is.na(mincross)|is.null(mincross) | missing(mincross)){
#     mincross <- npart
#   }
#
#   ## loop through each partition and gather results
#   # for(partition in seq_len(npart)){ ## lapply is better for now
#   results <- lapply(seq_len(npart), function(partition){
#     ## subset the full data according to the partion
#     subset <- (partition - 1)*n.p + (seq_len(n.p))
#     tmp <- fitGLS_R(X = X[subset, ], V = V[subset, subset], y = y[subset],
#                     X0 = X0[subset, ], nugget = nugget)
#
#     ## TBA: fit V matrix to individual partitions
#     ## TBA: allow for non-fixed nugget
#
#     out <- tmp[c("SSR", "SSE", "SSE0","betahat", "betahat0", "SE", "SE0",
#                  "Fstat", "pval.F", "logLik", "logLik0")]
#
#     ## include incvhol, xx, and xx0 for the first few subsets
#     if(!is.na(mincross) && partition <= mincross){
#       out$invcholV <- invert_cholR(V[subset, subset], nugget = nugget)
#       out$xx <- tmp$xx
#       out$xx0 <- tmp$xx0
#     } else{
#       out$invcholV <- NULL
#       out$xx <- NULL
#       out$xx0 <- NULL
#     }
#     return(out)})
#
#   ## Calculate pairwise cross-partition statistics
#   return(results)
#
# }
