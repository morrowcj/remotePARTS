
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
#' @rdname part_data
#'
#' @param part_i which partition to extract (integer)
#' @param part_form formula to make the model matrix
#' @param part_df object that contains all the data
#' @param part_mat a matrix whose columns contain indices for a partition.
#' @param part_locvars character vector of coordinate variable names. Default is c("lng", "lat)
#' @param part_form0 formula to make null model matrix (if NULL (default), "y ~ 1")
#'
#' @return
#' a list with the following elements:
#'
#' \describe{
#'     \item{\code{$y}}{the response vector of length \code{partsize}
#'     corresponding to parition i}
#'     \item{\code{$X}}{a model matrix for parition i}
#'     \item{\code{$X0}}{a null model matrix for parition i}
#'     \item{\code{$coords}}{a matrix or data frame of spatial coordinates for
#'     partition i. The first column are x coordinates and the second are y
#'     coordinates}
#' }
#'
#' @export
#'
#' @examples
#' ## using part_data and AK_ndvi_common-land.csv:
#'
#' n.pix = 30865 # pixels in AK_ndvi_common-land.csv
#' parts = sample_partitions(npix = n.pix, npart = 4, partsize = 1000)
#' data.file = system.file("extdata", "AK_ndvi_common-land.csv",
#'                         package = "remotePARTS")
#' df = data.table::fread(data.file, data.table = FALSE)
#'
#' A.out = part_data(1, part_form = cls.coef ~ 0 + land, part_df = df, part_mat = parts)
#'
#' ## look at first 6 values:
#' lapply(A.out, head)
#'
#' #' # dimensions of A.out elements:
#' dim(A.out$X)
#' dim(A.out$X0)
#' length(A.out$y)
#'
part_data <- function(part_i, part_form, part_df, part_mat, part_locvars = c("lng", "lat"), part_form0 = NULL){
  prt = part_mat[, part_i]
  df = as.data.frame(part_df)
  df_prt = df[prt, ]
  mf = stats::model.frame(formula(part_form), data = df_prt)
  resp = model.response(mf)

  X = stats::model.matrix(formula(part_form), data = df_prt)

  if(is.null(part_form0)){
    X0 = model.matrix(resp ~ 1)
  } else {
    X0 = model.matrix(formula(part_form0), data = df_prt)
  }

  return(list(y = as.vector(resp),
              X = as.matrix(X),
              X0 = as.matrix(X0),
              coords = as.matrix(df_prt[, part_locvars])))
}


#' @rdname part_data
#'
#' @param part_csv_path path to csv file containing data
#'
#' @details \code{part_csv()} is meant to be used with \code{fitGLS.partition()}
#' @seealso [fitGLS.partition()]
#'
#' @export
#'
#' @examples
#'
#' ## using part_csv:
#'
#' n.pix = 30865 # pixels in AK_ndvi_common-land.csv
#' parts = sample_partitions(npix = n.pix, npart = 4, partsize = 1000)
#'
#' data.file = system.file("extdata", "AK_ndvi_common-land.csv",
#'                         package = "remotePARTS")
#'
#' B.out = part_csv(part_i = 1, part_csv_path = data.file, part_mat = parts,
#'                  part_form = "cls.coef ~ 0 + land")
part_csv <- function(part_i, part_form, part_csv_path, part_mat, part_locvars = c("lng", "lat"), part_form0 = NULL){
  PartDF = data.frame(part = part_mat[, part_i]) #df with 1 column: current partition

  # use SQL syntax with read.csv.sql to only read specific rows (very fast)
  df_prt = sqldf::read.csv.sql(file = part_csv_path,
                             sql = "select file.* from file join PartDF on file.rowid = PartDF.part",
                             dbname = tempfile(),
                             header = TRUE)

  mf = stats::model.frame(formula(part_form), data = df_prt)
  resp = model.response(mf)
  X = stats::model.matrix(formula(part_form), data = df_prt)

  if(is.null(part_form0)){
    X0 = model.matrix(resp ~ 1)
  } else {
    X0 = model.matrix(formula(part_form0), data = df_prt)
  }

  return(list(y = as.vector(resp),
              X = as.matrix(X),
              X0 = as.matrix(X0),
              coords = as.matrix(df_prt[, part_locvars])))
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
#' @param part_f function to partition data. See details for more info
#' @param dist_f function to calculate distance. See details for more info
#' @param V.meth method passed to \code{fitV()}
#' @param spatcor spatial correlation used by \code{fitV()}
#' @param nug nugget, if NA (default), the nugget is estimated
#' @param partsize number of pixels in each partition
#' @param npart number of partitions
#' @param threads number of threads used by Eigen for matrix algebra
#' @param mincross number of partition pairs to calculate cross-partition
#' statistics from.
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
#' parts = sample_partitions(npix = n.pix, npart = 4, partsize = 200) # small partition matrix
#' data.file = system.file("extdata", "AK_ndvi_common-land.csv",
#'                         package = "remotePARTS")
#'
#' # Fit the partitioned GLS
#' GLS.part = fitGLS.partition(part_f = "part_csv", dist_f = "dist_km",
#'                             partsize = nrow(parts), npart = ncol(parts),
#'                             V.meth = "exponential", spatcor = .5,
#'                             # additional arguments passed to part_csv():
#'                             part_csv_path = data.file, part_mat = parts,
#'                             part_form = "cls.coef ~ 0 + land",
#'                             part_form0 = "cls.coef ~ 1")
fitGLS.partition <- function(part_f = "part_csv", dist_f = "dist_km",
                             V.meth = "exponential", spatcor, nug = NA,
                             partsize, npart, mincross = 6, threads = 1, ...){
  # Setup ----
  ## match the input functions
  func <- match.fun(part_f)
  D_func <- match.fun(dist_f)

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
      # stopifnot("part_f(i)$coords must be a matrix" = is.matrix(out.i$coords))

      ## calculate df
      dfs <- calc_dfpart(partsize, p = ncol(out.i$X), p0 = ncol(out.i$X0))
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
      if(is.na(nug)){
        gls.i <- GLS_worker(y = out.i$y, X = out.i$X, V = V.i, X0 = out.i$X0, save_xx = xx.print, threads = threads)
      } else {
        Xi = as.matrix(out.i$X)
        Vi = as.matrix(V.i)
        yi = as.matrix(out.i$y)
        X0i = as.matrix(out.i$X0)
        gls.i <- .Call(`_remotePARTS_fitGLS_cpp`, Xi, Vi, yi ,X0i ,nug, save_xx = xx.print, threads = threads)
        gls.i$nugget = nug
        gls.i$pval.t = sapply(gls.i$tstat, function(t){2*pt(abs(t), df=gls.i$dft)})
        gls.i$pval.F = stats::pf(gls.i$Fstat, gls.i$df.F[1], gls.i$df.F[2], lower.tail = FALSE)
        class(gls.i) <- append("remoteGLS", class(gls.i))
        attr(gls.i, "no_F") = FALSE
      }

      invchol_i <- invert_chol(V.i, gls.i$nugget)

    } else {
      ## copy out, V, and gls from previous j
      out.i <-  out.j
      V.i <- V.j
      gls.i <- gls.j
      invchol_i <- invchol_j
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
    if(is.na(nug)){
      gls.j <- GLS_worker(y = out.j$y, X = out.j$X, V = V.j, X0 = out.j$X0, save_xx = xx.print, threads = threads)
    } else {
      Xj = as.matrix(out.j$X)
      Vj = as.matrix(V.j)
      yj = as.matrix(out.j$y)
      X0j = as.matrix(out.j$X0)
      gls.j <- .Call(`_remotePARTS_fitGLS_cpp`, Xj, Vj, yj, X0j, nug, save_xx = xx.print, threads = threads)
      gls.j$nugget = nug
      gls.j$pval.t = sapply(gls.j$tstat, function(t){2*pt(abs(t), df=gls.j$dft)})
      gls.j$pval.F = stats::pf(gls.j$Fstat, gls.j$df.F[1], gls.j$df.F[2], lower.tail = FALSE)
      class(gls.j) <- append("remoteGLS", class(gls.j))
      attr(gls.j, "no_F") = FALSE
    }

    invchol_j <- invert_chol(V.j, gls.j$nugget)

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
                                  invChol_i = invchol_i,
                                  invChol_j = invchol_j,
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
  class(out.list) <- append("remoteGLS.parts", class(out.list))
  return(out.list)
}

#' @rdname fitGLS.partition
#'
#' @param ncores number of cores for parallel processing. Default is total cores - 1
#' @param export an optional character vector of names for any additional objects
#' needed for \code{part_f()}.
#' @param debug logical debug flag. If TRUE, prints the name of the step that is running
#'
#' @export
#'
#' @examples
#' ## now w. 2 cores:
#'
#' if(FALSE){ # change to TRUE to run multi-core version
#'   GLS.part.mc = fitGLS.partition.mc(part_f = "part_csv", dist_f = "dist_km",
#'                                     partsize = nrow(parts), npart = ncol(parts),
#'                                     V.meth = "exponential", spatcor = .5,
#'                                     part_csv_path = data.file, part_mat = parts,
#'                                     part_form = "cls.coef ~ 0 + land",
#'                                     part_form0 = "cls.coef ~ 1",
#'                                     ncores = 4)
#' }
fitGLS.partition.mc <- function(part_f = "part_csv", dist_f = "dist_km",
                                V.meth = "exponential", spatcor,
                                partsize, npart, mincross = 6,
                                ncores = parallel::detectCores() - 1,
                                export = NA, debug = TRUE,
                                ...){
  # Setup ----
  if(debug){print("1. setup")}
  ## setup cluster
  clst = parallel::makeCluster(ncores)
  ## After the function is run, close the cluster
  on.exit(parallel::stopCluster(clst))
  ## Register parallel backend
  doParallel::registerDoParallel(clst)

  ## match the input functions
  func <- match.fun(part_f)
  D_func <- match.fun(dist_f)

  ## GLS for each partition ----
  if(debug){print("2. part GLS")}
  part_out = foreach::foreach(i = 1:npart, .packages = "remotePARTS") %dopar% {
    out.i = func(i, ...) # get partition data
    # out.i = func(i, part_form, part_csv_path, part_mat, part_locvars, part_form0)
    if(i == 1){
      # checks
      stopifnot("part_f(i)$X must be a matrix" = is.matrix(out.i$X))
      stopifnot("part_f(i)$coords must be a matrix" = is.matrix(out.i$coords))
    }
    # out.i = func(i, csv.path = data.file, part.mat = parts)
    D.i <- D_func(out.i$coords) # calculate D
    V.i <- fitV(D.i, spatcor, V.meth) # calculate V
    gls.i <- GLS_worker(y = out.i$y, X = out.i$X, V = V.i, X0 = out.i$X0, save_xx = TRUE)
    # "return" statement

    list(data = out.i, GLS = gls.i)
  }

  ## calculate df
  dfs <- calc_dfpart(partsize, p = ncol(part_out[[1]]$data$X), p0 = ncol(part_out[[1]]$data$X0))


  ## Cross-partition setup ----
  ## all possible pairwise combinations
  mincross = ifelse(mincross > npart, npart, mincross)
  used.combs = t(utils::combn(mincross, 2))
  ncrosses = nrow(used.combs)
  # maxcross = nrow(combs)
  # ## adjust mincross if it is greater than maxcross
  # if (mincross > maxcross){
  #   ### NEED TO FIX THIS !!!!!
  #   mincross <- maxcross
  # }
  #
  # ## Calculate which combinations to use
  # if (mincross < maxcross){
  #   used.combs = combs[sample(nrow(combs), mincross), ]
  # } else {
  #   used.combs = combs[sample(maxcross, mincross), ]
  # }
  # used.combs = used.combs[order(used.combs[, 1]), ]

  ## Cross-partition GLS ----
  if(debug){print("3. cross-partition GLS")}
  cross = NULL
  cross_out = foreach::foreach(cross = iterators::iter(used.combs, by = "row"),
                               .packages = "remotePARTS") %dopar% {
                                 i = cross[1]; j = cross[2]

                                 ## recompute cholesky from distances
                                 Di <- D_func(part_out[[i]]$data$coords)
                                 Vi <- fitV(Di, spatcor, V.meth) # calculate V
                                 invchol_i = invert_chol(Vi, part_out[[i]]$GLS$nugget)

                                 Dj <- D_func(part_out[[j]]$data$coords)
                                 Vj <- fitV(Dj, spatcor, V.meth) # calculate V
                                 invchol_j = invert_chol(Vj, part_out[[j]]$GLS$nugget)

                                 ## cross-partition varcovar matrix
                                 D.ij <- D_func(part_out[[i]]$data$coords, part_out[[j]]$data$coords)
                                 V.ij <- fitV(D.ij, spatcor, V.meth)
                                 cross.ij = crosspart_worker(xxi = part_out[[i]]$GLS$xx, xxj = part_out[[j]]$GLS$xx,
                                                             xxi0 = part_out[[i]]$GLS$xx0, xxj0 = part_out[[j]]$GLS$xx0,
                                                             invChol_i = invchol_i,
                                                             invChol_j = invchol_j,
                                                             Vsub = V.ij,
                                                             nug_i = part_out[[i]]$GLS$nugget, nug_j = part_out[[j]]$GLS$nugget,
                                                             df1 = dfs[1], df2 = dfs[2])
                                 # Only return necessary output to save space
                                 list(rSSRij = cross.ij$rSSRij,
                                      rSSEij = cross.ij$rSSEij,
                                      rcoefij = cross.ij$rcoefij)
                               }


  # Results Collection ----
  if(debug){print("4. Results Collection")}
  f1 = part_out[[1]]$data
  ## setup output
  betas = matrix(NA, ncol = ncol(f1$X), nrow = npart,
                 dimnames = list(NULL, colnames(f1$X)))
  SEs = betas
  t.stats = betas
  pvals.t = betas
  mod.stats = matrix(NA, nrow = npart, ncol = 7,
                     dimnames = list(NULL, c("nugget","logLik", "SSE", "MSE",
                                             "MSR", "F.stat", "pval.F")))
  cross.SS = matrix(NA, nrow = ncrosses, ncol = 2,
                    dimnames = list(NULL, c("rSSR", "rSSE")))
  rcoefs = matrix(NA, nrow = ncrosses, ncol = ncol(f1$X),
                  dimnames = list(NULL, colnames(f1$X)))
  ## fill in the output
  # Partition stats
  for(i in 1:npart){
    betas[i, ] <- part_out[[i]]$GLS$betahat
    SEs[i, ] <- part_out[[i]]$GLS$SE
    t.stats[i, ] <- part_out[[i]]$GLS$tstat
    pvals.t[i, ] <- part_out[[i]]$GLS$pval.t
    mod.stats[i, ] <- c(part_out[[i]]$GLS$nugget, part_out[[i]]$GLS$logLik,
                        part_out[[i]]$GLS$SSE, part_out[[i]]$GLS$MSE,
                        part_out[[i]]$GLS$MSR,
                        part_out[[i]]$GLS$Fstat, part_out[[i]]$GLS$pval.F)
  }
  # Cross-partition stats
  for(ij in 1:length(cross_out)){
    cross.SS[ij, ] = c(cross_out[[ij]]$rSSRij, cross_out[[ij]]$rSSEij)
    rcoefs[ij, ] = as.vector(cross_out[[ij]]$rcoefij)
  }
  # Overall Statistics
  fmean = mean(mod.stats[, "F.stat"])
  rSSRmean = mean(cross.SS[, "rSSR"])
  rSSEmean = mean(cross.SS[, "rSSE"])
  rcoefmean = colMeans(rcoefs)
  coefmean = colMeans(betas)

  # Output ----
  if(debug){print("5. Output")}
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
  class(out.list) <- append("remoteGLS.parts", class(out.list))
  return(out.list)
}
