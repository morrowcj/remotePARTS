## fitGLS_partition ----
#' @title fit a partitioned GLS
#'
#' @description fit a GLS model to a large data set by partitioning the data
#' into smaller pieces (partitions) and processing these pieces individually and
#' summarizing output across partitions to conduct hypothesis tests.
#'
#' @rdname partGLS
#'
#' @family partitionedGLS
#'
#' @param formula a formula, passed to \code{part_FUN} and \code{fitGLS}
#' @param partmat a numeric partition matrix, with values containing indices of locations.
#' @param formula0 a formula for the null model, passed to \code{part_FUN} and \code{fitGLS}.
#' @param part_FUN a function to partition individual data. See details for more
#' information about requirements for this function.
#' @param distm_FUN a function to calculate distances from a coordinate matrix
#' @param covar_FUN a function to calculate covariances from a distance matrix
#' @param covar.pars a named list of parameters passed to \code{covar_FUN}
#' @param nugget a numeric fixed nugget component: if NA, the nugget is estimated for
#' each partition
#' @param ncross an integer indicating the number of partitions used to calculate
#' cross-partition statistics
#' @param save.GLS logical: should full GLS output be saved for each partition?
#' @param do.t.test logical: should a t-test of the GLS coefficients be conducted?
#' @param do.chisqr.test logical: should a correlated chi-squared test of the model
#' fit be conducted?
#' @param progressbar logical: should progress be tracked with a progress bar?
#' @param debug logical debug mode
#' @param ncores an optional integer indicating how many CPU threads to use for calculations.
#' @param parallel logical: should all calculations be done in parallel? See details for
#' more information
#' @param ... arguments passed to \code{part_FUN}
#'
#' @details
#'
#' The function specified by \code{part_FUN} is called internally to obtain
#' properly formatted subsets of the full data (i.e., partitions). Two functions
#' are provided in the \code{remotePARTs} package for this purpose: \code{part_data}
#' and \code{part_csv}. Both of these functions have required arguments that
#' must be specified through the call to \code{fitGLS_partition} (via \code{...}).
#' Check each function's argument list and see "\code{part_FUN} details" below
#' for more information.
#'
#' \code{partmat} is used to partition the data. \code{partmat} must be a complete
#' matrix, without any missing or non-finite values. Columns of \code{partmat} are
#' passed as the first argument \code{part_FUN} to obtain data, which is then
#' passed to \code{fitGLS}. Users are encouraged to use \code{sample_partitions()}
#' to obtain a valid \code{partmat}.
#'
#' The specific dimensions of \code{partmat} can have a substantial effect on the
#' efficiency of \code{fitGLS_partition}. For most systems, we do not recommend
#' fitting with partitions exceeding 3000 locations or pixels
#' (i.e., \code{partmat(partsize = 3000, ...)}). Any larger, and the covariance
#' matrix inversions may become quite slow (or impossible for some machines).
#' It may help performance to use smaller even partitions of around 1000-2000
#' locations.
#'
#' \code{ncross} determines how many partitions are used to estimate cross-partition
#' statistics. All partitions, up to \code{ncross} are compared with all others
#' in a pairwise fashion. There is no hard rule for setting \code{mincross}. More
#' crosses will ensure convergence, but we believe that the default of 6
#' (10 total comparisons) should be sufficient for most moderate-sized maps
#' if 1500-3000 pixel partitions are used. This may require testing with each
#' individual dataset to determine at what point convergence occurs.
#'
#' Covariance matrices for each partition are calculated with \code{covar_FUN}
#' from distances among points within the partition. Parameter values for
#' \code{covar_FUN} are given by \code{covar.pars}.
#'
#' The distances among points are calculated with \code{distm_FUN}.
#' \code{distm_FUN} can be any function, modeled after \code{geosphere::distm()},
#' that satisfies both: 1) returns a distance matrix among points when a single
#' coordinate matrix is given as first argument; and 2) returns a matrix
#' containing distances between two coordinate matrices if given as the first and
#' second arguments.
#'
#' If \code{nugget = NA}, a ML nugget is obtained for each partition. Otherwise,
#' a fixed nugget is used for all partitions.
#'
#' It is not required to use all partitions for cross-partition calculations, nor
#' is it recommended to do so for most large data sets.
#'
#' If \code{progressbar = TRUE} a text progress bar shows the current status
#' of the calculations in the console.
#'
#' @section parallel implementation:
#'
#' In order to be efficient and account for different user situations, parallel
#' processing is available natively in \code{fitGLS_partition}. There are a few
#' different specifications that will result in different behavior:
#'
#' When \code{parallel = TRUE} and \code{ncores > 1}, all calculations are done
#' completely in parallel (via \code{multicore_fitGLS_partition()}).
#' In this case, parallelization is implemented with the
#' \code{parallel}, \code{doParallel}, and \code{foreach} packages. In this version,
#' all matrix operations are serialized on each worker but multiple operations
#' can occur simultaneously..
#'
#' When \code{parallel = FALSE} and \code{ncores > 1}, then most calculations
#' are done on a single core but matrix opperations use multiple cores. In this
#' case, \code{ncores} is passed to fitGLS. In this option, it is suggested
#' to not exceed the number of physical cores (not threads).
#'
#' When \code{ncores <= 1}, then the calculations are completely serialized
#'
#' When \code{ncores = NA} (the default), only one core is used.
#'
#' In the parallel implementation of this function, a progress bar is not possible,
#' so \code{progressbar} is ignored.
#'
#' @section \code{part_FUN} details:
#'
#' \code{part_FUN} can be any function that satisfies the following criteria
#'
#' 1. the first argument of \code{part_FUN} must accept an index of pixels by which
#' to subset the data;
#'
#' 2. \code{part_FUN} must also accept \code{formula} and \code{formula0} from
#' \code{fitGLS_partition}; and
#'
#' 3. the output of \code{part_FUN} must be a list with at least the
#' following elements, which are passed to \code{fitGLS};
#'
#' \describe{
#'     \item{data}{a data frame containing all variables given by \code{formula}.
#'     Rows should correspond to pixels specified by the first argument}
#'     \item{coords}{a coordinate matrix or data frame. Rows should correspond to
#'     pixels specified by the first argument}
#' }
#'
#' Two functions that satisfy these criteria are provided by \code{remotePARTS}:
#' \code{part_data} and \code{part_csv}.
#'
#' \code{part_data} uses an in-memory data frame (\code{data})
#' as a data source. \code{part_csv}, instead reads data from a
#' csv file (\code{file}), one partition at a time, for efficient memory usage.
#' \code{part_csv} internally calls \code{sqldf::read.csv.sql()} for fast and
#' efficient row extraction.
#'
#' Both functions use \code{index} to subset rows of data and \code{formula} and
#' \code{formula0} (optional) to determine which variables to select.
#'
#' Both functions also use \code{coord.names} to indicate which variables contain
#' spatial coordinates. The name of the x-coordinate column should always preceed
#' the y-coordinate column: \code{c("x", "y")}.
#'
#' Users are encouraged to write their own \code{part_FUN} functions to meet their
#' needs. For example, one might be interested in using data stored in a raster
#' stack or any other file type. In this case, a user-defined \code{part_FUN}
#' function allows access to \code{fitGLS_partition} without saving reformatted
#' copies of data.
#'
#' @return \code{fitGLS_partition} returns a list object of class "partGLS" which
#' contains at least the following elements:
#'
#' \describe{
#'     \item{call}{the function call}
#'     \item{GLS}{an optional list of "remoteGLS" objects, one for each partition}
#'     \item{part}{statistics calculated from each partition: see below for further
#'     details}
#'     \item{cross}{statistics calculated from each pair of crossed partitions,
#'     determined by \code{ncross}: see below for further details}
#'     \item{overall}{summary statistics of the overall model: see below for further
#'     details}
#' }
#'
#' \code{part} is a sub-list containing the following elements
#'
#' \describe{
#'     \item{coefficients}{a numeric matrix of GLS coefficients for each partition}
#'     \item{SEs}{a numeric matrix of coefficient standard errors}
#'     \item{tstats}{a numeric matrix of coefficient t-statstitics}
#'     \item{pvals_t}{a numeric matrix of t-test pvalues}
#'     \item{nuggets}{a numeric vector of nuggets for each partition}
#'     \item{covar.pars}{\code{covar.pars} input vector}
#'     \item{modstats}{a numeric matrix with rows corresponding to partitions and
#'     columns corresponding to log-likelihoods (\code{logLik}),
#'     sum of square error (\code{SSE}), mean-squared error (\code{MSE}),
#'     regression mean-square (\code{MSR}), F-statistics (\code{Fstat}),
#'     and p-values from F-tests (\code{pval_F})}
#' }
#'
#' \code{cross} is a sub-list containing the following elements, which are use
#' to calculate the combined (across partitions) standard errors of the coefficient
#' estimates and statistical tests. See Ives et al. (2022).
#'
#' \describe{
#'     \item{rcoefs}{a numeric matrix of cross-partition correlations in the
#'     estimates of the coefficients}
#'     \item{rSSRs}{a numeric vector of cross-partition correlations in the
#'     regression sum of squares}
#'     \item{rSSEs}{a numeric vector of cross-partition correlations in the
#'     sum of squared errors}
#' }
#'
#' and \code{overall} is a sub-list containing the elements
#'
#' \describe{
#'     \item{coefficients}{a numeric vector of the average coefficient estimates
#'     across all partitions}
#'     \item{rcoefficients}{a numeric vector of the average cross-partition
#'     coefficient from across all crosses}
#'     \item{rSSR}{the average cross-partition correlation in the regression
#'     sum of squares}
#'     \item{rSSE}{the average cross-partition correlation in the sum of
#'     squared errors}
#'     \item{Fstat}{the average f-statistic across partitions}
#'     \item{dfs}{degrees of freedom to be used with partitioned GLS f-test}
#'     \item{partdims}{dimensions of \code{partmat}}
#'     \item{pval.chisqr}{if \code{chisqr.test = TRUE}, a p-value for the correlated
#'     chi-squared test}
#'     \item{t.test}{if \code{do.t.test = TRUE}, a table with t-test results, including
#'     the coefficient estimates, standard errors, t-statistics, and p-values}
#' }
#'
#' @references
#' Ives, A. R., L. Zhu, F. Wang, J. Zhu, C. J. Morrow, and V. C. Radeloff. in review.
#'     Statistical tests for non-independent partitions of large autocorrelated datasets.
#'     MethodsX.
#'
#' @examples
#' ## read data
#' data(ndvi_AK10000)
#' df = ndvi_AK10000[seq_len(1000), ] # first 1000 rows
#'
#' ## create partition matrix
#' pm = sample_partitions(nrow(df), npart = 3)
#'
#' ## fit GLS with fixed nugget
#' partGLS = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm,
#'                            data = df, nugget = 0, do.t.test = TRUE)
#'
#' ## hypothesis tests
#' chisqr(partGLS) # explanatory power of model
#' t.test(partGLS) # significance of predictors
#'
#' ## now with a numeric predictor
#' fitGLS_partition(formula = CLS_coef ~ lat, partmat = pm, data = df, nugget = 0)
#'
#' \donttest{
#' ## fit ML nugget for each partition (slow)
#' (partGLS.opt = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm,
#'                                 data = df, nugget = NA))
#' partGLS.opt$part$nuggets # ML nuggets
#'
#' # Certain model structures may not be useful:
#' ## 0 intercept with numeric predictor (produces NAs) and gives a warning in statistical tests
#' fitGLS_partition(formula = CLS_coef ~ 0 + lat, partmat = pm, data = df, nugget = 0)
#'
#' ## intercept-only, gives warning
#' fitGLS_partition(formula = CLS_coef ~ 1, partmat = pm, data = df, nugget = 0,
#'                  do.chisqr.test = FALSE)
#' }
#' @export
fitGLS_partition <- function(formula, partmat, formula0 = NULL,
                             part_FUN = "part_data",
                             distm_FUN = "distm_scaled", covar_FUN = "covar_exp",
                             covar.pars = c(range = .1), nugget = NA, ncross = 6,
                             save.GLS = FALSE, do.t.test = TRUE, do.chisqr.test = TRUE,
                             progressbar = TRUE, debug = FALSE, ncores = NA, parallel = TRUE,
                             ...){
  # Setup
  call = match.call()
  ## partmat dimensions
  partsize = nrow(partmat)
  npart = ncol(partmat)
  ## adjust number of partitions to cross
  ncross = ifelse(ncross > npart, npart, ncross) # can't be larger than npart
  cross.pairs = t(combn(seq_len(ncross), 2))
  npairs = nrow(cross.pairs)
  ## match functions
  part.f <- match.fun(part_FUN)
  dist.f <- match.fun(distm_FUN)
  covar.f <- match.fun(covar_FUN)
  if(is.null(formula0)){
    formula0 = update(as.formula(formula), . ~ 1)
  } else {
    formula0 = as.formula(formula0)
  }

  # Check for problematic or invalid formula objects ====
  ## formula are identical
  if(formula == formula0){
    warning("formula and formula0 are identical, which prevents model comparison calulations.",
            "\ndo.chisqr.test set to FALSE")
    do.chisqr.test = FALSE
  } else {
  }
  ## no right-hand side (i.e., formula(x ~ 0) or update(form, . ~ 0))
  if(formula[-2]=="~0" | formula0[-2]=="~0" | formula[-2]=="~1 - 1" | formula0[-2]=="~1 - 1"){
    stop("invalid formula: the right hand side may not be empty ('~0')")
  }
  # ====

  ## Decide if calculations should be parallelized
  if(!is.na(ncores) && parallel & ncores > 1){
    if(debug){cat("Conducting parallel paritioned GLS\n")}
    outlist = multicore_fitGLS_partition(formula = formula, partmat = partmat,
                                         formula0 = formula0, part_FUN = part_FUN,
                                         distm_FUN = distm_FUN, covar_FUN = covar_FUN,
                                         covar.pars = covar.pars, nugget = nugget,
                                         ncross = ncross, save.GLS = save.GLS,
                                         do.chisqr.test = do.chisqr.test, ncores = ncores,
                                         ...)
    outlist$call <- call # update call
  } else {
    if(debug){cat("Conducting partitioned GLS\n")}
    if(is.na(ncores)){
      ncores = 1L
    } else {
      ncores = as.integer(ncores)
    }

    # output setup
    partGLS = vector("list", npart)
    if(debug){crosspartGLS = vector("list", npairs)}
    if(progressbar){
      pb = txtProgressBar(min = 0, max = npart, style = 3)
    }
    for(i in 1:npart){
      # for (i in 1){
      if (debug) {cat("i =", i, "\n")}
      ## partition data
      idat <- part.f(partmat[, i], formula = formula, formula0 = formula0, ...)

      ## Calculate GLS, if not already done
      if (is.null(partGLS[[i]])){
        ## covariance of parition
        Vi = do.call(covar.f, args = append(list(d = dist.f(idat$coords)), as.list(covar.pars)))
        ## GLS of parition
        partGLS[[i]] <- fitGLS(formula = formula, data = idat$data, V = Vi,
                               nugget = nugget, formula0 = formula0, save.xx = (i <= ncross),
                               no.F = FALSE, save.invchol = (i <= ncross), logLik.only= FALSE,
                               ncores = ncores, suppress_compare_warning = TRUE)
        ## build some empty stat tables on first loop
        if (i == 1){
          p = ncol(partGLS[[1]]$xx)
          p0 = ncol(partGLS[[1]]$xx0)
          ## part stats
          coefs = SEs = tstats = tpvals =
            matrix(NA, nrow = npart, ncol = p,
                   dimnames = list(NULL, names(partGLS[[1]]$coefficients)))
          nuggets = LLs = SSEs = MSEs = MSRs = Fstats = Fpvals = rep(NA, times = npart)
          covar_coefs = array(NA, dim = c(p, p, npart),
                   dimnames = list(names(partGLS[[1]]$coefficients),
                                   names(partGLS[[1]]$coefficients),
                                   NULL))
          ## cross stats
          rSSRs = rSSEs = rep(NA, npairs)
          rcoefs = array(NA, dim = c(npairs, p, p),
                         dimnames = list(NULL,
                                         names(partGLS[[1]]$coefficients),
                                         names(partGLS[[1]]$coefficients)))
        }
        ## collect some stats
        coefs[i, ] = partGLS[[i]]$coefficients
        SEs[i, ] = partGLS[[i]]$SE
        covar_coefs[ , , i] = partGLS[[i]]$covar_coef
        tstats[i, ] = partGLS[[i]]$tstat
        tpvals[i, ] = partGLS[[i]]$pval_t
        nuggets[i] = partGLS[[i]]$nugget
        LLs[i] = partGLS[[i]]$logLik
        SSEs[i] = partGLS[[i]]$SSE
        MSEs[i] = partGLS[[i]]$MSE
        MSRs[i] = ifelse(is.null(partGLS[[i]]$MSR) | length(partGLS[[i]]$MSR) == 0, NA, partGLS[[i]]$MSR)
        Fstats[i] = partGLS[[i]]$Fstat
        Fpvals [i] = partGLS[[i]]$pval_F
      }
      ## calculate inverse cholesky it needed
      if (is.null(partGLS[[i]]$invcholV)){
        Vi = do.call(covar.f, args = append(list(d = dist.f(idat$coords)), as.list(covar.pars)))
        partGLS[[i]]$invcholV <- invert_chol(Vi)
      }
      if (i < ncross) for (j in (i+1):ncross) {
        if (debug) {cat("j =", j, "\n")}
        ## parition data
        jdat = part.f(partmat[, j], formula = formula, formula0 = formula0, ...)

        ## calculate GLS, if not already done
        if (is.null(partGLS[[j]])){
          Vj = do.call(covar.f, args = append(list(d = dist.f(jdat$coords)), as.list(covar.pars)))
          partGLS[[j]] <- fitGLS(formula = formula, data = jdat$data, V = Vj,
                                 nugget = nugget, formula0 = formula0, save.xx = TRUE,
                                 no.F = FALSE, save.invchol = TRUE, logLik.only= FALSE,
                                 ncores = ncores, suppress_compare_warning = TRUE)
          if(length(partGLS[[j]]$coefficients) != ncol(coefs)){
            stop("dimension mismatch: different number of coefficients between parts i and j. Filter data or try different partition matrix.")
          }
          ## collect some stats
          coefs[j, ] = partGLS[[j]]$coefficients
          SEs[j, ] = partGLS[[j]]$SE
          covar_coefs[ , , j] = partGLS[[j]]$covar_coef
          tstats[j, ] = partGLS[[j]]$tstat
          tpvals[j, ] = partGLS[[j]]$pval_t
          nuggets[j] = partGLS[[j]]$nugget
          LLs[j] = partGLS[[j]]$logLik
          SSEs[j] = partGLS[[j]]$SSE
          MSEs[j] = partGLS[[j]]$MSE
          MSRs[j] = partGLS[[j]]$MSR
          Fstats[j] = partGLS[[j]]$Fstat
          Fpvals[j] = partGLS[[j]]$pval_F
        }
        ## which cross are we on?
        cross = which( (cross.pairs[,1] == i) & (cross.pairs[, 2] == j) )
        if (debug) {cat( "cross #", cross, "\n")}
        ## calculate inverse cholesky if needed
        if (is.null(partGLS[[j]]$invcholV)){
          Vj = do.call(covar.f, args = append(list(d = dist.f(jdat$coords)), as.list(covar.pars)))
          partGLS[[j]]$invcholV <- invert_chol(Vj)
        }
        ## cross covariance
        if (debug) {cat("calulating Vij:\n")}
        Vij = do.call(covar.f, args = append(list(d = dist.f(idat$coords, jdat$coords)),
                                             as.list(covar.pars)))
        ## degrees of freedom
        dfs = calc_dfpart(partsize = partsize, p = ncol(partGLS[[j]]$xx), p0 = ncol(partGLS[[j]]$xx0))
        ## calculate cross-partition stats
        if (debug) {cat("calculating crossGLS\n")}
        rGLS <- crosspart_GLS(xxi = partGLS[[i]]$xx,
                              xxj = partGLS[[j]]$xx,
                              xxi0 = partGLS[[i]]$xx0,
                              xxj0 = partGLS[[j]]$xx0,
                              invChol_i = partGLS[[i]]$invcholV,
                              invChol_j = partGLS[[j]]$invcholV,
                              Vsub = Vij,
                              nug_i = partGLS[[i]]$nugget,
                              nug_j = partGLS[[j]]$nugget,
                              df1 = dfs[1], df2 = dfs[2],
                              small = FALSE, # return Vcoefij for use in correlations
                              ncores = ncores)

        ## delete large matrix for j
        partGLS[[j]]$invcholV = NULL
        if(debug){crosspartGLS[[cross]] = rGLS}
        ## collect stats
        rcoefs[cross, ,] <- rGLS$rcoefij
        # rcoefs[cross, ] <- as.vector(rGLS$rcoefij)
        rSSRs[cross] <- ifelse(is.na(rGLS$rSSRij) | is.infinite(rGLS$rSSRij),
                               NA,
                               rGLS$rSSRij)
        rSSEs[cross] <- rGLS$rSSEij
      }
      ## delete large matrix for i
      partGLS[[i]]$invcholV = NULL
      ## update progress bar
      if(progressbar){
        setTxtProgressBar(pb, i)
      }
    }

    ## make the coefficients
    rcoefficients = apply(rcoefs, MARGIN=c(2,3), FUN = function(x){mean(x, na.rm = TRUE)})
    # rcoefficients = apply(rcoefs, MARGIN = 2, FUN = function(x){mean(x, na.rm = TRUE)})

    ## collect and format output
    outlist = list(call = call,
                   GLS = if(save.GLS){partGLS}else{NULL},
                   part = list(coefficients = coefs, SEs = SEs,
                               covar_coefs = covar_coefs, tstats = tstats,
                               pvals_t = tpvals, nuggets = nuggets,
                               covar.pars = covar.pars,
                               modstats = cbind(LLs = LLs, SSEs = SSEs,
                                                MSEs = MSEs, MSRs = MSRs,
                                                Fstats = Fstats,
                                                pvals_F = Fpvals)),
                   cross = list(rcoefs = rcoefs, rSSRs = rSSRs, rSSEs = rSSEs),
                   overall = list(coefficients = colMeans(coefs, na.rm = TRUE),
                                  # rcoefficients = colMeans(rcoefs, na.rm = TRUE),
                                  rcoefficients =rcoefficients,
                                  rSSR = mean(rSSRs, na.rm = TRUE),
                                  rSSE = mean(rSSEs, na.rm = TRUE),
                                  Fstat = mean(Fstats, na.rm = TRUE),
                                  dfs = calc_dfpart(partsize, p, p0),
                                  partdims = c(npart = npart, partsize = partsize)))
    if(debug){outlist$crossGLS = crosspartGLS}
    class(outlist) <- append("partGLS", class(outlist))
    if(do.chisqr.test){
      if(as.formula(formula) == as.formula(formula0)){
        warning("chisqr test not valid when formula == formula0")
      } else {
      outlist$overall$pval.chisqr = tryCatch(chisqr(outlist),
                                             error = function(e){warning("error in chisqr()")},
                                             warning = function(w){warning("warning in chisqr()")})
      }
    }
    if(do.t.test){
      test.output = tryCatch(t.test(outlist),
                        error = function(e){warning("error in t.test()")},
                        warning = function(w){warning("warning in t.test()")})
      outlist$overall$t.test = test.output$p.t
      outlist$overall$covar_coef = test.output$covar_coef
    }
    close(pb)
  }
  ## final return statement
  return(outlist)
}

## calc_dfpart ----
#' calculate degrees of freedom for partitioned GLS
#'
#' @param partsize number of pixels in each partition
#' @param p number of predictors in alternate model
#' @param p0 number of parameters in null model
#'
#' @returns a named vector containing the first and second degrees of freedom ("df1" and "df2", respectively)
#'
calc_dfpart <- function(partsize, p, p0){
  stopifnot(length(partsize) == 1)
  df2 = partsize - (p - 1)
  df0 = partsize - (p0 - 1)
  df1 = df0 - df2
  return(c("df1" = df1, "df2" = df2))
}

## cross-partition worker ----
#' @title Calculate cross-partition statistics in a partitioned GLS
#'
#' @description Calculate cross-partition statistics between two GLS partitions
#'
#' @family partitionedGLS
#'
#' @param xxi numeric matrix xx from  partition i
#' @param xxj numeric matrix xx from  partition j
#' @param xxi0 numeric matrix xx0 from  partition i
#' @param xxj0 numeric matrix xx0 from  partition j
#' @param invChol_i numeric matrix invcholV from  partition i
#' @param invChol_j numeric matrix invcholV from  partition j
#' @param Vsub numeric variance matrix for Xij (upper block)
#' @param nug_i nugget from partition i
#' @param nug_j nugget from partition j
#' @param df1 first degree of freedom
#' @param df2 second degree of freedom
#' @param small logical: if \code{TRUE}, only return \code{rcoefij}, \code{rSSRij},
#' and \code{rSSEij}
#' @param ncores an optional integer indicating how many CPU threads to use for
#' matrix calculations.
#'
#' @return
#' \code{crosspart_GLS} returns a list of cross-partition statistics.
#'
#' If \code{small = FALSE}, the list contains the following elements
#'
#' \describe{
#'     \item{Rij}{}
#'     \item{Hi}{}
#'     \item{Hj}{}
#'     \item{Hi0}{}
#'     \item{Hj0}{}
#'     \item{SiR}{}
#'     \item{SjR}{}
#'     \item{rcoefij}{}
#'     \item{rSSRij}{}
#'     \item{rSSEij}{}
#'     \item{Vcoefij}{}
#' }
#'
#' If \code{small = FALSE}, the list only contains the necessary elements
#' \code{rcoefij}, \code{rSSRij}, and \code{rSSEij}.
#'
crosspart_GLS <- function(xxi, xxj, xxi0, xxj0, invChol_i, invChol_j, Vsub,
                          nug_i, nug_j, df1, df2, small = TRUE, ncores = NA){
  if(is.na(ncores)){
    ncores = 1L
  } else {
    ncores = as.integer(ncores)
  }
  # coerce input to matrices
  xxi = as.matrix(xxi)
  xxj = as.matrix(xxj)
  xxi0 = as.matrix(xxi0)
  xxj0 = as.matrix(xxj0)
  invChol_i = as.matrix(invChol_i)
  invChol_j = as.matrix(invChol_j)
  Vsub = as.matrix(Vsub)

  ## error handling
  stopifnot(all(is.double(xxi), is.double(xxj),
                is.double(xxi0), is.double(xxj0),
                is.double(invChol_i), is.double(invChol_j),
                is.double(Vsub)))
  stopifnot(all.equal(nrow(xxi), nrow(xxj), nrow(xxi0), nrow(xxj0)))
  stopifnot(all.equal(ncol(xxi), ncol(xxj)))
  stopifnot(all.equal(ncol(xxi0), ncol(xxj0)))
  # stopifnot(all(check_posdef(V)))

  vcoef = ifelse(small, FALSE, TRUE)

  outlist <- .Call(`_remotePARTS_crosspart_worker_cpp`, xxi, xxj, xxi0, xxj0,
                   invChol_i, invChol_j, Vsub, nug_i, nug_j,df1, df2, vcoef, ncores)
  if(small){
    return(list(rcoefij = outlist$rcoefij,
                rSSRij = outlist$rSSRij,
                rSSEij = outlist$rSSEij))
  } else {
    return(outlist)
  }
}

## part_data ----
#' @rdname partGLS
#' @family partitionedGLS
#'
#' @param index a vector of pixels with which to subset the data
#' @param formula a formula for the GLS model
#' @param data a data frame
#' @param formula0 an optional formula for the null GLS model
#' @param coord.names a vector containing names of spatial coordinate variables
#' (x and y, respectively)
#'
#' @return \code{part_data} and \code{part_csv} both return a list with two elements:
#'
#' \describe{
#'     \item{data}{a dataframe, containing the data subset}
#'     \item{coords}{a coordinate matrix for the subset}
#' }
#'
#'
#' @examples
#' ## part_data examples
#' part_data(1:20, CLS_coef ~ 0 + land, data = ndvi_AK10000)
#'
#' @export
part_data <- function(index, formula, data, formula0 = NULL, coord.names = c("lng", "lat")){
  if(missing(data)){
    stop("missing argument 'data'. 'data' required for function part_data()")
  }
  stopifnot(is.data.frame(data) | (is.matrix(data) & is.numeric(data)))

  if(is.null(formula0)){
    formula0 = update(as.formula(formula), . ~ 1)
  } else {
    formula0 = as.formula(formula0)
  }

  cols = unique(rownames(attr(terms.formula(as.formula(formula)), "factors")),
                rownames(attr(terms.formula(as.formula(formula0)), "factors")))

  if(length(cols) < 1){
    cols = as.character(as.formula(formula)[[2]])
  }

  df = as.data.frame(data[index, cols])
  names(df) = cols

  return(list(data = df,
              coords = data[index, coord.names]))
}

## part_csv ----
#' @rdname partGLS
#' @family partitionedGLS
#'
#' @param file a text string indicating the csv file from which to read data
#'
#' @examples
#' \donttest{
#' ## part_csv examples - ## CAUTION: examples for part_csv() include manipulation side-effects:
#' # first, create a .csv file from ndviAK
#' data(ndvi_AK10000)
#' file.path = file.path(tempdir(), "ndviAK10000-remotePARTS.csv")
#' write.csv(ndvi_AK10000, file = file.path)
#'
#' # build a partition from the first 30 pixels in the file
#' part_csv(1:20, formula = CLS_coef ~ 0 + land, file = file.path)
#'
#' # now with a random 20 pixels
#' part_csv(sample(3000, 20), formula = CLS_coef ~ 0 + land, file = file.path)
#'
#' # remove the example csv file from disk
#' file.remove(file.path)
#' }
#'
#' @export
part_csv <- function(index, formula, file, formula0 = NULL, coord.names = c("lng", "lat")){

  if(is.null(formula0)){
    formula0 = update(as.formula(formula), . ~ 1)
  } else {
    formula0 = as.formula(formula0)
  }

  cols = unique(rownames(attr(terms.formula(as.formula(formula)), "factors")),
                rownames(attr(terms.formula(as.formula(formula0)), "factors")))

  if(length(cols) < 1){
    cols = as.character(as.formula(formula)[[2]])
  }

  PartDF = data.frame(part = index)
  df = sqldf::read.csv.sql(file = file,
                           sql = "select file.* from file join PartDF on file.rowid = PartDF.part",
                           dbname = tempfile(),
                           header = TRUE)

  df.cols <- as.data.frame(df[, cols])
  names(df.cols) <- cols

  return(list(data =  df.cols,
              coords = df[, coord.names]))
}

## sample_partitions ----
#' @title Randomly sample a partition matrix for partitioned GLS
#'
#' @family partitionedGLS
#'
#' @description Create a matrix whose columns contain indices of non-overlapping
#' random samples.
#'
#' @param npix number of pixels in full dataset
#' @param npart number of partitions to create
#' @param partsize size of each partition
#' @param pixels vector of pixel indexes to sample from
#' @param verbose logical: TRUE prints additional info
#'
#' @details
#'
#' If both \code{npart} and \code{partsize} is specified, a partition matrix with
#' these dimensions is returned. If only \code{npart}, is specified,
#' \code{partsize} is selected as the largest integer possible  that creates
#' equal sized partitions. Similarly, if only \code{npart = NA}, then \code{npart}
#' is selected to obtain as many partitions as possible.
#'
#' @return \code{sample_partitions} returns a matrix with \code{partsize}
#' rows and \code{npart} columns. Columns contain random, non-overlapping samples
#' from \code{1:npix}
#'
#' @examples
#' # dummy data with 100 pixels and 20 time points
#' dat.M <- matrix(rnorm(100*20), ncol = 20)
#'
#' # 4 partitions (exhaustive)
#' sample_partitions(npix = nrow(dat.M), npart = 4)
#'
#' # partitions with 10 pixels each (exhaustive)
#' sample_partitions(npix = nrow(dat.M), partsize = 10)
#'
#' # 4 partitions each with 10 pixels (non-exhaustive, produces warning)
#' sample_partitions(npix = nrow(dat.M), npart = 4, partsize = 10)
#'
#' # index of 50 pixels to use as subset
#' sub.indx <- c(1:10, 21:25, 30:62, 70:71)
#'
#' # 5 partitions (exhaustive) from only the specified pixel subset
#' sample_partitions(npix = nrow(dat.M), npart = 5, pixels = sub.indx)
#'
#' @export
sample_partitions <- function(npix, npart = 10, partsize = NA,
                              pixels = NA, verbose = FALSE){

  if(all(!is.na(pixels)) & (length(pixels) > 1)){
    npix = length(pixels)
    from = pixels
  } else {
    from = 1:npix
  }

  ## check which npart of partsize was given
  no.partsize <- (is.na(partsize) | is.null(partsize))
  no.npart <- (is.na(npart) | is.null(npart))

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

  # remainder = npix %% partsize
  # samp_size = npix - remainder
  samp_size = npart * partsize
  samp <- sample(from, size = samp_size, replace = FALSE)

  part.mat <- matrix(samp, ncol = npart, nrow = partsize)
  colnames(part.mat) <- paste("part",1:npart, sep = ".")

  return(part.mat)
}
