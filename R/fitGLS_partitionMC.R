#' @title fit a parallel partitioned GLS
#'
#' @rdname partGLS
#'
#' @return a "MC_partGLS", which is a precursor to a "partGLS" object
#'
MC_GLSpart <- function(formula, partmat, formula0 = NULL, part_FUN = "part_data",
                       distm_FUN = "distm_scaled", covar_FUN = "covar_exp",
                       covar.pars = c(range = .1), nugget = NA, ncross = 6,
                       save.GLS = FALSE,
                       ncores = parallel::detectCores(logical = FALSE) - 1,
                       debug = FALSE,
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


  # Parallel setup
  ## check for required packages
  requireNamespace(c("parallel", "doParallel", "foreach"))
  # setup cluster
  clst = parallel::makeCluster(ncores)
  ## After the function is run, close the cluster
  on.exit(parallel::stopCluster(clst))
  ## Register parallel backend
  doParallel::registerDoParallel(clst)


  # Run parallel computations, collect all output in a list
  i <- NULL # variable declaration for iterator
  out_list <- foreach::foreach(i = 1:npart, .packages = "remotePARTS") %dopar% {
    # Partition i
    ## partition data
    idat <- part.f(partmat[, i], formula = formula, formula0 = formula0, ...)
    # idat <- part.f(partmat[, i], formula = formula, formula0 = formula0, data = data)

    ## covariance of parition i
    Vi = do.call(covar.f, args = append(list(d = dist.f(idat$coords)), as.list(covar.pars)))
    ## GLS of parition i
    partGLS.i <- fitGLS(formula = formula, data = idat$data, V = Vi,
                           nugget = nugget, formula0 = formula0, save.xx = (i <= ncross),
                           no.F = FALSE, save.invchol = (i <= ncross), logLik.only= FALSE,
                           ncores = 1)

    # Partition J
    rGLS.list <- list() # empty list
    if (i < ncross) for (j in (i+1):ncross){
      j.lab <- paste0("j.",j)
      ## parition data
      jdat = part.f(partmat[, j], formula = formula, formula0 = formula0, ...)
      # jdat = part.f(partmat[, j], formula = formula, formula0 = formula0, data = data)

      ## Covariance of parititon j
      Vj = do.call(covar.f, args = append(list(d = dist.f(jdat$coords)), as.list(covar.pars)))
      ## GLS of partition j
      partGLS.j <- fitGLS(formula = formula, data = jdat$data, V = Vj,
                             nugget = nugget, formula0 = formula0, save.xx = TRUE,
                             no.F = FALSE, save.invchol = TRUE, logLik.only= FALSE,
                             ncores = 1)
      ## check for dimension mismatch
      dim.mismatch <- ifelse(length(partGLS.j$coefficients) != length(partGLS.i$coefficients), TRUE, FALSE)

      # Cross partition
      if(dim.mismatch){
        ## output text string if dimension mistmatch
        rGLS.list[[j.lab]] <- "Error: Dimension mistmatch between parition i and j."
      } else {
        ## cross-covariance
        Vij <- do.call(covar.f, args = append(list(d = dist.f(idat$coords, jdat$coords)),
                                              as.list(covar.pars)))
        ## degrees of freedom
        dfs = calc_dfpart(partsize = partsize, p = ncol(partGLS.j$xx), p0 = ncol(partGLS.j$xx0))
        ## cross-GLS
        rGLS <- crosspart_GLS(xxi = partGLS.i$xx,
                              xxj = partGLS.j$xx,
                              xxi0 = partGLS.i$xx0,
                              xxj0 = partGLS.j$xx0,
                              invChol_i = partGLS.i$invcholV,
                              invChol_j = partGLS.j$invcholV,
                              Vsub = Vij,
                              nug_i = partGLS.i$nugget,
                              nug_j = partGLS.j$nugget,
                              df1 = dfs[1], df2 = dfs[2],
                              ncores = ncores)
        if(debug){
          ## output full cross, if debug mode
          rGLS.list[[j.lab]] <- rGLS
        } else {
          ## otherwise, only important stats
          rGLS.list[[j.lab]] <- list(rcoefij = rGLS$rcoefij,
                                               rSSRij = rGLS$rSSRij,
                                               rSSEij = rGLS$rSSEij)
        }
      }
    }

    ## remove large variables that are unnecessary, unless in debug mode or if save.GLS specified
    if(!debug & !save.GLS){
      partGLS.i$xx <- partGLS.i$xx0 <- partGLS.i$invcholV <- NULL
    }

    ## combine partitioned results and cross-partitioned output
    out.i <- list(partGLS = partGLS.i, crossGLS = rGLS.list)

    out.i
  }

  ## change output class and print output
  class(out_list) <- append("MC_partGLS", class(out_list))
  out_list
}

#' Convert raw MC_GLSpart() output to partGLS format
#'
#' @param MCpartGLS object resulting from MC_partGLS()
#' @param partsize number of locations per partition
#'
#' @rdname partGLS
#'
#' @return a "partGLS" object
#'
MCGLS_partsummary <- function(MCpartGLS, covar.pars = c(range = .1),
                              save.GLS = FALSE, partsize){
  # collect info
  nparts = length(MCpartGLS)  # check number of partitions
  first_coefs = MCpartGLS[[1]]$partGLS$coefficients
  ncoefs = length(first_coefs)  # number of independent variables
  coef_names = names(first_coefs)  # names of independent variables
  npairs = sum(sapply(MCpartGLS, function(x)length(x$crossGLS)))
  # null model info
  first_null_coefs = MCpartGLS[[1]]$partGLS$coefficients0
  ncoefs_null = length(first_null_coefs)
  null_coef_names = names(first_null_coefs)

  # setup empty data structures
  ## matrices
  matrix_skeleton = matrix(NA, nrow = nparts, ncol = ncoefs,
                           dimnames = list(NULL, coef_names))
  coefs = SEs = tstats = tpvals = matrix_skeleton
  ## vectors (known length)
  vector_skeleton = rep(NA, times = nparts)
  nuggets = LLs = SSEs = MSEs = MSRs = Fstats = Fpvals = vector_skeleton
  ## vectors (unknown length)
  rSSEs = NULL
  rSSRs = NULL
  if(save.GLS){ # conditionally create GLS list
    GLS_list = list()
  } else {
    GLS_list = NULL
  }
  ## arrays
  rcoefs = array(NA, dim = c(npairs, ncoefs, ncoefs),
                 dimnames = list(NULL, coef_names, coef_names))
  covar_coefs = array(NA, dim = c(ncoefs, ncoefs, nparts),
                      dimnames = list(coef_names, coef_names, NULL))

  # loop through each partition
  cross_counter = 0  # to keep track of which cross we're on
  for(i in seq_len(nparts) ){ # loop through paritions
    this_part = MCpartGLS[[i]]  # current parition
    partGLS = this_part$partGLS  # GLS of current parition
    # conditionally save GLS to the list
    if(save.GLS){
      GLS_list[[i]] <- partGLS
    }
    # place output in their proper place
    ## matrices
    coefs[i, ] <- rbind(partGLS$coefficients)
    SEs[i, ] <- rbind(partGLS$SE)
    tstats[i, ] <- rbind(partGLS$tstat)
    tpvals[i, ] <- rbind(partGLS$pval_t)
    ## vectors
    nuggets[i] <- partGLS$nugget
    LLs[i] <- partGLS$logLik
    SSEs[i] <- partGLS$SSE
    MSEs[i] <- partGLS$MSE
    MSRs[i] <- partGLS$MSR
    Fstats[i] <- partGLS$Fstat
    Fpvals[i] <- partGLS$pval_F
    ## arrays
    covar_coefs[,,i] <- partGLS$covar_coef


    # cross-partition
    cross_part = this_part$crossGLS
    if(length(cross_part) >= 1){ # loop through crosses with this partition
      for(j in seq_len(length(cross_part))){
        cross_counter = cross_counter + 1
        # collect cross-parititon statistics
        rSSRs = c(rSSRs, cross_part[[j]]$rSSRij)
        rSSEs = c(rSSEs, cross_part[[j]]$rSSEij)
        rcoefs[cross_counter, ,] <- cross_part[[j]]$rcoefij
      } # end cross loop
    }
  } # end part loop

  # Calculate summary statistics
  coefmeans = colMeans(coefs, na.rm = TRUE)
  rSSRmean = mean(rSSRs, na.rm = TRUE)
  rSSEmean = mean(rSSEs, na.rm = TRUE)
  Fmean = mean(Fstats, na.rm = TRUE)

  # arrange output
  ## collect model statistics
  modstats = cbind(LLs = LLs, SSEs = SSEs, MSEs = MSEs, MSRs = MSRs,
                   Fstats = Fstats, Fpvals = Fpvals)
  ## parition summary
  part = list(coefficients = coefs, SEs = SEs, tstats = tstats, tpvals = tpvals,
              nuggets = nuggets, covar.pars = covar.pars, modstats = modstats,
              covar_coefs = covar_coefs)
  ## cross-parition summary
  rcoefficients = apply(rcoefs, MARGIN = c(2, 3),
                        FUN = function(x){mean(x, na.rm = TRUE)})
  cross = list(rcoefs = rcoefs, rSSRs = rSSRs, rSSEs = rSSEs)
  ## overall summary
  overall = list(coefficients = coefmeans, rcoefficients = rcoefficients,
                 rSSR = rSSRmean, rSSE = rSSEmean, Fstat = Fmean,
                 dfs = calc_dfpart(partsize, ncoefs, ncoefs_null),
                 partdims = c(npart = nparts, partsize = partsize)
  )

  ## return formatted output (partGLS)
  return_list = list(call = match.call(), GLS = GLS_list,
                     part = part, cross = cross, overall = overall)
  class(return_list) <- append("partGLS", class(return_list))

  return(return_list)
}

#' Conduct multi-core partitioned GLS
#'
#' @return "partGLS" object
#' @rdname partGLS
#'
#' @export
multicore_fitGLS_partition <- function(formula, partmat, formula0 = NULL,
                                       part_FUN = "part_data",
                                       distm_FUN = "distm_scaled",
                                       covar_FUN = "covar_exp",
                                       covar.pars = c(range = 0.1),
                                       nugget = NA, ncross = 6, save.GLS = FALSE,
                                       ncores = parallel::detectCores(logical = FALSE) - 1,
                                       do.t.test = TRUE,
                                       do.chisqr.test = TRUE,
                                       debug = FALSE,
                                       ...){
  if(is.null(formula0)){
    formula0 = update(as.formula(formula), . ~ 1)
  } else {
    formula0 = as.formula(formula0)
  }
  # run GLS
  if(debug){cat("\nconducting GLS\n")}
  MCGLS_raw = MC_GLSpart(formula = formula, partmat = partmat, formula0 = formula0,
                         part_FUN = part_FUN, distm_FUN = distm_FUN,
                         covar_FUN = covar_FUN, covar.pars = covar.pars,
                         nugget = nugget, ncross = ncross, save.GLS = save.GLS,
                         ncores = ncores, debug = debug, ...)
  # reformat
  if(debug){cat("\nreformatting GLS\n")}
  MCGLS = MCGLS_partsummary(MCGLS_raw, covar.pars = covar.pars,
                            save.GLS = save.GLS, partsize=nrow(partmat))
  if(debug){cat("\nupdating call\n")}
  MCGLS$call = match.call() # update call

  # add optional tests
  if(do.chisqr.test & as.formula(formula) != as.formula(formula0)){
    if(debug){"\nperforming chisqr test\n"}
    MCGLS$overall$pval.chisqr = tryCatch(
      chisqr(MCGLS),
      error = function(e){base::warning("error in chisqr()")},
      warning = function(w){base::warning("warning in chisqr()")})
  }

  if(do.t.test){
    if(debug){"\nperforming t-test\n"}
    test_results = t.test(MCGLS)
    MCGLS$overall$t.test = test_results$p.t
  }

  return(MCGLS)
}
