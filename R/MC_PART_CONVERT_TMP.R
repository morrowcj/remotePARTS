# Function to convert raw MC_GLSpart() output to partGLS format
MCGLS_partsummary <- function(MCpartGLS, covar.pars = NA, save.GLS = FALSE,
                              partsize){
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

## full multicore GLS function
multicore_fitGLS <- function(formula, partmat, formula0 = NULL, part_FUN,
                             distm_FUN, covar_FUN, covar.pars ,nugget, ncross,
                             save.GLS, ncores, do.t.test=TRUE,
                             do.chisqr.test = TRUE, ...){
  # run GLS
  MCGLS_raw = MC_GLSpart(formula = formula, partmat = partmat, formula0 = formula0,
                     part_FUN = part_FUN, distm_FUN = distm_FUN,
                     covar_FUN = covar_FUN, covar.pars = covar.pars,
                     nugget = nugget, ncross = ncross, save.GLS = save.GLS,
                     ncores = ncores, debug = debug, ...)
  # reformat
  MCGLS = MCGLS_partsummary(MCGLS_raw, covar.pars = covar.pars,
                            save.GLS = save.GLS, partsize=nrow(partmat))
  MCGLS$call = match.call() # update call

  # add optional tests
  if(do.chisqr.test & as.formula(formula) != as.formula(formula0)){
    MCGLS$overall$pval.chisqr = tryCatch(
      chisqr(MCGLS),
      error = function(e){base::warning("error in chisqr()")},
      warning = function(w){base::warning("warning in chisqr()")})
  }

  if(do.t.test){
    test_results = t.test(MCGLS)
    MCGLS$overall$t.test = test_results$p.t
    MCGLS$overall$covar_coef = test_results$covar_coef
  }

  return(MCGLS)
}
