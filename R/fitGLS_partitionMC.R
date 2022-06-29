#' @title fit a parallel parititoned GLS
#'
#' @rdname partGLS
#'
#' @return
#'
#' @examples
#' \dontrun{
#' ## read data
#' data(ndvi_AK10000)
#' df = ndvi_AK10000[seq_len(2000), ] # first 1000 rows
#'
#' ## create partition matrix
#' pm = sample_partitions(nrow(df), npart = 3)
#'
#' ## fit GLS with fixed nugget
#' MCpartGLS = remotePARTS:::MC_GLSpart(formula = CLS_coef ~ 0 + land, partmat = pm,
#'                                      data = df, nugget = 0, ncores = 2L)
#'
#' (partGLS.mc = remotePARTS:::MC_to_partGLS(MCpartGLS, partsize = nrow(pm), npart = ncol(pm)))
#'
#' MCpartGLS2 = remotePARTS:::MC_GLSpart(formula = CLS_coef ~ lat, partmat = pm,
#'                                      data = df, nugget = 0, ncores = 2L)
#'
#' (partGLS2.mc = remotePARTS:::MC_to_partGLS(MCpartGLS2, partsize = nrow(pm), npart = ncol(pm)))
#' }
MC_GLSpart <- function(formula, partmat, formula0 = NULL,
                               part_FUN = "part_data",
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

    ## combine paritioned results and cross-partitioned output
    out.i <- list(partGLS = partGLS.i, crossGLS = rGLS.list)

    out.i
  }

  ## change output class and print output
  class(out_list) <- append("MC_partGLS", class(out_list))
  out_list
}


#' Convert parallel parititoned GLS into partGLS
#'
#' @param object "MC_partGLS" object
#' @param partsize number of pixels per partition
#' @param npart number of partitions
#'
#' @rdname partGLS
MC_to_partGLS <- function(object, covar.pars = c(range = .1), partsize, npart,
                          save.GLS = FALSE,
                          do.t.test = TRUE, do.chisqr.test = TRUE, debug = FALSE){
  stopifnot("MC_partGLS" %in% class(object))

  ## Setup output to be like fitGLS_parititon
  part = list(coefficients = t(sapply(object, function(x)x$partGLS$coefficients)),
              SEs = t(sapply(object, function(x)x$partGLS$SE)),
              tstats = t(sapply(object, function(x)x$partGLS$tstat)),
              tpvals = t(sapply(object, function(x)x$partGLS$pval_t)),
              nuggets = sapply(object, function(x)x$partGLS$nugget),
              covar.pars = covar.pars,
              modstats = cbind(LLs = sapply(object, function(x)x$partGLS$logLik),
                               SSEs = sapply(object, function(x)x$partGLS$SSE),
                               MSEs = sapply(object, function(x)x$partGLS$MSE),
                               MSRs = sapply(object, function(x)x$partGLS$MSR),
                               Fstats = sapply(object, function(x)x$partGLS$Fstat),
                               Fpvals = sapply(object, function(x)x$partGLS$pval_F)))

  p = ncol(part$coefficients)
  p0 = length(object[[1]]$partGLS$coefficients0)

  rSSRs <- NULL
  rSSEs <- NULL
  # rcoefs <- NULL
  npairs = sum(sapply(object, function(x)length(x$crossGLS)))

  rcoefs = array(NA, dim = c(npairs, p, p), dimnames = list(NULL, names(object[[1]]$partGLS$coefficients), names(object[[1]]$partGLS$coefficients)))
  covar_coefs = array(NA, dim = c(p, p, npart), dimnames = list(names(object[[1]]$partGLS$coefficients), names(object[[1]]$partGLS$coefficients), NULL))

  tmp <- lapply(object, function(x)x$partGLS$covar_coef)

  m = 0
  for(i in seq_len(length(object))){
    x = object[[i]]
    # covar_coefs[, , i] = x$partGLS$covar_coef
    covar_coefs[,,i] = tmp[[i]]

    if (!length(x$crossGLS) < 1){
      for(j in seq_len(length(x$crossGLS))){
        rSSRs = c(rSSRs, x$crossGLS[[j]]$rSSRij)
        rSSEs = c(rSSEs, x$crossGLS[[j]]$rSSEij)
        m = m + 1
        rcoefs[m, ,] = x$crossGLS[[j]]$rcoefij
      }
    }
  }


  colnames(rcoefs) = colnames(part$coefficients)
  cross = list(rcoefs = rcoefs, rSSRs = rSSRs, rSSEs = rSSEs)

  rcoefficients = apply(cross$rcoefs, MARGIN = c(2, 3), FUN = function(x){mean(x, na.rm = TRUE)})
  dimnames(rcoefficients)[[1]] = dimnames(rcoefficients)[[2]] = dimnames(rcoefs)[[2]]

  dimnames(covar_coefs)[[1]] = dimnames(covar_coefs)[[2]] = dimnames(rcoefs)[[2]]
  part$covar_coefs = covar_coefs

  overall = list(coefficients = colMeans(part$coefficients, na.rm = TRUE),
                 rcoefficients = rcoefficients,
                 rSSR = mean(rSSRs, na.rm = TRUE),
                 rSSE = mean(rSSEs, na.rm = TRUE),
                 Fstat = mean(part$modstats[, "Fstats"], na.rm = TRUE),
                 dfs = calc_dfpart(partsize, p, p0),
                 partdims = c(npart = length(object), partsize = partsize))

  outlist <- list(call = match.call(),
                  GLS = if(save.GLS){lapply(object, function(x)x$partGLS)}else{NULL},
                  part = part,
                  cross = cross,
                  overall = overall)

  if(debug){outlist$crossGLS = lapply(object, function(x)x$crossGLS)}
  class(outlist) <- append("partGLS", class(outlist))
  if(do.chisqr.test){
    if(as.formula(object[[1]]$partGLS$formula) == as.formula(object[[1]]$partGLS$formula0)){
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

  return(outlist)
}
