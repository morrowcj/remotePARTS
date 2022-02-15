# #' @rdname fitGLS.partition
# #'
# #' @param ncores number of cores for parallel processing. Default is total cores - 1
# #' @param export an optional character vector of names for any additional objects
# #' needed for \code{part_f()}.
# #' @param debug logical debug flag. If TRUE, prints the name of the step that is running
# #'
# #' @export
# #'
# #' @examples
# #' ## now w. 2 cores:
# #'
# #' if(FALSE){ # change to TRUE to run multi-core version
# #'   GLS.part.mc = fitGLS.partition.mc(part_f = "part_csv", dist_f = "dist_km",
# #'                                     partsize = nrow(parts), npart = ncol(parts),
# #'                                     V.meth = "exponential", spatcor = .5,
# #'                                     part_csv_path = data.file, part_mat = parts,
# #'                                     part_form = "cls.coef ~ 0 + land",
# #'                                     part_form0 = "cls.coef ~ 1",
# #'                                     ncores = 4)
# #' }
# fitGLS.partition.mc <- function(part_f = "part_csv", dist_f = "dist_km",
#                                 V.meth = "exponential", spatcor,
#                                 partsize, npart, mincross = 6,
#                                 ncores = parallel::detectCores() - 1,
#                                 export = NA, debug = TRUE,
#                                 ...){
#   # Setup ----
#   if(debug){print("1. setup")}
#   ## setup cluster
#   clst = parallel::makeCluster(ncores)
#   ## After the function is run, close the cluster
#   on.exit(parallel::stopCluster(clst))
#   ## Register parallel backend
#   doParallel::registerDoParallel(clst)
#
#   ## match the input functions
#   func <- match.fun(part_f)
#   D_func <- match.fun(dist_f)
#
#   ## GLS for each partition ----
#   if(debug){print("2. part GLS")}
#   part_out = foreach::foreach(i = 1:npart, .packages = "remotePARTS") %dopar% {
#     out.i = func(i, ...) # get partition data
#     # out.i = func(i, part_form, part_csv_path, part_mat, part_locvars, part_form0)
#     if(i == 1){
#       # checks
#       stopifnot("part_f(i)$X must be a matrix" = is.matrix(out.i$X))
#       stopifnot("part_f(i)$coords must be a matrix" = is.matrix(out.i$coords))
#     }
#     # out.i = func(i, csv.path = data.file, part.mat = parts)
#     D.i <- D_func(out.i$coords) # calculate D
#     V.i <- fitV(D.i, spatcor, V.meth) # calculate V
#     gls.i <- GLS_worker(y = out.i$y, X = out.i$X, V = V.i, X0 = out.i$X0, save_xx = TRUE)
#     # "return" statement
#
#     list(data = out.i, GLS = gls.i)
#   }
#
#   ## calculate df
#   dfs <- calc_dfpart(partsize, p = ncol(part_out[[1]]$data$X), p0 = ncol(part_out[[1]]$data$X0))
#
#
#   ## Cross-partition setup ----
#   ## all possible pairwise combinations
#   mincross = ifelse(mincross > npart, npart, mincross)
#   used.combs = t(utils::combn(mincross, 2))
#   ncrosses = nrow(used.combs)
#   # maxcross = nrow(combs)
#   # ## adjust mincross if it is greater than maxcross
#   # if (mincross > maxcross){
#   #   ### NEED TO FIX THIS !!!!!
#   #   mincross <- maxcross
#   # }
#   #
#   # ## Calculate which combinations to use
#   # if (mincross < maxcross){
#   #   used.combs = combs[sample(nrow(combs), mincross), ]
#   # } else {
#   #   used.combs = combs[sample(maxcross, mincross), ]
#   # }
#   # used.combs = used.combs[order(used.combs[, 1]), ]
#
#   ## Cross-partition GLS ----
#   if(debug){print("3. cross-partition GLS")}
#   cross = NULL
#   cross_out = foreach::foreach(cross = iterators::iter(used.combs, by = "row"),
#                                .packages = "remotePARTS") %dopar% {
#                                  i = cross[1]; j = cross[2]
#
#                                  ## recompute cholesky from distances
#                                  Di <- D_func(part_out[[i]]$data$coords)
#                                  Vi <- fitV(Di, spatcor, V.meth) # calculate V
#                                  invchol_i = invert_chol(Vi, part_out[[i]]$GLS$nugget)
#
#                                  Dj <- D_func(part_out[[j]]$data$coords)
#                                  Vj <- fitV(Dj, spatcor, V.meth) # calculate V
#                                  invchol_j = invert_chol(Vj, part_out[[j]]$GLS$nugget)
#
#                                  ## cross-partition varcovar matrix
#                                  D.ij <- D_func(part_out[[i]]$data$coords, part_out[[j]]$data$coords)
#                                  V.ij <- fitV(D.ij, spatcor, V.meth)
#                                  cross.ij = crosspart_worker(xxi = part_out[[i]]$GLS$xx, xxj = part_out[[j]]$GLS$xx,
#                                                              xxi0 = part_out[[i]]$GLS$xx0, xxj0 = part_out[[j]]$GLS$xx0,
#                                                              invChol_i = invchol_i,
#                                                              invChol_j = invchol_j,
#                                                              Vsub = V.ij,
#                                                              nug_i = part_out[[i]]$GLS$nugget, nug_j = part_out[[j]]$GLS$nugget,
#                                                              df1 = dfs[1], df2 = dfs[2])
#                                  # Only return necessary output to save space
#                                  list(rSSRij = cross.ij$rSSRij,
#                                       rSSEij = cross.ij$rSSEij,
#                                       rcoefij = cross.ij$rcoefij)
#                                }
#
#
#   # Results Collection ----
#   if(debug){print("4. Results Collection")}
#   f1 = part_out[[1]]$data
#   ## setup output
#   betas = matrix(NA, ncol = ncol(f1$X), nrow = npart,
#                  dimnames = list(NULL, colnames(f1$X)))
#   SEs = betas
#   t.stats = betas
#   pvals.t = betas
#   mod.stats = matrix(NA, nrow = npart, ncol = 7,
#                      dimnames = list(NULL, c("nugget","logLik", "SSE", "MSE",
#                                              "MSR", "F.stat", "pval.F")))
#   cross.SS = matrix(NA, nrow = ncrosses, ncol = 2,
#                     dimnames = list(NULL, c("rSSR", "rSSE")))
#   rcoefs = matrix(NA, nrow = ncrosses, ncol = ncol(f1$X),
#                   dimnames = list(NULL, colnames(f1$X)))
#   ## fill in the output
#   # Partition stats
#   for(i in 1:npart){
#     betas[i, ] <- part_out[[i]]$GLS$betahat
#     SEs[i, ] <- part_out[[i]]$GLS$SE
#     t.stats[i, ] <- part_out[[i]]$GLS$tstat
#     pvals.t[i, ] <- part_out[[i]]$GLS$pval.t
#     mod.stats[i, ] <- c(part_out[[i]]$GLS$nugget, part_out[[i]]$GLS$logLik,
#                         part_out[[i]]$GLS$SSE, part_out[[i]]$GLS$MSE,
#                         part_out[[i]]$GLS$MSR,
#                         part_out[[i]]$GLS$Fstat, part_out[[i]]$GLS$pval.F)
#   }
#   # Cross-partition stats
#   for(ij in 1:length(cross_out)){
#     cross.SS[ij, ] = c(cross_out[[ij]]$rSSRij, cross_out[[ij]]$rSSEij)
#     rcoefs[ij, ] = as.vector(cross_out[[ij]]$rcoefij)
#   }
#   # Overall Statistics
#   fmean = mean(mod.stats[, "F.stat"])
#   rSSRmean = mean(cross.SS[, "rSSR"])
#   rSSEmean = mean(cross.SS[, "rSSE"])
#   rcoefmean = colMeans(rcoefs)
#   coefmean = colMeans(betas)
#
#   # Output ----
#   if(debug){print("5. Output")}
#   out.list <- list(call = match.call(),
#                    part.stats = list("coefficients" = betas,
#                                      "SEs" = SEs,
#                                      "t.stats" = t.stats,
#                                      "pvals.t" = pvals.t,
#                                      "mod.stats" = mod.stats),
#                    cross.stats = list(cross.SS = cross.SS, rcoefs = rcoefs),
#                    overall.stats = list("dfs" = dfs,
#                                         "coefmean" = coefmean,
#                                         "rcoefmean" = rcoefmean,
#                                         "meanstats" = c("fmean" = fmean,
#                                                         "rSSRmean" = rSSRmean,
#                                                         "rSSEmean" = rSSEmean))
#   )
#   class(out.list) <- append("remoteGLS.parts", class(out.list))
#   return(out.list)
# }
