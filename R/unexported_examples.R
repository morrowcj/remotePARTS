# # @examples for unexported functions.
#
# # ---- AR_fun ----
#
# ## simulate dummy data
# x = rnorm(31)
# time = 1:30
# x = x[2:31] + x[1:30] + 0.3*time #AR(1) process + time trend
# U = stats::model.matrix(formula(x ~ time))
#
# ## fit an AR
# remotePARTS:::AR_fun(par = .2, y = x, X = U, logLik.only = FALSE)
#
# ## get the partial logLik of the AR parameter, given the data.
# remotePARTS:::AR_fun(par = .2, y = x, X = U, logLik.only = FALSE)
#
# ## show that minimizing the partial logLik maximizes the true logLik (NOT RUN)
# n = 100
# out.mat = matrix(NA, nrow = n, ncol = 3,
#                  dimnames = list(NULL, c("par", "logLik", "partialLL")))
# out.mat[, "par"] = seq(-10, 10, length.out = n)
# for (i in seq_len(n) ) {
#    p = out.mat[i, "par"]
#    out.mat[i, "logLik"] = remotePARTS:::AR_fun(par = p, y = x, X = U, logLik.only = TRUE)
#    out.mat[i, "partialLL"] = remotePARTS:::AR_fun(par = p, y = x, X = U,
#                                                   logLik.only = FALSE)$logLik
# }
# plot(x = out.mat[, "partialLL"], y = out.mat[, "logLik"])
#
# # ---- fitGLS_opt_FUN ----
# ## read data
# data(ndvi_AK10000)
# df = ndvi_AK10000[seq_len(200), ] # first 500 rows
# coords = df[, c("lng", "lat")]
# remotePARTS:::fitGLS_opt_FUN(op = c(range = .1, nugget = .2),
#                              formula = CLS_coef ~ 0 + land, data = df,
#                              coords = coords)
# remotePARTS:::fitGLS_opt_FUN(op = c(range = .1), fp = c(nugget = 0),
#                              formula = CLS_coef ~ 0 + land, data = df,
#                              coords = coords)
#
# logit <- function(p) {log(p / (1 - p))}
# inv_logit <- function(l) {1 / (1 + exp(-l))}
#
# ## input logit-transformed range parameters
# remotePARTS:::fitGLS_opt_FUN(op = c(range = .1, nugget = logit(.2)),
#                              formula = CLS_coef ~ 0 + land, data = df,
#                              coords = coords, is.trans = TRUE,
#                              backtrans = list(nugget = inv_logit))
# ## transformed range and nugget
# remotePARTS:::fitGLS_opt_FUN(op = c(range = logit(.1), nugget = logit(.2)),
#                              formula = CLS_coef ~ 0 + land, data = df,
#                              coords = coords, is.trans = TRUE,
#                              backtrans = list(nugget = inv_logit, range = inv_logit))
#
#
# # ---- MC_GLSpart ----
# ## read data
# data(ndvi_AK10000)
# df = ndvi_AK10000[seq_len(1000), ] # first 1000 rows
#
# ## create partition matrix
# pm = sample_partitions(nrow(df), npart = 3)
#
# ## fit GLS with fixed nugget
# MCpartGLS = remotePARTS:::MC_GLSpart(formula = CLS_coef ~ 0 + land, partmat = pm,
#                                      data = df, nugget = 0, ncores = 2L)
#
# (partGLS.mc = remotePARTS:::MCGLS_partsummary(MCpartGLS, partsize = nrow(pm)))
#
# MCpartGLS2 = remotePARTS:::MC_GLSpart(formula = CLS_coef ~ lat, partmat = pm,
#                                      data = df, nugget = 0, ncores = 2L)
#
# (partGLS2.mc = remotePARTS:::MCGLS_partsummary(MCpartGLS2, partsize = nrow(pm)))
#
# MCpartGLS_int = remotePARTS:::MC_GLSpart(formula = CLS_coef ~ 1, partmat = pm,
#                                          data = df, nugget = 0, ncores = 2L)
#
# (partGLS_int.mc = remotePARTS:::MCGLS_partsummary(MCpartGLS_int, partsize = nrow(pm)))
#
# # --- optimize_nugget ----
# ## read data
# data(ndvi_AK10000)
# df = ndvi_AK10000[seq_len(200), ] # first 200 rows
#
# ## format data
# X = stats::model.matrix(CLS_coef ~ 0 + land, data = df)
#
# ## fit covariance matrix
# V = covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = .01)
#
# ## find the ML nugget
# remotePARTS:::optimize_nugget(X = X, V = V, y = df$CLS_coef, debug = TRUE)
#
# # ---- part_chisqr ----
# remotePARTS:::part_chisqr(Fmean = 3.6, rSSR = .021, df1 = 2, npart = 5)
#
# # ---- calc_dfpart ----
# calc_dfpart(partsize = 2000, p = 4, p0 = 1)
#
# # ---- crosspart_GLS ----
# ## read data
# data(ndvi_AK10000)
# df = ndvi_AK10000[seq_len(1000), ] # first rows
#
# # partition matrix
# pm = sample_partitions(nrow(df), npart = 2, partsize = 500)
#
# ## partition data
# data.i = df[pm[, 1], ]
# data.j = df[pm[, 2], ]
#
# ## partition coordinates
# coords.i = data.i[, c("lng", "lat")]
# coords.j = data.j[, c("lng", "lat")]
#
# ## partition covariance
# V.i = covar_exp(distm_scaled(coords.i), range = .01)
# V.j = covar_exp(distm_scaled(coords.j), range = .01)
#
# ## partition GLS
# GLS.i = fitGLS(CLS_coef ~ 0 + land, data.i, V.i, nugget = 0, save.xx = TRUE,
#                save.invchol = TRUE, no.F = FALSE)
# GLS.j = fitGLS(CLS_coef ~ 0 + land, data.j, V.j, nugget = 0, save.xx = TRUE,
#                save.invchol = TRUE, no.F = FALSE)
#
# ## cross-covariance
# V.ij = covar_exp(distm_scaled(coords.i, coords.j), range = .01)
#
# ## degrees of freedom
# dfs = remotePARTS:::calc_dfpart(partsize = nrow(pm), p = ncol(GLS.i$xx),
#                                 p0 = ncol(GLS.i$xx0))
#
# # Calculate cross-partition statistics
# (crossGLS = remotePARTS:::crosspart_GLS(xxi = GLS.i$xx,
#                                        xxj = GLS.j$xx,
#                                        xxi0 = GLS.i$xx0,
#                                        xxj0 = GLS.j$xx0,
#                                        invChol_i = GLS.i$invcholV,
#                                        invChol_j = GLS.j$invcholV,
#                                        Vsub = V.ij,
#                                        nug_i = GLS.i$nugget,
#                                        nug_j = GLS.j$nugget,
#                                        df1 = dfs[1], df2 = dfs[2]))
#
# # ---- max_dist ----
# coords <- matrix(stats::rnorm(20e6), ncol = 2)  # cloud of 20 million pixels
# remotePARTS:::max_dist(coords)
#
# remotePARTS:::max_dist(coords, dist_FUN = "distm_scaled")
#
#
# # ---- test_covar_fun ----
# # distance vector
# d = seq(0, 1, length.out = 10)
# # named parameter list
# test_covar_fun(d = d, covar_FUN = "covar_exppow", covar.pars = list(range = .5))
# test_covar_fun(d = d, covar_FUN = "covar_exppow", covar.pars = list(range = .5, shape = .5))
# # unnamed parameter vector
# test_covar_fun(d = d, covar_FUN = "covar_exppow", covar.pars = .5)
# test_covar_fun(d = d, covar_FUN = "covar_exppow", covar.pars = c(.5, .5))
# # different function
# test_covar_fun(d = d, covar_FUN = "covar_taper", covar.pars = list(theta = .5))
# # user-defined function, with no extra parameters
# test_covar_fun(d = d, covar_FUN = function(d){return(d)}, covar.pars = NULL)
# test_covar_fun(d = d, covar_FUN = function(d){return(d)}, covar.pars = list())
#
# map.width = 5 # square map width
# coords = expand.grid(x = 1:map.width, y = 1:map.width) # coordinate matrix
#
# # calculate distance
# D = geosphere::distm(coords) # distance matrix
#
# test_covar_fun(D, covar.pars = list(range = .1*max(D)))
#
#
