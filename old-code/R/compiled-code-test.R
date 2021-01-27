# Test R vs. compiled versions.
n.samples <- 100

## spatialcor.fit ----

spatialcor.fit.data.c <- compiler::cmpfun(spatialcor.fit.data)

R.time <- matrix(ncol = 5, nrow = n.samples)
C.time <- matrix(ncol = 5, nrow = n.samples)

for(i in 1:n.samples){
R.time[i, ] <- system.time(r.est <- spatialcor.fit.data(X, t.scale, data = dat,
                                         fit.n.sample = fit.n.sample,
                                         plot.fig = F,
                                         FUN = "exponential-power"))


C.time[i, ] <- system.time(r.est.c <- spatialcor.fit.data.c(X, t.scale, data = dat,
                                             fit.n.sample = fit.n.sample,
                                             plot.fig = F,
                                             FUN = "exponential-power"))
}

tmp <- rbind(R = colMeans(R.time), C = colMeans(C.time))[ ,1:3]
colnames(tmp) <- c("user", "system", "elapsed")
tmp

# user system elapsed
# R 0.2061 0.0027  0.2127
# C 0.2066 0.0014  0.2725

# compiled is slightly slower

## V.fit ----
V.fit.c <- compiler::cmpfun(V.fit)

V.time.R <- matrix(ncol = 5, nrow = n.samples)
V.time.C <- matrix(ncol = 5, nrow = n.samples)
for(i in 1:n.samples){
  V.time.R[i, ] <- system.time(V.r <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential-power"))
  V.time.C[i, ] <- system.time(V.c <- V.fit.c(Dist, spatialcor=r.est$spatialcor, FUN="exponential-power"))
}

tmp <- rbind(R = colMeans(V.time.R), C = colMeans(V.time.C))[, 1:3]
colnames(tmp) <- c("user", "system", "elapsed")
tmp

# user system elapsed
# R 0.9533  1e-04  0.9532
# C 0.9533  2e-04  0.9533

# compiled is effectively identical.

# ## inverse cholesky ----
#
# # I don't think these functions are doing the same thing.
#
# inv.chol <- function(M){t(backsolve( chol(M), diag(dim(M)[1]) ))}
# inv.chol.C <- compiler::cmpfun(inv.chol)
# inv.cpp <- Rcpp::cppFunction('
#   arma::mat funcCpp (arma::mat M){
#     arma::mat C = inv(trimatu( chol(M) ));
#     return C;
#   }', depends = "RcppArmadillo")
#
# invchol.A <- matrix(ncol = 5, nrow = n.samples)
# invchol.B <- matrix(ncol = 5, nrow = n.samples)
# invchol.cpp <- matrix(ncol = 5, nrow = n.samples)
# invchol.C <- matrix(ncol = 5, nrow = n.samples)
#
# for (i in 1:n.samples){
#   invchol.A[i, ] <- system.time(inv.chol(V))
#   invchol.B[i, ] <- system.time(chol2inv(chol(V)))
#   invchol.C[i, ] <- system.time(inv.chol.C(V))
#   invchol.cpp[i, ] <- system.time(inv.cpp(V))
# }
#
# rbind(R = colMeans(invchol.A), R.alt = colMeans(invchol.B),
#       C = colMeans(invchol.C), Cpp = colMeans(invchol.cpp))


