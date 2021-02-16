# Compare with old code

# setup ----
source("original-code.R")

## starting parameters
run.GLS = TRUE
n.pix = 100 # number of pixels to test with

## get data
set.seed(916)
data("ndvi_AK3000")
subsamp = sample.int(n = nrow(ndvi_AK3000), size = n.pix)
df = ndvi_AK3000[subsamp, ] # subset the data

X = as.matrix(df[, grep(pattern = "ndvi", names(df))])
n.t = ncol(X)
time = seq_len(ncol(X))
time = ( time - min(time) )/( max(time) - min(time) )
land = df[, c("land")]
modmat = model.matrix(~ 0 + land, data = df)
null.modmat = model.matrix(~ 1, data = df)
coords = df[, c("lng", "lat")]
D = geosphere::distm(coords)/1000 # long time
form = "y ~ 0 + land"

# CLS functions ----
old.CLS <-  CLS.fit(X, time)
new.CLS <- fitCLS.map(X, time, TRUE, TRUE, TRUE, TRUE)
expect_equivalent(old.CLS[, "mean"], new.CLS$mean, info = "check mean")
expect_equivalent(old.CLS[, c("c", "t", "p")],
                  new.CLS$time.coef[, c("Est", "t", "p.t")],
                  info = "check coef table")
expect_equivalent(old.CLS[, c("b", "MSE")],
                  data.frame(new.CLS$xi.coef[,"Est"], new.CLS[["MSE"]]),
                  info = "check beta and MSE")
## scale response for GLS: y
y = new.CLS$time.coef$Est/new.CLS$mean
df$y = y

# spatial correlation functions ----
old.r <- spatialcor.fit(X, time, D, fit.n.sample = n.pix)
new.r <- fit_spatialcor(X, time, dist = D, location = coords,
                        fit.n = n.pix)
## compare output
expect_equivalent(old.r$spatialcor, new.r$spatialcor,
                  info = "compare parameters")
expect_equivalent(old.r$spatialcor.sigma, sigma(new.r$mod),
                  info = "compare variance")
## assign r
r = new.r$spatialcor

# Variance fitting ----
expect_equivalent(og.V.fit(D, r), (V <- fitV(D, r)),
                  info = "compare varcovar matrix")

# nugget fitting functions ----
old.nug <- nugget.fit(formula = form, data = df, V = V)
new.nug_C <- optimize_nugget(X = modmat, V = V, y = y)
## compare output
expect_equivalent(old.nug, new.nug_C, info = "compare nugget fitting")

# GLS ----
if(run.GLS){
old.GLS <- GLS.fit(form, data = df,  V = V,
                   invcholV = invert_chol(V, new.nug_C))
GLS.R <- fitGLS_R(X = modmat, V = V, y = y, X0 = null.modmat, nugget = new.nug_C)
new.GLS <- fitGLS(X = modmat, y = y,V = V, new.nug_C, X0 = null.modmat,
                  save_xx = TRUE, threads = 1)
GLS.wrk = GLS_worker(y, modmat, V, null.modmat)
## compare output
expect_equivalent(old.GLS$coef, new.GLS$betahat, GLS.wrk$betahat)
expect_equivalent(old.GLS$t, new.GLS$tstat, GLS.wrk$tstat)
expect_equivalent(old.GLS$SSE, new.GLS$SSE, GLS.wrk$SSE)
expect_equivalent(old.GLS$SSE0, new.GLS$SSE0, GLS.wrk$SSE0)
expect_equivalent(old.GLS$xx, new.GLS$xx, GLS.wrk$xx)
expect_equivalent(old.GLS$xx0, new.GLS$xx0, GLS.wrk$xx0)
expect_equivalent(old.GLS$F, new.GLS$Fstat, GLS.wrk$Fstat)

## Now just check the R function
expect_equivalent(GLS.R$betahat, new.GLS$betahat)
expect_equivalent(GLS.R$tstat, new.GLS$tstat)
expect_equivalent(GLS.R$SSE, new.GLS$SSE)
expect_equivalent(GLS.R$SSE0, new.GLS$SSE0)
expect_equivalent(GLS.R$Fstat, new.GLS$Fstat)
}

# Partitioned GLS ----
# TBA
