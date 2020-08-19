## packages used
# data.table
# microbenchmark
set.seed(58)


# prepare data ----
## This section is not needed for the user, as *.rda files will be saved

## Full data
# read in the data
dtab <- data.table::fread(file = "data-raw/north_america_checked.csv")

# convert $land collumn to a factor
land.classes = c("Evergr needle","Evergr broad","Decid needle","Decid broad",
                 "Mixed forest","Shrubland","Savanna","Grassland","Cropland",
                 "Cropland mosaics") # all land classes (in order of lvs)

dtab$land <- factor(dtab$land, labels = land.classes)

# save the rdata file
if(!file.exists("data/north-america_ndvi.rda")){
  saveRDS(dtab, "data/north-america_ndvi.rda")
}
## Alaska subset
dataAK <- dtab[dtab$lng < -141, ] #North Am. west of -141 is approx. AK

# relative contirbution of each land class to AK
tmp <- aggregate(
  x = dataAK$land,
  by = list(landclass = dataAK$land),
  FUN = function(x)length(x)/nrow(dataAK)
)

# any land class that occurs less than 2% is rare
dataAK$rare.land <- tmp[match(dataAK$land, tmp$landclass), "x"] <= 0.02

# reorder columns so ndvi are all at the end
ord <- c(grep("ndvi", names(dataAK), invert = TRUE), #columns without 'ndvi'
         grep("ndvi", names(dataAK))) # cols with ndvi
dataAK <- dataAK[, ..ord] #..x is data.table notation

# save AK data file
if(!file.exists("data/alaska_ndvi.rda")){
  saveRDS(dataAK, "data/alaska_ndvi.rda")
}

## subset AK further
dataAK_small <- dataAK[!dataAK$rare.land, ] # only common land classes
dataAK_small <- dataAK_small[sample.int(n = nrow(dataAK_small), size = 3000), ] # 3000 pts

# save small AK data file
if(!file.exists("data/alaska_ndvi_3000pts.rda")){
  saveRDS(dataAK_small, "data/alaska_ndvi_3000pts.rda")
}

# load data ----
##

rm(list = ls()) # clear all objects from memory
source("R/lsos.R") # load memory tracking functions

dat <- readRDS("data/alaska_ndvi_3000pts.rda")

library(microbenchmark)


# CLS ----
##

source("R/fitCLS.R") # source the functions

dat <- dat[!dat$rare.land, ] # remove rare land classes
# X data
Xmat <-  as.matrix(dat[, -c(1:6)]) #rows = location, cols = time
n = dim(Xmat)[1] # sites
t.n = dim(Xmat)[2] # time points
t.scale = ((1:t.n) - min(1))/t.n # scaled time

# test speed and check accuracy:
if (FALSE) {
  fit <- star <- matrix(NA, 5, 5)
  for (i in 1:5) {
    fit[i, ] <- system.time(tmp1 <- CLS.fit(Xmat, t.scale))
    star[i, ] <- system.time(tmp2 <- cls_star(Xmat, t.scale))
  }
  rbind("fit" = colMeans(fit, na.rm = TRUE),
        "star" = colMeans(star, na.rm = TRUE))[, 1:3]
  tmp1[1:5, ]
  tmp2[1:5, ]
#       user  sys elapsed
# #fit  2.862 0   2.864
# #star 2.794 0   2.796
# note: cls_star() does not loose speed and can handle covariates:
}


# example with covars:
if (FALSE) {
  # generate some random covariates (just for structure):
  # one covariate
  covar.mat1 <- matrix(rnorm(n*t.n), n, t.n)
  colnames(covar.mat1) = gsub("ndvi", "Z", colnames(Xmat))

  cls_star(Xmat[1:5, ], t.scale, covar.list = list(covar.mat1[1:5, ]))

  ## now 2 covaritates
  covar.mat2 <- matrix(rnorm(n*t.n), n, t.n)
  colnames(covar.mat2) = gsub("ndvi", "U", colnames(Xmat))

  cls_star(X = Xmat[1:5, ], t = t.scale,
           covar.list = list(covar.mat1[1:5, ], covar.mat2[1:5, ]))
}

# fit the initial model and check for outliers
tmp <- cls_star(Xmat, t.scale)
tmp[, c("lat", "lng")] <- dat[, c("lat", "lng")]
tmp$rel.stat <-  (tmp$Est/(1 - tmp$x_t0.EST))/tmp$mean # relative CLS statistic

outl = abs(scale(tmp$rel.stat)) > -qnorm(p = 1/nrow(tmp)/10) # outliers
# # 11 outliers

# plot with outliers
## ggplot way
if (FALSE) {
  library(ggplot2)
  ggplot(tmp[!outl, ], aes(x = lng, y = lat, col = rel.stat)) +
    geom_tile(size = 2) +
    scale_color_gradient2(low = "red", mid = "grey", high = "green",
                          midpoint = 0) +
    geom_tile(data = tmp[outl, ], col = "black", size = 1) # add outliers
  detach("package:ggplot2", unload = TRUE)
}

## base R way
base.col = ifelse(test = outl, yes = "black",
                  no = colorRampPalette(
                    c("orange", "grey", "darkgreen")
                    )(nrow(tmp))[ordered(tmp$rel.stat)]
)
plot(dat$lat ~ dat$lng, pch = 15, cex = .7, col = base.col,
     xlab = "longitude", ylab = "latitude")
legend(x = "bottomright", fill = c("darkgreen", "grey","orange"),
       legend = c("high", "med", "low"), title = "NDVI")
# # Note that the outliers often occur on borders where data is missing.

## relative to overall distribution
hist(tmp$mean[outl], freq = FALSE, breaks = 0:20, ylim = c(0,.5), col = "grey",
     lty = 3, xlab = "site average NDVI", main = "Histogram of NDVI")
hist(tmp$mean[!outl], freq = FALSE, breaks = 0:20, ylim = c(0,.5), col = NULL,
     add = TRUE)
legend("topright", legend = c("TRUE","FALSE"), fill = c("grey", "white"),
       title = "outlier")

# remove outliers and plot
# # Do this later...
# tmp <- tmp[!outl, ]

# GLS ----

# load in functions
source("R/fitSpatCor.R")
source("R/fitGLS.R")

# Distance matrix
location <- dat[, c("lng","lat")]
Dist <- geosphere::distm(location)/1000

# Spatial Correlation
r.est <- fit_spatialcor(X = Xmat, t = t.scale, fit.n = 200 ,
                        location = location,
                        fun = "exp-pwr", scale.dist = TRUE,
                        dist.scl = 1000, plot.fig = FALSE)
## compare speed
if (FALSE){ # compare speed
Old <- NULL
New <- NULL
for (j in 1:10){
old <- system.time({
  r.est.0 <- spatialcor.fit.data(X = Xmat, data = dat, t.scale = t.scale,
                             fit.n.sample = 200,
                             FUN = "taper-spherical", plot.fig = TRUE)
})
new <- system.time({
  r.est <- fit_spatialcor(X = Xmat, t = t.scale, fit.n = 200 ,
                         location = dat[, c("lng","lat")],
                         fun = "taper", scale.dist = TRUE,
                         dist.scl = 1000, plot.fig = TRUE,
                         )
})
Old <- rbind(Old, old)
New <- rbind(New, new)
}
rbind(old = colMeans(Old)[1:3], new = colMeans(New)[1:3])
#     user.self sys.self elapsed
# old     0.353    2.281   2.683
# new     0.320    1.947   2.272
## new code is a bit faster with every formula
}

# fit variance matrix

(bench <- microbenchmark(
  new.R = (V <- fitDistVar_R(Dist, r.est$spatialcor, fun = "exp-pwr")),
  old.R = (tmp <- V.fit(Dist, r.est$spatialcor, FUN = "exponential-power")),
  times = 1L))
stopifnot(all.equal(V, tmp))
# ggplot2::autoplot(bench)

# decomp V
(bench <- microbenchmark(
  inv.R = (tmp <- t(backsolve(chol(V), diag(nrow(V))))),
  inv.Rfunc = (invcholV <- invert_cholR(V)),times = 1L))
stopifnot(all.equal(invcholV, tmp))
# ggplot2::autoplot(bench)

(AA <- crossprod(matrix(1:6, ncol = 2)))
cM <- chol(AA) # cholesky matrix (Upper)
crossprod(cM) # t(cm) %*% cm = original matrix


inv.A <- t(backsolve(cM, diag(2))) # this is the inverse OF the chol matrix
t(inv.A) %*% AA
t(inv.A) %*% cM

inv.B <- chol2inv(cM) # this is the inverse of the original matrix through chol
round(inv.B %*% AA, digits = 4)

inv.C <- solve(AA, diag(2))
inv.C %*% AA

if(FALSE){
  solveone <- function(X, y, V){
    invcholV <- invert_cholR(V)
    xx <- invcholV %*% X
    yy <- invcholV %*% y
    coef <- as.numeric(solve(crossprod(xx), crossprod(xx,yy)))
  }

  solvetwo <- function(X, y, V){
    invV <- chol2inv(chol(V))
    XVX <- t(X) %*% invV %*% X
    XVY <- t(X) %*% invV %*% y
    coef <- as.numeric(solve(XVX, XVY))
  }

  y <- rnorm(nrow(Xmat))

  (microbenchmark(
    OG = (beta <- solveone(Xmat, y, V)),
    new = (beta.new <- solvetwo(Xmat, y , V)),
    times = 1L))
  ## solveone is WAY faster.
  stopifnot(all.equal(beta, beta.new))
}


# fit GLS partiion

source('R/fitGLS.R', echo=TRUE)
source("R/fitSpatCor.R")
load(file = "R/vignettes-and-examples/test-gls.rda")
data <- readRDS("data/alaska_ndvi_3000pts.rda")

n.p = 50; npart = 2
partition <- matrix(sample(1:nrow(data), size = n.p * npart), ncol = npart)


## For now, I'm going to explore the distributed computing option

### partition 1
y1 <- rnorm(n.p)
X1 <- as.matrix(data[partition[,1], -c(1:6)])
loc1 <- data[partition[, 1], c("lat","lng")]
V1 <- V.fit(Dist[partition[,1], partition[,1]],
            spatialcor = r.est$spatialcor, FUN = "exponential-power")
### partition 2
y2 <- rnorm(n.p)
X2 <- as.matrix(data[partition[,2], -c(1:6)])
loc2 <- data[partition[, 2], c("lat","lng")]
V2 <- V.fit(Dist[partition[,1], partition[,1]],
            spatialcor = r.est$spatialcor, FUN = "exponential-power")

V12 <- V.fit(Dist[as.vector(partition), as.vector(partition)],
             spatialcor = r.est$spatialcor, FUN = "exponential-power")
Xnull <- matrix(1, nrow = nrow(X1))

df2 <- n.p - (ncol(X1) - 1)
df0 <- n.p - (ncol(Xnull) - 1)
df1 <- df0 - df2

out1 <- GLS_worker_cpp(y1, X1, V1, Xnull, save_xx = TRUE)
out2 <- GLS_worker_cpp(y2, X2, V2, Xnull, save_xx = TRUE)

out.cross <- crosspart_worker_cpp(out1, out2, V12, df1, df2) # not working yet

## Now I'll do the whole thing..
fitGLS.partition_rcpp(X = X, y = rnorm(nrow(X)), X0 = matrix(1, nrow = nrow(X)), Dist = Dist, spatcor = r.est$spatialcor, nugget = 0)
# fitGLS.partition(X = Xmat, V = V, y = rnorm(nrow(Xmat)),
#                  X0 = matrix(1, nrow = nrow(Xmat)), nugget = 0,
#                  npart = 2, mincross = 2)

# fit nugget (re-fit GLS)

