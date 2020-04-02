library(Matrix)
library(LatticeKrig)
library(geosphere)
library(colorspace)

# source('remote_sensing_tools_24Mar20.R')
source("R/remote-sensing-functions.R")

#### Base analysis for a subset for WI data

## Setup ----
# data <- read.csv(file="north_america_checked.csv")
# summary(data)

# # to create a play dataset, isolate WI
# wi.data <- data[(data$lng > -92) & (data$lng < -87) & (data$lat > 42) &
#                   (data$lat < 47),]
# write.table(file="wisconsin_checked.csv", wi.data, sep=",", row.names=F)

# data <- read.csv(file="wisconsin_checked.csv")
load("data/wisconsin.example.rda", verbose = TRUE)
data <- wisconsin.example
summary(data)

sort(unique(data$land))
n.obs <- 32
t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

## Construct landcover classes ----

## landclasses: 1 Evergreen needleleaf forests 2 evergreen broadleaf forests,
##   3 deciduous needleleaf forests, 4 deciduous broadleaf forests, 5 mixed
##   forests, 6 shrublands, 8 savannas, 10 grasslands, 12 croplands, 14
##   croplands/natural vegetation mosaics.

# Do some basic data construction to sort landcover classes
landclasses <- c("Evergr needle","Evergr broad","Decid needle","Decid broad",
                 "Mixed forest","Shrubland","Savanna","Grassland","Cropland",
                 "Cropland mosaics")
data$landclass <- landclasses[1]
count <- 0
for(i in sort(unique(data$land))) {
  count <- count+1
  data$landclass[data$land==i] <- landclasses[count]
}
data$landclass <- as.factor(data$landclass)
n.classes <- aggregate(data$landclass, by=list(landscape=data$landclass),
                       FUN=length)
rare.class.threshold <- 0.005 * nrow(data)
rare.classes <- n.classes$landscape[n.classes$x <= rare.class.threshold]
data <- data[!is.element(data$landclass, rare.classes),]
data$landclass <- droplevels(data$landclass)
levels(data$landclass)

landclasses <- c("Evergr needle","Evergr broad","Decid broad","Mixed forest",
                 "Shrubland","Savanna","Grassland","Cropland")

data$landclass.num <- 1
for(i in 1:length(landclasses)) {
  data$landclass.num[data$landclass == landclasses[i]] <- i
}

n.classes <- aggregate(data$landclass, by=list(landscape=data$landclass),
                       FUN=length)
n.classes
## landscape    x
## 1  Decid needle  755
## 2  Evergr broad   78
## 3 Evergr needle  781
## 4  Mixed forest 1055
## 5     Shrubland   51

## plot map of data ----
minlng <- min(data$lng)
maxlng <- max(data$lng)

col.pal <- terrain.colors(nlevels(data$landclass))
plot(lat ~ lng, data=data, pch=20, cex=.1, col=col.pal[landclass.num],
     xlim=c(minlng, maxlng), xlab="", ylab="")
legend(x=-180, y=45, legend=landclasses, col=col.pal, pch=20, bty="n")

## perform CLS ----

# create distance matrix in kilometers
location <- data[,c('lng','lat')]
Dist <- distm(location, fun=distGeo)/1000

# fit CLS
X <- as.matrix(data[,6:37])
dat.map <- CLS.fit(X, t.scale)
n <- nrow(dat.map)

# this adds the landscape variable to dat.map
dat.map$landclass <- data$landclass
dat.map$landclass <- droplevels(dat.map$landclass)
dat.map$lat <- data$lat
dat.map$lng <- data$lng
dat.map$c.cls <- dat.map$c
dat.map$c <- dat.map$c.cls/(1-dat.map$b)
dat.map$rel.c.cls <- dat.map$c.cls/dat.map$mean
dat.map$abslat <- abs(dat.map$lat)

# fit the spatial correlation for a random subset of fit.n.sample points
fit.n.sample <- 2000
r.est <- spatialcor.fit.data(X, t.scale, data=data, fit.n.sample=fit.n.sample,
                             r.start=.1, plot.fig=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## perform the full analysis using GLS.fit ----
# NOTE: GLS.fit is a base-level function that requires fitting the nugget
#   separately and constructing the covariance matrix
V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")

# construct the GLS correlation matrix and fit the GLS
# absolute c.cls
nugget <- nugget.fit(formula='rel.c.cls ~ 0 + landclass', dat.map, V,
                     nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.GLS.fit <- GLS.fit(rel.c.cls ~ 0 + landclass, data=dat.map,
                     invcholV=invcholVn)
names(z.GLS.fit)
as.data.frame(z.GLS.fit[c("coef","se","t","df.t","p.t")])
as.data.frame(z.GLS.fit[c("F","df1.F","df2.F","p.F","logLik","logLik0")])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## perform the analysis by partitioning data using GLS.partition.data ----

# number of partitions is npart
npart <- 4
z.GLS.partition <- GLS.partition.data(rel.c.cls ~ 0 + landclass,
                                      formula0=rel.c.cls ~ 1, data=dat.map,
                                      r=r.est$spatialcor, est.nugget=T,
                                      npart=npart)
names(z.GLS.partition)
z.pvalue <- GLS.partition.pvalue(z.GLS.partition, nboot = 10^5)
as.data.frame(z.GLS.partition[c("coef","Fmean", "df1", "df2")])
as.data.frame(z.GLS.partition[c("F.part", "p.F.part")])
z.pvalue
