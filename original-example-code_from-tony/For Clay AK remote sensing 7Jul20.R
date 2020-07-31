library(Matrix)
library(geosphere)
library(colorspace)

source('../original-example-code_from-tony/remote_sensing_tools_8Jul20.R')

###########################################################
# just AK(ish)
###########################################################

data <- read.csv(file="data-raw/north_america_checked.csv")
summary(data)
sort(unique(data$land))
n.obs <- 32

dataAK <- data[data$lng < -141,]

# make dataset smaller
dataAK <- dataAK[sample.int(size=3000, n = nrow(dataAK)),]

summary(dataAK)
nrow(data)
nrow(dataAK)
sort(unique(dataAK$land))

# look for outliers
X <- as.matrix(dataAK[,6:37])
hist(X)

# landclasses: 1 Evergreen needleleaf forests 2 evergreen broadleaf forests, 3 deciduous needleleaf forests, 4 deciduous broadleaf forests, 5 mixed forests, 6 shrublands, 8 savannas, 10 grasslands, 12 croplands, 14 croplands/natural vegetation mosaics.

land.df <- data.frame(class = c("Evergr needle","Evergr broad","Decid needle","Decid broad","Mixed forest","Shrubland","Savanna","Grassland","Cropland","Cropland mosaics"), num = c(1,2,3,4,5,6,8,10,12,14))

dataAK$landclass <- land.df$class[match(dataAK$land, land.df$num)]
dataAK$landclass <- as.factor(dataAK$landclass)
n.classes <- aggregate(dataAK$landclass, by=list(landscape=dataAK$landclass), FUN=length)
rare.class.threshold <- 0.02 * nrow(dataAK)
rare.classes <- n.classes$landscape[n.classes$x <= rare.class.threshold]
dataAK <- dataAK[!is.element(dataAK$landclass, rare.classes),]
dataAK$landclass <- droplevels(dataAK$landclass)
levels(dataAK$landclass)
dim(dataAK)

# remove
landclasses <- c("Shrubland","Savanna","Grassland")

dataAK$landclass.num <- 1
for(i in 1:length(landclasses)) dataAK$landclass.num[dataAK$landclass == landclasses[i]] <- i

n.classes <- aggregate(dataAK$landclass, by=list(landscape=dataAK$landclass), FUN=length)
n.classes
  # landscape     x
# 1 Grassland  6046
# 2   Savanna 12305
# 3 Shrubland 12514

###########################################################
# AK outliers
###########################################################
dat <- dataAK
dim(dat)
# [1] 31486  37

# set up a scaled time variable
t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

# fit CLS
X <- as.matrix(dat[,6:37])
dat.map <- CLS.fit(X, t.scale)

# this adds the landscape variable to dat.map
dat.map$landclass <- dat$landclass
dat.map$landclass <- droplevels(dat.map$landclass)
dat.map$lat <- dat$lat
dat.map$lng <- dat$lng
dat.map$c.cls <- dat.map$c
dat.map$c <- dat.map$c.cls/(1-dat.map$b)
dat.map$rel.c.cls <- dat.map$c.cls/dat.map$mean
dat.map$rel.c <- dat.map$c/dat.map$mean
dat.map$rel.MSE <- dat.map$MSE/dat.map$mean^2
dat.map$rel.MSE.25 <- dat.map$rel.MSE^.25

# look for outliers in rel.c.cls
hist(scale(dat.map$rel.c.cls))
min(scale(dat.map$rel.c.cls))
max(scale(dat.map$rel.c.cls))

# look for outliers in c.cls
hist(scale(dat.map$c.cls))
min(scale(dat.map$c.cls))
max(scale(dat.map$c.cls))

nrow(dat.map)
q <- -qnorm(p = 1/nrow(dat.map)/10)
dat.map$outlier <- F
dat.map$outlier[abs(scale(dat.map$rel.c.cls)) > q] <- T
sum(dat.map$outlier)
mean(dat.map$outlier)

# plot with outliers
col.pal.trend <- hcl.colors(100, "RdBu", rev = TRUE)
palette(col.pal.trend)
dat.map$col <- dat.map$rel.c.cls/(max(dat.map$rel.c.cls) - min(dat.map$rel.c.cls))
minmax <- max(abs(min(dat.map$col)), max(dat.map$col))
dat.map$col <- .5*(1 + dat.map$col/minmax)

par(mfrow = c(1,1))
plot(lat ~ lng, data=dat.map, pch=15, cex=.7, col=col.pal.trend[100*dat.map$col], xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n")
points(lat ~ lng, data=dat.map[dat.map$outlier,], pch=15, cex=.7, col="black")

dat.map$rel.c.cls[dat.map$outlier]

hist(dat.map$mean[dat.map$outlier], freq=F, ylim=c(0,.5), breaks=0:20)
hist(dat.map$mean, freq=F, add=T, col=NULL, breaks=0:20)

##########
# remove outliers
dim(dat.map)
dat.map <- dat.map[!dat.map$outlier,]
dim(dat.map)

###########################################################
# maps
if(F){
	# plot map of data
	minlng <- min(dataAK$lng)
	maxlng <- max(dataAK$lng)

	col.pal.land <- terrain.colors(8)
	col.pal.trend <- hcl.colors(100, "RdBu", rev = TRUE)

	dat.map$cls.sig <- F
	dat.map$cls.sig[dat.map$p < 0.05] <- T


	par(mfrow = c(1,2), mai = c(0,0,.5,0))

	# A
	palette(col.pal.land)
	plot(lat ~ lng, data=dataAK, pch=15, cex=.4, col=col.pal.land[landclass.num+2], xlim=c(minlng, maxlng), xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n")
	legend(x=-180, y=60, legend=landclasses, col=col.pal.land[-(1:2)], pch=20, bty="n", cex = 1.2)

	dat.map$cls.sig <- F
	dat.map$cls.sig[dat.map$p < 0.05 & dat.map$rel.c.cls > 0] <- T
	points(lat ~ lng, data=dat.map[dat.map$cls.sig == T,], pch=20, cex=.15, col="red")

	dat.map$cls.sig <- F
	dat.map$cls.sig[dat.map$p < 0.05 & dat.map$rel.c.cls < 0] <- T
	points(lat ~ lng, data=dat.map[dat.map$cls.sig == T,], pch=20, cex=.15, col="blue")

	mtext("A", side = 3, adj=1, cex = 2)

	# B
	palette(col.pal.trend)
	dat.map$col <- dat.map$rel.c.cls/(max(dat.map$rel.c.cls) - min(dat.map$rel.c.cls))
	minmax <- max(abs(min(dat.map$col)), max(dat.map$col))
	dat.map$col <- .5*(1 + dat.map$col/minmax)

	plot(lat ~ lng, data=dat.map, pch=15, cex=.4, col=col.pal.trend[100*dat.map$col], xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n")
	mtext("B", side = 3, adj=1, cex = 2)

}

###########################################################
# Analysis of full map
###########################################################
fit.n.sample <- 200
par(mfrow=c(1,1))
r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=fit.n.sample,
                             plot.fig=T, FUN="exponential-power")
# r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=fit.n.sample,
#                              plot.fig=T, FUN="exponential")
r.est

location <- dat.map[,c('lng','lat')]
Dist <- distm(location, fun=distGeo)/1000
V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential-power")
lsos()

# estimating the nugget
system.time(nugget <- nugget.fit(formula='rel.c.cls ~ 0 + landclass', dat.map, V, nugget.tol = 0.0001, verbose = T))
Vn <- (1 - nugget) * V + nugget * diag(nrow(V))

system.time(invcholVn <- t(backsolve(chol(Vn), diag(nrow(Vn)))))

rm(Vn)
system.time(z.rel.cls <- GLS.fit(rel.c.cls ~ 0 + landclass, data=dat.map, invcholV=invcholVn))
dim(dat.map)
nugget
z.rel.cls$coef
z.rel.cls$p.F
z.rel.cls$logLik


system.time(z.rel.cls <- GLS.fit(rel.c.cls ~ lat, data=dat.map, invcholV=invcholVn))
dim(dat.map)
nugget
z.rel.cls$coef
z.rel.cls$p.F
z.rel.cls$logLik


###########################################################
# Partition analysis npart
###########################################################
r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=200, plot.fig=T, FUN="exponential-power")
r.est

for(npart in c(10, 5)){
	min.num.rSS <- 10

	z.lat <- GLS.partition.data(rel.c.cls ~ 1 + lat, formula0=rel.c.cls ~ 1, data=dat.map, spatial.autocor.FUN = "exponential-power", spatialcor=r.est$spatialcor, est.nugget=T, npart=npart, min.num.rSS = min.num.rSS, verbose = F)
	pvalue.lat <- GLS.partition.pvalue(z.lat, nboot = 10^5)
	show(z.lat$coef)
	show(z.lat$nugget)
	show(z.lat$nugget.part)
	show(unlist(pvalue.lat))

	formula = list(rel.c.cls ~ 0 + landclass, rel.c.cls ~ 1 + lat, rel.c.cls ~ 0 + lat*landclass, rel.c.cls ~ 0 + lat*landclass, rel.c.cls ~ 0 + lat*landclass)
	formula0 = list(rel.c.cls ~ 1, rel.c.cls ~ 1, rel.c.cls ~ 1, rel.c.cls ~ 0 + lat + landclass, rel.c.cls ~ 0 + landclass)
	z <- GLS.partition.data.multiformula(formula = formula, formula0 = formula0, data=dat.map, spatial.autocor.FUN = "exponential-power", spatialcor=r.est$spatialcor, est.nugget=T, npart=npart, min.num.rSS = min.num.rSS, verbose = T, nugget.interval = c(.1,.5))
	pvalue <- list()
	for(i in 1:length(formula)){
		pvalue[[i]] <- GLS.partition.pvalue(z[[i]], nboot = 10^5)
	}
	for(i in 1:length(formula)){
		show(formula[i])
		show(formula0[i])
		show(z[[i]]$coef)
		show(z[[i]]$nugget)
		show(z[[i]]$nugget.part)
		print(unlist(pvalue[[i]]))
	}

}


