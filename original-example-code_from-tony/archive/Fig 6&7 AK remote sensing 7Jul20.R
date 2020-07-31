setwd("D:beautydata/arives/RSE ms/")

library(Matrix)
library(geosphere)
library(colorspace)

source('remote_sensing_tools_6Jul20.R')

###########################################################
# just AK(ish)
###########################################################

data <- read.csv(file="north_america_checked.csv")
summary(data)
sort(unique(data$land))
n.obs <- 32

dataAK <- data[data$lng < -141,]

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

outliers <- dat.map[dat.map$outlier, c("lng","lat")]
dat$outlier <- F
for(i.outlier in 1:nrow(outliers)) dat$outlier[is.element(dat$lng, outlier$lng[i.outlier]) & is.element(dat$lat, outlier$lat[i.outlier])] <- T
sum(dat$outlier)

palette(rainbow(100))
matplot(t(dat[dat$outlier,6:37]), typ="l", lty=1, col=1:sum(dat.map$outlier), main="Outliers", lwd=2)

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

	pdf(file="Fig 6 AK 7Jul20.pdf", width=12, height=7)
	
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
		
	dev.off()
}

###########################################################
# Analysis of full map
###########################################################
fit.n.sample <- 2000
par(mfrow=c(1,1))
r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=fit.n.sample, plot.fig=T, FUN="exponential-power")
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
saveRDS(file=paste0("z full land 7Jul20.RDS"), z.rel.cls)
dim(dat.map)
nugget
z.rel.cls$coef
z.rel.cls$p.F
z.rel.cls$logLik


system.time(z.rel.cls <- GLS.fit(rel.c.cls ~ lat, data=dat.map, invcholV=invcholVn))
saveRDS(file=paste0("z full lat 7Jul20.RDS"), z.rel.cls)
dim(dat.map)
nugget
z.rel.cls$coef
z.rel.cls$p.F
z.rel.cls$logLik


###########################################################
# Partition analysis npart
###########################################################
r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=2000, plot.fig=T, FUN="exponential-power")
r.est

for(npart in c(20, 15, 10, 5)){
	min.num.rSS <- 10
	
	formula = list(rel.c.cls ~ 0 + landclass, rel.c.cls ~ 1 + lat, rel.c.cls ~ 0 + lat*landclass, rel.c.cls ~ 0 + lat*landclass, rel.c.cls ~ 0 + lat*landclass)
	formula0 = list(rel.c.cls ~ 1, rel.c.cls ~ 1, rel.c.cls ~ 1, rel.c.cls ~ 0 + lat + landclass, rel.c.cls ~ 0 + landclass)
	z <- GLS.partition.data.multiformula(formula = formula, formula0 = formula0, data=dat.map, spatial.autocor.FUN = "exponential-power", spatialcor=r.est$spatialcor, est.nugget=T, npart=npart, min.num.rSS = min.num.rSS, verbose = T, nugget.interval = c(.1,.5))
	pvalue <- list()
	for(i in 1:length(formula)){
		pvalue[[i]] <- GLS.partition.pvalue(z[[i]], nboot = 10^5)
	}
	saveRDS(file=paste0("z npart=",npart," 7Jul20.RDS"), z)
	saveRDS(file=paste0("pvalue npart=",npart," 7Jul20.RDS"), pvalue)
	for(i in 1:length(formula)){
		show(formula[i])
		show(formula0[i])
		show(z[[i]]$coef)
		show(z[[i]]$nugget)
		show(z[[i]]$nugget.part)
		print(unlist(pvalue[[i]]))
	}
	
	z.lat <- GLS.partition.data(rel.c.cls ~ 1 + lat, formula0=rel.c.cls ~ 1, data=dat.map, spatial.autocor.FUN = "exponential-power", spatialcor=r.est$spatialcor, est.nugget=T, npart=npart, min.num.rSS = min.num.rSS, verbose = F)	
	pvalue.lat <- GLS.partition.pvalue(z.lat, nboot = 10^5)
	saveRDS(file=paste0("z.lat npart=",npart," 7Jul20.RDS"), z.lat)
	saveRDS(file=paste0("pvalue.lat npart=",npart," 7Jul20.RDS"), pvalue.lat)
	show(z.lat$coef)
	show(z.lat$nugget)
	show(z.lat$nugget.part)
	show(unlist(pvalue.lat))
}


###########################################################
###########################################################
# Upload files
###########################################################
###########################################################

z <- readRDS(file="z 2Jul20.RDS")

pvalue <- list()
for(i in 1:3){
	pvalue[[i]] <- GLS.partition.pvalue(z[[i]], nboot = 10^5)
}
for(i in 1:3){
	show(z[[i]]$coef)
	show(z[[i]]$nugget)
	show(z[[i]]$nugget.part)
	print(unlist(pvalue[[i]]))
}


###########################################################
# non-spatial analysis
###########################################################
lm.rel.c.cls <- lm(rel.c.cls ~ 0 + landclass, data=dat.map)
summary(lm.rel.c.cls)
lm.rel.c.cls0 <- lm(rel.c.cls ~ 1, data=dat.map)
summary(lm.rel.c.cls0)
anova(lm.rel.c.cls, lm.rel.c.cls0)


# Proportion of trends that are positive
mean(dat.map$c.cls[dat.map$p < 0.05] > 0)

# Histograms of c, c.cls and rel.c.cls
colMeans(dat.map[,2:7])
     # mean         c         t         p         b       MSE 
# 9.2988495 0.1177657 0.4006831 0.2593884 0.1890124 0.2012181 

# Values of b among land classes
summary(lm(b ~ 0 + landclass, data=dat.map))

par(mfrow=c(1,2))
hist(dat.map$c[dat.map$c>-10 & dat.map$c < 10], breaks=.5*(-20:20), main="", xlab="c and c.cls", ylim = c(0,5000))
hist(dat.map$c.cls[dat.map$c.cls>-10 & dat.map$c.cls < 10], add=T, col="red", breaks=.5*(-20:20))
hist(dat.map$rel.c.cls[dat.map$rel.c.cls > -10 & dat.map$rel.c.cls < 10], add=T, col="green", breaks=.5*(-20:20))

# check distribution of MSE
hist(dat.map$rel.MSE.25, breaks=40)

###########################################################
# Fig. 7: plot results
###########################################################

# Fig. 7: plot results
pdf("Fig 7 remote sensing rel.cls 19Nov19.pdf", height=6, width=6)
z.cls <- readRDS(file="Fig7 north_america rel.cls 19Nov19.RDS")

par(mfrow=c(2,1), mai=c(.1,1,.2,.4))

# relative c.cls
S <- z.cls$varcov
p <- length(z.cls$coef)
se.cond <- array(p)
for(i in 1:p){
	S11 <- S[i,i]
	S22 <- S[-i,-i]
	S12 <- S[i,-i]
	se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
}


arg <- 1:length(z.cls$coef)
ymax <- max(z.cls$coef + z.cls$se+.005)
ymin <- min(z.cls$coef - z.cls$se-.005)
plot(arg+.05, z.cls$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", pch=15, ylab=expression(paste("Slope ", italic(c[cls]),"/mean NDVI")))
arrows(x0=arg, y0=z.cls$coef - z.cls$se, y1=z.cls$coef + z.cls$se, angle=90, code=3, length=.05)
arrows(x0=arg+.1, y0=z.cls$coef - se.cond, y1=z.cls$coef + se.cond, angle=90, code=3, length=.05, col="blue")
lines(c(0,20), c(0,0), lty=2)
mtext("A", at=9, cex=1.2)

n.data <- aggregate(dat.map$rel.c.cls, by=list(landscape=dat.map$landclass), FUN=length)[,2]
lm.est <- lm(rel.c.cls ~ 0 + landclass, data=dat.map)
c.data <- summary(lm.est)$coef[,1]
se.data <- summary(lm.est)$coef[,2]
points(arg + .3, c.data, col="red")
arrows(x0=arg + .3, y0=c.data - se.data, y1=c.data + se.data, angle=90, code=3, length=.05, col="red")

par(mai=c(1.55,1,.1,.3))

p.LS.est <- aggregate(dat.map$p.LS, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.LS.est, ylim=c(0,1), xlab="", xlim=c(0, 10), ylab=expression(paste(italic(P), "-values")), col="white")
p.est <- aggregate(dat.map$p, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.est, ylim=c(0,1), xlab="", xlim=c(0, 11), add=T, col="black")

at.min <- .8
at.max <- max(arg) + 2.35
at <- at.min + (at.max - at.min)*(0:(length(arg)-1))/length(arg)
axis(side=1, at=at, las=3, levels(dat.map$landclass))
lines(c(0,11), c(.05,.05), lty=1, col="red")
text(x=at, y=p.LS.est+.15, n.data, srt=90, cex=.9)
mtext("B", at=at.max, cex=1.2)

dev.off()

######################################
# Diagnostics
# Fig. S2: plot results
pdf("Fig 7 (S2) remote sensing b 19Nov19.pdf", height=6, width=6)
z.cls <- readRDS(file="Fig7 north_america b 19Nov19.RDS")

par(mfrow=c(1,1), mai=c(.1,1,.2,.4))

# relative c.cls
S <- z.cls$varcov
p <- length(z.cls$coef)
se.cond <- array(p)
for(i in 1:p){
	S11 <- S[i,i]
	S22 <- S[-i,-i]
	S12 <- S[i,-i]
	se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
}

arg <- 1:length(z.cls$coef)
ymax <- max(z.cls$coef + z.cls$se+.005)
ymin <- min(z.cls$coef - z.cls$se-.005)
plot(arg+.05, z.cls$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", pch=15, ylab=expression(paste("Slope ", italic(c[cls]),"/mean NDVI")))
arrows(x0=arg, y0=z.cls$coef - z.cls$se, y1=z.cls$coef + z.cls$se, angle=90, code=3, length=.05)
arrows(x0=arg+.1, y0=z.cls$coef - se.cond, y1=z.cls$coef + se.cond, angle=90, code=3, length=.05, col="blue")
lines(c(0,20), c(0,0), lty=2)

n.data <- aggregate(dat.map$rel.c.cls, by=list(landscape=dat.map$landclass), FUN=length)[,2]
lm.est <- lm(rel.c.cls ~ 0 + landclass, data=dat.map)
c.data <- summary(lm.est)$coef[,1]
se.data <- summary(lm.est)$coef[,2]
points(arg + .3, c.data, col="red")
arrows(x0=arg + .3, y0=c.data - se.data, y1=c.data + se.data, angle=90, code=3, length=.05, col="red")

dev.off()



##############################
# Diagnostic maps of b and MSE

plot(rel.MSE^.5 ~ as.numeric(landclass), data=dat.map)
plot(order(rel.MSE) ~ landclass, data=dat.map)

par(mfrow=c(2,4))
for(i in levels(dat.map$landclass)) hist(dat.map$MSE[dat.map$landclass == i]^.5, main=i, breaks=40)

pdf(file="FigS2(7) b and cls 19Nov19.pdf", width=10, height=10)

par(mfcol=c(2,2), mai=c(.8,.8,.4,.1))

minlng <- min(dat.map$lng)
maxlng <- -50
col.pal <- rainbow(80)

z <- dat.map$b
col.x <- 40 + 40*(z - min(z))/(max(z) - min(z))
plot(lat ~ lng, data=dat.map, pch=15, cex=.35, col=col.pal[col.x], xlim=c(minlng, maxlng), xlab="", ylab="", main=expression(italic(b)), cex.main=2)
legend(x=-170, y=25, legend=c(paste0("min = ",round(min(z), digits=2)), paste0("max = ",round(max(z), digits=2))), col=col.pal[c(41,80)], bty="n", pch=15, cex=1.5)

par(mai=c(1.4,.8,.1,.1))
zz <- aggregate(b ~ landclass, data=dat.map, FUN=mean)
zz$sd <- aggregate(b ~ landclass, data=dat.map, FUN=sd)[,2]
arg <- 1:nrow(zz)
plot(b ~ landclass, data=zz, xlab="", xaxt="n", ylim=c(-.2,.6), ylab=expression(italic(b)), cex.lab=1.5)
arrows(x0=arg, y0=zz[,2] - zz[,3], y1=zz[,2] + zz[,3], angle=90, code=3, length=.05)
at.min <- 1
at.max <- max(arg)+1
at <- at.min + (at.max - at.min)*(0:(length(arg)-1))/length(arg)
axis(side=1, at=at, las=3, levels(dat.map$landclass), cex=2)


par(mai=c(.8,.8,.4,.1))
plot(rel.c ~ rel.c.cls, data=dat.map, pch=19, xlim=c(-2,4.5), ylim=c(-2,4.5), xlab=expression(italic(c[cls.rel])), ylab=expression(italic(c[rel])), cex.lab=1.5)

par(mai=c(1.4,.8,.1,.1))
zz <- aggregate(rel.c.cls ~ landclass, data=dat.map, FUN=mean)
zz$sd <- aggregate(rel.c.cls ~ landclass, data=dat.map, FUN=sd)[,2]
plot(rel.c.cls ~ landclass, data=zz, xlab="", xaxt="n", ylim=c(-.2,.2), ylab=expression(italic(c[cls.rel])), cex.lab=1.5)
arrows(x0=arg, y0=zz[,2] - zz[,3], y1=zz[,2] + zz[,3], angle=90, code=3, length=.05)
axis(side=1, at=at, las=3, levels(dat.map$landclass), cex=2)

dev.off()


# pdf(file="Fig7 b and MSE 19Nov19.pdf", width=10, height=10)
# par(mfcol=c(2,2))

# minlng <- min(dat.map$lng)
# maxlng <- -50

# col.pal <- rainbow(80)
# z <- dat.map$b
# col.x <- 40 + 40*(z - min(z))/(max(z) - min(z))
# plot(lat ~ lng, data=dat.map, pch=15, cex=.35, col=col.pal[col.x], xlim=c(minlng, maxlng), xlab="", ylab="", main="Autocorrelation b")
# legend(x=-170, y=25, legend=c(paste0("min = ",round(min(z), digits=2)), paste0("max = ",round(max(z), digits=2))), col=col.pal[c(41,80)], bty="n", pch=15, cex=1.5)

# hist(z, breaks=40, xlab="b", main="")


# # standardize MSE for rel.c.cls
# #z <- dat.map$MSE^.5/dat.map$mean
# z <- dat.map$MSE
# col.x <- 40 + 40*(z^.5 - min(z^.5))/(max(z^.5) - min(z^.5))
# plot(lat ~ lng, data=dat.map, pch=15, cex=.35, col=col.pal[col.x], xlim=c(minlng, maxlng), xlab="", ylab="", main="MSE")
# legend(x=-180, y=25, legend=c(paste0("min = ",round(min(z), digits=2)), paste0("max = ",round(max(z), digits=2))), col=col.pal[c(41,80)], bty="n", pch=15, cex=1.5)

# hist(z[z<1.5], breaks=40, xlab="MSE", main="", xlim=c(0,1.5))

# dev.off()


# This estimates r. It is set up now to subsample pixels (r.fit.n.smaple) and it re-computes the CLS values. 
fit.n.sample <- 2000
r.est <- spatialcor.fit.data(X, t.scale, data=dat, fit.n.sample=fit.n.sample, r.start=10, plot.fig=T)
r.est
# $spatialcor
       # r 
# 550.4929 

# $spatialcor.sigma
# [1] 0.2195705
V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")

# construct the GLS correlation matrix and fit the GLS
# absolute c.cls
nugget <- nugget.fit(formula='c.cls ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget
#[1] 0.3129171

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.cls <- GLS.fit(c.cls ~ 0 + landclass, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america cls 19Nov19.RDS", z.cls)


# relative c.cls
nugget <- nugget.fit(formula='rel.c.cls ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget
# [1] 0.4136469

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.rel.cls <- GLS.fit(rel.c.cls ~ 0 + landclass, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america rel.cls 19Nov19.RDS", z.rel.cls)

# b
nugget <- nugget.fit(formula='b ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.b <- GLS.fit(b ~ 0 + landclass, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america b 19Nov19.RDS", z.b)

# MSE
nugget <- nugget.fit(formula='MSE ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget
# [1] 0.1252042

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.MSE <- GLS.fit(MSE ~ 0 + landclass, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america MSE 19Nov19.RDS", z.MSE)

# rel.MSE.25
nugget <- nugget.fit(formula='rel.MSE.25 ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget
#[1] 0.1394387

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.rel.MSE.25 <- GLS.fit(rel.MSE.25 ~ 0 + landclass, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america rel.MSE.25 19Nov19.RDS", z.rel.MSE.25)




# contrast test for Ever needle vs. grassland
dat.map$landclass0 <- dat.map$landclass
dat.map$landclass0[dat.map$landclass0 == "Grassland"] <- "Evergr needle"
dat.map$landclass0 <- droplevels(dat.map$landclass0)
nugget <- nugget.fit(formula='rel.c.cls ~ 0 + landclass', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget
#[1] 0.4136469

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.rel.cls0 <- GLS.fit(rel.c.cls ~ 0 + landclass, formula0=rel.c.cls ~ 0 + landclass0, data=dat.map, invcholV=invcholVn)
saveRDS(file="Fig7 north_america rel.cls.Evergr needle.Grassland 19Nov19.RDS", z.rel.cls0)


################################################################
z.cls <- readRDS(file="Fig7 north_america cls 19Nov19.RDS")
z.rel.cls <- readRDS(file="Fig7 north_america rel.cls 19Nov19.RDS")
z.b <- readRDS(file="Fig7 north_america b 19Nov19.RDS")
z.MSE <- readRDS(file="Fig7 north_america MSE 19Nov19.RDS")
z.rel.MSE.25 <- readRDS(file="Fig7 north_america rel.MSE.25 19Nov19.RDS")
z.rel.cls.Evergrneedle.Grassland <- readRDS(file="Fig7 north_america rel.cls.Evergr needle.Grassland 19Nov19.RDS")

z.cls$p.t
z.cls$F
z.cls$p.F
z.cls$coef0/z.cls$se0
2*pt(abs(z.cls$coef0/z.cls$se0), df=10000, lower.tail=F)
     # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest       landclassSavanna 
            # 0.11359198             0.10246940             0.45051146             0.47412972             0.11531029             0.06196635             0.34542443 
    # landclassShrubland 
            # 0.04413186 
# [1] 5.533894
# [1] 2.257141e-06
# [1,] 1.326737
# [1,] 0.184626

z.rel.cls$p.t
z.rel.cls$F
z.rel.cls$p.F
z.rel.cls$coef0/z.rel.cls$se0
2*pt(abs(z.rel.cls$coef0/z.rel.cls$se0), df=10000, lower.tail=F)


    # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest       landclassSavanna 
            # 0.10642832             0.24608554             0.53077754             0.52099238             0.06986644             0.22262794             0.43976299 
    # landclassShrubland 
            # 0.03451549 
# [1] 5.079996
# [1] 8.98405e-06
# [1,] 1.317407
# [1,] 0.1877323

z.b$coef
z.b$p.t
z.b$F
z.b$p.F
z.b$coef0/z.b$se0
2*pt(abs(z.b$coef0/z.b$se0), df=10000, lower.tail=F)
     # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # 0.1902902              0.1960579              0.2366679              0.2078597              0.2042857              0.2409578 
      # landclassSavanna     landclassShrubland 
             # 0.2164642              0.2123562 

# [1] 1.730487
# [1] 0.09702765
# [1,] 4.848346
# [1,] 1.263699e-06

z.rel.MSE.25$coef
z.rel.MSE.25$p.t
z.rel.MSE.25$F
z.rel.MSE.25$p.F
z.rel.MSE.25$coef0/z.b$se0
2*pt(abs(z.rel.MSE.25$coef0/z.rel.MSE.25$se0), df=10000, lower.tail=F)

     # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # 0.2545895              0.2434556              0.2441044              0.2410735              0.2694370              0.2416177 
      # landclassSavanna     landclassShrubland 
             # 0.2446680              0.2585908 
# [1] 31.3673
# [1] 2.532423e-43
# [1,] 5.838165
# [1,] 2.180285e-82

z.MSE$coef
z.MSE$p.t
z.MSE$F
z.MSE$p.F
z.MSE$coef0/z.b$se0
2*pt(abs(z.MSE$coef0/z.MSE$se0), df=10000, lower.tail=F)
     # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # 0.2727545              0.2980761              0.3324454              0.3287236              0.3014445              0.3073751 
      # landclassSavanna     landclassShrubland 
             # 0.2986511              0.2969421 

# [1] 5.115238
# [1] 8.075039e-06
# [1,] 6.93892
# [1,] 1.475574e-10


z.rel.cls.Evergrneedle.Grassland$p.t
z.rel.cls.Evergrneedle.Grassland$p.F
     # landclassCropland   landclassDecid broad  landclassEvergr broad landclassEvergr needle     landclassGrassland  landclassMixed forest       landclassSavanna 
            # 0.10642832             0.24608554             0.53077754             0.52099238             0.06986644             0.22262794             0.43976299 
    # landclassShrubland 
            # 0.03451549 
# [1] 0.0009614116



# Fig. 7: plot results
pdf("Fig 7 remote sensing cls 19Nov19.pdf", height=6, width=6)
z.cls <- readRDS(file="Fig7 north_america cls 19Nov19.RDS")

par(mfrow=c(2,1), mai=c(.1,1,.2,.2))

# absolute c.cls
S <- z.cls$varcov
p <- length(z.cls$coef)
se.cond <- array(p)
for(i in 1:p){
	S11 <- S[i,i]
	S22 <- S[-i,-i]
	S12 <- S[i,-i]
	se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
}


arg <- 1:length(z.cls$coef)
ymax <- max(z.cls$coef + z.cls$se)
ymin <- min(z.cls$coef - z.cls$se)
plot(arg+.05, z.cls$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", pch=15, ylab=expression(paste("Slope ", italic(c[cls]))))
arrows(x0=arg, y0=z.cls$coef - z.cls$se, y1=z.cls$coef + z.cls$se, angle=90, code=3, length=.05)
arrows(x0=arg+.1, y0=z.cls$coef - se.cond, y1=z.cls$coef + se.cond, angle=90, code=3, length=.05, col="blue")
lines(c(0,20), c(0,0), lty=2)
mtext("A", at=10, cex=1.2)

n.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landclass), FUN=length)[,2]
lm.est <- lm(c.cls ~ 0 + landclass, data=dat.map)
c.data <- summary(lm.est)$coef[,1]
se.data <- summary(lm.est)$coef[,2]
points(arg + .3, c.data, col="red")
arrows(x0=arg + .3, y0=c.data - se.data, y1=c.data + se.data, angle=90, code=3, length=.05, col="red")


par(mai=c(1.55,1,.1,.5))

p.LS.est <- aggregate(dat.map$p.LS, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.LS.est, ylim=c(0,1), xlab="", xlim=c(0, 11), ylab=expression(paste(italic(P), "-values")), col="white")
p.est <- aggregate(dat.map$p, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.est, ylim=c(0,1), xlab="", xlim=c(0, 11), add=T, col="black")

at.min <- .8
at.max <- max(arg) + 2.5
at <- at.min + (at.max - at.min)*(0:(length(arg)-1))/length(arg)
axis(side=1, at=at, las=3, levels(dat.map$landclass))
lines(c(0,11), c(.05,.05), lty=1, col="red")
text(x=at, y=p.LS.est+.08, n.data)
mtext("B", at=at.max, cex=1.2)

dev.off()


# Fig. 7: plot results
pdf("Fig 7 remote sensing rel.cls 19Nov19.pdf", height=6, width=6)
z.cls <- readRDS(file="Fig7 north_america rel.cls 19Nov19.RDS")

par(mfrow=c(2,1), mai=c(.1,1,.2,.4))

# relative c.cls
S <- z.cls$varcov
p <- length(z.cls$coef)
se.cond <- array(p)
for(i in 1:p){
	S11 <- S[i,i]
	S22 <- S[-i,-i]
	S12 <- S[i,-i]
	se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
}


arg <- 1:length(z.cls$coef)
ymax <- max(z.cls$coef + z.cls$se+.005)
ymin <- min(z.cls$coef - z.cls$se-.005)
plot(arg+.05, z.cls$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", pch=15, ylab=expression(paste("Slope ", italic(c[cls]),"/mean NDVI")))
arrows(x0=arg, y0=z.cls$coef - z.cls$se, y1=z.cls$coef + z.cls$se, angle=90, code=3, length=.05)
arrows(x0=arg+.1, y0=z.cls$coef - se.cond, y1=z.cls$coef + se.cond, angle=90, code=3, length=.05, col="blue")
lines(c(0,20), c(0,0), lty=2)
mtext("A", at=9, cex=1.2)

n.data <- aggregate(dat.map$rel.c.cls, by=list(landscape=dat.map$landclass), FUN=length)[,2]
lm.est <- lm(rel.c.cls ~ 0 + landclass, data=dat.map)
c.data <- summary(lm.est)$coef[,1]
se.data <- summary(lm.est)$coef[,2]
points(arg + .3, c.data, col="red")
arrows(x0=arg + .3, y0=c.data - se.data, y1=c.data + se.data, angle=90, code=3, length=.05, col="red")

par(mai=c(1.55,1,.1,.3))

p.LS.est <- aggregate(dat.map$p.LS, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.LS.est, ylim=c(0,1), xlab="", xlim=c(0, 10), ylab=expression(paste(italic(P), "-values")), col="white")
p.est <- aggregate(dat.map$p, by=list(landscape=dat.map$landclass), FUN=function(x) mean(x<0.05))[,2]
barplot(p.est, ylim=c(0,1), xlab="", xlim=c(0, 11), add=T, col="black")

at.min <- .8
at.max <- max(arg) + 2.35
at <- at.min + (at.max - at.min)*(0:(length(arg)-1))/length(arg)
axis(side=1, at=at, las=3, levels(dat.map$landclass))
lines(c(0,11), c(.05,.05), lty=1, col="red")
text(x=at, y=p.LS.est+.15, n.data, srt=90, cex=.9)
mtext("B", at=at.max, cex=1.2)

dev.off()

######################################
# Diagnostics
# Fig. S2: plot results
pdf("Fig 7 (S2) remote sensing b 19Nov19.pdf", height=6, width=6)
z.cls <- readRDS(file="Fig7 north_america b 19Nov19.RDS")

par(mfrow=c(1,1), mai=c(.1,1,.2,.4))

# relative c.cls
S <- z.cls$varcov
p <- length(z.cls$coef)
se.cond <- array(p)
for(i in 1:p){
	S11 <- S[i,i]
	S22 <- S[-i,-i]
	S12 <- S[i,-i]
	se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
}

arg <- 1:length(z.cls$coef)
ymax <- max(z.cls$coef + z.cls$se+.005)
ymin <- min(z.cls$coef - z.cls$se-.005)
plot(arg+.05, z.cls$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", pch=15, ylab=expression(paste("Slope ", italic(c[cls]),"/mean NDVI")))
arrows(x0=arg, y0=z.cls$coef - z.cls$se, y1=z.cls$coef + z.cls$se, angle=90, code=3, length=.05)
arrows(x0=arg+.1, y0=z.cls$coef - se.cond, y1=z.cls$coef + se.cond, angle=90, code=3, length=.05, col="blue")
lines(c(0,20), c(0,0), lty=2)

n.data <- aggregate(dat.map$rel.c.cls, by=list(landscape=dat.map$landclass), FUN=length)[,2]
lm.est <- lm(rel.c.cls ~ 0 + landclass, data=dat.map)
c.data <- summary(lm.est)$coef[,1]
se.data <- summary(lm.est)$coef[,2]
points(arg + .3, c.data, col="red")
arrows(x0=arg + .3, y0=c.data - se.data, y1=c.data + se.data, angle=90, code=3, length=.05, col="red")

dev.off()



# check for homogeneity of variances among landclasses
sd.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landclass), FUN=sd)[,2]
         # landscape         x
# 1         Cropland 0.6606138
# 2 Cropland mosaics 0.8347937
# 3      Decid broad 0.5433243
# 4     Evergr broad 0.7300473
# 5    Evergr needle 0.8489938
# 6        Grassland 0.5320982
# 7     Mixed forest 0.8624745
# 8          Savanna 0.7012845
# 9        Shrubland 0.5923595

par(mfrow=c(1,1), mai=c(1.6,1,.1,.1))
plot(sd.data ~ arg, xlab="", xaxt="n", ylab="SD", ylim=c(min(sd.data), max(sd.data)+.1))
axis(side=1, at=arg, las=3, levels(dat.map$landclass))
text(x=arg, y=sd.data+.02, n.data)

# check for differences in b among landclasses
z.b <- GLS.fit(b ~ 0 + landclass, data=dat.map, invcholV=invcholV)
# $coef
        # landclassCropland landclassCropland mosaics      landclassDecid broad     landclassEvergr broad    landclassEvergr needle        landclassGrassland     landclassMixed forest 
                # 0.2108406                 0.1810689                 0.2094276                 0.2513201                 0.2096670                 0.2195387                 0.2465843 
         # landclassSavanna        landclassShrubland 
                # 0.2304908                 0.2251309 

# $se
        # landclassCropland landclassCropland mosaics      landclassDecid broad     landclassEvergr broad    landclassEvergr needle        landclassGrassland     landclassMixed forest 
                # 0.1829649                 0.1843195                 0.1830824                 0.1836050                 0.1829523                 0.1828266                 0.1830738 
         # landclassSavanna        landclassShrubland 
                # 0.1828174                 0.1829184 

# $t
        # landclassCropland landclassCropland mosaics      landclassDecid broad     landclassEvergr broad    landclassEvergr needle        landclassGrassland     landclassMixed forest 
                 # 1.152355                  0.982364                  1.143898                  1.368809                  1.146020                  1.200803                  1.346912 
         # landclassSavanna        landclassShrubland 
                 # 1.260770                  1.230772 

# $p.t
        # landclassCropland landclassCropland mosaics      landclassDecid broad     landclassEvergr broad    landclassEvergr needle        landclassGrassland     landclassMixed forest 
                # 0.2491860                 0.3259301                 0.2526768                 0.1710713                 0.2517980                 0.2298389                 0.1780209 
         # landclassSavanna        landclassShrubland 
                # 0.2074035                 0.2184197 

# $FF
# [1] 2.974835

# $p.FF
# [1] 0.002484895
z.b <- GLS.fit(b ~ 1, data=dat.map, invcholV=invcholV)

# $coef
# (Intercept) 
  # 0.2251349 

# $se
# (Intercept) 
  # 0.1828304 

# $t
# (Intercept) 
   # 1.231386 

# $p.t
# (Intercept) 
  # 0.2181901 


par(mfrow=c(1,1), mai=c(1.6,1,.2,.1))
arg <- 1:length(z.b$coef)
ymax <- max(z.b$coef + z.cls$se)
ymin <- min(z.b$coef - z.b$se)
plot(arg, z.b$coef, xlim=c(min(arg)-.5, max(arg)+.5), ylim=c(ymin,ymax), xlab="", xaxt="n", ylab=expression(paste("Autocorrelation ", italic(b))))
arrows(x0=arg, y0=z.b$coef - z.b$se, y1=z.b$coef + z.b$se, angle=90, code=3, length=.05)
lines(c(0,20), c(0,0), lty=2)
axis(side=1, at=arg, las=3, levels(dat.map$landclass))



# Variance in the estmates of r
r.ests <- NULL
for(rep in 1:10) r.ests <- c(r.ests, spatialcor.fit(X, t.scale, Dist, fit.n.sample, plot.fig=F)$spatialcor)
r.ests
c(mean(r.ests), sd(r.ests))

######################################
# OLD results for full map
######################################


# create distance matrix in kilometers
location <- dat[,c('lng','lat')]
Dist <- distm(location, fun=distGeo)/1000
V <- V.fit(Dist, spatialcor=r.est$spatialcor, FUN="exponential")
lsos()
rm(Dist)

# skip estimating the nugget
nugget <- 0.11

Vn <- (1 - nugget) * V
rm(V)
Vn <- Vn + nugget * diag(nrow(Vn))

system.time(invcholVn <- t(backsolve(chol(Vn), diag(nrow(Vn)))))
#    user   system  elapsed 
#16889.60    22.63 19437.09 

rm(Vn)
system.time(z.rel.cls <- GLS.fit(rel.c.cls ~ 0 + landclass, data=dat.map, invcholV=invcholVn))
rm(invcholVn)
cbind(z.rel.cls$coef, z.rel.cls$se, z.rel.cls$p.t)
#                             [,1]       [,2]      [,3]
#landclassEvergr needle 0.02973995 0.04983386 0.5506572
#landclassGrassland     0.03578251 0.04965064 0.4711083
#landclassMixed forest  0.03430228 0.04984387 0.4913353
#landclassSavanna       0.03329040 0.04967334 0.5027438
#landclassShrubland     0.03756878 0.04966850 0.4494208

z.rel.cls$p.F
#[1] 0.207918

saveRDS(file="Fig6 AK z.rel.cls.RDS", z.rel.cls)
 z.rel.cls[c(1:22)]
# $coef
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
            # 0.02973995             0.03578251             0.03430228 
      # landclassSavanna     landclassShrubland 
            # 0.03329040             0.03756878 

# $se
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
            # 0.04983386             0.04965064             0.04984387 
      # landclassSavanna     landclassShrubland 
            # 0.04967334             0.04966850 

# $t
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # 0.5967821              0.7206857              0.6881946 
      # landclassSavanna     landclassShrubland 
             # 0.6701865              0.7563905 

# $df.t
# [1] 31451

# $p.t
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # 0.5506572              0.4711083              0.4913353 
      # landclassSavanna     landclassShrubland 
             # 0.5027438              0.4494208 

# $F
# [1] 1.471292

# $df1.F
# [1] 4

# $df2.F
# [1] 31451

# $p.F
# [1] 0.207918

# $logLik
# [1] 41680.2

# $logLik0
# [1] 41675.26

# $MSE
# [1] 0.02785214

# $MSE0
# [1] 0.02785735

# $MSR
# [1] 0.04097862

# $SSE
# [1] 875.9777

# $SSE0
# [1] 876.1416

# $SSR
# [1] 0.1639145

# $coef0
           # [,1]
# [1,] 0.03523029

# $se0
# [1] 0.04965107

# $varX
                       # landclassEvergr needle landclassGrassland
# landclassEvergr needle              1580.6152          -195.2835
# landclassGrassland                  -195.2835          9334.3889
# landclassMixed forest               -120.3249          -222.0805
# landclassSavanna                   -1253.1481         -2271.0390
# landclassShrubland                   -11.6351         -6638.2821
                       # landclassMixed forest landclassSavanna
# landclassEvergr needle             -120.3249       -1253.1481
# landclassGrassland                 -222.0805       -2271.0390
# landclassMixed forest              1402.5767        -936.7923
# landclassSavanna                   -936.7923        9910.2017
# landclassShrubland                 -122.9239       -5447.0190
                       # landclassShrubland
# landclassEvergr needle           -11.6351
# landclassGrassland             -6638.2821
# landclassMixed forest           -122.9239
# landclassSavanna               -5447.0190
# landclassShrubland             12220.5744

# $varcov
                       # landclassEvergr needle landclassGrassland
# landclassEvergr needle            0.002483413        0.002463620
# landclassGrassland                0.002463620        0.002465186
# landclassMixed forest             0.002466551        0.002463286
# landclassSavanna                  0.002466507        0.002463840
# landclassShrubland                0.002464811        0.002464422
                       # landclassMixed forest landclassSavanna
# landclassEvergr needle           0.002466551      0.002466507
# landclassGrassland               0.002463286      0.002463840
# landclassMixed forest            0.002484412      0.002465805
# landclassSavanna                 0.002465805      0.002467441
# landclassShrubland               0.002464480      0.002465323
                       # landclassShrubland
# landclassEvergr needle        0.002464811
# landclassGrassland            0.002464422
# landclassMixed forest         0.002464480
# landclassSavanna              0.002465323
# landclassShrubland            0.002466960

# $varcov0
            # [,1]
# [1,] 0.002465229

system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, data=dat.map, invcholV=z.rel.cls$invcholV))
cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)

z.rel.cls.latland$p.F

z.rel.cls.latland[c(1:22)]

# lat                       -0.0021990084 0.007695004 0.775055477
# landclassEvergr needle     0.1646706954 0.474883182 0.728773629
# landclassGrassland         0.1323650866 0.459076411 0.773096737
# landclassMixed forest      0.6660412357 0.490525304 0.174533030
# landclassSavanna           0.1327046203 0.460840743 0.773377727
# landclassShrubland         0.0796329429 0.460280203 0.862644911
# lat:landclassGrassland     0.0006121685 0.001936263 0.751883296
# lat:landclassMixed forest -0.0079664751 0.003072207 0.009516494
# lat:landclassSavanna       0.0005454071 0.001894469 0.773429399
# lat:landclassShrubland     0.0014422624 0.001974136 0.465041859
# > 
# > z.rel.cls.latland$p.F
# [1] 0.03133547
# > 
# > z.rel.cls.latland[c(1:22)]
# $coef
                      # lat    landclassEvergr needle        landclassGrassland 
            # -0.0021990084              0.1646706954              0.1323650866 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
             # 0.6660412357              0.1327046203              0.0796329429 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
             # 0.0006121685             -0.0079664751              0.0005454071 
   # lat:landclassShrubland 
             # 0.0014422624 

# $se
                      # lat    landclassEvergr needle        landclassGrassland 
              # 0.007695004               0.474883182               0.459076411 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
              # 0.490525304               0.460840743               0.460280203 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # 0.001936263               0.003072207               0.001894469 
   # lat:landclassShrubland 
              # 0.001974136 

# $t
                      # lat    landclassEvergr needle        landclassGrassland 
               # -0.2857709                 0.3467604                 0.2883291 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
                # 1.3578122                 0.2879620                 0.1730097 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
                # 0.3161598                -2.5930791                 0.2878945 
   # lat:landclassShrubland 
                # 0.7305789 

# $df.t
# [1] 31446

# $p.t
                      # lat    landclassEvergr needle        landclassGrassland 
              # 0.775055477               0.728773629               0.773096737 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
              # 0.174533030               0.773377727               0.862644911 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # 0.751883296               0.009516494               0.773429399 
   # lat:landclassShrubland 
              # 0.465041859 

# $F
# [1] 2.039127

# $df1.F
# [1] 9

# $df2.F
# [1] 31446

# $p.F
# [1] 0.03133547

# $logLik
# [1] 41686.44

# $logLik0
# [1] 41672.76

# $MSE
# [1] 0.02784553

# $MSE0
# [1] 0.02786178

# $MSR
# [1] 0.05678058

# $SSE
# [1] 875.6306

# $SSE0
# [1] 876.1416

# $SSR
# [1] 0.5110253

# $coef0
           # [,1]
# [1,] 0.03523029

# $se0
# [1] 0.04965502

# $varX
                                  # lat landclassEvergr needle landclassGrassland
# lat                       42838.25621               13.30892           473.4153
# landclassEvergr needle       13.30892             1580.61519          -195.2835
# landclassGrassland          473.41531             -195.28347          9334.3889
# landclassMixed forest        26.10946             -120.32492          -222.0805
# landclassSavanna            134.38443            -1253.14806         -2271.0390
# landclassShrubland           44.45825              -11.63510         -6638.2821
# lat:landclassGrassland    29296.83993           -11768.72863        601669.5827
# lat:landclassMixed forest  1564.67707            -7447.32709        -13574.2216
# lat:landclassSavanna       8308.80805           -81266.06965       -139632.8786
# lat:landclassShrubland     2852.21899             -763.14054       -436223.9625
                          # landclassMixed forest landclassSavanna
# lat                                    26.10946         134.3844
# landclassEvergr needle               -120.32492       -1253.1481
# landclassGrassland                   -222.08048       -2271.0390
# landclassMixed forest                1402.57669        -936.7923
# landclassSavanna                     -936.79233        9910.2017
# landclassShrubland                   -122.92392       -5447.0190
# lat:landclassGrassland             -13575.01602     -139634.1119
# lat:landclassMixed forest           87684.40686      -59019.9212
# lat:landclassSavanna               -59021.27847      629056.5270
# lat:landclassShrubland              -7614.87607     -348999.5751
                          # landclassShrubland lat:landclassGrassland
# lat                                 44.45825               29296.84
# landclassEvergr needle             -11.63510              -11768.73
# landclassGrassland               -6638.28214              601669.58
# landclassMixed forest             -122.92392              -13575.02
# landclassSavanna                 -5447.01902             -139634.11
# landclassShrubland               12220.57439             -436222.43
# lat:landclassGrassland         -436222.42668            38909977.37
# lat:landclassMixed forest        -7615.67626             -829988.12
# lat:landclassSavanna           -349000.71091            -8603945.06
# lat:landclassShrubland          793647.41511           -28737553.53
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                        1564.677             8308.808
# landclassEvergr needle                    -7447.327           -81266.070
# landclassGrassland                       -13574.222          -139632.879
# landclassMixed forest                     87684.407           -59021.278
# landclassSavanna                         -59019.921           629056.527
# landclassShrubland                        -7615.676          -349000.711
# lat:landclassGrassland                  -829988.122         -8603945.060
# lat:landclassMixed forest               5485493.846         -3720351.492
# lat:landclassSavanna                   -3720351.492         40001681.086
# lat:landclassShrubland                  -471934.073        -22395709.373
                          # lat:landclassShrubland
# lat                                 2.852219e+03
# landclassEvergr needle             -7.631405e+02
# landclassGrassland                 -4.362240e+05
# landclassMixed forest              -7.614876e+03
# landclassSavanna                   -3.489996e+05
# landclassShrubland                  7.936474e+05
# lat:landclassGrassland             -2.873755e+07
# lat:landclassMixed forest          -4.719341e+05
# lat:landclassSavanna               -2.239571e+07
# lat:landclassShrubland              5.165821e+07

# $varcov
                                    # lat landclassEvergr needle landclassGrassland
# lat                        5.921309e-05          -0.0036338593      -3.398586e-03
# landclassEvergr needle    -3.633859e-03           0.2255140363       2.104818e-01
# landclassGrassland        -3.398586e-03           0.2104817553       2.107512e-01
# landclassMixed forest     -3.460620e-03           0.2143910927       2.102677e-01
# landclassSavanna          -3.416009e-03           0.2115549344       2.105658e-01
# landclassShrubland        -3.401925e-03           0.2106926361       2.106996e-01
# lat:landclassGrassland    -3.685140e-06           0.0002351284      -4.212895e-06
# lat:landclassMixed forest -2.695515e-06           0.0001728165       3.403102e-06
# lat:landclassSavanna      -3.399942e-06           0.0002175882      -1.304658e-06
# lat:landclassShrubland    -3.627959e-06           0.0002315361      -3.412479e-06
                          # landclassMixed forest landclassSavanna
# lat                               -3.460620e-03    -3.416009e-03
# landclassEvergr needle             2.143911e-01     2.115549e-01
# landclassGrassland                 2.102677e-01     2.105658e-01
# landclassMixed forest              2.406151e-01     2.113035e-01
# landclassSavanna                   2.113035e-01     2.123742e-01
# landclassShrubland                 2.105436e-01     2.110549e-01
# lat:landclassGrassland             6.483991e-05     1.602203e-05
# lat:landclassMixed forest         -4.205778e-04     4.226203e-06
# lat:landclassSavanna               4.785231e-05    -1.257547e-05
# lat:landclassShrubland             6.027155e-05     8.328038e-06
                          # landclassShrubland lat:landclassGrassland
# lat                            -3.401925e-03          -3.685140e-06
# landclassEvergr needle          2.106926e-01           2.351284e-04
# landclassGrassland              2.106996e-01          -4.212895e-06
# landclassMixed forest           2.105436e-01           6.483991e-05
# landclassSavanna                2.110549e-01           1.602203e-05
# landclassShrubland              2.118579e-01          -1.741599e-07
# lat:landclassGrassland         -1.741599e-07           3.749113e-06
# lat:landclassMixed forest       2.499185e-06           2.648619e-06
# lat:landclassSavanna           -5.561948e-06           3.419470e-06
# lat:landclassShrubland         -1.790259e-05           3.680879e-06
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                   -2.695515e-06        -3.399942e-06
# landclassEvergr needle                 1.728165e-04         2.175882e-04
# landclassGrassland                     3.403102e-06        -1.304658e-06
# landclassMixed forest                 -4.205778e-04         4.785231e-05
# landclassSavanna                       4.226203e-06        -1.257547e-05
# landclassShrubland                     2.499185e-06        -5.561948e-06
# lat:landclassGrassland                 2.648619e-06         3.419470e-06
# lat:landclassMixed forest              9.438454e-06         2.636315e-06
# lat:landclassSavanna                   2.636315e-06         3.589012e-06
# lat:landclassShrubland                 2.662195e-06         3.482239e-06
                          # lat:landclassShrubland
# lat                                -3.627959e-06
# landclassEvergr needle              2.315361e-04
# landclassGrassland                 -3.412479e-06
# landclassMixed forest               6.027155e-05
# landclassSavanna                    8.328038e-06
# landclassShrubland                 -1.790259e-05
# lat:landclassGrassland              3.680879e-06
# lat:landclassMixed forest           2.662195e-06
# lat:landclassSavanna                3.482239e-06
# lat:landclassShrubland              3.897215e-06

# $varcov0
            # [,1]
# [1,] 0.002465621

system.time(z.rel.cls.lat <- GLS.fit(rel.c.cls ~ 1 + lat, data=dat.map, invcholV=z.rel.cls$invcholV))
cbind(z.rel.cls.lat$coef, z.rel.cls.lat$se, z.rel.cls.lat$p.t)

z.rel.cls.lat$p.F
                    # [,1]        [,2]      [,3]
# (Intercept)  0.137341105 0.459124705 0.7648376
# lat         -0.001668212 0.007456861 0.8229804
# > 
# > z.rel.cls.lat$p.F
# [1] 0.8229804

system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, formula0 = ~0 + landclass, data=dat.map, invcholV=z.rel.cls$invcholV))
cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)

z.rel.cls.latland$p.F

z.rel.cls.latland[c(1:22)]

# # > cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)
                                   # [,1]        [,2]        [,3]
# lat                       -0.0021990084 0.007695004 0.775055477
# landclassEvergr needle     0.1646706954 0.474883182 0.728773629
# landclassGrassland         0.1323650866 0.459076411 0.773096737
# landclassMixed forest      0.6660412357 0.490525304 0.174533030
# landclassSavanna           0.1327046203 0.460840743 0.773377727
# landclassShrubland         0.0796329429 0.460280203 0.862644911
# lat:landclassGrassland     0.0006121685 0.001936263 0.751883296
# lat:landclassMixed forest -0.0079664751 0.003072207 0.009516494
# lat:landclassSavanna       0.0005454071 0.001894469 0.773429399
# lat:landclassShrubland     0.0014422624 0.001974136 0.465041859
# > 
# > z.rel.cls.latland$p.F
# [1] 0.02895773
# > 
# > z.rel.cls.latland[c(1:22)]
# $coef
                      # lat    landclassEvergr needle        landclassGrassland 
            # -0.0021990084              0.1646706954              0.1323650866 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
             # 0.6660412357              0.1327046203              0.0796329429 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
             # 0.0006121685             -0.0079664751              0.0005454071 
   # lat:landclassShrubland 
             # 0.0014422624 

# $se
                      # lat    landclassEvergr needle        landclassGrassland 
              # 0.007695004               0.474883182               0.459076411 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
              # 0.490525304               0.460840743               0.460280203 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # 0.001936263               0.003072207               0.001894469 
   # lat:landclassShrubland 
              # 0.001974136 

# $t
                      # lat    landclassEvergr needle        landclassGrassland 
               # -0.2857709                 0.3467604                 0.2883291 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
                # 1.3578122                 0.2879620                 0.1730097 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
                # 0.3161598                -2.5930791                 0.2878945 
   # lat:landclassShrubland 
                # 0.7305789 

# $df.t
# [1] 31446

# $p.t
                      # lat    landclassEvergr needle        landclassGrassland 
              # 0.775055477               0.728773629               0.773096737 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
              # 0.174533030               0.773377727               0.862644911 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # 0.751883296               0.009516494               0.773429399 
   # lat:landclassShrubland 
              # 0.465041859 

# $F
# [1] 2.493117

# $df1.F
# [1] 5

# $df2.F
# [1] 31446

# $p.F
# [1] 0.02895773

# $logLik
# [1] 41686.44

# $logLik0
# [1] 41675.7

# $MSE
# [1] 0.02784553

# $MSE0
# [1] 0.02785657

# $MSR
# [1] 0.06942215

# $SSE
# [1] 875.6306

# $SSE0
# [1] 875.9777

# $SSR
# [1] 0.3471108

# $coef0
                             # [,1]
# landclassEvergr needle 0.02973995
# landclassGrassland     0.03578251
# landclassMixed forest  0.03430228
# landclassSavanna       0.03329040
# landclassShrubland     0.03756878

# $se0
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
            # 0.04983782             0.04965459             0.04984783 
      # landclassSavanna     landclassShrubland 
            # 0.04967729             0.04967245 

# $varX
                                  # lat landclassEvergr needle landclassGrassland
# lat                       42838.25621               13.30892           473.4153
# landclassEvergr needle       13.30892             1580.61519          -195.2835
# landclassGrassland          473.41531             -195.28347          9334.3889
# landclassMixed forest        26.10946             -120.32492          -222.0805
# landclassSavanna            134.38443            -1253.14806         -2271.0390
# landclassShrubland           44.45825              -11.63510         -6638.2821
# lat:landclassGrassland    29296.83993           -11768.72863        601669.5827
# lat:landclassMixed forest  1564.67707            -7447.32709        -13574.2216
# lat:landclassSavanna       8308.80805           -81266.06965       -139632.8786
# lat:landclassShrubland     2852.21899             -763.14054       -436223.9625
                          # landclassMixed forest landclassSavanna
# lat                                    26.10946         134.3844
# landclassEvergr needle               -120.32492       -1253.1481
# landclassGrassland                   -222.08048       -2271.0390
# landclassMixed forest                1402.57669        -936.7923
# landclassSavanna                     -936.79233        9910.2017
# landclassShrubland                   -122.92392       -5447.0190
# lat:landclassGrassland             -13575.01602     -139634.1119
# lat:landclassMixed forest           87684.40686      -59019.9212
# lat:landclassSavanna               -59021.27847      629056.5270
# lat:landclassShrubland              -7614.87607     -348999.5751
                          # landclassShrubland lat:landclassGrassland
# lat                                 44.45825               29296.84
# landclassEvergr needle             -11.63510              -11768.73
# landclassGrassland               -6638.28214              601669.58
# landclassMixed forest             -122.92392              -13575.02
# landclassSavanna                 -5447.01902             -139634.11
# landclassShrubland               12220.57439             -436222.43
# lat:landclassGrassland         -436222.42668            38909977.37
# lat:landclassMixed forest        -7615.67626             -829988.12
# lat:landclassSavanna           -349000.71091            -8603945.06
# lat:landclassShrubland          793647.41511           -28737553.53
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                        1564.677             8308.808
# landclassEvergr needle                    -7447.327           -81266.070
# landclassGrassland                       -13574.222          -139632.879
# landclassMixed forest                     87684.407           -59021.278
# landclassSavanna                         -59019.921           629056.527
# landclassShrubland                        -7615.676          -349000.711
# lat:landclassGrassland                  -829988.122         -8603945.060
# lat:landclassMixed forest               5485493.846         -3720351.492
# lat:landclassSavanna                   -3720351.492         40001681.086
# lat:landclassShrubland                  -471934.073        -22395709.373
                          # lat:landclassShrubland
# lat                                 2.852219e+03
# landclassEvergr needle             -7.631405e+02
# landclassGrassland                 -4.362240e+05
# landclassMixed forest              -7.614876e+03
# landclassSavanna                   -3.489996e+05
# landclassShrubland                  7.936474e+05
# lat:landclassGrassland             -2.873755e+07
# lat:landclassMixed forest          -4.719341e+05
# lat:landclassSavanna               -2.239571e+07
# lat:landclassShrubland              5.165821e+07

# $varcov
                                    # lat landclassEvergr needle landclassGrassland
# lat                        5.921309e-05          -0.0036338593      -3.398586e-03
# landclassEvergr needle    -3.633859e-03           0.2255140363       2.104818e-01
# landclassGrassland        -3.398586e-03           0.2104817553       2.107512e-01
# landclassMixed forest     -3.460620e-03           0.2143910927       2.102677e-01
# landclassSavanna          -3.416009e-03           0.2115549344       2.105658e-01
# landclassShrubland        -3.401925e-03           0.2106926361       2.106996e-01
# lat:landclassGrassland    -3.685140e-06           0.0002351284      -4.212895e-06
# lat:landclassMixed forest -2.695515e-06           0.0001728165       3.403102e-06
# lat:landclassSavanna      -3.399942e-06           0.0002175882      -1.304658e-06
# lat:landclassShrubland    -3.627959e-06           0.0002315361      -3.412479e-06
                          # landclassMixed forest landclassSavanna
# lat                               -3.460620e-03    -3.416009e-03
# landclassEvergr needle             2.143911e-01     2.115549e-01
# landclassGrassland                 2.102677e-01     2.105658e-01
# landclassMixed forest              2.406151e-01     2.113035e-01
# landclassSavanna                   2.113035e-01     2.123742e-01
# landclassShrubland                 2.105436e-01     2.110549e-01
# lat:landclassGrassland             6.483991e-05     1.602203e-05
# lat:landclassMixed forest         -4.205778e-04     4.226203e-06
# lat:landclassSavanna               4.785231e-05    -1.257547e-05
# lat:landclassShrubland             6.027155e-05     8.328038e-06
                          # landclassShrubland lat:landclassGrassland
# lat                            -3.401925e-03          -3.685140e-06
# landclassEvergr needle          2.106926e-01           2.351284e-04
# landclassGrassland              2.106996e-01          -4.212895e-06
# landclassMixed forest           2.105436e-01           6.483991e-05
# landclassSavanna                2.110549e-01           1.602203e-05
# landclassShrubland              2.118579e-01          -1.741599e-07
# lat:landclassGrassland         -1.741599e-07           3.749113e-06
# lat:landclassMixed forest       2.499185e-06           2.648619e-06
# lat:landclassSavanna           -5.561948e-06           3.419470e-06
# lat:landclassShrubland         -1.790259e-05           3.680879e-06
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                   -2.695515e-06        -3.399942e-06
# landclassEvergr needle                 1.728165e-04         2.175882e-04
# landclassGrassland                     3.403102e-06        -1.304658e-06
# landclassMixed forest                 -4.205778e-04         4.785231e-05
# landclassSavanna                       4.226203e-06        -1.257547e-05
# landclassShrubland                     2.499185e-06        -5.561948e-06
# lat:landclassGrassland                 2.648619e-06         3.419470e-06
# lat:landclassMixed forest              9.438454e-06         2.636315e-06
# lat:landclassSavanna                   2.636315e-06         3.589012e-06
# lat:landclassShrubland                 2.662195e-06         3.482239e-06
                          # lat:landclassShrubland
# lat                                -3.627959e-06
# landclassEvergr needle              2.315361e-04
# landclassGrassland                 -3.412479e-06
# landclassMixed forest               6.027155e-05
# landclassSavanna                    8.328038e-06
# landclassShrubland                 -1.790259e-05
# lat:landclassGrassland              3.680879e-06
# lat:landclassMixed forest           2.662195e-06
# lat:landclassSavanna                3.482239e-06
# lat:landclassShrubland              3.897215e-06

# $varcov0
                       # landclassEvergr needle landclassGrassland
# landclassEvergr needle            0.002483808        0.002464012
# landclassGrassland                0.002464012        0.002465578
# landclassMixed forest             0.002466943        0.002463678
# landclassSavanna                  0.002466899        0.002464232
# landclassShrubland                0.002465203        0.002464813
                       # landclassMixed forest landclassSavanna landclassShrubland
# landclassEvergr needle           0.002466943      0.002466899        0.002465203
# landclassGrassland               0.002463678      0.002464232        0.002464813
# landclassMixed forest            0.002484807      0.002466197        0.002464872
# landclassSavanna                 0.002466197      0.002467833        0.002465715
# landclassShrubland     


system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, formula0 = ~0 + lat, data=dat.map, invcholV=z.rel.cls$invcholV))
cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)

z.rel.cls.latland$p.F

# > cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)
                                   # [,1]        [,2]        [,3]
# lat                       -0.0021990084 0.007695004 0.775055477
# landclassEvergr needle     0.1646706954 0.474883182 0.728773629
# landclassGrassland         0.1323650866 0.459076411 0.773096737
# landclassMixed forest      0.6660412357 0.490525304 0.174533030
# landclassSavanna           0.1327046203 0.460840743 0.773377727
# landclassShrubland         0.0796329429 0.460280203 0.862644911
# lat:landclassGrassland     0.0006121685 0.001936263 0.751883296
# lat:landclassMixed forest -0.0079664751 0.003072207 0.009516494
# lat:landclassSavanna       0.0005454071 0.001894469 0.773429399
# lat:landclassShrubland     0.0014422624 0.001974136 0.465041859
# > 
# > z.rel.cls.latland$p.F
# [1] 0.03092704

system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, formula0 = ~0 + lat + landclass, data=dat.map, invcholV=z.rel.cls$invcholV))
cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)

z.rel.cls.latland$p.F
                                   [,1]        [,2]        [,3]
# lat                       -0.0021990084 0.007695004 0.775055477
# landclassEvergr needle     0.1646706954 0.474883182 0.728773629
# landclassGrassland         0.1323650866 0.459076411 0.773096737
# landclassMixed forest      0.6660412357 0.490525304 0.174533030
# landclassSavanna           0.1327046203 0.460840743 0.773377727
# landclassShrubland         0.0796329429 0.460280203 0.862644911
# lat:landclassGrassland     0.0006121685 0.001936263 0.751883296
# lat:landclassMixed forest -0.0079664751 0.003072207 0.009516494
# lat:landclassSavanna       0.0005454071 0.001894469 0.773429399
# lat:landclassShrubland     0.0014422624 0.001974136 0.465041859
# > 
# > z.rel.cls.latland$p.F
# [1] 0.01453302


# system.time(z.rel.cls.1 <- GLS.fit(rel.c.cls ~ 1, data=dat.map, invcholV=z.rel.cls$invcholV))
# cbind(z.rel.cls.1$coef, z.rel.cls.1$se, z.rel.cls.1$p.t)

# z.rel.cls.1$p.F
                  # [,1]       [,2]      [,3]
# (Intercept) 0.03523029 0.04964792 0.4779559
# > 
# > z.rel.cls.1$p.F
# [1] 1

#########
# PROBLEM
> system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, formula0 = ~1, data=dat.map, invcholV=z.rel.cls$invcholV))
   user  system elapsed 
  20.66    0.01   20.69 
> z.rel.cls.latland$p.F
[1] 0.03133547
> 

###########################################################
# CLS ANOVA
###########################################################
system.time(z.rel.cls.1 <- GLS.fit(rel.c.cls ~ 1, data=dat.map))
cbind(z.rel.cls.1$coef, z.rel.cls.1$se, z.rel.cls.1$p.t)

z.rel.cls.1$p.F

                   # [,1]         [,2]         [,3]
# (Intercept) 0.009111838 0.0005270967 1.202155e-66
# > 
# > z.rel.cls.1$p.F
# [1] 1
# > 
# > 


system.time(z.rel.cls.land <- GLS.fit(rel.c.cls ~ 0 + landclass, data=dat.map))
cbind(z.rel.cls.land$coef, z.rel.cls.land$se, z.rel.cls.land$p.t)

z.rel.cls.land$p.F

z.rel.cls.land[c(1:22)]

# #    user  system elapsed 
  # 18.39    0.79   19.20 
# > cbind(z.rel.cls.land$coef, z.rel.cls.land$se, z.rel.cls.land$p.t)
                              # [,1]         [,2]          [,3]
# landclassEvergr needle -0.04565036 0.0050483972  1.616153e-19
# landclassGrassland      0.01932420 0.0011740682  1.292808e-60
# landclassMixed forest  -0.04137005 0.0056185653  1.840887e-13
# landclassSavanna       -0.01282482 0.0008229749  1.508534e-54
# landclassShrubland      0.02824410 0.0008160736 1.244230e-257
# > 
# > z.rel.cls.land$p.F
# [1] 3.952525e-323
# > 
# > z.rel.cls.land[c(1:22)]
# $coef
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
           # -0.04565036             0.01932420            -0.04137005 
      # landclassSavanna     landclassShrubland 
           # -0.01282482             0.02824410 

# $se
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
          # 0.0050483972           0.0011740682           0.0056185653 
      # landclassSavanna     landclassShrubland 
          # 0.0008229749           0.0008160736 

# $t
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
             # -9.042544              16.459178              -7.363099 
      # landclassSavanna     landclassShrubland 
            # -15.583491              34.609748 

# $df.t
# [1] 31451

# $p.t
# landclassEvergr needle     landclassGrassland  landclassMixed forest 
          # 1.616153e-19           1.292808e-60           1.840887e-13 
      # landclassSavanna     landclassShrubland 
          # 1.508534e-54          1.244230e-257 

# $F
# [1] 383.5487

# $df1.F
# [1] 4

# $df2.F
# [1] 31451

# $p.F
# [1] 3.952525e-323

# $logLik
# [1] 30664.73

# $logLik0
# [1] 29913.64

# $MSE
# [1] 0.008334025

# $MSE0
# [1] 0.008740562

# $MSR
# [1] 3.196504

# $SSE
# [1] 262.1134

# $SSE0
# [1] 274.8994

# $SSR
# [1] 12.78602

# $coef0
            # [,1]
# [1,] 0.009111838

# $se0
# [1] 0.0005271303

# $varX
                       # landclassEvergr needle landclassGrassland
# landclassEvergr needle                    327                  0
# landclassGrassland                          0               6046
# landclassMixed forest                       0                  0
# landclassSavanna                            0                  0
# landclassShrubland                          0                  0
                       # landclassMixed forest landclassSavanna landclassShrubland
# landclassEvergr needle                     0                0                  0
# landclassGrassland                         0                0                  0
# landclassMixed forest                    264                0                  0
# landclassSavanna                           0            12305                  0
# landclassShrubland                         0                0              12514

# $varcov
                       # landclassEvergr needle landclassGrassland
# landclassEvergr needle           2.548631e-05       0.000000e+00
# landclassGrassland               0.000000e+00       1.378436e-06
# landclassMixed forest            0.000000e+00       0.000000e+00
# landclassSavanna                 0.000000e+00       0.000000e+00
# landclassShrubland               0.000000e+00       0.000000e+00
                       # landclassMixed forest landclassSavanna landclassShrubland
# landclassEvergr needle          0.000000e+00     0.000000e+00       0.000000e+00
# landclassGrassland              0.000000e+00     0.000000e+00       0.000000e+00
# landclassMixed forest           3.156828e-05     0.000000e+00       0.000000e+00
# landclassSavanna                0.000000e+00     6.772877e-07       0.000000e+00
# landclassShrubland              0.000000e+00     0.000000e+00       6.659761e-07

# $varcov0
             # [,1]
# [1,] 2.778663e-07


system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, data=dat.map))
cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)

z.rel.cls.latland$p.F

z.rel.cls.latland[c(1:22)]

# > system.time(z.rel.cls.latland <- GLS.fit(rel.c.cls ~ 0 + lat*landclass, data=dat.map))
   # user  system elapsed 
  # 24.33    0.98   25.33 
# > cbind(z.rel.cls.latland$coef, z.rel.cls.latland$se, z.rel.cls.latland$p.t)
                                   # [,1]        [,2]          [,3]
# lat                        1.089694e-02 0.001862487  4.942045e-09
# landclassEvergr needle    -7.402164e-01 0.118799182  4.698727e-10
# landclassGrassland        -6.799880e-01 0.014698494  0.000000e+00
# landclassMixed forest      8.148888e-02 0.191623364  6.706534e-01
# landclassSavanna          -6.081987e-01 0.022064156 2.710720e-165
# landclassShrubland        -1.033308e+00 0.014812744  0.000000e+00
# lat:landclassGrassland    -3.101612e-05 0.001876367  9.868118e-01
# lat:landclassMixed forest -1.286559e-02 0.003590317  3.396466e-04
# lat:landclassSavanna      -1.609342e-03 0.001893989  3.954926e-01
# lat:landclassShrubland     5.090661e-03 0.001875768  6.653099e-03
# > 
# > z.rel.cls.latland$p.F
# [1] 0
# > 
# > z.rel.cls.latland[c(1:22)]
# $coef
                      # lat    landclassEvergr needle        landclassGrassland 
             # 1.089694e-02             -7.402164e-01             -6.799880e-01 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
             # 8.148888e-02             -6.081987e-01             -1.033308e+00 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
            # -3.101612e-05             -1.286559e-02             -1.609342e-03 
   # lat:landclassShrubland 
             # 5.090661e-03 

# $se
                      # lat    landclassEvergr needle        landclassGrassland 
              # 0.001862487               0.118799182               0.014698494 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
              # 0.191623364               0.022064156               0.014812744 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # 0.001876367               0.003590317               0.001893989 
   # lat:landclassShrubland 
              # 0.001875768 

# $t
                      # lat    landclassEvergr needle        landclassGrassland 
               # 5.85074970               -6.23082092              -46.26242730 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
               # 0.42525547              -27.56501135              -69.75807147 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
              # -0.01652988               -3.58341272               -0.84971038 
   # lat:landclassShrubland 
               # 2.71390734 

# $df.t
# [1] 31446

# $p.t
                      # lat    landclassEvergr needle        landclassGrassland 
             # 4.942045e-09              4.698727e-10              0.000000e+00 
    # landclassMixed forest          landclassSavanna        landclassShrubland 
             # 6.706534e-01             2.710720e-165              0.000000e+00 
   # lat:landclassGrassland lat:landclassMixed forest      lat:landclassSavanna 
             # 9.868118e-01              3.396466e-04              3.954926e-01 
   # lat:landclassShrubland 
             # 6.653099e-03 

# $F
# [1] 1124.47

# $df1.F
# [1] 9

# $df2.F
# [1] 31446

# $p.F
# [1] 0

# $logLik
# [1] 34304.01

# $logLik0
# [1] 29911.14

# $MSE
# [1] 0.006613528

# $MSE0
# [1] 0.008741952

# $MSR
# [1] 7.436715

# $SSE
# [1] 207.969

# $SSE0
# [1] 274.8994

# $SSR
# [1] 66.93044

# $coef0
            # [,1]
# [1,] 0.009111838

# $se0
# [1] 0.0005271722

# $varX
                                   # lat landclassEvergr needle landclassGrassland
# lat                       133454966.45               20842.83             389110
# landclassEvergr needle        20842.83                 327.00                  0
# landclassGrassland           389110.04                   0.00               6046
# landclassMixed forest         16475.70                   0.00                  0
# landclassSavanna             788801.65                   0.00                  0
# landclassShrubland           830910.51                   0.00                  0
# lat:landclassGrassland     25169883.06                   0.00             389110
# lat:landclassMixed forest   1028916.17                   0.00                  0
# lat:landclassSavanna       50621351.44                   0.00                  0
# lat:landclassShrubland     55304397.33                   0.00                  0
                          # landclassMixed forest landclassSavanna
# lat                                     16475.7         788801.6
# landclassEvergr needle                      0.0              0.0
# landclassGrassland                          0.0              0.0
# landclassMixed forest                     264.0              0.0
# landclassSavanna                            0.0          12305.0
# landclassShrubland                          0.0              0.0
# lat:landclassGrassland                      0.0              0.0
# lat:landclassMixed forest               16475.7              0.0
# lat:landclassSavanna                        0.0         788801.6
# lat:landclassShrubland                      0.0              0.0
                          # landclassShrubland lat:landclassGrassland
# lat                                 830910.5               25169883
# landclassEvergr needle                   0.0                      0
# landclassGrassland                       0.0                 389110
# landclassMixed forest                    0.0                      0
# landclassSavanna                         0.0                      0
# landclassShrubland                   12514.0                      0
# lat:landclassGrassland                   0.0               25169883
# lat:landclassMixed forest                0.0                      0
# lat:landclassSavanna                     0.0                      0
# lat:landclassShrubland              830910.5                      0
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                       1028916.2           50621351.4
# landclassEvergr needle                          0.0                  0.0
# landclassGrassland                              0.0                  0.0
# landclassMixed forest                       16475.7                  0.0
# landclassSavanna                                0.0             788801.6
# landclassShrubland                              0.0                  0.0
# lat:landclassGrassland                          0.0                  0.0
# lat:landclassMixed forest                 1028916.2                  0.0
# lat:landclassSavanna                            0.0           50621351.4
# lat:landclassShrubland                          0.0                  0.0
                          # lat:landclassShrubland
# lat                                   55304397.3
# landclassEvergr needle                       0.0
# landclassGrassland                           0.0
# landclassMixed forest                        0.0
# landclassSavanna                             0.0
# landclassShrubland                      830910.5
# lat:landclassGrassland                       0.0
# lat:landclassMixed forest                    0.0
# lat:landclassSavanna                         0.0
# lat:landclassShrubland                55304397.3

# $varcov
                                    # lat landclassEvergr needle landclassGrassland
# lat                        3.468857e-06          -2.211033e-04      -3.185368e-18
# landclassEvergr needle    -2.211033e-04           1.411325e-02       2.026492e-16
# landclassGrassland        -8.975390e-21           2.526587e-18       2.160457e-04
# landclassMixed forest      9.022695e-18          -5.716384e-16       8.642006e-18
# landclassSavanna           7.394075e-18          -4.732493e-16       7.135798e-18
# landclassShrubland        -1.046705e-19           6.671650e-18       9.636569e-32
# lat:landclassGrassland    -3.468857e-06           2.211033e-04      -3.339926e-06
# lat:landclassMixed forest -3.468857e-06           2.211033e-04       3.047003e-18
# lat:landclassSavanna      -3.468857e-06           2.211033e-04       3.074163e-18
# lat:landclassShrubland    -3.468857e-06           2.211033e-04       3.185479e-18
                          # landclassMixed forest landclassSavanna
# lat                                1.795154e-17     0.000000e+00
# landclassEvergr needle            -1.140576e-15    -2.500290e-32
# landclassGrassland                -8.734607e-18     8.031876e-34
# landclassMixed forest              3.671951e-02     1.266124e-32
# landclassSavanna                  -1.073476e-19     4.868270e-04
# landclassShrubland                -4.295344e-19     4.005923e-20
# lat:landclassGrassland            -1.781619e-17     5.429206e-34
# lat:landclassMixed forest         -5.879775e-04     3.526304e-34
# lat:landclassSavanna              -1.795023e-17    -7.585928e-06
# lat:landclassShrubland            -1.794477e-17    -6.018624e-22
                          # landclassShrubland lat:landclassGrassland
# lat                            -3.041397e-18          -3.468857e-06
# landclassEvergr needle          1.944842e-16           2.211033e-04
# landclassGrassland              4.390149e-18          -3.339926e-06
# landclassMixed forest           6.936545e-17          -9.113874e-18
# landclassSavanna               -3.583962e-17          -7.546461e-18
# landclassShrubland              2.194174e-04           1.046705e-19
# lat:landclassGrassland          2.973496e-18           3.520753e-06
# lat:landclassMixed forest       1.930309e-18           3.468857e-06
# lat:landclassSavanna            3.600014e-18           3.468857e-06
# lat:landclassShrubland         -3.296595e-06           3.468857e-06
                          # lat:landclassMixed forest lat:landclassSavanna
# lat                                   -3.468857e-06        -3.468857e-06
# landclassEvergr needle                 2.211033e-04         2.211033e-04
# landclassGrassland                     1.395145e-19         8.613938e-21
# landclassMixed forest                 -5.879775e-04        -9.023044e-18
# landclassSavanna                      -7.386561e-18        -7.585928e-06
# landclassShrubland                     1.117425e-19         1.042986e-19
# lat:landclassGrassland                 3.468857e-06         3.468857e-06
# lat:landclassMixed forest              1.289037e-05         3.468857e-06
# lat:landclassSavanna                   3.468857e-06         3.587194e-06
# lat:landclassShrubland                 3.468857e-06         3.468857e-06
                          # lat:landclassShrubland
# lat                                -3.468857e-06
# landclassEvergr needle              2.211033e-04
# landclassGrassland                 -3.429804e-20
# landclassMixed forest              -9.964509e-18
# landclassSavanna                   -6.845272e-18
# landclassShrubland                 -3.296595e-06
# lat:landclassGrassland              3.468857e-06
# lat:landclassMixed forest           3.468857e-06
# lat:landclassSavanna                3.468857e-06
# lat:landclassShrubland              3.518506e-06

# $varcov0
             # [,1]
# [1,] 2.779105e-07




