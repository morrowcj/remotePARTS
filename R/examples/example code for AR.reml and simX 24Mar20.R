#This is the base code used to produce figures 3-5 in the RSE manuscript. It is not the full code, though. Here I just present it to illustrate use of AR.reml() and simX()

library(mvtnorm)
library(geosphere)
library(mblm)
library(fields)

source('remote_sensing_tools_24Mar20.R')

###############################
# Set up spatial covariance matrix
###############################

# set up variables
nSpace <- 12
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
Dist <- Dist/(max(Dist)/2^.5)

# distribution of 4 land classes
n.size <- nSpace/4
n.cluster <- 2
landscape <- kronecker(kronecker(matrix(1,ncol=n.cluster,nrow=n.cluster), matrix(1:4, nrow=2, ncol=2)), matrix(1, nrow=n.size, ncol=n.size))
dat <- data.frame(landclass=matrix(landscape, ncol=1) - 1)
dat$landclass <- as.factor(dat$landclass)

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- 0
seed <- 0

slopes <- 0:3

t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

nugget.fit.flag <- T

r <- .1

if(r > 0) {
	Dr <- t(chol(exp(-Dist/r)))
} else {
	Dr <- diag(nrow(Dist))
}

alpha <- .1

slopes <- 0:3
cc <- 0
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

##################################
# simulate data with simX
##################################
X <- simX('~0 + c0', data=dat, coef=c(1), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)

z <- data.frame(site=1:n)
for(counter in 1:dim(X)[1]){
	x <- X[counter,]
	z$r[counter] <- r
	z$b[counter] <- b
	z$c[counter] <- cc

	z.LS <- lm(x ~ t.scale)
	z$cc.LS[counter] <- summary(z.LS)$coef[2,1]
	z$t.LS[counter] <- summary(z.LS)$coef[2,3]
	z$p.LS[counter] <- summary(z.LS)$coef[2,4]

	if(r > 0){
		z.Kendall <- cor.test(x, t.scale, method="kendall")
		z$tau.Kendall[counter] <- z.Kendall$estimate
		z$p.Kendall[counter] <- z.Kendall$p.value
		z$cc.Kendall[counter] <- mblm(x ~ t.scale, repeated = F)$coef[2]

		z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x)-1)] + t.scale[2:length(x)])
		z$cc.CLS[counter] <- summary(z.CLS)$coef[3,1]/(1 - summary(z.CLS)$coef[2,1])
		z$t.CLS[counter] <- summary(z.CLS)$coef[3,3]
		z$p.CLS[counter] <- summary(z.CLS)$coef[3,4]
		z$b.CLS[counter] <- summary(z.CLS)$coef[2,1]

##################################
# analyze data with AR.reml
##################################
		z.AR.reml <- AR.reml(x ~ t.scale)
		z$cc.AR.reml[counter] <- z.AR.reml$beta[2]
		z$b.AR.reml[counter] <- z.AR.reml$b
		z$p.AR.reml[counter] <- z.AR.reml$Pr[2]
	}
}

mean(z$p.LS < alpha)
mean(z$p.Kendall < alpha)
mean(z$p.CLS < alpha)
mean(z$p.AR.reml < alpha)

#################
# Fig 3 LS only
par(mfcol = c(1,1), mai=c(.1,.1,.1,.1))
col.pal <- hcl.colors(21, "RdBu", rev = TRUE)

zmax <- max(abs(z$cc.CLS))

matrix.fig <- matrix(z$cc.LS, ncol=ydim)
mask <- matrix(z$p.LS, ncol=ydim)
matrix.fig[mask > alpha] <- NA
image.plot(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
text("LS", x=.5, y=.95, cex=2.5)

#################
# Fig 4 All methods
par(mfcol = c(2,2), mai=c(.1,.1,.1,.5))
col.pal <- hcl.colors(21, "RdBu", rev = TRUE)

zmax <- max(abs(z$cc.CLS))

matrix.fig <- matrix(z$cc.LS, ncol=ydim)
mask <- matrix(z$p.LS, ncol=ydim)
matrix.fig[mask > alpha] <- NA
image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
text("LS", x=.5, y=.95, cex=1.5)

matrix.fig <- matrix(z$cc.CLS, ncol=ydim)
mask <- matrix(z$p.CLS, ncol=ydim)
matrix.fig[mask > alpha] <- NA
image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
text("CLS", x=.5, y=.95, cex=1.5)

matrix.fig <- matrix(z$cc.Kendall, ncol=ydim)
mask <- matrix(z$p.Kendall, ncol=ydim)
matrix.fig[mask > alpha] <- NA
image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
text("Mann-Kendall", x=.5, y=.95, cex=1.5)

par(mai=c(.1,.1,.1,.1))
matrix.fig <- matrix(z$cc.AR.reml, ncol=ydim)
mask <- matrix(z$p.AR.reml, ncol=ydim)
matrix.fig[mask > alpha] <- NA
image.plot(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, zlim=c(-zmax,zmax))
text("AR-REML", x=.5, y=.95, cex=1.5)








###############################
#Fig5: land classes
###############################

Fig5.plot <- function(M, lab){
	col.pal <- hcl.colors(20, "RdBu", rev = TRUE)

	matrix.fig <- M$matrix.fig
	mask <- M$mask
	c.data <- M$c.data
	se.data <- M$se.data
	c.est <- M$c.est
	se.est <- M$se.est

	par(mai=c(.22,.22,.22,.22))

	matrix.fig[mask > alpha] <- NA
	image(landscape, xaxt="n", yaxt="n")
	mtext(lab[1], side=4, at=.98, adj=-.5, las=1, cex=1.2)
	image(matrix.fig, xlab="", ylab="", main="", xaxt="n", yaxt="n", col=col.pal, add=T, zlim=c(-zmax,zmax))

	par(mai=c(.8,.8,.3,.3))

	arg <- 1:4
	plot(arg+.15, c.data, xlim=c(.5,4.5), ylim=c(-.4,.5), ylab="Slope", xlab="Land class", main="", pch=1, cex.lab=1.5, col="red")
	mtext(lab[2], side=4, at=.5, adj=-.5, las=1, cex=1.2)
	arrows(x0=arg+.15, y0=c.data-se.data, y1=c.data+se.data, angle=90, code=3, length=.05, col="red")

	# coef
	points(arg+.05, M$coef, col="black", pch=15)
	arrows(x0=arg, y0=M$coef-M$se, y1=M$coef+M$se, angle=90, length=.05, code=3)

	# for conditional confidence intervals
	S <- M$varcov
	p <- length(M$coef)
	se.cond <- array(p)
	for(i in 1:p){
		S11 <- S[i,i]
		S22 <- S[-i,-i]
		S12 <- S[i,-i]
		se.cond[i] <- (S11 - t(S12) %*% solve(S22) %*% S12)^.5
	}
	arrows(x0=arg+.1, y0=M$coef - se.cond, y1=M$coef + se.cond, angle=90, code=3, length=.05, col="blue")

}

#################################################
# set up variables
nSpace <- 100
xdim <- nSpace
ydim <- nSpace
n <- nSpace^2

location <- cbind(rep(1:nSpace,times=nSpace),rep(1:nSpace,each=nSpace)) * 10^-3
colnames(location) <- c("lng", "lat")
Dist <- distm(location)
Dist <- Dist/(max(Dist)/2^.5)

# distribution of 4 land classes
n.size <- nSpace/4
n.cluster <- 2
landscape <- kronecker(kronecker(matrix(1,ncol=n.cluster,nrow=n.cluster), matrix(1:4, nrow=2, ncol=2)), matrix(1, nrow=n.size, ncol=n.size))
dat <- data.frame(landclass=matrix(landscape, ncol=1) - 1)
dat$landclass <- as.factor(dat$landclass)

n.obs <- 30
n.burn <- 10

b0 <- 0
b <- .2
s <- 1
c0.sd <- .25

slopes <- 0:3

t.scale <- 1:n.obs
t.scale <- (t.scale-min(t.scale))/max(t.scale)

nugget.fit.flag <- T

r <- .1

if(r > 0) {
	Dr <- t(chol(exp(-Dist/r)))
} else {
	Dr <- diag(nrow(Dist))
}

##################################
# simulation

# for cc = 0
seed <- 10

slopes <- 0:3
cc <- 0
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

X <- simX('~0 + c0 + landclass', data=dat, coef=c(1, cc * slopes), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)

dat.map <- CLS.fit(X,t.scale)
dat.map$landscape <- as.factor(dat$landclass)
dat.map$lng <- location[,1]
dat.map$lat <- location[,2]
dat.map$c.cls <- dat.map$c
dat.map$c.cor <- dat.map$c/(1-dat.map$b)

# look for correlation between c.cls and b
par(mfrow=c(3,2))
hist(dat.map$b)
hist(dat.map$c.cls)
hist(dat.map$c.cor)
plot(c.cls ~ b, data=dat.map)
plot(c.cls ~ c.cor, data=dat.map)

sd(dat.map$c.cls)
sd(dat.map$c.cor)

fit.n.sample <- 2000
r.est <- spatialcor.fit(X, t.scale, Dist=Dist, fit.n.sample=fit.n.sample, plot.fig=T)
r.est

V <- exp(-Dist/r.est$spatialcor)
nugget <- nugget.fit(formula='c.cls ~ 0 + landscape', dat.map, V, nugget.tol = 0.00001, verbose = T)
nugget

Vn <- (1 - nugget) * V + nugget * diag(n)

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.cls <- GLS.fit(c.cls ~ 0 + landscape, data=dat.map, invcholV=invcholVn)

# save output
c.cls.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landscape), FUN=mean)[,2]
se.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landscape), FUN=mean)[,2]/(n/4)^.5
matrix.fig <- matrix(dat.map$c.cls, ncol=ydim)
mask <- matrix(dat.map$p, ncol=ydim)

M <- list(c.data=c.cls.data, se.data=se.data, matrix.fig=matrix.fig, mask=mask, coef=z.cls$coef, se=z.cls$se, FF=z.cls$F, p.FF=z.cls$p.F, varcov=z.cls$varcov, r=r.est$spatialcor, nugget=nugget, coef0=z.cls$coef0, se0=z.cls$se0)

saveRDS(file=paste0("Fig5 rel.cls cc=",cc," seed=",seed," 24Mar20.RDS"), M)

# for cc = 0.1
seed <- 10

slopes <- 0:3
cc <- .1
dat$c0 <- rnorm(n=n, 0, sd=c0.sd)

X <- simX('~0 + c0 + landclass', data=dat, coef=c(1, cc * slopes), b=b, s=s, Dr=Dr, t.scale=t.scale, n=n, n.obs=n.obs, n.burn=n.burn, seed=seed)

dat.map <- CLS.fit(X,t.scale)
dat.map$landscape <- as.factor(dat$landclass)
dat.map$lng <- location[,1]
dat.map$lat <- location[,2]

dat.map$c.cls <- dat.map$c

fit.n.sample <- 2000
r.est <- spatialcor.fit(X, t.scale, Dist=Dist, fit.n.sample=fit.n.sample, plot.fig=T)
r.est


V <- exp(-Dist/r.est$spatialcor)
nugget <- nugget.fit(formula='c.cls ~ 0 + landscape', dat.map, V, nugget.tol = 0.00001, verbose = T)
Vn <- (1 - nugget) * V + nugget * diag(n)
nugget

invcholVn <- t(backsolve(chol(Vn), diag(n)))
z.cls <- GLS.fit(c.cls ~ 0 + landscape, data=dat.map, invcholV=invcholVn)


# save output
c.cls.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landscape), FUN=mean)[,2]
se.data <- aggregate(dat.map$c.cls, by=list(landscape=dat.map$landscape), FUN=mean)[,2]/(n/4)^.5
matrix.fig <- matrix(dat.map$c.cls, ncol=ydim)
mask <- matrix(dat.map$p, ncol=ydim)

M <- list(c.data=c.cls.data, se.data=se.data, matrix.fig=matrix.fig, mask=mask, coef=z.cls$coef, se=z.cls$se, FF=z.cls$F, p.FF=z.cls$p.F, varcov=z.cls$varcov, r=r.est$spatialcor, nugget=nugget, coef0=z.cls$coef0, se0=z.cls$se0)

saveRDS(file=paste0("Fig5 rel.cls cc=",cc," seed=",seed," 24Mar20.RDS"), M)


##################################
# plotting
cc <- 0
M1 <- readRDS(file=paste0("Fig5 rel.cls cc=",cc," seed=",seed," 24Mar20.RDS"))

cc <- .1
M2 <- readRDS(file=paste0("Fig5 rel.cls cc=",cc," seed=",seed," 24Mar20.RDS"))


par(mfcol=c(2,2))
zmax <- max(abs(M1$matrix.fig), abs(M2$matrix.fig))

alpha <- .1
Fig5.plot(M1, lab=c("A","B"))
Fig5.plot(M2, lab=c("C","D"))


pdf(file=paste0("Fig5 landscape r=",r, " cc=", cc, " b=", b, " n.obs=", n.obs, " alpha=", alpha, " seed=",seed," 22Nov19.pdf"), height=8, width=8)
par(mfcol=c(2,2))

zmax <- max(abs(M1$matrix.fig), abs(M2$matrix.fig))
alpha <- .1

Fig5.plot(M1, lab=c("A","B"))
Fig5.plot(M2, lab=c("C","D"))

dev.off()

M1$FF
M1$p.FF
M1$coef0/M1$se0
2*min(pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=T),pt(abs(M1$coef0/M1$se0), df=10000, lower.tail=F))

anova(lm(matrix(M1$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
summary(lm(matrix(M1$matrix.fig,ncol=1) ~ 1))


M2$FF
M2$p.FF
M2$coef0/M2$se0
2*min(pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=T),pt(abs(M2$coef0/M2$se0), df=10000, lower.tail=F))

anova(lm(matrix(M2$matrix.fig,ncol=1) ~ as.factor(matrix(landscape,ncol=1))))
summary(lm(matrix(M2$matrix.fig,ncol=1) ~ 1))
