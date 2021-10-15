#setwd("D:beautydata/arives/RSE ms/")
# Tony's Code ----

## Chisqr
correlated.chisq <- function(Fmean.obs, rSSR, df1, npart) {

  require(CompQuadForm)

  rZ <- (rSSR/df1)^0.5
  v.MSR <- diag(df1) - rZ
  V.MSR <- kronecker(diag(npart), v.MSR) + rZ
  D.MSR <- t(chol(V.MSR, pivot = T))
  if (attr(D.MSR, "rank") < npart*df1) {
    rank.deficient.MSR <- npart*df1-attr(D.MSR, "rank")
    v.MSR <- diag(df1) - (0.99/df1)
    V.MSR <- kronecker(diag(npart), v.MSR) + (0.99/df1)
  } else {
    rank.deficient.MSR <- 0
  }

  lambda <- eigen(V.MSR)$values
  pvalue <- imhof(q = npart * df1 * Fmean.obs, lambda = lambda)$Qq
  return(list(pvalue=pvalue,  rank.deficient.MSR = rank.deficient.MSR))
}

## T test
correlated.t <- function(coef, se.part, rcoef, df2, npart) {

  secoef <- matrix(NA, length(coef), 1)
  for (i.coef in 1:length(coef)) {
    R <- (1 - rcoef[i.coef]) * diag(npart) + rcoef[i.coef] * matrix(1, npart, npart)
    secoef[i.coef, ] <- (se.part[i.coef, ] %*% R %*% se.part[i.coef, ])^0.5/npart
  }

  tscore <- coef/secoef
  # add 25Nov20
  pvalue <- 2 * pt(abs(tscore), df = sum(df2), lower.tail = F)

  ttest <- cbind(coef, secoef, tscore, pvalue)
  colnames(ttest) <- c("coef", "se", "tscore", "P")

  return(p.t = ttest)
}

## Clay's Code ----
# Chisqr
cor_chisq <- function(Fmean, rSSR, df1, npart){
  rZ <- rSSR^.5/df1
  v.MSR <- diag(df1) - rZ
  V.MSR <- kronecker(diag(npart),v.MSR) + rZ
  lambda <- eigen(V.MSR)$values
  pvalue <- suppressWarnings(CompQuadForm::imhof(q = npart * df1 * Fmean, lambda = lambda)$Qq)
  pvalue = ifelse(pvalue <= 1e-06, 1e-06, pvalue) # prevent from being negative/too low
  return(c("pval.chisqr" = pvalue))
}

# T test
cor_t <- function(coefs, part.SEs, rcoef, df2, npart){

  secoef <- matrix(NA, length(coefs), 1)
  for(i in seq_len(length(coefs))){
    R <- (1 - rcoef[i]) * diag(npart) + rcoef[i] * matrix(1, npart, npart)
    secoef[i, ] <- (part.SEs[,i] %*% R %*% part.SEs[, i])^.5/npart
  }

  tscore <- coefs/secoef
  pvalue <- 2 * pt(abs(tscore), df=df2 * npart, lower.tail = FALSE)

  ttest <- cbind(coefs, secoef, tscore, pvalue)
  colnames(ttest) <- c("Est", "SE", "t.stat", "pval.t")

  return(p.t = ttest)
}
###########################################################
###########################################################
###########################################################
###########################################################

library(Matrix)
library(geosphere)
library(colorspace)

###########################################################
# just AK(ish) but with all classes
###########################################################
land.df <- data.frame(class = c("Evergr needle","Evergr broad","Decid needle","Decid broad","Mixed forest","Shrubland","Savanna","Grassland","Cropland","Cropland mosaics", "Mixed"), num = c(1,2,3,4,5,6,8,10,12,14,16))

data <- read.csv(file="north_america_checked.csv")

# remove mixed/unstable
data <- data[data$land != 16,]
summary(data)
dim(data)
n.obs <- 34

data$landclass <- land.df$class[match(data$land, land.df$num)]

sort(unique(data$landclass))
n.classes <- aggregate(data$landclass, by=list(landscape=data$landclass), FUN=length)
n.classes

dataAK <- data[data$lng < -141,]
n.classes <- aggregate(dataAK$landclass, by=list(landscape=dataAK$landclass), FUN=length)
n.classes
      # landscape     x
# 1   Decid broad     1
# 2 Evergr needle    83
# 3     Grassland  2924
# 4  Mixed forest    36
# 5       Savanna  7361
# 6     Shrubland 10409

summary(dataAK)
nrow(data)
nrow(dataAK)
sort(unique(dataAK$land))

# landclasses: 1 Evergreen needleleaf forests 2 evergreen broadleaf forests, 3 deciduous needleleaf forests, 4 deciduous broadleaf forests, 5 mixed forests, 6 shrublands, 8 savannas, 10 grasslands, 12 croplands, 14 croplands/natural vegetation mosaics.

dim(dataAK)
dataAK$landclass <- land.df$class[match(dataAK$land, land.df$num)]
dataAK$landclass <- as.factor(dataAK$landclass)
n.classes <- aggregate(dataAK$landclass, by=list(landscape=dataAK$landclass), FUN=length)

dataAK$meanNDVI <- rowMeans(dataAK[,6:(ncol(dataAK)-2)])

write.table(dataAK, file="dataAK.csv", sep=",", row.names=F)


###########################################################
# maps

minlng <- min(dataAK$lng)
maxlng <- max(dataAK$lng)

col.pal.NDVI <- hcl.colors(100, "Greens", rev = TRUE)
cex.pt <- .4

pdf(file="Prediction Fig 1 AK 12Oct21.pdf", width=7, height=7)
#png(file="Fig 6 AK 29Dec20.png", width=1050, height=1200)
	par(mfrow = c(1,1), mai = c(0,0,.5,0))
	# # A
	# palette(col.pal.land)
	# plot(lat ~ lng, data=dataAK, pch=15, cex=.5, col=col.pal.land[landclass.num+2], xlim=c(minlng, maxlng), xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n")
	# legend(x=-152, y=60, legend=landclasses, col=col.pal.land[-(1:2)], pch=20, bty="n", cex = 1.5)
	# mtext("A", side = 3, adj=1, cex = 2)

	# B
	palette(col.pal.NDVI)
	dataAK$col <- dataAK$meanNDVI/(max(dataAK$meanNDVI) - min(dataAK$meanNDVI))
	minmax <- max(abs(min(dataAK$col)), max(dataAK$col))
	# dataAK$col <- .5*(1 + dataAK$col/minmax)
	dataAK$col <- dataAK$col/minmax
	plot(lat ~ lng, data=dataAK, pch=15, cex=cex.pt, col=col.pal.NDVI[100*dataAK$col], xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n")

	#D
	# par(mai = c(2,1,2,1))
	# hist(dat.map$b, freq=T, main="", xlab=expression(italic(b)), col="white",cex.lab=2)
	# mtext("D", side = 3, adj=1.25, padj=-6.5, cex = 2)
	# mtext("Temporal autocorrelation", side = 3, padj=-2,cex = 1.5)

dev.off()



###########################################################
# Analysis of full map
###########################################################

library(remotePARTS)
# function
optimize_GLS_TI <- function (formula, D, V.meth = "exponential-power", nugget = NA,
    spcor = NA, verbose = FALSE, contrasts = NULL, data, pars.start = c(r = 0.5, a = 1, nug = 0))
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    if (is.matrix(y)) {
        stop("response is a matrix: must be a vector")
    }
    ny <- length(y)
    modmat <- model.matrix(mt, mf, contrasts)
    rm(mf)

    pars.fixed = pars.start
    pars.fixed["r"] = spcor[1]
    if (V.meth == "exponential-power") {
        pars.fixed["a"] = spcor[2]
    }
    pars.fixed["nug"] = nugget
    pars.full = pars.start
    pars.full[!is.na(pars.fixed)] = pars.fixed[!is.na(pars.fixed)]
    pars = pars.full[is.na(pars.fixed)]
    if (verbose) {
        cat("optimizing parameters:\n")
    }
    opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat,
        D = D, verbose = verbose, control=list(reltol=1e-6), method = "Nelder-Mead")
    # opt.out = optim(pars, fn = remotePARTS:::optim_GLS_func, y = y, V.meth = V.meth, modmat = modmat,
        # D = D, verbose = verbose, control=list(reltol=1e-6), method = "L-BFGS-B", lower=c(1e-3,1e-3), upper=c(1,1))
    if (verbose) {
        cat("\n")
    }
    r.ml = opt.out$par["r"]
    a.ml = opt.out$par["a"]
    nug.ml = opt.out$par["nug"]
    if (V.meth == "exponential-power") {
        spcor.ml = c(r.ml, a.ml)
    }
    else {
        spcor.ml = r.ml
    }
    V = fitV(Dist = D/max(D), spatialcor = spcor.ml, method = V.meth)
#browser()

    GLS.out = fitGLS2(formula = y ~ 0 + modmat, V = V, nugget = nug.ml, save_xx=T)
    GLS.out$model.info$call = match.call()
    return(list(GLS.out=GLS.out, spatial.pars=list(r=r.ml, a=a.ml, nug=nug.ml)))
}

# partition data
pick <- sample.int(size=2000, n=nrow(dataAK))
data <- dataAK[pick,]

# Analysis
location <- data[,c('lng','lat')]
Dist.km <- distm(location, fun=distGeo)/1000
Dist <- Dist.km/max(Dist.km)

# mod <- optimize_GLS_TI(meanNDVI ~ lat, data=data, V.meth="exponential", D=Dist, verbose=T, pars.start = c(r = .01, a = 1, nug = .15))
# mod

r.est <- max(Dist.km) * mean(c(0.03308253, 0.03316746, 0.02512171, 0.03836761))

data.file = system.file("dataAK.csv")

parts = sample_partitions(npix = nrow(dataAK), partsize = 2000)
mod.parts <- fitGLS.partition.mc(
  part_f = "part_csv",
  dist_f = "dist_km",
  V.meth = "exponential",
  spatcor = r.est,
  part_csv_path = "dataAK.csv",
  part_mat = parts,
  part_form = "meanNDVI ~ lat",
  part_form0 = "meanNDVI ~ 1",
  partsize = nrow(parts), npart = ncol(parts),
  ncores = 3
)

# mod.parts2 <- fitGLS.partition(part_f = "part_csv",
#                                dist_f = "dist_km",
#                                V.meth = "exponential",
#                                spatcor = r.est,
#                                part_csv_path = "dataAK.csv",
#                                part_mat = parts,
#                                part_form = "meanNDVI ~ lat",
#                                part_form0 = "meanNDVI ~ 1",
#                                partsize = nrow(parts), npart = ncol(parts))

cor_chisq.test(mod.parts)
# pval.chisqr
# 0.003845231

cor_t.test(mod.parts)
                   # Est         SE     t.stat       pval.t
# (Intercept) 19.2256437 6.95215861   2.765421 5.688418e-03
# lat         -0.1860874 0.01848737 -10.065650 8.549760e-24

cor_t(coefs = mod.parts$overall.stats$coefmean,
            part.SEs = mod.parts$part.stats$SEs,
            rcoef = mod.parts$overall.stats$rcoefmean,
            df2 = mod.parts$overall.stats$dfs[2],
            npart = nrow(mod.parts$part.stats$coefficients))
                        # Est         SE     t.stat       pval.t
# (Intercept) 19.2256437 6.95215861   2.765421 5.688418e-03
# lat         -0.1860874 0.01848737 -10.065650 8.549760e-24

correlated.t(coef = mod.parts$overall.stats$coefmean,
                   se.part = t(mod.parts$part.stats$SEs),
                   rcoef = mod.parts$overall.stats$rcoefmean,
                   df2 = mod.parts$overall.stats$dfs[2],
                   npart = nrow(mod.parts$part.stats$coefficients))
                  # coef         se     tscore            P
# (Intercept) 19.2256437 6.95215861   2.765421 5.737182e-03
# lat         -0.1860874 0.01848737 -10.065650 2.781104e-23


source("remote_sensing_tools_24May21.R")

parts.list <- list(parts[,1])
for(i in 2:ncol(parts)) parts.list[[i]] <- parts[,i]

z <- GLS.partition.data(formula = "meanNDVI ~ lat", formula0 = "meanNDVI ~ 1", data = dataAK, spatial.autocor.FUN = "exponential", spatialcor = r.est, est.nugget = T, partition = parts.list, nugget.interval = c(0.0001, .9999), verbose = T)

GLS.partition.pvalue(z, doFtest = F, nboot = 1e+05)
# $p.chisq
# [1] 0.004140434

# $rank.deficient.MSR
# [1] 0

# $p.Fhochberg
# [1] 0.002117635

# $p.Fhommel
# [1] 0.002117635

# $p.Ffdr
# [1] 0.002117635

# $p.t
                  # coef         se    tscore            P
# (Intercept) 19.2256436 4.11286617  4.674512 2.959343e-06
# lat         -0.1860874 0.06533853 -2.848050 4.401794e-03

mod.parts$overall.stats
z$coef
z$rcoef
z$rSSR
z$rSSE

# > mod.parts$overall.stats
# $dfs
 # df1  df2
   # 1 1999

# $coefmean
# (Intercept)         lat
 # 19.2256437  -0.1860874

# $rcoefmean
# (Intercept)         lat
 # 2.50878150  0.00145467

# $meanstats
    # fmean  rSSRmean  rSSEmean
# 7.1559064 0.6807561 0.2242049

# > z$coef
# (Intercept)         lat
 # 19.2256436  -0.1860874
# > z$rcoef
# [1] 0.8314795 0.8410280
# > z$rSSR
# [1] 0.708616
# > z$rSSE
# [1] 0.2294176


## To correct old code: ----
# Original partition matrix:
og.parts = sample_partitions(npix = nrow(data), partsize = 2000)
# Original GLS - using mincross = 6:
old.gls <- fitGLS.partition.mc(..., part_mat = og.parts, npart = ncol(og.parts),
                               partsize = nrow(og.parts), mincross = 6)
# New GLS to get rcoef, using npart = mincross:
new.gls = fitGLS.partition.mc(..., part_mat = og.parts[1:6], npart = 6,
                              partsize = nrow(og.parts), mincross = 6)
# Update relevant statistics in original gls:
old.gls$cross.stats <- new.gls$cross.stats
old.gls$overall.stats$rcoefmean <- new.gls$overall.stats$rcoefmean
old.gls$overall.stats$meanstats["rSSRmean"] = new.gls$overall.stats["rSSRmean"]
old.gls$overall.stats$meanstats["rSSEmean"] = new.gls$overall.stats["rSSEmean"]
# re-run t-test:
cor_t.test(old.gls)
