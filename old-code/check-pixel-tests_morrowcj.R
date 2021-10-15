library(remotePARTS)

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

# Run on same dataset ----
## make reproducible
set.seed(156)
## load data
dataAK <- read.csv(system.file("extdata", "AK_ndvi_common-land.csv", package = "remotePARTS"))
## create partitions
parts <-  sample_partitions(npix = nrow(dataAK), npart = 5, partsize = 2000)
## estimate spatial paramters (from 1 partition)
r.est <- fit_spatialcor(as.matrix(dataAK[parts[,1], -c(1:7)]),
                        t = 1:(ncol(dataAK) - 7),
                        location = dataAK[parts[,1], c("lng", "lat")])$spatialcor
         #300.8817 (ignores seed for some reason)
GLS <- fitGLS.partition.mc(part_f = "part_data", V.meth = "exponential", spatcor = 300.88,
                        partsize = nrow(parts), npart = ncol(parts), ncores = 4,
                        part_df = dataAK, part_mat = parts,
                        part_form = "cls.coef ~ lat")
## Tests
chi1 <- cor_chisq.test(GLS) # remotePARTS
  # 0.148542
chi2 <- cor_chisq(Fmean = GLS$overall.stats$meanstats["fmean"],
                  rSSR = GLS$overall.stats$meanstats["rSSRmean"],
                  df1 = GLS$overall.stats$dfs[1],
                  npart  = nrow(GLS$part.stats$coefficients))
  # 0.148542
chi3 <- correlated.chisq(Fmean.obs = GLS$overall.stats$meanstats["fmean"],
                         rSSR = GLS$overall.stats$meanstats["rSSRmean"],
                         df1 = GLS$overall.stats$dfs[1],
                         npart = nrow(GLS$part.stats$coefficients))
  # 0.0148542
t1 <- cor_t.test(GLS)  # remotePARTS
#                      Est           SE     t.stat     pval.t
# (Intercept) -0.091747680 0.1813240797 -0.5059873 0.61287673
# lat          0.001445664 0.0006038048  2.3942577 0.01667236
t2 <- cor_t(coefs = GLS$overall.stats$coefmean,
            part.SEs = GLS$part.stats$SEs,
            rcoef = GLS$overall.stats$rcoefmean,
            df2 = GLS$overall.stats$dfs[2],
            npart = nrow(GLS$part.stats$coefficients))
#                      Est           SE     t.stat     pval.t
# (Intercept) -0.091747680 0.1813240797 -0.5059873 0.61287673
# lat          0.001445664 0.0006038048  2.3942577 0.01667236
t3 <- correlated.t(coef = GLS$overall.stats$coefmean,
                   se.part = t(GLS$part.stats$SEs),
                   rcoef = GLS$overall.stats$rcoefmean,
                   df2 = GLS$overall.stats$dfs[2],
                   npart = nrow(GLS$part.stats$coefficients))
#                     coef           se     tscore          P
# (Intercept) -0.091747680 0.1813240797 -0.5059873 0.61292136
# lat          0.001445664 0.0006038048  2.3942577 0.01674567


# The differences between my code and tony's code are in the calculation of rSSE
## setup some variables to share
mincross = 6
npart = 5
## Within Tony's GLS.partition(), cross-partitions are handled by:
crosses = ifelse(mincross <= npart, mincross, npart)
combos.T = NULL
for(i in 1:(crosses - 1)) for (j in (i + 1):crosses){
  combos.T = rbind(combos.T, c(i, j))
}
combos.T
#       [,1] [,2]
# [1,]     1    2
# [2,]     1    3
# [3,]     1    4
# [4,]     1    5
# [5,]     2    3
# [6,]     2    4
# [7,]     2    5
# [8,]     3    4
# [9,]     3    5
# [10,]    4    5

## Whereas my code handles it by:
possible.combos <- t(utils::combn(npart, 2))
maxcross = nrow(possible.combos)
if(mincross < maxcross){
  used.combos = possible.combos[sample(nrow(possible.combos), mincross), ]
} else {
  used.combos = possible.combos[sample(maxcross, mincross), ]
}
used.combos
#      [,1] [,2]
# [1,]    2    4
# [2,]    1    3
# [3,]    1    4
# [4,]    1    2
# [5,]    2    3
# [6,]    3    5

## Because we use different partitions (and, sometimes, a different number of partitions),
  ## we end up with different rSSE.
