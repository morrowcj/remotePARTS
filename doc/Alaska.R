## ----chunk_setup, include = FALSE---------------------------------------------
library(knitr)
knitr::opts_chunk$set(collapse = TRUE, comment = "##", dpi = 50)

## ----public_library-----------------------------------------------------------
library(remotePARTS)

## ----private_library----------------------------------------------------------
library(dplyr)
library(ggplot2)

## ---- echo = FALSE------------------------------------------------------------
## set default ggplot theme
theme_set(theme(panel.grid = element_blank(), 
                strip.background = element_blank(),
                panel.background = element_blank(), 
                text = element_text(size = 10),
                legend.text = element_text(size = 10), 
                axis.text = element_text(size = 10)
                ))

## -----------------------------------------------------------------------------
data("ndvi_AK")
data("ndvi_AK3000")

## -----------------------------------------------------------------------------
str(ndvi_AK)

## ---- fig.width = 5, fig.asp = .4---------------------------------------------
reshape2::melt(ndvi_AK, measure = c("ndvi1982", "ndvi1998", "ndvi2013")) %>% 
  ggplot(aes(x = lng, y = lat, col = value )) + 
  geom_tile() +
  labs(col = "ndvi") +
  facet_wrap(~ gsub("ndvi", "", variable), ncol = 3) +
  scale_color_viridis_c(option = "magma") +
  labs(x = "Longitude", y = "Latitude")

## ----fig.width = 3, fig.asp = .9----------------------------------------------
ndvi_AK %>% filter(rare.land == FALSE) %>% 
ggplot(aes(x = lng, y = lat, fill = land, col = land)) + 
  geom_tile() + 
  scale_fill_viridis_d(direction = -1, end = .9) + 
  scale_color_viridis_d(direction = -1, end = .9) +
  labs(y = "Latitude", x = "Longitude", col = "Land cover", fill = "Land cover")

## -----------------------------------------------------------------------------
ndvi.cols = grep("ndvi", names(ndvi_AK3000), value = TRUE)
Y = as.matrix(ndvi_AK3000[, ndvi.cols])

## -----------------------------------------------------------------------------
coords = as.matrix(ndvi_AK3000[, c("lng", "lat")])

## -----------------------------------------------------------------------------
ARfit = fitAR_map(Y = Y, coords = coords)

## -----------------------------------------------------------------------------
head(coefficients(ARfit))
ndvi_AK3000$AR_coef = coefficients(ARfit)[, "t"] # save time trend coefficient

## ----fig.width = 3, fig.asp = .9----------------------------------------------
ndvi_AK %>% filter(rare.land == FALSE) %>% 
  ggplot(aes(x = lng, y = lat, fill = AR_coef, col = AR_coef)) + 
  geom_tile() + 
  scale_color_gradient2(high = "forestgreen", low = "darkorchid", 
                        mid = "grey90", midpoint = 0) +
  guides(fill = "none") + 
  labs(y = "Latitude", x = "Longitude", col = expression(beta[1]))

## -----------------------------------------------------------------------------
D = distm_scaled(coords)

## ---- fig.asp = 1, fig.width = 3, include = FALSE-----------------------------
curve(covar_exp(x, r = .1), xlab = "distance", ylab = "covar_exp(d, r)")
curve(covar_exp(x, r = .2), add = TRUE, col = "red")
legend("topright", legend = c("0.1", "0.2"), title = "r", 
       col = c("black", "red"), lty = 1)

## ---- fig.asp = 1, fig.width = 3, include = FALSE-----------------------------
nugget = .3
curve(covar_exp(x, r = .1), xlab = "distance", ylab = "covar_exp(d, r)")
curve((1-nugget)*covar_exp(x, r = .1), add = TRUE, col = "red")
legend("topright", legend = c(0, .2), title = expression(eta), 
       col = c("black", "red"), lty = 1)

## -----------------------------------------------------------------------------
r = 0.1
V = covar_exp(D, r)

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  image(V)

## -----------------------------------------------------------------------------
nugget = 0.2 
I = diag(nrow(V)) # identity matrix
Sigma = nugget*I + (1-nugget)*V

## -----------------------------------------------------------------------------
GLS.0 <- fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = V, nugget = nugget)

## ---- eval = FALSE------------------------------------------------------------
#  fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = Sigma, nugget = 0) # equivalent

## -----------------------------------------------------------------------------
coefficients(GLS.0)

## -----------------------------------------------------------------------------
corfit = fitCor(resids = residuals(ARfit), coords = coords, 
                covar_FUN = "covar_exp", start = list(range = 0.1), 
                fit.n = 3000)
(range.opt = corfit$spcor)

## ---- eval = FALSE------------------------------------------------------------
#  max.dist <- max(distm_km(coords))
#  corfit.km = fitCor(resids = residuals(ARfit), coords = coords,
#                     covar_FUN = "covar_exp", start = list(range = max.dist*0.1),
#                     distm_FUN = "distm_km", fit.n = 3000)

## -----------------------------------------------------------------------------
V.opt = covar_exp(D, range.opt)

## -----------------------------------------------------------------------------
GLS.opt = fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = V.opt, nugget = NA, 
                 no_F = FALSE)
(nug.opt = GLS.opt$nugget)
coefficients(GLS.opt)

## -----------------------------------------------------------------------------
rbind(GLS.0 = c(range = r, nugget = GLS.0$nugget, LL = GLS.0$LL, MSE = GLS.0$MSE),
      GLS.opt = c(range = range.opt, nugget = GLS.opt$nugget, LL = GLS.opt$LL, MSE = GLS.opt$MSE))

## ---- eval = FALSE------------------------------------------------------------
#  fitopt <- optimize_GLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000,
#                         coords = ndvi_AK3000[, c("lng", "lat")],
#                         covar_FUN = "covar_exp",
#                         start = c(range = .1, nugget = .2))
#  fitopt$opt$par
#  ##      range     nugget
#  ## 0.02907116 0.30272255
#  fitopt$GLS$LL
#  ## [1] 6941.524
#  fitopt$GLS$MSE
#  ## [1] 0.001033084

## -----------------------------------------------------------------------------
(GLS.int <- fitGLS(AR_coef ~ 1, data = ndvi_AK3000,V = V.opt, nugget = nug.opt))

## -----------------------------------------------------------------------------
GLS.opt

## -----------------------------------------------------------------------------
(GLS.lat <- fitGLS(AR_coef ~ 1 + lat, data = ndvi_AK3000, V = V.opt, nugget = nug.opt,
                   no_F = FALSE))

## -----------------------------------------------------------------------------
df = ndvi_AK[!ndvi_AK$rare.land, ]

## -----------------------------------------------------------------------------
pm <- sample_partitions(npix = nrow(df), partsize = 1500)
dim(pm)

## ---- eval = FALSE------------------------------------------------------------
#  partGLS_ndviAK <- fitGLS_partition(formula = AR_coef ~ 0 + land, data = df,
#                                     partmat = pm, covar_FUN = "covar_exp",
#                                     covar_pars = list(range = range.opt),
#                                     nugget = nug.opt)

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  save(partGLS_ndviAK, file = "data/partGLS_ndviAK.rda", compress = "xz")

## -----------------------------------------------------------------------------
data(partGLS_ndviAK)

## -----------------------------------------------------------------------------
partGLS_ndviAK$overall$t_test

## -----------------------------------------------------------------------------
partGLS_ndviAK$overall$pval.chisqr

## -----------------------------------------------------------------------------
(stats = with(partGLS_ndviAK$overall, c(Fmean = Fstat, rSSR = rSSR, dfs[1], partdims["npart"])))

## use function
remotePARTS:::part_chisqr(stats["Fmean"], stats["rSSR"], stats["df1"], stats["npart"])

## test manually
(rZ = (stats["rSSR"]/stats["df1"])^0.5)
(v.MSR = diag(stats["df1"]) - rZ)
V.MSR = kronecker(diag(stats["npart"]), v.MSR) + rZ
(lambda = eigen(V.MSR)$values)
(pval = CompQuadForm::imhof(q = stats["npart"] * stats["df1"] * stats["Fmean"], lambda = lambda)$Qq)
(pval.print = ifelse(pval < 1e-6, 1e-6, pval))

