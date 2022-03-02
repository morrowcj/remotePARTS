## ----load_remotePARTS---------------------------------------------------------
library(remotePARTS)

## ----load_tidyverse-----------------------------------------------------------
library(dplyr)
library(ggplot2)

## ----set_ggtheme, echo = FALSE------------------------------------------------
## set default ggplot theme
theme_set(theme(panel.grid = element_blank(), 
                strip.background = element_blank(),
                panel.background = element_blank(), 
                text = element_text(size = 10),
                legend.text = element_text(size = 10), 
                axis.text = element_text(size = 10)
                ))

## ----load_data----------------------------------------------------------------
data("ndvi_AK10000")
ndvi_AK3000 <- ndvi_AK10000[seq_len(3000),] # first 3000 pixels from the random 10K

## ----data_structure-----------------------------------------------------------
str(ndvi_AK10000)

## ----plot_ndvi_time, fig.width = 6.5, fig.asp = .4----------------------------
reshape2::melt(ndvi_AK10000, measure = c("ndvi1982", "ndvi1998", "ndvi2013")) %>% 
  ggplot(aes(x = lng, y = lat, col = value )) + 
  geom_point(size = .1) +
  labs(col = "ndvi") +
  facet_wrap(~ gsub("ndvi", "", variable), ncol = 3) +
  scale_color_viridis_c(option = "magma") +
  labs(x = "Longitude", y = "Latitude")

## ----plot_land, fig.width = 4.5, fig.asp = .8---------------------------------
ndvi_AK10000 %>%  
ggplot(aes(x = lng, y = lat, col = land)) + 
  geom_point(size = .1) + 
  scale_color_viridis_d(direction = -1, end = .9) +
  labs(y = "Latitude", x = "Longitude", col = "Land cover", fill = "Land cover")

## ----extract_Y----------------------------------------------------------------
ndvi.cols <- grep("ndvi", names(ndvi_AK3000), value = TRUE)
Y <- as.matrix(ndvi_AK3000[, ndvi.cols])

## ----extract_coords-----------------------------------------------------------
coords <- as.matrix(ndvi_AK3000[, c("lng", "lat")])

## ----time_series_regression---------------------------------------------------
ARfit <- fitAR_map(Y = Y, coords = coords)

## ----coefficient_slice--------------------------------------------------------
head(coefficients(ARfit))

## ----standardize_coefficients-------------------------------------------------
ARfit$coefficients[, "t"] <- ARfit$coefficients[,"t"]/rowMeans(ndvi_AK3000[, ndvi.cols])
ndvi_AK3000$AR_coef <- coefficients(ARfit)[, "t"] # save time trend coefficient

## ----plot_ndvi_trend, fig.width = 4.5, fig.asp = .8---------------------------
ndvi_AK10000 %>% 
  ggplot(aes(x = lng, y = lat, col = AR_coef)) + 
  geom_point(size = .1) + 
  scale_color_gradient2(high = "red", low = "blue", 
                        mid = "grey90", midpoint = 0) + 
  guides(fill = "none") + 
  labs(y = "Latitude", x = "Longitude", col = expression(beta[1]))

## ----calc_distance------------------------------------------------------------
D <- distm_scaled(coords)

## ----visualize_range, fig.asp = 1, fig.width = 4.5, include = FALSE-----------
curve(covar_exp(x, r = .1), xlab = "distance", ylab = "covar_exp(d, r)")
curve(covar_exp(x, r = .2), add = TRUE, col = "red")
legend("topright", legend = c("0.1", "0.2"), title = "r", 
       col = c("black", "red"), lty = 1)

## ----visualize_nugget, fig.asp = 1, fig.width = 4.5, include = FALSE----------
nugget <- .3
curve(covar_exp(x, r = .1), xlab = "distance", ylab = "covar_exp(d, r)")
curve((1-nugget)*covar_exp(x, r = .1), add = TRUE, col = "red")
legend("topright", legend = c(0, .2), title = expression(eta), 
       col = c("black", "red"), lty = 1)

## ----calc_covariance----------------------------------------------------------
r <- 0.1
V <- covar_exp(D, r)

## ----visualize_V, eval = FALSE, echo = FALSE----------------------------------
#  image(V)

## ----add_nugget---------------------------------------------------------------
nugget <- 0.2 
I <- diag(nrow(V)) # identity matrix
Sigma <- nugget*I + (1-nugget)*V

## ----fit_land_GLS-------------------------------------------------------------
GLS.0 <- fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = V, nugget = nugget)

## ----alt_fit_land_GLS, eval = FALSE-------------------------------------------
#  fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = Sigma, nugget = 0) # equivalent

## ----extract_GLS_coefs--------------------------------------------------------
coefficients(GLS.0)

## ----print_GLS----------------------------------------------------------------
GLS.0

## ----fit_range_parameter------------------------------------------------------
corfit <- fitCor(resids = residuals(ARfit), coords = coords, 
                covar_FUN = "covar_exp", start = list(range = 0.1), 
                fit.n = 3000)
(range.opt = corfit$spcor)

## ----fit_range_km, eval = FALSE-----------------------------------------------
#  max.dist <- max(distm_km(coords))
#  corfit.km <- fitCor(resids = residuals(ARfit), coords = coords,
#                     covar_FUN = "covar_exp", start = list(range = max.dist*0.1),
#                     distm_FUN = "distm_km", fit.n = 3000)

## ----optimized_covariance-----------------------------------------------------
V.opt <- covar_exp(D, range.opt)

## ----optimized_nugget---------------------------------------------------------
GLS.opt <- fitGLS(formula = AR_coef ~ 0 + land, data = ndvi_AK3000, V = V.opt, nugget = NA, 
                 no.F = FALSE)
(nug.opt = GLS.opt$nugget)
coefficients(GLS.opt)

## ----compare_GLS--------------------------------------------------------------
rbind(GLS.0 = c(range = r, nugget = GLS.0$nugget, logLik = GLS.0$logLik, MSE = GLS.0$MSE),
      GLS.opt = c(range = range.opt, nugget = GLS.opt$nugget, logLik = GLS.opt$logLik, MSE = GLS.opt$MSE))

## ----optimized_GLS, eval = FALSE----------------------------------------------
#  fitopt <- fitGLS_opt(formula = AR_coef ~ 0 + land, data = ndvi_AK3000,
#                         coords = ndvi_AK3000[, c("lng", "lat")],
#                         covar_FUN = "covar_exp",
#                         start = c(range = .1, nugget = .2))
#  fitopt$opt$par
#  ##      range     nugget
#  ## 0.03660171 0.26859993
#  fitopt$GLS$logLik
#  ## [1] 13276.15
#  fitopt$GLS$MSE
#  ## [1] 1.720775e-05

## ----GLS_intercept_only-------------------------------------------------------
(GLS.int <- fitGLS(AR_coef ~ 1, data = ndvi_AK3000, V = V.opt, nugget = nug.opt, no.F = TRUE))

## ----GLS_land_optimized-------------------------------------------------------
GLS.opt

## ----fit_GLS_latitude---------------------------------------------------------
(GLS.lat <- fitGLS(AR_coef ~ 1 + lat, data = ndvi_AK3000, V = V.opt, nugget = nug.opt,
                   no.F = FALSE))

## ----remove_rare_land---------------------------------------------------------
df = ndvi_AK10000

## ----partition_matrix---------------------------------------------------------
pm <- sample_partitions(npix = nrow(df), partsize = 1500, npart = NA)
dim(pm)

## ----partitioned_GLS, eval = FALSE--------------------------------------------
#  partGLS_ndviAK <- fitGLS_partition(formula = AR_coef ~ 0 + land, data = df,
#                                     partmat = pm, covar_FUN = "covar_exp",
#                                     covar.pars = list(range = range.opt),
#                                     nugget = nug.opt)

## ----load_partitioned_GLS-----------------------------------------------------
data(partGLS_ndviAK)

## ----partitioned_t_test-------------------------------------------------------
partGLS_ndviAK$overall$t_test

## ----chisqr_test--------------------------------------------------------------
partGLS_ndviAK$overall$pval.chisqr

