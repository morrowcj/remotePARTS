# Test AR functions

## load and prepare test data
data("ndvi_AK3000")
X = as.matrix(ndvi_AK3000[1:5, -c(1:6)]) # time series for 20 pixels
t = 1:ncol(X); npix = nrow(X)
x = unlist(X[1, ]) # 1 pixel
U = stats::model.matrix(formula(x ~ t))

# single-pixel ----
fm.pix = fitAR(x ~ t)
## classes
expect_equal(class(fm.pix), c("remoteAR", "pixel"), info = "check class")
## Values
pix.coef = rbind(c(7.02561924, 0.2162664293, 15.107411,  1.450210e-15),
                 c(-0.04416808, 0.0005669092,  1.855033, 7.344275e-02))
expect_equivalent(as.matrix(fm.pix$coef), pix.coef, tolerance = 1e-7,
                  info = "check regression output")
expect_equivalent(fm.pix$b, 0.681149, tolerance = 1e-7,
                  info = "check AR parameter")

# multi-pixel ----
AR.map = fitAR.map(X, t)
## class
expect_equal(class(AR.map), c("remoteAR", "map"), info = "check class")
## time coef
time.coef = matrix(data = c(-0.044168076, 5.669092e-04, 1.8550334, 0.073442747,
                             0.002635761, 2.005760e-05, 0.5885272, 0.560583757,
                             0.016828823, 9.220295e-05, 1.7525941, 0.089887730,
                            -0.091918154, 9.229887e-04, 3.0255414, 0.005054021,
                            -0.010345184, 2.316411e-05, 2.1494652, 0.039784780),
                   nrow = npix, byrow = TRUE)
expect_equivalent(as.matrix(AR.map$time.coef), time.coef, tolerance = 1e-7,
                  info = "check time coefficients")
AR.par = c(0.681148974, 0.178839083, 0.028116052, 0.291800648, -0.002310766)
expect_equivalent(AR.map$AR.par, AR.par, tolerance = 1e-7,
                  info = "check AR parameters")


# With covariates ----

## TBA
