# test CLS functions

## load and prepare test data
data("ndvi_AK3000")
X = as.matrix(ndvi_AK3000[1:5, -c(1:6)]) # time series for 20 pixels
t = 1:ncol(X); npix = nrow(X)
x = unlist(X[1, ]) # 1 pixel

# single-pixel ----
fm.pix = fitCLS(x, t)
## classes
expect_equal(class(fm.pix), c("remoteCLS", "pixel"), info = "check class")
expect_true(class(fm.pix$fm) == "lm", info = "check for lm object")
## Values
pix.coef = c("(Intercept)" = 3.11658604, "x.ti" = 0.57096874, "t.j" = -0.02443475)
expect_equal(fm.pix$fm$coefficients, pix.coef, tolerance = 1e-7,
             info = "check output")

# multi-pixel ----
CLS.map = fitCLS.map(X, t, TRUE, TRUE, TRUE, TRUE)
## class
expect_equal(class(CLS.map), c("remoteCLS", "map"), info = "check class")
## time coef
time.coef = matrix(data = c(-0.024434751, 0.011586119, -2.1089677, 0.044019823,
                             0.002639177, 0.004065018,  0.6492412, 0.521473005,
                             0.016588435, 0.010456328,  1.5864494, 0.123867675,
                            -0.084174430, 0.029977911, -2.8078818, 0.008982851,
                            -0.012529777, 0.005399981, -2.3203370, 0.027830405),
                   nrow = npix, byrow = TRUE)
expect_equivalent(as.matrix(CLS.map$time.coef), time.coef, tolerance = 1e-7,
                  info = "check time coefficients")
## AR coef (xi)
xi.coef = matrix(data = c( 0.57096874, 0.1547736,  3.6890568, 0.0009608687,
                           0.10798651, 0.1914750,  0.5639718, 0.5772631213,
                          -0.03857352, 0.1926047, -0.2002730, 0.8427145926,
                           0.20350091, 0.1820448,  1.1178616, 0.2731259808,
                          -0.06431173, 0.1874393, -0.3431070, 0.7340807702),
                 nrow = npix, byrow = TRUE)
expect_equivalent(as.matrix(CLS.map$xi.coef), xi.coef, tolerance = 1e-7,
                  info = "check AR coefficients")
## intercept
int.coef = matrix(data = c( 3.1165860, 1.1140162, 2.797613, 9.207937e-03,
                            0.9976204, 0.2319537, 4.300946, 1.866231e-04,
                           13.0221331, 2.4242401, 5.371635, 1.004959e-05,
                           11.9671851, 2.7659574, 4.326598, 1.740831e-04,
                           12.3161785, 2.1640913, 5.691155, 4.214759e-06),
                 nrow = npix, byrow = TRUE)
expect_equivalent(as.matrix(CLS.map$int.coef), int.coef, tolerance = 1e-7,
                  info = "check intercepts")

# With covariates ----

## TBA

