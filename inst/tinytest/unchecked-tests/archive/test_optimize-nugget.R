# test optimize_nugget()

# Setup ----
n.pix = 100 # number of pixels to test with
set.seed(916)
data("ndvi_AK3000")
subsamp = sample.int(n = nrow(ndvi_AK3000), size = n.pix)
df = ndvi_AK3000[subsamp, ] # subset the data

modmat = model.matrix(~ 0 + land, data = df)


nug.ml <- optimize_nugget(X = modmat, V = V, y = y)

# expect_true(FALSE, info = paste("need to come up with and write appropriate",
#                                 "tests of optimize_nugget()"))
#
# expect_true(FALSE, info = "need to depricate fitNugget()")
# expect_true(FALSE, info = "need to depricate fitNugget_Rcpp()")

