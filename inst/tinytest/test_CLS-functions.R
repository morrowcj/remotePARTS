# test fitCLS and cls_star

## load test data
load("data-for-tests.rda")
load("cls-output.rda")

## single-pixel
pixel.cls <- fitCLS(X.small[1, ], 1:ncol(X.small))
expect_equal(pixel.cls, cls.out,
             info = paste("test if the output has changed"))

## multi-pixel
full.CLS <- cls_star(X.small, t = 1:ncol(X.small))
expect_equal(full.CLS, clsstar.out,
             info = paste("test if the output has changed"))

## TBA
expect_true(FALSE, info = "need to better implement and test use of covariates")
expect_true(FALSE, info = "need re-write: switch to using fastLM()")
expect_true(FALSE, info = "need to implement a starmod.cls class and methods")
expect_true(FALSE, info = "need to add error handling")
