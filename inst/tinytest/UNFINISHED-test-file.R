


bench <- microbenchmark(
  new.full = remoteSTAR::fitGLS2(formula = test.y ~ 0 + test.land,
                                 data = test.df, V =  test.V,
                                 nugget = 0, form.0 = NULL, save_xx = FALSE,
                                 threads = 1, contrasts = NULL),
  old.full = {
    modmat = model.matrix(~ 0 + test.land)
    remoteSTAR::fitGLS(X = modmat, V = test.V, y = test.y,
           X0 = cbind(rep(1, nrow(test.modmat))),
           nugget = 0, save_xx = FALSE,
           threads = 1)
  },
  old.R = {
    modmat = model.matrix(~ 0 + test.land)
    remoteSTAR::fitGLS_R(X = modmat, V = test.V, y = test.y,
             X0 = cbind(rep(1, nrow(test.modmat))),
             nugget = 0)
  }, times = 10L)

bench
ggplot2::autoplot(bench)

X0 = model.matrix(test.y ~ 1)
tmp.GLS = remoteGLS()

bench2 <- microbenchmark(
  new.cpp = {
    # tmp.GLS = remoteGLS()
    remoteSTAR:::.fitGLS2_cpp(L = tmp.GLS, X = test.modmat, V = test.V, y = test.y, X0 = X0,
                 nugget = 0, save_xx = FALSE, threads = 1)
  },
  old.cpp = {
    remoteSTAR:::.fitGLS_cpp(X = test.modmat, V = test.V, y = test.y, X0 = X0,
                nugget = 0, save_xx = FALSE, threads = 1)
  },
  times = 10L)

bench2
ggplot2::autoplot(bench2)


