library(remotePARTS)
data(ndvi_AK10000)
df = ndvi_AK10000[seq_len(1000), ]
pm = sample_partitions(nrow(df), npart = 3)

proc_ratio <- function(expr){
  tab = system.time(expr)
  ratio = (tab["user.self"] + tab["sys.self"]) / tab["elapsed"]
  return(ratio)
}

test_1 <- proc_ratio(fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df,
                                      nugget = 0))
test_2 <- proc_ratio(fitGLS_partition(formula = CLS_coef ~ lat, partmat = pm, data = df,
                                      nugget = 0))
test_3 <- proc_ratio(fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df,
                                      nugget = NA))
test_4 <- proc_ratio(multicore_fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm,
                                                data = df, nugget = 0, ncores = 2))
test_5 <- proc_ratio(fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df,
                                      nugget = 0, ncores = 2, parallel = TRUE))
test_6 <- proc_ratio(fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df,
                                      nugget = 0, ncores = 2, parallel = FALSE))

all_tests <- c(test_1, test_2, test_3, test_4, test_5, test_6)

all_tests
any(all_tests >= 2.5)

