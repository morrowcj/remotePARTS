time_fun <- function(expr){
  tab = system.time(expr)
  print(tab)
  ratio = (tab["user.self"] + tab["sys.self"]) / tab["elapsed"]
  names(ratio) = "ratio"
  print(ratio)
}

## PartGLS funs

library(remotePARTS)
## read data
data(ndvi_AK10000)
df = ndvi_AK10000[seq_len(1000), ] # first 1000 rows

## create partition matrix
pm = sample_partitions(nrow(df), npart = 3)

## fit GLS with fixed nugget
time_fun({partGLS = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df, nugget = 0, do.t.test = TRUE)})

## hypothesis tests
chisqr(partGLS) # explanatory power of model
t.test(partGLS) # significance of predictors

## now with a numeric predictor
time_fun(fitGLS_partition(formula = CLS_coef ~ lat, partmat = pm, data = df, nugget = 0))

## fit ML nugget for each partition (slow)
time_fun({partGLS.opt = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df, nugget = NA)})
partGLS.opt$part$nuggets # ML nuggets

## Explicitly use multicore_fitGLS_partition()
time_fun(multicore_fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df, nugget = 0, ncores = 2L))

time_fun(multicore_fitGLS_partition(formula = CLS_coef ~ 1, partmat = pm, ncores = 2L, data = df, nugget = 0, do.chisqr.test = FALSE))

## fully parallel, using 2 cores
time_fun({MC_GLSpart = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm, data = df, nugget = 0, ncores = 2, parallel = TRUE, debug = FALSE)})

time_fun(fitGLS_partition(formula = CLS_coef ~ lat, partmat = pm, data = df, nugget = 0, ncores = 2, parallel = TRUE, debug = FALSE))

time_fun(fitGLS_partition(formula = CLS_coef ~ 1, partmat = pm, data = df, nugget = 0, parallel = TRUE, ncores = 2, do.chisqr.test = FALSE))

# Certain model structures may not be useful:
## 0 intercept with numeric predictor (produces NAs) and gives a warning in statistical tests
time_fun(fitGLS_partition(formula = CLS_coef ~ 0 + lat, partmat = pm, data = df, nugget = 0))

## intercept-only, gives warning
time_fun(fitGLS_partition(formula = CLS_coef ~ 1, partmat = pm, data = df, nugget = 0, do.chisqr.test = FALSE))

data(ndvi_AK10000)
file.path = "ndviAK10000-remotePARTS.csv"
time_fun(write.csv(ndvi_AK10000, file = file.path))
time_fun(part_csv(1:20, formula = CLS_coef ~ 0 + land, file = file.path))
time_fun(part_csv(sample(3000, 20), formula = CLS_coef ~ 0 + land, file = file.path))
time_fun(file.remove(file.path))
