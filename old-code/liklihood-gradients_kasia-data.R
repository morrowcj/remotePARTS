# Log-liklihood gradient for Kasia's parititons

## Load libraries and data
library(data.table);library(ggplot2);library(remotePARTS);library(dplyr);library(geosphere)
load("old-code/from_Kasia/soil_3rd_GLS_subset4Clay.RData", verbose = TRUE)

## Setup paramter table
each.par = 10
nug.vec = seq(0, 1, length.out = each.par)
r.vec = seq(1e-9, 1, length.out = each.par)
a.vec = seq(1e-9, 2, length.out = ceiling(each.par/2))
par.list <- expand.grid(r = r.vec, a = a.vec, nug = nug.vec)

## Setup output
LL.mat <- matrix(NA, nrow = nrow(par.list), ncol = ncol(parts_sub),
                 dimnames = list(rownames(par.list), colnames(parts_sub)))
## Setup progress bar and iterator
pb <- txtProgressBar(min = 0, max = nrow(par.list)*ncol(parts_sub), style = 3)
iter = 0

# Main Loop
for(i in seq_len(ncol(parts_sub))){
  ## partition
  part <- parts_sub[, i]
  indx = seq(from = (2000*(i - 1)) + 1, to = (2000*(i - 1)) + 2000)
  ## data
  df.part <- dat.df0_sub[indx, ]
  loc <- df.part[, c("lng", "lat")]
  D <- geosphere::distm(loc)

  for(j in seq_len(nrow(par.list))){
    ## parameters
    r = par.list[j, "r"]
    a = par.list[j, "a"]
    nug = par.list[j, "nug"]
    ## variance matrix
    V <- fitV(Dist = D, spatialcor = c("r" = r, "a" = a),
              method = "exponential-power")
    ## Log-Liklihood
    LL.mat[j, i] <- fitGLS2(formula = relEst ~ 0 + Grassland,
                            V = V, nugget = nug, data = df.part, LL_only = TRUE)
    ## increment iterator and progress bar
    iter = iter + 1
    setTxtProgressBar(pb, iter)
  }

}

LL.dat <- cbind(par.list, LL = LL.mat)

write.csv(LL.dat, "liklihood-gradient_results.csv")

ggplot(LL.dat, aes(x = r, y = nug, z = LL.part.3, col = LL.part.3)) +
  geom_point() +
  facet_wrap(~ round(a, 2)) +
  scale_color_gradient2(low = "blue", high = 'red', midpoint = 5500.5)

 # ## Partition
# part.num = c(3,4,6,10,11)
# i = 1; j = 1
# ## Correlation
# resids <- get(paste0("resids_p", part.num[i]))
# # cor.resid <- cor(t(resids))
# # vec_cor = cor.resid[upper.tri(cor.resid, diag = TRUE)]
#
# ## Distance
#
# # vec_dist <- scales::rescale(D[upper.tri(D, diag = TRUE)], to = c(0, 1))

# # Spatial correlation function
# spcor.func = function(d, r, a, nug){
#   return((1 - nug)*exp(-(d/r)^a))
# }
#
# par.list$z.test <- with(par.list, .2*nug^4 + .3*r^2 + .1*a^3)
#
# ## z as function of r and a
# ggplot(par.list, aes(x = r, y = a)) +
#   geom_contour_filled(aes(z = z.test))
# ## z as a function of r and nug
# ggplot(par.list, aes(x = r, y = nug)) +
#   geom_contour_filled(aes(z = z.test))
# ## z as a function of a and nug
# ggplot(par.list, aes(x = nug, y = a)) +
#   geom_contour_filled(aes(z = z.test))
