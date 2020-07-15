## packages used
# data.table


# prepare data ----
## This section is not needed for the user, as *.rda files will be saved

## Full data
# read in the data
dtab <- data.table::fread(file = "data-raw/north_america_checked.csv")

# convert $land collumn to a factor
land.classes = c("Evergr needle","Evergr broad","Decid needle","Decid broad",
                 "Mixed forest","Shrubland","Savanna","Grassland","Cropland",
                 "Cropland mosaics") # all land classes (in order of lvs)

dtab$land <- factor(dtab$land, labels = land.classes)

# save the rdata file
saveRDS(dtab, "data/north-america_ndvi.rda")

## Alaska subset
dataAK <- dtab[dtab$lng < -141, ] #North Am. west of -141 is approx. AK

# relative contirbution of each land class to AK
tmp <- aggregate(
  x = dataAK$land,
  by = list(landclass = dataAK$land),
  FUN = function(x)length(x)/nrow(dataAK)
)

# any land class that occurs less than 2% is rare
dataAK$rare.land <- tmp[match(dataAK$land, tmp$landclass), "x"] <= 0.02

# reorder columns so ndvi are all at the end
ord <- c(grep("ndvi", names(dataAK), invert = TRUE), #columns without 'ndvi'
         grep("ndvi", names(dataAK))) # cols with ndvi
dataAK <- dataAK[, ..ord] #..x is data.table notation

# save AK data file
saveRDS(dataAK, "data/alaska_ndvi.rda")

## subset AK further
dataAK_small <- dataAK[!dataAK$rare.land, ] # only common land classes
dataAK_small <- dataAK_small[sample.int(n = nrow(dataAK_small), size = 3000), ] # 3000 pts

# save small AK data file
saveRDS(dataAK_small, "data/alaska_ndvi_3000pts.rda")

# load data ----
##

rm(list = ls()) # clear all objects from memory
source("R/lsos.R") # load memory tracking functions

dat <- readRDS("data/alaska_ndvi_3000pts.rda")

# CLS ----
##

source("R/fitCLS.R") # source the functions

dat <- dat[!dat$rare.land, ] # remove rare land classes
# X data
Xmat <-  as.matrix(dat[, -c(1:6)]) #rows = location, cols = time
n = dim(Xmat)[1] # sites
t.n = dim(Xmat)[2] # time points
t.scale = ((1:t.n) - min(1))/t.n # scaled time

# test speed and check accuracy
if (FALSE) {
  fit <- star <- matrix(NA, 5, 5)
  for (i in 1:5) {
    fit[i, ] <- system.time(tmp1 <- CLS.fit(Xmat, t.scale))
    star[i, ] <- system.time(tmp2 <- cls_star(Xmat, t.scale))
  }
  rbind("fit" = colMeans(fit, na.rm = TRUE),
        "star" = colMeans(star, na.rm = TRUE))[, 1:3]
  tmp1[1:5, ]
  tmp2[1:5, ]
}
#       user  sys elapsed
# #fit  2.862 0   2.864
# #star 2.794 0   2.796

# note: cls_star() does not loose speed and can handle covariates:

if (FALSE) {
  # generate some random covariates (just for structure):
  # one covariate
  covar.mat1 <- matrix(rnorm(n*t.n), n, t.n)
  colnames(covar.mat1) = gsub("ndvi", "Z", colnames(Xmat))

  cls_star(Xmat[1:5, ], t.scale, covar.list = list(covar.mat1[1:5, ]))

  ## now 2 covaritates
  covar.mat2 <- matrix(rnorm(n*t.n), n, t.n)
  colnames(covar.mat2) = gsub("ndvi", "U", colnames(Xmat))

  cls_star(X = Xmat[1:5, ], t = t.scale,
           covar.list = list(covar.mat1[1:5, ], covar.mat2[1:5, ]))
}

# fit the initial model and check for outliers
tmp <- cls_star(Xmat, t.scale)
tmp[, c("lat", "lng")] <- dat[, c("lat", "lng")]
tmp$rel.stat <-  (tmp$Est/(1 - tmp$x_t0.EST))/tmp$mean # relative CLS statistic

outl = abs(scale(tmp$rel.stat)) > -qnorm(p = 1/nrow(tmp)/10) # outliers

# plot with outliers
library(ggplot2)
ggplot(tmp[!outl, ], aes(x = lng, y = lat, col = rel.stat)) +
  geom_tile(size = 2) +
  scale_color_gradient2(low = "red", mid = "grey", high = "green",
                        midpoint = 0) +
  geom_tile(data = tmp[outl, ], col = "black", size = 1) # add outliers
detach("package:ggplot2", unload = TRUE)

base.col = ifelse(
  test = outl,
  yes = "black",
  no = colorRampPalette(
    c("orange", "grey", "darkgreen")
  )(nrow(tmp))[ordered(tmp$rel.stat)]
)
plot(dat$lat ~ dat$lng, pch = 15, cex = .7, col = base.col)
