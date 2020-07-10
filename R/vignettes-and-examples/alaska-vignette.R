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

Xmat <-  as.matrix(dat[, -c(1:6)]) #rows = location, cols = time

n = dim(Xmat)[1]
t = dim(Xmat)[2]
t.scale = ((1:t) - min(1))/t

fit.timemod <- function(x, t){
  x_t = x[2:length(x)] #X_{t}
  x_t0 = x[1:(length(x) - 1)] #X_{t-1}
  time = t[2:length(x)] #T_{t}
  fm <- lm(x_t ~ x_t0 + time)
  return(coef(fm))
}

for.vers <- app.vers <- matrix(ncol = 5, nrow = 100)
for (tests in 1:50){
  app.vers[tests, ] <- system.time(t(apply(Xmat, MARGIN = 1, FUN = fit.timemod, t = t.scale)))


  for.vers[tests, ] <- system.time(
  for (i in 1:n) {
    fit.timemod(Xmat[i, ], t.scale)
  })
}
rbind(colMeans(for.vers, na.rm = TRUE), colMeans(app.vers, na.rm = TRUE))
rbind(colSums(for.vers, na.rm = TRUE), colSums(app.vers, na.rm = TRUE))
