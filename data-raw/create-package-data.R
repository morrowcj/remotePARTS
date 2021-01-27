# Create data for remotePARTS package
  ## All paths are relative to the romotePARTS root directory
library(data.table)


## Full data ----
# read in the data
ndvi <- data.table::fread(file = "data-raw/north_america_checked.csv")

# convert $land column to a factor
land.classes = c("Evergr needle","Evergr broad","Decid needle","Decid broad",
                 "Mixed forest","Shrubland","Savanna","Grassland","Cropland",
                 "Cropland mosaics") # all land classes (in order of lvs)

ndvi$land <- factor(ndvi$land, labels = land.classes)

# save the rdata file
if(!file.exists("data-raw/ndvi.rda")){
  save(ndvi, file =  "data-raw/ndvi.rda", compress = "xz")
}

## Alaska subset ----
ndvi_AK <- ndvi[ndvi$lng < -141, ] #North Am. west of -141 is approx. AK
ndvi_AK$land <- droplevels(ndvi_AK$land)

# relative contirbution of each land class to AK
tmp <- aggregate(
  x = ndvi_AK$land,
  by = list(landclass = ndvi_AK$land),
  FUN = function(x)length(x)/nrow(ndvi_AK)
)

# any land class that occurs less than 2% is rare
ndvi_AK$rare.land <- tmp[match(ndvi_AK$land, tmp$landclass), "x"] <= 0.02

# reorder columns so ndvi are all at the end
ord <- c(grep("ndvi", names(ndvi_AK), invert = TRUE), #columns without 'ndvi'
         grep("ndvi", names(ndvi_AK))) # cols with ndvi
ndvi_AK <- ndvi_AK[, ..ord] #..x is data.table notation

# save AK data file
if(!file.exists("data/ndvi_AK.rda")){
  save(ndvi_AK, file = "data/ndvi_AK.rda", compress = "xz")
}

## subset AK further
ndvi_AK3000 <- ndvi_AK[!ndvi_AK$rare.land, ] # only common land classes
ndvi_AK3000 <- ndvi_AK3000[sample.int(n = nrow(ndvi_AK3000), size = 3000), ] # 3000 pts

# save small AK data file
if(!file.exists("data/ndvi_AK3000.rda")){
  save(ndvi_AK3000, file = "data/ndvi_AK3000.rda", compress = "xz")
}

## Save ndvi_AK as .csv file
## added 27-Jan-2021
if (!file.exists("data-raw/AK_ndvi_common-land.csv")) {
  load("data/ndvi_AK.rda")
  AK = ndvi_AK[!ndvi_AK$rare.land, ] # remove rare land classes
  AK$land = droplevels(AK$land) # drop unused land classes
  write.csv(AK, "data-raw/AK_ndvi_common-land.csv", row.names = FALSE)
}
