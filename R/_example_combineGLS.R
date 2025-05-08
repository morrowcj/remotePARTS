# ---- Load packages ----

library(tidyverse)
library(remotePARTS)

# ---- Setup and parameters ----

# load the example data
data(ndvi_AK10000)

# create partition matrix, with 20 partitions
pm = sample_partitions(npix = nrow(ndvi_AK10000), npart = 20)

# specify the number of pairs that will be considered for the cross-partition analysis
ncross = 10

# specify where the RDS objects will get saved
save_dir = "data-raw/GLS_part_obs"

# should results be overwritten if they already exist?
overwrite = TRUE

# number of cores to use
ncores = 6

# ---- Fit and save the GLS ----

# create the directory if it doesn't exist
if(!dir.exists(save_dir)) {
  dir.create(save_dir)
}

# create a counter for the number of GLS actually fit
n_GLS_fitted = 0

# fit a GLS object for each partition
for (i in seq_len(ncol(pm))) {
  # convert i into a fixed-length string, with padded zeros
  i_str = str_pad(i, 2, side = "left", pad = "0")

  # should the xx and invchol be saved for this iteration?
  save_extra = (i <= ncross)

  # convert to a string to identify the files with these
  xx_str = ifelse(save_extra, "_with-xx", "")

  # build the file name for this partition
  name_i = paste0("AK10000_GLS_results_partition", i_str, xx_str, ".rds")

  # and the full path
  file_i = file.path(save_dir, name_i)

  # skip to partition if this file already exists (unless overwrite=TRUE)
  if (file.exists(file_i) & !overwrite) {
    next
  }


  # get the column of partition matrix to use
  part_i = pm[, i]

  # subset the data by the partition
  data_i <- ndvi_AK10000[part_i, ]

  # get/set parameter values (set to constants for ease here)
  nug = 0
  rng = 0.01

  # get the coordinates for this subset
  coords_i = data_i[, c("lng", "lat")]

  # calculate the varcov matrix
  V = covar_exp(distm_km(coords_i), range = rng)

  # fit the GLS - note that only the first 10 (ncross) will have xx and invchol
  GLS_i <- fitGLS(
    CLS_coef ~ lat + land, data = data_i, V = V, nugget = nug,
    save.xx = save_extra, save.invchol = save_extra, ncores = ncores
  )

  # Add coordinates to the GLS for the partitions to cross
  if (save_extra){
    GLS_i$coords <- coords_i
  }
  ## TODO: this should be added to any function that automates this.

  # save the RDS to a file
  saveRDS(GLS_i, file = file_i)

  # increment the iterator
  n_GLS_fitted = n_GLS_fitted + 1
}

# ---- Calculate cross-partition statistics ----

# get a list of all the partitioned results that have xx
file_list = list.files(save_dir, pattern = "with-xx", full.names = TRUE)

# calculate the number of partitions to cross
ncross = length(file_list)

# get a table of all the possible crosses
cross.pairs = t(combn(seq_len(ncross), 2))

# indicate when a new GLS needs to be loaded for a given pair
load_table = rbind(c(TRUE, TRUE), diff(cross.pairs) != 0)

# set up a progress bar to track progress
pb = txtProgressBar(min = 0, max = nrow(cross.pairs), style = 3)

# build empty arrays
rSSRs = rSSEs = rep(NA, nrow(cross.pairs))
rcoefs = array(
  NA,
  dim = c(
    nrow(cross.pairs), length(GLS_i$coefficients), length(GLS_i$coefficients)
  ),
  dimnames = list(
    NULL,names(GLS_i$coefficients),names(GLS_i$coefficients)
  )
)

# loop through the cross pairs table
for (pair in seq_len(nrow(cross.pairs))) {
  # get the index of the pairs and associated files
  i = cross.pairs[pair, 1]
  path_i = file_list[i]
  j = cross.pairs[pair, 2]
  path_j = file_list[j]

  # load the GLS from disk, when it is needed
  ## i
  if (load_table[pair, 1]) {
    GLS_i = readRDS(path_i)
  }
  ## j
  if (load_table[pair, 2]) {
    GLS_j = readRDS(path_j)
  }

  # calculate the joint varcovar matrix
  dist_ij = distm_km(coords = GLS_i$coords, coords2 = GLS_j$coords)
  V_ij = covar_exp(dist_ij, range = rng)

  # calculate the degrees of freedom
  degrees_freedom = remotePARTS:::calc_dfpart(
    partsize = nrow(GLS_i$coords),
    p = length(GLS_i$coefficients),
    p0 = length(GLS_i$coefficients0)
  )

  # perform the cross
  pair_cross = remotePARTS:::crosspart_GLS(
    xxi = GLS_i$xx, xxj = GLS_j$xx,
    xxi0 = GLS_i$xx0, xxj0 = GLS_j$xx0,
    invChol_i = GLS_i$invcholV, invChol_j = GLS_j$invcholV,
    Vsub = V_ij,
    nug_i = GLS_i$nugget, nug_j = GLS_j$nugget,
    df1 = degrees_freedom[1], df2 = degrees_freedom[2],
    ncores = ncores
  )

  ## collect stats
  # TODO this could be done better
  rcoefs[pair, ,] <- pair_cross$rcoefij
  # rcoefs[cross, ] <- as.vector(rGLS$rcoefij)
  rSSRs[pair] <- ifelse(
    is.na(pair_cross$rSSRij) | is.infinite(pair_cross$rSSRij),
    NA, pair_cross$rSSRij
  )
  rSSEs[pair] <- pair_cross$rSSEij
}

## make the coefficients
rcoefficients = apply(rcoefs, MARGIN=c(2,3), FUN = function(x){mean(x, na.rm = TRUE)})
