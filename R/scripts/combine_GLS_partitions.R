# ---- Introduction ----

## This script walks a user through how to
## 1) run GLS on individual partitions and save the results to disk.
  # Note that it is important to use specific options for *at least* all of the
  # partitions that you want to cross, including the coordinates
  # (see use_coord_files and save_extra options),
## 2) calculate cross-partition statistics by loading these GLS objects two at a
  # time, and collecting the associated statistics,
## 3) collect important statistics from *all* partitions, including those not
  # used in a cross, and
## 4) perform the overall tests of the full model (t-test and chi-squared test)

# ---- Load packages ----

# tidyverse for data manipulation
library(tidyverse)
# remote parts for key functions
library(remotePARTS)

# ---- Setup and parameters ----

# load the example data to demonstrate the process
data(ndvi_AK10000)

# create partition matrix, with 20 partitions
pm = sample_partitions(npix = nrow(ndvi_AK10000), npart = 20)

# specify the number of pairs that will be considered for the cross-partition
  # analysis. This number should be set according to how many partitions you
  # want to cross. This probably should not be every single partition.
ncross = 10

# Note that the number of total "crossings" that happen are given
# by choose(ncross, 2). Use this knowledge, and the computation required here
# to inform the decision of ncross values.
choose(ncross, 2)

# specify where the RDS objects will get saved
save_dir = "data-raw/GLS_part_obs"

# option: Should coordinates be saved be loaded from separate files
use_coord_files = TRUE

# specify a directory where coordinates for each partition will be saved
coords_dir = "data-raw/GLS_part_coords"

# option: should results be overwritten if they already exist?
overwrite = TRUE

# number of cores to use for matrix math
ncores = 6

# create the GLS save directory if it doesn't exist
if(!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}
# create the coordinate directory if it doesn't exist (and option used)
if (!dir.exists(coords_dir) & use_coord_files) {
  dir.create(coords_dir, recursive = TRUE)
}

# ---- 1. Fit and save the GLS ----

# fit a GLS to each partition, by looping through the partition matrix columns
for (i in seq_len(ncol(pm))) {
  # convert i into a fixed-length string, with padded zeros
  i_str = str_pad(i, 3, side = "left", pad = "0")

  # should the xx and invchol be saved for this iteration?
    # this method makes it so only the partitions that will be crossed will
    # have the extra information saved during the fit process
  save_extra = (i <= ncross)

  # convert to a string that identifies the files that should be crossed
  xx_str = ifelse(save_extra, "_to-cross", "_dont-cross")

  # build the file name for this partition
  name_i = paste0("AK10000_GLS_results_partition", i_str, xx_str, ".rds")

  # build the the full path to this file
  file_i = file.path(save_dir, name_i)

  # skip to partition if this file already exists (unless overwrite=TRUE)
    # When overwrite = FALSE, this prevents the script from rerunning/saving
    # GLS objects for partitions that already exist in the folder to save time.
  if (file.exists(file_i) & !overwrite) {
    next
  }

  # get the column of partition matrix to use
  part_i = pm[, i]

  # subset the data by the partition (you may have a different process for this
    # step if the full dataset does not fit in memory)
  data_i <- ndvi_AK10000[part_i, ]

  # get/set parameter values (set to constants for ease here, but use or
    # calculate values according to your needs)
  nug = 0 # nugget
  rng = 0.01 # range of spatial autocorrelation

  # get the coordinates for this subset (partition)
  coords_i = data_i[, c("lng", "lat")]

  # calculate the varcov matrix
  V = covar_exp(distm_km(coords_i), range = rng)

  # fit the GLS - note that, because of save_extra, only the first ncross (10)
    # partitions will have xx and invchol saved
  GLS_i <- fitGLS(
    CLS_coef ~ lat + land, data = data_i, V = V, nugget = nug,
    save.xx = save_extra, save.invchol = save_extra, ncores = ncores
  )

  # Save coordinates to for the partitions that will be crossed
  if (save_extra){
    if (use_coord_files) {
      # save the coordinates to separate files
      coords_file_name = gsub("GLS_results", "coordinates", name_i)
      write.csv(
        coords_i, file = file.path(coords_dir, coords_file_name),
        row.names = FALSE
      )
    } else {
      # alternatively, save coordinates to the GLS objects directly
      GLS_i$coords <- coords_i
    }
  }

  # save the GLS to an RDS file
  saveRDS(GLS_i, file = file_i)
}

# ---- 2. Calculate and collect cross-partition statistics ----

# get a list of all the partitioned results that will be crossed
cross_files = list.files(save_dir, pattern = "to-cross", full.names = TRUE)

# get the coordinates from the files, if the option is used
if (use_coord_files){
  cross_coords_files = list.files(
    coords_dir, pattern = "to-cross", full.names = TRUE
  )
}

# get a list of the remaining partitions (that won't be crossed but still used)
uncrossed_files = list.files(
  save_dir, pattern = "dont-cross", full.names = TRUE
)

# (re)calculate the number of partitions to cross
ncross = min(ncross, length(cross_files))

# re(calculate) the number of total partitions
npart = length(uncrossed_files) + ncross

# get a table of all the possible crosses - this determines how crossing happens
cross.pairs = t(combn(seq_len(ncross), 2))

# Create a table that indicates the iterations on which a new GLS needs to be
  # loaded for a given pair of GLS. This will simply prevent the same GLS object
  # from getting loaded over and over again.
load_table = rbind(c(TRUE, TRUE), diff(cross.pairs) != 0)

# set up a progress bar to track progress
pb = txtProgressBar(min = 0, max = nrow(cross.pairs), style = 3)

# build empty vectors to store correlated sum of squares values
rSSRs <- rSSEs <- rep(NA, nrow(cross.pairs))

# build empty list to store correlated coefficient estimates
rcoefs <- list()

# read the first GLS object to facilitate calculation of some parameters
GLS_i <- readRDS(cross_files[1])
# number of coefficients:
n_coefs <- length(GLS_i$coefficients)
# number of null model coeficients:
n_null_coefs <- length(GLS_i$coefficients0)
# names of the coefficients:
coef_names <- names(GLS_i$coefficients)
# size of each partition:
partsize <- nrow(GLS_i$xx)
# degrees of freedom:
dfs <- remotePARTS:::calc_dfpart(partsize, n_coefs, n_null_coefs)

# create an empty table to hold coefficients from all partitions
coefs <- matrix(
  NA, nrow = npart, ncol = n_coefs,
  dimnames = list(NULL, coef_names)
)

# create an empty list to hold covariance of coefficients
covar_coefs <- vector("list", npart)

# create empty vectors to hold F-statistics
F_stats <- rep(NA, times = npart)

# loop through the cross pairs table, to perform the crossing
for (pair in seq_len(nrow(cross.pairs))) {
  # update the progress bar, to indicate what step we're on
  setTxtProgressBar(pb, value = pair)

  # get the number of the first partition, from our cross table
  i = cross.pairs[pair, 1]
  # get the file path to the first partition's GLS object
  path_i = cross_files[i]
  # get the number of the second partition, from our cross table
  j = cross.pairs[pair, 2]
  # get the file path to the second partition's GLS object
  path_j = cross_files[j]

  # load the first GLS from its file, unless it doesn't need to be loaded
  if (load_table[pair, 1]) {
    GLS_i = readRDS(path_i)
  }
  # load the second GLS from its file, unless it doesn't need to be loaded
  if (load_table[pair, 2]) {
    GLS_j = readRDS(path_j)
  }

  # fill the coefficient table for both partitions
  coefs[i, ] <- GLS_i$coefficients
  coefs[j, ] <- GLS_j$coefficients
  # do the same for the covariance among the coefficients
  covar_coefs[[i]] <- GLS_i$covar_coef
  covar_coefs[[j]] <- GLS_j$covar_coef
  # and F-statistics
  F_stats[i] <- GLS_i$Fstat
  F_stats[j] <- GLS_j$Fstat

  # calculate the pairwise distance between the two partitions
  if (use_coord_files) {
    # either based on coordinates saved in files
    coords_i = read.csv(cross_coords_files[i])
    coords_j = read.csv(cross_coords_files[j])
    dist_ij = distm_km(coords = coords_i, coords2 = coords_j)
  } else {
    # or based on coordinates saved in the GLS objects
    dist_ij = distm_km(coords = GLS_i$coords, coords2 = GLS_j$coords)
  }

  # calculate the joint varcovar matrix
  V_ij = covar_exp(dist_ij, range = rng)

  # perform the cross, using an unlisted remotePARTS function
    # (accessed with 3 colons ":::").
    # Not that we use the "xx" and "xx0" values of the GLS object,
    # which were conditionally saved in step 1 (according to save_extra)
    # with the "save.xx" argument in fitGLS. We also use the invcholV matrix
    # saved with "save.invchol" argument.
  pair_cross = remotePARTS:::crosspart_GLS(
    xxi = GLS_i$xx, xxj = GLS_j$xx,
    xxi0 = GLS_i$xx0, xxj0 = GLS_j$xx0,
    invChol_i = GLS_i$invcholV, invChol_j = GLS_j$invcholV,
    Vsub = V_ij,
    nug_i = GLS_i$nugget, nug_j = GLS_j$nugget,
    df1 = dfs[1], df2 = dfs[2],
    ncores = ncores
  )

  # rename the cross-correlation coefficients' dimensions
  rownames(pair_cross$rcoefij) <- colnames(pair_cross$rcoefij) <- coef_names

  # collect the cross-correlations, from this cross, into our vector
  rcoefs[[pair]] <- pair_cross$rcoefij

  # collect the regression sum of squares (NA will replace infinite values)
  rSSRs[pair] <- ifelse(
    is.na(pair_cross$rSSRij) | is.infinite(pair_cross$rSSRij),
    NA, pair_cross$rSSRij
  )

  # collect the sum of squares errors
  rSSEs[pair] <- pair_cross$rSSEij

}

# ---- 3. Collect GLS statistics from remaining non-crossed partitions ----

# loop through the remaining GLS objects for partitions not being crossed
for (i in seq_len(length(uncrossed_files))) {

  # identify this partition number (out of *all* partitions)
  part_i = ncross + i

  # load the GLS object for this partition
  GLS_i <- readRDS(uncrossed_files[i])

  # collect the regression coefficients for this partition, into the table
  coefs[part_i, ] <- GLS_i$coefficients

  # collect the covariance among coefficients
  covar_coefs[[part_i]] <- GLS_i$covar_coef

  # collect the F statistics
  F_stats[part_i] <- GLS_i$Fstat
}

# ---- process some of the stats ----

# average the coefficients across all partitions
mean_coefs = colMeans(coefs, na.rm = TRUE)

# array of the coefficient cross-correlations among crossed partition pairs
  # this step simply converts the list into an array in the proper arrangement
  # for the next step
rcoefs_array <- simplify2array(rcoefs) |> aperm(c(3, 1, 2))

# average the coefficient cross-coefficients across all partitions
rcoef_means = apply(rcoefs_array, c(2,3), function(x)mean(x, na.rm = TRUE))

# array of the covariance among coefficients for each partition
  # again, this is converted to an array for computational purposes
covar_coef_array <- simplify2array(covar_coefs) |> aperm(c(1, 2, 3))

# average the F statistics across all partitions
F_mean <- mean(F_stats, na.rm = TRUE)

# average the correlated regression sum of squares (SSR) values
mean_rSSR = mean(rSSRs, na.rm = TRUE)

# ---- 4.1 T test ----

# run the correlated t-test of the coefficients
t_test_results <- remotePARTS:::part_ttest(
  mean_coefs, covar_coef_array, rcoef_means, dfs[2], npart
)

# print the t-test results
  # pval.t < 0.05 indicates that the coefficient significantly affects the
  # response
t_test_results

# ---- 4.2 Chi-squared test ----

# run the chi-squared test of overall model fit
chisqr_results <- remotePARTS:::part_chisqr(F_mean, mean_rSSR, dfs[1], npart)

# print the chi-square test results
  # pval.chisqr < 0.05 indicates this model explains the variation in the
  # data well.
chisqr_results
