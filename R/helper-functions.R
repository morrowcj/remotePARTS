## invert_chol() helpers ----

#' check if matrix is positive definite
#'
#' @details check if a matrix is 1) square, 2) symmetric, and 3) positive
#' definite
#'
#' @param M numeric matrix
#'
#' @export
check_posdef <- function(M){

  res <- c(sqr = FALSE, sym = FALSE, posdef = FALSE)

  res["sqr"] = nrow(M) == ncol(M)
  # if(!res["sqr"]){
  #   stop("M not square")
  # } # not square

  if(res["sqr"]){

    res["sym"] = isSymmetric(M)
    # if(!res["sym"]){
    #   stop("M not symmetric")
    # } # not symetric

    if(res["sym"]){

      res["posdef"] = all(eigen(M, only.values = TRUE)$values > 1e-8)
      # if(!res["posdef"]){
      #   stop("M not positive definite")
      # }
    }
  }
  return(res)
}

## GLS Helpers ----

## GLS Partition Helpers ----

#' function to calculate partition size or number of partitions
#'
#' @param npix number of pixels in full dataset
#' @param pixels vector of pixel indexes to sample from
#' @export
#'
#' @examples
#' # setup data
#' dat.M <- matrix(rnorm(3000*20), ncol = 20)
#' # 4 partitions (exhaustive)
#' sample_partitions(npix = nrow(dat.M), npart = 4)
#' # partitions with 500 pixels each (exhaustive)
#' sample_partitions(npix = nrow(dat.M), partsize = 500)
#' # 4 partitions each with 500 pixels (non-exhaustive)
#' sample_partitions(npix = nrow(dat.M), npart = 4, partsize = 500)
#'
#' # index of pixels to subset
#' sub.indx <- 1:1000
#' # 4 partitions (exhaustive) using only the specified pixels
#' sample_partitions(npix = nrow(dat.M), npart = 4, pixels = sub.indx)
sample_partitions <- function(npix, npart = 10, partsize = NA,
                              pixels = NA, verbose = TRUE){

  if(!is.na(pixels) && length(pixels) > 1){
    npix = length(pixels)
    from = pixels
  } else {
    from = 1:npix
  }

  ## check which npart of partsize was given
  no.partsize <- (missing(partsize) || is.na(partsize) | is.null(partsize))
  no.npart <- (missing(npart) || is.na(npart) | is.null(npart))

  ## caclulate partition size
  if(no.partsize){
    if(verbose){print("calculating partsize")}
    partsize = (npix - (npix%%npart)) / npart
  }

  ## OR calculate number of partitions
  if(no.npart){
    if(verbose){print("calculating npart")}
    npart = floor(npix/partsize)
  }

  if(npart * partsize > npix){
    stop("npart * partsize may not be greater than npix")
  }

  remainder = npix %% partsize
  samp <- sample(from, size = npix - remainder, replace = FALSE)

  part.mat <- matrix(samp, ncol = npart, nrow = partsize)
  colnames(part.mat) <- paste("part",1:npart, sep = ".")

  return(part.mat)
}

#' calculate degrees of freedom for partitioned GLS
#'
#' @param part.size number of pixels in each partition
#' @param p number of predictors in alternate model
#' @param p0 number of parameters in null model
#'
#' @export
#'
#' @examples
#' calc_df(partsize = 2000, p = 4, p0 = 1)
calc_dfpart <- function(partsize, p, p0){
  stopifnot(length(partsize) == 1)
  df2 = partsize - (p - 1)
  df0 = partsize - (p0 - 1)
  df1 = df0 - df2
  return(c("df1" = df1, "df2" = df2))
}
