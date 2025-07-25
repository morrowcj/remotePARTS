## TODO This file should contain funcitions that facilitate the partitioned GLS
# from individual GLS objects (calculated on individual partitions and saved
# separately). See R/scripts/combine_GLS_partitions.R for a general process.

# ---- functions to aid cross-GLS ----
#'
#'
#' #' get the partition statistics from a GLS object
#' #'
#' #' @param GLS
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' get_part_stats <- function(GLS = NULL){
#'   # initialize output list
#'   out = list(
#'     coefficients = numeric(), SEs = numeric(), covar_coefs = numeric(),
#'     tstats = numeric(), pvals_t = numeric(), nuggets = numeric(),
#'     covar.pars = numeric(),
#'     modstats = data.frame(
#'       LLs = numeric(), SSEs = numeric(), MSEs = numeric(), MSRs = numeric(),
#'       Fstats = numeric(), pvals_F = numeric()
#'     )
#'   )
#'   class(out) <- unique(c("part_stats", class(out)))
#'
#'   # return an empty object if no GLS given
#'   if (is.null(GLS)) {
#'     return(out)
#'   }
#'
#'   # out$coefficients = GLS$coefficients
#'   # out$SEs = GLS$SE
#'   # out$covar_coefs = GLS$covar_coef
#'   # out$tstats = GLS$tstat
#'   # out$pvals_t = GLS$pval_t
#'   # out$nuggets = GLS$nugget
#'   # out$covar.pars = GLS$covar.pars
#'   # out$modstats = data.frame(LL = GLS$logLik)
#'
#'   # TODO: CODE FOR ACTUAL GLS obj ....
#' }
#'
#' get_cross_stats <- function(GLS1, GLS2){
#'
#' }
#'
#' #' Check that the part_stats object has all the correct parts
#' #'
#' #' @param GLS the part_stats object to check
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' check_part_stats <- function(GLS){
#'
#'   stopifnot(is.list(GLS))
#'
#'   element_names = c(
#'     "coefficients", "SEs", "covar_coefs", "tstats", "pvals_t",
#'     "nuggets", "covar.pars", "modstats"
#'   )
#'
#'   stopifnot(all(names(GLS) %in% element_names))
#'
#'   modstats_names = c("LLs" , "SSEs" , "MSEs" , "MSRs" , "Fstats" , "pvals_F")
#'
#'   stopifnot(all(names(GLS[["modstats"]]) %in% modstats_names))
#'
#'   if (!"part_stats" %in% class(GLS)) {
#'     warning(paste(substitute(GLS), "does not have the 'part_stats' class"))
#'   }
#' }
#'
#' #' append the partitions stats from one part_stats object to another
#' #'
#' #' @param GLS1 the first part_stats object
#' #' @param GLS2 the second part_stats object
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' append_part_stats <- function(GLS1, GLS2) {
#'   # check that both GLS1 and GLS2 are GLS objects
#'   # TODO: improve this
#'   check_part_stats(GLS1)
#'   check_part_stats(GLS2)
#'
#'   out <- GLS1
#'
#'   append_cols <- c(
#'     "coefficients", "SEs", "covar_coefs", "tstats", "pvals_t",
#'     "nuggets", "covar.pars"
#'   )
#'
#'   for (elem in append_cols) {
#'     out[[elem]] <- append(GLS1[[elem]], GLS2[[elem]])
#'   }
#'
#'   out[["modstats"]] <- rbind(GLS1[["modstats"]], GLS2[["modstats"]])
#'
#'   return(out)
#'
#'   # TODO generalize to more than 2 objects? ...
#'   # # convert the objects to a list
#'   # objects <- list(...)
#'   #
#'   # # how many objects are there?
#'   # len = length(objects)
#'   #
#'   # if (len < 2) {stop("at least 2 GLS objects needed")}
#'
#' }
#'
