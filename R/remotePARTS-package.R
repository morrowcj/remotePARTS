
#' @keywords internal
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @importFrom stats aggregate coef cor formula lm logLik model.matrix
#' model.response nls nls.control optimize p.adjust pf pt residuals rnorm sigma
#' optim
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#' @importFrom stats t.test as.formula update anova terms terms.formula quantile
#' @importFrom utils txtProgressBar setTxtProgressBar combn
#' @docType package
## usethis namespace: start
#' @useDynLib remotePARTS, .registration = TRUE
## usethis namespace: end
NULL

# # get rid of warning for "i" iterator in foreach
# utils::globalVariables("i") # alternative to i <- NULL in fitGLS_partitionMC
