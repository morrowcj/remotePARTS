
#' @title Fit a PARTS GLS model.
#'
#' @details conduct generalized least-squares regression of
#' spatiotemporal trends
#'
#' @param formula a model formula
#' @param data an optional data frame environment in which to search for
#' variables given by \code{formula}
#' @param V a covariance matrix, which must be positive definitive. This argument
#' is optional if \code{coords}, \code{distm_FUN}, \code{covar_FUN}, and
#' \code{covar.pars} are given instead.
#' @param nugget an optional numeric nugget, must be positive
#' @param formula0 an optional formula for the null model to be compared with
#' \code{formula} by an F-test
#' @param save.xx logical: should information needed for cross-partition
#' comparisons be returned?
#' @param save.invchol logical: should the inverse of the Cholesky matrix be
#' returned?
#' @param logLik.only logical: should calculations stop after calculating parital
#' log-likelihood?
#' @param no.F logical: should F-test calculations be made?
#' @param coords optional coordinate matrix for calculating \code{V} internally
#' @param distm_FUN optional function for calculating a distance matrix from
#' \code{coords}, when calculating \code{V} internally
#' @param covar_FUN optional distance-based covariance function for calculating
#' \code{V} internally
#' @param covar.pars an optional named list of parameters passed to \code{covar_FUN}
#' when calculating \code{V} internally
#' @param invCholV optional pre-calculated inverse cholesky matrix to use in place
#' of \code{V}
#' @param ... additional arguments passed to \code{optimize_nugget}, which are
#' only used if if \code{nugget = NA}
#'
#' @details \code{fitGLS} fits a GLS model, using terms specified in \code{formula}.
#' In the PARTS method, generally the left side of \code{formula} should be
#' pixel-level trend estimates and the right side should be spatial predictors.
#' The errors of the GLS are correlated according to covariance matrix \code{V}.
#'
#' If \code{nugget = NA}, an ML nugget is estimated from the data using the
#' \code{optimize_nugget()} function. Arguments additional arguments (\code{...})
#' are passed to \code{optimize_nugget} in this case. \code{V} must be provided
#' for nugget optimization.
#'
#' If \code{formula0} is not specified, the default is to fit an intercept-only
#' null model.
#'
#' \code{save.xx} is included to allow for manually conducting a partitioned
#' GLS analyses. Because most users will not need this feature, opting instead
#' to use \code{fitGLS_parition()}, \code{save.xx = FALSE} by default.
#'
#' Similarly, \code{save.invchol} is included to allow for recycling of the
#' inverse cholesky matrix. Often, inverting the large cholesky matrix
#' (i.e., \code{invert_chol(V)}) is the slowest part of GLS. This argument exists
#' to allow users to recycle this process, though no \code{remotePARTS} function
#' currently exists that can use \code{invert_chol(V)} to fit the GLS.
#'
#' \code{logLik.only = TRUE} will return only the partial log-likelihood, which can
#' minimized to obtain the maximum likelihood for a given set of data.
#'
#' If \code{no.F = TRUE}, then the model given by \code{formula} is not compared
#' to the model given by \code{formula0}.
#'
#' If \code{V} is not provided, it can be fit internally by specifying all of
#' \code{coords}, \code{distm_FUN}, \code{covar_FUN}, and \code{covar.pars}.
#' The function given by \code{distm_FUN} will calculate a distance matrix from
#' \code{coords}, which is then transformed into a distance-based covariance
#' matrix with \code{covar_FUN} and parameters given by \code{covar.pars}.
#'
#' This function uses C++ code that uses the Eigen matrix library (RcppEigen
#' package) to fit models as efficiently as possible. As such, all available
#' CPU cores are used for matrix calculations on systems with OpenMP
#' support.
#'
#' @return \code{fitGLS} returns a list object of class "remoteGLS", if
#' \code{logLik.only = FALSE}. The list contains at least the following elements:
#'
#' \describe{
#'     \item{coefficients}{coefficient estimates for predictor variables}
#'     \item{SSE}{sum of squares error}
#'     \item{MSE}{mean squared error}
#'     \item{SE}{standard errors}
#'     \item{df_t}{degrees of freedom for the t-test}
#'     \item{logDetV}{log-determinant of V}
#'     \item{tstat}{t-test statistic}
#'     \item{pval_t}{p-value of the t-statistic}
#'     \item{logLik}{the Log-likelihood of the model}
#'     \item{nugget}{the nugget used in fitting}
#' }
#'
#' If \code{no.F = FALSE}, the following elements, corresponding to the null
#' model and F-test are also calculated:
#'
#' \describe{
#'     \item{coefficients0}{coefficient estimates for the null model}
#'     \item{SSE0}{sum of squares error for the null model}
#'     \item{MSE0}{mean squared error for the null model}
#'     \item{SE0}{the standard errors for null coefficients}
#'     \item{MSR}{the regression mean square}
#'     \item{df0}{the null model F-test degrees of freedom}
#'     \item{LL0}{the log-likelihood of the null model}
#'     \item{df_F}{the F-test degrees of freedom, for the main model}
#'     \item{Fstat}{the F-statistic}
#'     \item{pval_F}{the F-test p-value}
#'     \item{formula}{the alternate formula used}
#'     \item{formula0}{the null formula used}
#' }
#'
#' An attribute called also set to \code{"no.F"} is set to the value of
#' argument \code{no.F}, which signals to generic methods how to handle the output.
#'
#' If \code{logLik.only = TRUE}, a single numeric output containing the partial
#' log-likelihood is returned. This value is primarily for ML estimation. In
#' parameter space for a given dataset, the minimum partial likelihood corresponds
#' to the maximum true likelihood.
#'
#' @examples
#' ## read data
#' data(ndvi_AK3000)
#' df = ndvi_AK3000[seq_len(1000), ] # first 1000 rows
#'
#' ## fit covariance matrix
#' V = covar_exp(distm_scaled(cbind(df$lng, df$lat)), range = .01)
#'
#' ## run GLS
#' (GLS = fitGLS(CLS_coef ~ 0 + land, data = df, V = V))
#'
#' ## with F-test calculations to compare with the NULL model
#' (GLS.F = fitGLS(CLS_coef ~ 0 + land, data = df, V = V, no.F = FALSE))
#'
#' ## find ML nugget
#' fitGLS(CLS_coef ~ 0 + land, data = df, V = V, no.F = FALSE, nugget = NA)
#'
#' ## calculate V internally
#' coords = cbind(df$lng, df$lat)
#' fitGLS(CLS_coef ~ 0 + land, data = df, logLik.only = FALSE, coords = coords,
#'        distm_FUN = "distm_scaled", covar_FUN = "covar_exp", covar.pars = list(range = .01))
#'
#' ## use inverse cholesky
#' fitGLS(CLS_coef ~ 0 + land, data = df, invCholV = invert_chol(V))
#'
#' ## save inverse cholesky matrix
#' invchol = fitGLS(CLS_coef ~ 0 + land, data = df, V = V, save.invchol = TRUE)$invcholV
#'
#' ## re-use inverse cholesky instead of V
#' fitGLS(CLS_coef ~ 0 + land, data = df, invCholV = invchol)
#'
#' ## Log-likelihood (fast)
#' fitGLS(CLS_coef ~ 0 + land, data = df, V = V, logLik.only = TRUE)
#'
#' @export
fitGLS <- function(formula, data, V, nugget = 0, formula0 = NULL, save.xx = FALSE,
                   save.invchol = FALSE, logLik.only = FALSE, no.F = TRUE,
                   coords, distm_FUN ,covar_FUN, covar.pars, invCholV,
                   ...){

  # Parse formula arguments to make model matrix
  call <- match.call() # function call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame) # rename the function call
  mf$drop.unused.levels = TRUE
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- as.double(stats::model.response(mf, "numeric"))
  if (is.matrix(y)) {stop("response is a matrix: must be a vector")}
  X <- stats::model.matrix(mt, mf, contrasts = NULL)
  rm(mf) # delete the large model frame from memory

  # Use invCholV if provided
  if (missing(invCholV) || is.null(invCholV)){
    # print("no inverse cholesky provided")
    invCholV = diag(1) # default matrix to pass to C++
    use_invCholV = FALSE
  } else {
    if (is.na(nugget) & missing(V)){
      stop("nugget cannot be optimized with pre-calculated invCholV, use V instead")
    }
    # print("inverse cholesky WAS provided")
    invCholV = invCholV
    use_invCholV = TRUE
  }

  # calculate V, if needed
  if (missing(V) & !use_invCholV) {
    if (any(missing(coords), missing(distm_FUN), missing(covar_FUN))) {
      stop("If V or invCholV are not provided, then coords, distm_FUN, covar_FUN, and covar.pars are needed to calculate V")
    }
    if (nrow(coords) != nrow(X)){
      stop("Rows of coords do not match rows in data")
    }
    dist.f = match.fun(distm_FUN)
    V <- dist.f(coords) # distance
    pars = append(list(V), covar.pars) # pars list
    covar.f = match.fun(covar_FUN)
    V <- do.call(covar.f, pars) # covariance matrix
    nv = nrow(V)
  } else if (missing(V) & use_invCholV){
    V = diag(1) # default matrix to pass to C++
    nv <- nrow(invCholV)
  } else {nv = nrow(V)}

  # Build null model
  if (is.null(formula0)){ # formula
    formula0 = update(as.formula(formula), . ~ 1)
  } else {
    formula0 = as.formula(formula0)
  }
  X0 <- if (missing(data) || is.null(data)) { # conditionally assign X0
    model.matrix(formula0)
  } else {
    model.matrix(formula0, data)
  }

  # Checks
  ## check positive definitive
  if (!all(check_posdef(V))) {
    if (logLik.only){
      warning("V is not positive definitive")
      return(NA) # return NA for logLik function
    }
    stop("V is not positive definitive")
  }
  # if (!all(check_posdef(invCholV))){
  #   stop("invCholV is not positive definitive")
  # }
  ## check for correct dimensions
  if (!all.equal(length(y), nrow(X), nrow(X0), nv)) {
    stop("Input dimension mismatch")
  }
  ## check that all variables are numeric
  if (!all(is.double(y), is.double(X), is.double(V), is.double(invCholV), is.double(X0))) {
    stop("All inputs must be numeric (double precision)")
  }

  # Handle missing nugget
  if (is.na(nugget)) {
    nugget = optimize_nugget(X = X, y = y, V = V, ...)
  }

  # Run GLS
  GLS <- .Call(`_remotePARTS_fitGLS_cpp`, X, V, y, X0,
               nugget, save.xx, save.invchol, logLik.only, no.F,
               optimize_nugget = FALSE, nug_l = 0, nug_u = 1, nug_tol = 1e-5,
               invCholV = invCholV,  use_invCholV = use_invCholV)

  # return only log-likelihood, if prompted
  if(logLik.only){
    return(unlist(GLS$logLik))
  }

  # pvalues
  GLS$pval_t <- 2 * pt(abs(GLS$tstat), df = GLS$df_t, lower.tail = F)
  if(!no.F){GLS$pval_F <- pf(GLS$Fstat, df1 = GLS$df_F[1], df2 = GLS$df_F[2], lower.tail = F)}

  # Update list elements
  GLS <- append(list(call = call), GLS)
  names(GLS$coefficients) = names(GLS$SE) = names(GLS$tstat) = names(GLS$pval_t) = colnames(X)
  # GLS$predictors = colnames(X)
  # GLS$nugget = nugget
  GLS$formula = deparse(as.formula(formula))
  GLS$formula0 = deparse(as.formula(formula0))

  ## class and attributes
  class(GLS) <- append("remoteGLS", class(GLS))
  attr(GLS, "no.F") <- no.F

  ## Return ----
  return(GLS)
}
