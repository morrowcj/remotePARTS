## remotePARTS objects (S3) and methods

#' generic method to extract fitted model from a complex object
#'
#' @param object complex S3 object
#' @param ... additional arguments passed to methods
#'
#' @seealso [fitCLS()] for example usage
#'
#' @export
get_fm = function(object, ...){ ## get the model object from complex objects
  UseMethod("get_fm")
}

# CLS ----

#' Build AR data frame
#'
#' @param x length p time series vector
#' @param t length p vector temporal vector
#' @param Z optional vector or design matrix of covariates. if only one
#' covariate, \code{length(Z)} = p. If more than one, \code{nrow(Z)} = p
#'
#' @return an AR data frame with columns corresponding to the value of x
#' at time t (\code{$x.tj}), the value of x at time t-1 (\code{$x.ti}),
#' the temporal variable at time t (\code{$t.j}), and the value of Z
#' at time t {\code{$Z.*}}.
#'
#'
#' @export
#'
#' @examples
#' AR_df(x = rnorm(20), t = 1:20)
#' AR_df(x = rnorm(20), t = 1:20, Z = cbind(rnorm(20), rnorm(20)))
AR_df <- function(x, t, Z = NULL){
  # variables
  t_n = length(t)

  # covar handling
  if (is.null(Z)) {
    Zs = NA
  } else {
    Zs = as.matrix(Z)
    stopifnot(dim(Z)[1] == t_n)
    Zs = Zs[2:t_n, ]
  }

  # data frame for the model
  AR.df <- data.frame(
    x.tj = x[2:t_n], #X_{t}
    x.ti = x[1:(t_n - 1)], #X_{t-1}
    t.j = t[2:t_n] #T_{t}
  )

  # fit the model
  if (!is.null(Z)) AR.df$z = Zs

  return(AR.df)
}

#' get linear model object from remote CLS
#'
#' @param obj remoteCLS object
#'
#' @return lm class object from stats::lm
#' @export
get_fm.remoteCLS <- function(obj){
  stopifnot("remoteCLS" %in% class(obj))
  if ("pixel" %in% class(obj)) {
    return(obj$fm)
  } else if ("map" %in% class(obj)) {

  } else
    NULL
}

#' printmethod for remoteCLS object
#'
#' @param obj remoteCLS object
#' @export
print.remoteCLS <- function(obj){
  stopifnot("remoteCLS" %in% class(obj))

  ## Pixel version
  if ("pixel" %in% class(obj)) {
    # call
    cat("Call: ",deparse(obj$call), "\n",
        "CLS model: ", deparse(obj$fm$call$formula), "\n\n",
        "Coefficeints:", "\n", sep = "")
    # then lm
    print(obj$fm$coefficients)

  } else if ("map" %in% class(obj)) {
    ## Map version
    cat("Call:",deparse(obj$call), "\n")
    cat("\n","Time Coefficients:","\n", sep = "")
    print(format(as.data.frame(obj$time.coef)), digits = 3, nsmall = 2, scientific = 1)
  }
}


#' Extract coefficeints from remoteCLS object
#'
#' @param obj remoteCLS object
#'
#' @return coefficient data frame
#' @export
coef.remoteCLS <- function(obj){
  stopifnot("remoteCLS" %in% class(obj))

  if ("pixel" %in% class(obj)) {
   return(as.data.frame(obj$fm$coefficeints))
  } else if ("map" %in% class(obj)) {
    return(as.data.frame(obj$time.coef))
  }
}

#' summary of remoteCLS object
#'
#' @param obj remoteCLS object
#'
#' @export
summary.remoteCLS <- function(obj){
  stopifnot("remoteCLS" %in% class(obj))

  ## Pixel version
  if ("pixel" %in% class(obj)) {
    fm <- get_fm(obj)
    smry <- summary(fm)
    MSE = sigma(fm)^2

    # coefficient table
    cat("coefficients:\n")
    print(format(as.data.frame(smry$coefficients),
           digits = 2, nsmall = 2, scientific = 1))
    # MSE
    cat("\nMSE =", format(MSE, digits = 3))


  } else if ("map" %in% class(obj)) {
    ## Map version
    cat("Summary\n")
    cat("\neffect of time:\n")
    print(summary(as.data.frame(obj$time.coef)))
    if(!is.null(obj$xi.coef)){
      cat("\neffect of previous x:\n")
      print(summary(as.data.frame(obj$xi.coef)))
    }
    if(!is.null(obj$int.coef)){
      cat("\nintercept:\n")
      print(summary(as.data.frame(obj$int.coef)))
    }
    if(!is.null(obj$MSE)){
      cat("\nMSE:\n")
      print(summary(as.data.frame(obj$MSE)))
    }
  }
}

#' resisduals of remoteCLS object
#'
#' @param obj remoteCLS object
#'
#' @return model residuals
#' @export
residuals.remoteCLS <- function(obj){
  stopifnot("remoteCLS" %in% class(obj))
  ## Pixel version
  if ("pixel" %in% class(obj)) {
    fm = get_fm(obj)
    return(fm$residuals)
  } else if ("map" %in% class(obj)) {
    obj$residuals
  }
}


# AR ----
#' print method for remoteAR object
#'
#' @param obj remoteAR object
#'
#' @return
#' @export
print.remoteAR <- function(obj){
  stopifnot("remoteAR" %in% class(obj))

  ## Pixel version
  if ("pixel" %in% class(obj)) {
    # call
    cat("Call:",deparse(obj$call), "\n")
    # then lm
    cat("\n","Coefficients:","\n", sep = "")
    print(format(obj$coef, nsmall = 2, digits = 3, scientific = 1))
    cat("\n",
        "AR parameter estimate: ",
        format(obj$b, digits = 3, nsmall = 2, scientific = 1),
        "\n",
        "MSE: ", obj$MSE, "\n",
        "log-likelihood: ", obj$logLik, sep = "")


  } else if ("map" %in% class(obj)) {
    ## Map version
    cat("Call: ", deparse(obj$call), "\n\n",
        "Time Coefficeints: ", "\n", sep = "")
    print(format(as.data.frame(obj$time.coef)), digits = 3, nsmall = 2, scientific = 1)
  }
}

#' Extract coefficients from remoteAR object
#'
#' @param obj remoteAR object
#'
#' @return coefficeint data frame
#' @export
coef.remoteAR <- function(obj){
  stopifnot("remoteAR" %in% class(obj))

  ## Pixel version
  if ("pixel" %in% class(obj)) {
    return(obj$coef)
  } else if ("map" %in% class(obj)) {
    return(as.data.frame(obj$time.coef))
  }
}

#' Extract residuals from remoteAR object
#'
#' @param obj remoteAR object
#'
#' @return residuals
#' @export
residuals.remoteAR <- function(obj){
  stopifnot("remoteAR" %in% class(obj))

  ## Pixel version
  if ("pixel" %in% class(obj)) {
    return(obj$resids)
  } else if ("map" %in% class(obj)) {
    return(obj$resids)
  }
}

# GLS ----

#' remoteGLS constructor (S3)
#'
#' @return an empty S3 object of class "remoteGLS"
#' @export
#'
#' @examples
#' tmp <- remoteGLS() #empty remoteGLS object
remoteGLS <- function(form){
  # create empty list
  elements = c("model.info", # model info
               "betahat", "SSE", "MSE", "SE",
               "dft", "logDetV",
               "tstat", "pval.t", "logLik",
               "betahat0", "SSE0", "MSE0", "SE0",
               "MSR", "df0", "logLik0", "df.F",
               "Fstat","pval.F",
               "xx", "xx0", "invcholV", "nugget")
  GLS.obj = vector("list", length(elements))
  names(GLS.obj) = elements

  GLS.obj$model.info = list("call" = match.call(),
                            "formula" = NULL,
                            "response" = NULL,
                            "predictors" = NULL,
                            "coef.names" = NULL)
  if(!missing(form)){
    GLS.obj$model.info$formula = formula(form)
    GLS.obj$model.info$response = all.vars(formula(form))[1]
    GLS.obj$model.info$predictors = all.vars(formula(form))[-1]
  }

  # GLS.obj$call = match.call()
  # print(match.call(expand.dots = FALSE))

  class(GLS.obj) <- c("remoteGLS")

  return(GLS.obj) # return the empty list
}

#' print method for remoteGLS
#'
#' @param obj remoteGLS object
#'
#' @return formatted output for remoteGLS object
#' @export
print.remoteGLS <- function(obj, print.call = TRUE){
  stopifnot("remoteGLS" %in% class(obj))
# print.remoteGLS <- function(obj){
  if(is.null(unlist(obj[-1]))){
    cat("empty remoteGLS object\ncall: ")
    if(!is.null(unlist(obj[1]))){
      print(obj$model.info)
    } else {print(obj$model.info$call)}
  } else {

    ## Coefficient table
    coefs = data.frame("Est" = obj$betahat, "t stat" = obj$tstat, "pval t" = obj$pval.t)

    ## Model stats
    mod.stats = data.frame("df" = obj$df.F,
                           "SSE" = c(obj$SSE, obj$SSE0),
                           "MSE" = c(obj$MSE, obj$MSE0),
                           "logLik" = c(obj$logLik, obj$logLik0),
                           "F stat" = c(obj$Fstat, NA),
                           "pval F" = c(obj$pval.F, NA))
    rownames(mod.stats) = c("mod_A", "mod_0")

    ## return (needs cleaning)
    if (print.call){
      cat("call: ");print(obj$model.info$call)
    }
    cat("response:", obj$model.info$response,"\n\n")
    cat("t tests:\n")
    print(format(coefs, digits = 2, nsmall = 2))
    cat("\nF test:\n")
    print(format(mod.stats, digits = 2, nsmall = 2, scientific = -2))

  }
}

## summary function
# summary.remoteGLS <- function(obj){


# Other ----

## Test a small C++ function to alter list elements ----
if (FALSE) {
  tmp <- remoteGLS()
  ## Create a test funciton that replaces
  sourceCpp(
  code =
    '#include <Rcpp.h>

    using namespace Rcpp;

    // [[Rcpp::export(.test_remoteGLS)]]
    void test_remoteGLS(List L){

      // by index
      L[0] = 20.1;

      // by name
      L["dft"] = 5.;

      // character
      L["SSE"] = CharacterVector::create("boom!");

      // vector
      L["SE"] = NumericVector::create(1,3,5,7,9);

      // using a variable
      NumericVector v = {1,2,5.1};
      L["MSE"] = v;
    }'
  )

  print(tmp)
  .test_remoteGLS(tmp) # changes the values within tmp
  print(tmp)
}



## ----

## More tests
# df = as.data.frame(matrix(rnorm(50), ncol = 5))
# df$fact <- factor(ceiling(runif(10, 0, 3)))
# names(df) = LETTERS[1:5]

# mod.frame = test.func(A ~ B + D, data = df, method = "model.frame")
# mod.mat = model.matrix(A ~ B + D, mod.frame)
# tmp <- test.func(A ~ B + D, data = df, method = "qr")
#
# model.matrix(tmp$mt, mod.frame)



# build.mod(A ~ B + D + fact, data = df)

