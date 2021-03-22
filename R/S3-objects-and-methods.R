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
#' @param object remoteCLS object
#' @param ... additional arguments
#'
#' @return lm class object from stats::lm
#' @export
get_fm.remoteCLS <- function(object, ...){
  stopifnot("remoteCLS" %in% class(object))
  if ("pixel" %in% class(object)) {
    return(object$fm)
  } else if ("map" %in% class(object)) {
    stop("no applicable method for remoteCLS.map")
  } else
    NULL
}

#' printmethod for remoteCLS object
#'
#' @param x remoteCLS object
#' @param ... additional arguments
#' @export
print.remoteCLS <- function(x, ...){
  stopifnot("remoteCLS" %in% class(x))

  ## Pixel version
  if ("pixel" %in% class(x)) {
    # call
    cat("Call: ",deparse(x$call), "\n",
        "CLS model: ", deparse(x$fm$call$formula), "\n\n",
        "Coefficeints:", "\n", sep = "")
    # then lm
    print(x$fm$coefficients)

  } else if ("map" %in% class(x)) {
    ## Map version
    cat("Call:",deparse(x$call), "\n")
    cat("\n","Time Coefficients:","\n", sep = "")
    print(format(as.data.frame(x$time.coef)), digits = 3, nsmall = 2, scientific = 1)
  }
}


#' Extract coefficeints from remoteCLS object
#'
#' @param object remoteCLS object
#' @param var desired coefficients from remoteCLS.map object. One of "time"
#' (default), "xi", or "intercept".
#' @param ... additional arguments
#'
#' @return coefficient data frame
#' @export
coef.remoteCLS <- function(object, var, ...){
  stopifnot("remoteCLS" %in% class(object))

  if ("pixel" %in% class(object)) {
   return(as.data.frame(object$fm$coefficeints))
  } else if ("map" %in% class(object)) {
    if(missing(var) || var == "time"){
      return(as.data.frame(object$time.coef))
    } else if (var == "xi") {
      return(as.data.frame(object$xi.coef))
    } else if (var == "intercept") {
      return(as.data.frame(object$int.coef))
    }
  }
}

#' summary of remoteCLS object
#'
#' @param object remoteCLS object
#' @param ... additional arguments
#'
#' @export
summary.remoteCLS <- function(object, ...){
  stopifnot("remoteCLS" %in% class(object))

  ## Pixel version
  if ("pixel" %in% class(object)) {
    fm <- get_fm(object)
    smry <- summary(fm)
    MSE = sigma(fm)^2

    # coefficient table
    cat("coefficients:\n")
    print(format(as.data.frame(smry$coefficients),
           digits = 2, nsmall = 2, scientific = 1))
    # MSE
    cat("\nMSE =", format(MSE, digits = 3))


  } else if ("map" %in% class(object)) {
    ## Map version
    cat("Effect of time (time.coef):\n")
    print(summary(as.data.frame(object$time.coef)))
    if(!is.null(object$xi.coef)){
      cat("\nEffect of previous x (xi.coef):\n")
      print(summary(as.data.frame(object$xi.coef)))
    }
    if(!is.null(object$int.coef)){
      cat("\nIntercept (int.coef):\n")
      print(summary(as.data.frame(object$int.coef)))
    }
    if(!is.null(object$MSE)){
      cat("\nModel MSE:\n")
      print(summary(object$MSE))
    }
  }
}

#' resisduals of remoteCLS object
#'
#' @param object remoteCLS object
#' @param ... additional arguments
#'
#' @return model residuals
#' @export
residuals.remoteCLS <- function(object, ...){
  stopifnot("remoteCLS" %in% class(object))
  ## Pixel version
  if ("pixel" %in% class(object)) {
    fm = get_fm(object)
    return(fm$residuals)
  } else if ("map" %in% class(object)) {
    object$residuals
  }
}


# AR ----
#' print method for remoteAR object
#'
#' @param x remoteAR object
#' @param ... additional arguments
#'
#' @return print formatted remoteAR
#' @export
print.remoteAR <- function(x, ...){
  stopifnot("remoteAR" %in% class(x))

  ## Pixel version
  if ("pixel" %in% class(x)) {
    # call
    cat("Call:",deparse(x$call), "\n")
    # then lm
    cat("\n","Coefficients:","\n", sep = "")
    print(format(x$coef, nsmall = 2, digits = 3, scientific = 1))
    cat("\n",
        "AR parameter estimate: ",
        format(x$b, digits = 3, nsmall = 2, scientific = 1),
        "\n",
        "MSE: ", x$MSE, "\n",
        "log-likelihood: ", x$logLik, sep = "")


  } else if ("map" %in% class(x)) {
    ## Map version
    cat("Call: ", deparse(x$call), "\n\n",
        "Time Coefficeints: ", "\n", sep = "")
    print(format(as.data.frame(x$time.coef)), digits = 3, nsmall = 2, scientific = 1)
  }
}

#' Extract coefficients from remoteAR object
#'
#' @param object remoteAR object
#' @param ... additional arguments
#'
#' @return coefficeint data frame
#' @export
coef.remoteAR <- function(object, ...){
  stopifnot("remoteAR" %in% class(object))

  ## Pixel version
  if ("pixel" %in% class(object)) {
    return(object$coef)
  } else if ("map" %in% class(object)) {
    return(as.data.frame(object$time.coef))
  }
}

#' Extract residuals from remoteAR object
#'
#' @param object remoteAR object
#' @param ... additional arguments
#'
#' @return residuals
#' @export
residuals.remoteAR <- function(object, ...){
  stopifnot("remoteAR" %in% class(object))

  ## Pixel version
  if ("pixel" %in% class(object)) {
    return(object$resids)
  } else if ("map" %in% class(object)) {
    return(object$resids)
  }
}

# GLS ----

#' remoteGLS constructor (S3)
#'
#' @param form optional argument specifying the GLS formula
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
                            "predictors" = NULL
                            # "coef.names" = NULL
                            )
  if(!missing(form)){
    GLS.obj$model.info$formula = formula(form)
    GLS.obj$model.info$response = all.vars(formula(form))[1]
    GLS.obj$model.info$predictors = all.vars(formula(form))[-1]
  }

  # GLS.obj$call = match.call()
  # print(match.call(expand.dots = FALSE))

  class(GLS.obj) <- c("remoteGLS")
  attr(GLS.obj, "no_F") = FALSE

  return(GLS.obj) # return the empty list
}

#' print method for remoteGLS
#'
#' @param x remoteGLS object
#' @param ret_call should the function call be returned?
#' @param ... additional arguments
#'
#' @return formatted output for remoteGLS object
#' @export
print.remoteGLS <- function(x, ret_call = TRUE, ...){
  stopifnot("remoteGLS" %in% class(x))
  if(is.null(unlist(x[-1]))){
    cat("empty remoteGLS object\ncall: ")
    if(!is.null(unlist(x[1]))){
      print(x$model.info)
    } else {print(x$model.info$call)}
  } else {

    ## Coefficient table
    coefs = data.frame("Est" = x$betahat, "t stat" = x$tstat, "pval t" = x$pval.t)

    ## Model stats
    if(!attr(x, "no_F")){
      mod.stats = data.frame("df.F" = x$df.F,
                             "SSE" = c(x$SSE, x$SSE0),
                             "MSE" = c(x$MSE, x$MSE0),
                             "logLik" = c(x$logLik, x$logLik0),
                             "F stat" = c(x$Fstat, NA),
                             "pval F" = c(x$pval.F, NA))
      rownames(mod.stats) = c("mod_A", "mod_0")
    } else {
      mod.stats = data.frame("SSE" = x$SSE,
                             "MSE" = x$MSE,
                             "logLik" = x$logLik)
    }
    ## return (needs cleaning)
    if (ret_call){
      cat("call: ");print(x$model.info$call);cat("\n")
    }
    if (!is.null(x$model.info$response)){
      cat("response:", x$model.info$response,"\n\n")
    }
    cat("t tests:\n")
    print(format(coefs, digits = 2, nsmall = 2))

    if(!attr(x, "no_F")){
      cat("\nF test:\n")
    } else {
      cat("\nmodel stats:\n")
    }
    print(format(mod.stats, digits = 2, nsmall = 2, scientific = -2))
  }
}

## summary function
# summary.remoteGLS <- function(object){



## GLS.partition

#' T test of partitioned GLS
#'
#' @param object remoteGLS.parts object
#'
#' @return t-table
#' @export
cor_t.test <- function(object){
  stopifnot("remoteGLS.parts" %in% class(object))

  ## Correlated t-test
  return(cor_t(coefs = object$overall.stats$coefmean,
               part.SEs = object$part.stats$SEs,
               rcoef = object$overall.stats$rcoefmean,
               df2 = object$overall.stats$dfs[2],
               npart = nrow(object$part.stats$coefficients)
               ))
}

#' correlated chi-squared test of partitioned GLS
#'
#' @param object remoteGLS.parts object
#'
#' @return p value
#' @export
cor_chisq.test <- function(object){
  stopifnot("remoteGLS.parts" %in% class(object))
  return(cor_chisq(Fmean = object$overall.stats$meanstats["fmean"],
                   rSSR = object$overall.stats$meanstats["rSSRmean"],
                   df1 = object$overall.stats$dfs[1],
                   npart = nrow(object$part.stats$coefficients))
  )
}

#' correlated F test of partitioned GLS
#'
#' @param object remoteGLS.parts object
#' @param nboot number of bootstrap iterations
#'
#' @return P value
#' @export
cor_F.test <- function(object, nboot = 1000){
  stopifnot("remoteGLS.parts" %in% class(object))
  test = boot_corF(Fmean.obs = object$overall.stats$meanstats["fmean"],
            rSSR = object$overall.stats$meanstats["rSSRmean"],
            rSSE = object$overall.stats$meanstats["rSSEmean"],
            df1 = object$overall.stats$dfs[1],
            df2 = object$overall.stats$dfs[2],
            npart = nrow(object$part.stats$coefficients),
            nboot = nboot)
  return(c("pval.F" = test$pvalue))
}

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

