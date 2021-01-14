## remotePARTS objects (S3) and methods

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
# summary.remoteGLS <- function(obj){}

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

