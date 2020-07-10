## Original ----
# original CLS with no scaled time variable
CLS.fit.U <- function(X, U=NULL) {

  n <- dim(X)[1]

  # CLS for entire map
  d <- data.frame(site = 1:n)
  for (i in 1:dim(X)[1]) {
    x <- X[i, ]

    d$mean[i] <- mean(x)

    if(is.null(U)){
      z.CLS <- lm(x[-1] ~ x[-length(x)])
      d$b0[i] <- summary(z.CLS)$coef[1, 1]
      d$b[i] <- summary(z.CLS)$coef[2, 1]
      d$MSE[i] <- summary(z.CLS)$sigma^2
    }else{
      u <- U[i, ]
      z.CLS <- lm(x[-1] ~ x[-length(x)] + u[-length(x)])
      d$b0[i] <- summary(z.CLS)$coef[1, 1]
      d$b[i] <- summary(z.CLS)$coef[2, 1]
      d$c[i] <- summary(z.CLS)$coef[3, 1]
      d$MSE[i] <- summary(z.CLS)$sigma^2
    }
  }
  return(d)
}

CLS.fit <- function(X, t.scale) {

  n <- dim(X)[1]

  # CLS for entire map
  d <- data.frame(site = 1:n)
  for (i in 1:dim(X)[1]) {
    x <- X[i, ]

    d$mean[i] <- mean(x)

    z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
    d$c[i] <- summary(z.CLS)$coef[3, 1]
    d$t[i] <- summary(z.CLS)$coef[3, 3]
    d$p[i] <- summary(z.CLS)$coef[3, 4]
    d$b[i] <- summary(z.CLS)$coef[2, 1]
    d$MSE[i] <- summary(z.CLS)$sigma^2
  }
  return(d)
}

LS.fit <- function(X, t.scale) {

  n <- dim(X)[1]

  # LS for entire map
  d <- data.frame(site = 1:n)
  for (i in 1:dim(X)[1]) {
    x <- X[i, ]

    d$mean[i] <- mean(x)

    z.LS <- lm(x ~ t.scale)
    d$c.LS[i] <- summary(z.LS)$coef[2, 1]
    d$p.LS[i] <- summary(z.LS)$coef[2, 4]
  }
  return(d)
}

