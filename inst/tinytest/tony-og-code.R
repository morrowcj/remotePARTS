# from GarrettMooney/moonmisc: Personal Utility Functions

# memory reporting function
.ls.objects <- function(pos = 1, pattern, order.by, decreasing = FALSE, head = FALSE, n = 5) {
	napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.prettysize <- napply(names, function(x) {
		format(utils::object.size(x), units = "auto")
	})
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
	names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing = decreasing), ]
	if (head)
		out <- head(out, n)
	out
}

# shorthand
lsos <- function(..., n = 10) {
	.ls.objects(..., order.by = "Size", decreasing = TRUE, head = TRUE, n = n)
}



AR.reml <- function(formula, data = list()) {

	AR.reml.funct <- function(par, x, u) {
		b <- par
		n.obs <- length(x)
		q <- dim(u)[2]
		B <- diag(n.obs)
		diag(B[-1, ]) <- -b

		iS <- diag(n.obs)
		iS[1, 1] <- (1 - b^2)
		iV <- t(B) %*% iS %*% B
		logdetV <- -determinant(iV)$modulus[1]

		beta <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
		H <- x - u %*% beta

		s2 <- (t(H) %*% iV %*% H)/(n.obs - q)
		LL <- 0.5 * ((n.obs - q) * log(s2) + logdetV + determinant(t(u) %*% iV %*% u)$modulus[1] + (n.obs - q))
		#show(c(LL,b))
		return(LL)
	}


	mf <- model.frame(formula = formula, data = data)
	u <- model.matrix(attr(mf, "terms"), data = mf)
	x <- model.response(mf)

	q <- dim(u)[2]

	opt <- optim(fn = AR.reml.funct, par = 0.2, method = "Brent", upper = 1, lower = -1, control = list(maxit = 10^4), x = x, u = u)
	b <- opt$par

	n.obs <- length(x)
	q <- dim(u)[2]
	B <- diag(n.obs)
	diag(B[-1, ]) <- -b

	iS <- diag(n.obs)
	iS[1, 1] <- (1 - b^2)
	iV <- t(B) %*% iS %*% B
	logdetV <- -determinant(iV)$modulus[1]

	beta <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
	H <- x - u %*% beta

	MSE <- as.numeric((t(H) %*% iV %*% H)/(n.obs - q))
	s2beta <- MSE * solve(t(u) %*% iV %*% u)
	Pr <- 1:q
	for (i in 1:q) Pr[i] <- 2 * pt(abs(beta[i])/s2beta[i, i]^0.5, df = n.obs - q, lower.tail = F)

	logLik <- 0.5 * (n.obs - q) * log(2 * pi) + determinant(t(u) %*% u)$modulus[1] - opt$value

	return(list(beta = beta, b = b, MSE = MSE, s2beta = s2beta, Pr = Pr, logLik = logLik))
}


#################################################
simX <- function(formula, data = data.frame(rep(1, n)), coef, b, s, Dr = NULL, t.scale, n, n.obs, n.burn, seed = 0) {
	set.seed(seed = seed)
	mf <- model.frame(formula = formula, data = data)
	u <- model.matrix(attr(mf, "terms"), data = mf)
	if (!is.matrix(coef))
		coef <- matrix(coef, ncol = 1)

	if (nrow(coef) != ncol(u)) {
		stop("Length of coef must equal the number of independent variables (including the intercept).")
	}

	XX <- matrix(0, nrow = n, ncol = n.obs)
	x <- matrix(0, nrow = n, ncol = 1)
	d <- 0
	for (t in 1:(n.burn + n.obs)) {
		if (is.null(Dr)) {
			e <- rnorm(n, sd = s)
		} else {
			e <- Dr %*% rnorm(n, sd = s)
		}
		if (t <= n.burn) {
			d <- b * d + e
			x <- d
		} else {
			d <- b * d + e
			x <- t.scale[t - n.burn] * as.numeric(u %*% coef) + d
		}
		if (t > n.burn)
			XX[, t - n.burn] <- as.matrix(x)
	}
	XX <- XX - rowMeans(XX)
	return(XX)
}

#################################################
simX_shock <- function(formula, data = data.frame(rep(1, n)), coef, b, s, Dr = NULL, t.scale, n, n.obs, n.burn, seed = 0, shock = 0, n.shock = 0) {
	set.seed(seed = seed)
	mf <- model.frame(formula = formula, data = data)
	u <- model.matrix(attr(mf, "terms"), data = mf)
	if (!is.matrix(coef))
		coef <- matrix(coef, ncol = 1)

	if (nrow(coef) != ncol(u)) {
		stop("Length of coef must equal the number of independent variables (including the intercept).")
	}

	XX <- matrix(0, nrow = n, ncol = n.obs)
	x <- matrix(0, nrow = n, ncol = 1)
	d <- 0
	for (t in 1:(n.burn + n.obs)) {
		if (is.null(Dr)) {
			e <- rnorm(n, sd = s)
		} else {
			e <- Dr %*% rnorm(n, sd = s)
		}

		if (t <= n.burn) {
			d <- b * d + e
			x <- d
		} else {
			d <- b * d + e
			if (t == (n.burn + n.shock)) {
				d <- d + shock
			}
			x <- t.scale[t - n.burn] * as.numeric(u %*% coef) + d
		}
		if (t > n.burn)
			XX[, t - n.burn] <- as.matrix(x)
	}
	return(XX)
}


#################################################
#CLS with no scaled time variable
CLS.fit.U <- function(X, U = NULL) {

	n <- dim(X)[1]

	# CLS for entire map
	d <- data.frame(site = 1:n)
	for (i in 1:dim(X)[1]) {
		x <- X[i, ]

		d$mean[i] <- mean(x)

		if (is.null(U)) {
			z.CLS <- lm(x[-1] ~ x[-length(x)])
			d$b0[i] <- summary(z.CLS)$coef[1, 1]
			d$b[i] <- summary(z.CLS)$coef[2, 1]
			d$MSE[i] <- summary(z.CLS)$sigma^2
		} else {
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
		d$c.ls[i] <- summary(z.LS)$coef[2, 1]
		d$p.ls[i] <- summary(z.LS)$coef[2, 4]
	}
	return(d)
}

MK.fit <- function(X, t.scale) {

	n <- dim(X)[1]

	# MK for entire map
	d <- data.frame(site = 1:n)
	for (i in 1:dim(X)[1]) {
		x <- X[i, ]

		z.MK <- cor.test(x, t.scale, method = "kendall")

		d$c.mk[i] <- z.MK$estimate
		d$p.mk[i] <- z.MK$p.value
	}
	return(d)
}


#################################################
taper.spherical <- function(d, beta) {
	x <- d
	x[d > beta] <- 0
	x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
	x
}

taper.spherical.dif <- function(d, cor, b) {
	beta <- exp(-b)
	x <- d
	x[d > beta] <- 0
	x[d <= beta] <- ((1 - d[d <= beta]/beta)^2) * (1 + d[d <= beta]/(2 * beta))
	cor - x
}

spatialcor.fit <- function(X, t.scale, Dist, r.start = 0.1, fit.n.sample, FUN = "exponential", plot.fig = F, col.plot = NULL) {

	n <- nrow(X)

	# subsample for r.fit
	fit.pick <- sample.int(n = n, size = fit.n.sample)

	resid <- matrix(0, nrow = fit.n.sample, ncol = length(t.scale) - 1) # corrected
	for (i in 1:fit.n.sample) {
		x <- X[fit.pick[i], ]

		z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
		resid[i, ] <- z.CLS$resid
	}
	cor.resid <- cor(t(resid))
	dist <- Dist[fit.pick, fit.pick]/max(Dist)

	# colors for plotting
	if (is.null(col.plot)) {
		col.plot <- "black"
	} else {
		col.plot <- col.plot[fit.pick]
	}


	cor.resid[lower.tri(cor.resid)] <- NA
	dist[lower.tri(dist)] <- NA

	v.cor.resid <- matrix(cor.resid, ncol = 1)
	v.dist <- matrix(dist, ncol = 1)
	v.cor.resid <- v.cor.resid[!is.na(v.cor.resid)]
	v.dist <- v.dist[!is.na(v.dist)]

	w <- as.data.frame(cbind(v.dist, v.cor.resid))
	names(w) <- c("dist", "cor")

	if (FUN == "exponential") {
		fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start), nls.control(maxiter = 500))
		spatialcor <- coef(fit) * max(Dist)
	}
	if (FUN == "exponential-power") {
		fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = 1), nls.control(maxiter = 500))
		spatialcor <- coef(fit) * max(Dist)
	}
	if (FUN == "taper-spherical") {
		fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5), nls.control(maxiter = 500))
		spatialcor <- exp(-coef(fit)) * max(Dist)
	}
	if (plot.fig) {
		plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
		x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
		if (FUN == "exponential")
			lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 1, lwd = 2)
		if (FUN == "exponential-power")
			lines(x.dist, exp(-(x.dist/spatialcor[1])^spatialcor[2]), col = "red", lty = 1, lwd = 2)
		if (FUN == "taper-spherical")
			lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 1, lwd = 2)
	}
	return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}

spatialcor.fit.U <- function(X, U = NULL, Dist, r.start = 0.1, fit.n.sample, FUN = "exponential", plot.fig = F, col.plot = NULL) {

	n <- nrow(X)

	# subsample for r.fit
	fit.pick <- sample.int(n = n, size = fit.n.sample)

	resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
	for (i in 1:fit.n.sample) {
		x <- X[fit.pick[i], ]

		if (is.null(U)) {
			z.CLS <- lm(x[-1] ~ x[-length(x)])
		} else {
			u <- U[i, ]
			z.CLS <- lm(x[-1] ~ x[-length(x)] + u[-length(x)])
		}
		resid[i, ] <- z.CLS$resid
	}
	cor.resid <- cor(t(resid))
	dist <- Dist[fit.pick, fit.pick]/max(Dist)

	# colors for plotting
	if (is.null(col.plot)) {
		col.plot <- "black"
	} else {
		col.plot <- col.plot[fit.pick]
	}

	cor.resid[lower.tri(cor.resid)] <- NA
	dist[lower.tri(dist)] <- NA

	v.cor.resid <- matrix(cor.resid, ncol = 1)
	v.dist <- matrix(dist, ncol = 1)
	v.cor.resid <- v.cor.resid[!is.na(v.cor.resid)]
	v.dist <- v.dist[!is.na(v.dist)]

	w <- as.data.frame(cbind(v.dist, v.cor.resid))
	names(w) <- c("dist", "cor")

	if (FUN == "exponential") {
		fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start), nls.control(maxiter = 500))
		spatialcor <- coef(fit) * max(Dist)
	}
	if (FUN == "taper-spherical") {
		fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5), nls.control(maxiter = 500))
		spatialcor <- exp(-coef(fit)) * max(Dist)
	}
	if (plot.fig) {
		plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
		x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
		if (FUN == "exponential")
			lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 2)
		if (FUN == "taper-spherical")
			lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 2)
	}
	return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma))
}


# NOTE: This returns the spatial correlation scaled to the input distance matrix in km (from distm)
spatialcor.fit.data <- function(X, t.scale, data, r.start = 0.1, a.start = 1, fit.n.sample, FUN = "exponential", plot.fig = F, col.plot = NULL) {

	n <- nrow(X)

	# subsample for r.fit
	fit.pick <- sample.int(n = n, size = fit.n.sample)

	resid <- matrix(0, nrow = fit.n.sample, ncol = n.obs - 1)
	for (i in 1:fit.n.sample) {
		x <- X[fit.pick[i], ]

		z.CLS <- lm(x[2:length(x)] ~ x[1:(length(x) - 1)] + t.scale[2:length(x)])
		resid[i, ] <- z.CLS$resid
	}
	cor.resid <- cor(t(resid))

	# create distance matrix in kilometers
	location <- data[fit.pick, c("lng", "lat")]
	Dist <- geosphere::distm(location, fun = distGeo)/1000
	dist <- Dist/max(Dist)

	# colors for plotting
	if (is.null(col.plot)) {
		col.plot <- "black"
	} else {
		col.plot <- col.plot[fit.pick]
	}

	cor.resid[lower.tri(cor.resid)] <- NA
	dist[lower.tri(dist)] <- NA

	v.cor.resid <- matrix(cor.resid, ncol = 1)
	v.dist <- matrix(dist, ncol = 1)
	v.cor.resid <- v.cor.resid[!is.na(v.cor.resid)]
	v.dist <- v.dist[!is.na(v.dist)]

	w <- as.data.frame(cbind(v.dist, v.cor.resid))
	names(w) <- c("dist", "cor")

	if (FUN == "exponential") {
		fit <- nls(cor ~ exp(-dist/r), data = w, start = list(r = r.start))
		spatialcor <- coef(fit) * max(Dist)
	}
	if (FUN == "exponential-power") {
		#fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
		#browser()
# f <- function(x) {
# r <- x[1]
# a <- x[2]
# return(sum(exp((-w$dist/r)^a) - w$cor)^2)
# }
# fit <- optim(f, par = c(r.start, a.start))

		fit <- nls(cor ~ exp(-(dist/r)^a), data = w, start = list(r = r.start, a = a.start), nls.control(maxiter = 500))
		spatialcor <- coef(fit)
		spatialcor[1] <- spatialcor[1] * max(Dist)
	}
	if (FUN == "taper-spherical") {
		fit <- nls(~taper.spherical.dif(d = dist, cor = cor, b = b), data = w, start = list(b = 0.5))
		spatialcor <- exp(-coef(fit)) * max(Dist)
	}
	if (plot.fig) {
		plot(dist * max(Dist), cor.resid, pch = 20, cex = 0.5, col = col.plot)
		x.dist <- (1:fit.n.sample)/fit.n.sample * max(Dist)
		if (FUN == "exponential")
			lines(x.dist, exp(-x.dist/spatialcor), col = "red", lty = 1, lwd = 2)
		if (FUN == "exponential-power")
			lines(x.dist, exp(-(x.dist/spatialcor[1])^spatialcor[2]), col = "red", lty = 1, lwd = 2)
		if (FUN == "taper-spherical")
			lines(x.dist, taper.spherical(d = x.dist, beta = spatialcor), col = "red", lty = 1, lwd = 2)
	}
	return(list(spatialcor = spatialcor, spatialcor.sigma = summary(fit)$sigma, logLik = logLik(fit)))
}


#################################################
V.fit <- function(Dist, spatialcor, FUN = "exponential") {

	if (FUN == "exponential")
		return(exp(-Dist/spatialcor))

	if (FUN == "exponential-power")
		return(exp(-(Dist/spatialcor[1])^spatialcor[2]))

	if (FUN == "taper-spherical")
		return(taper.spherical(Dist, spatialcor))

}

#################################################
nugget.fit.funct <- function(nugget, formula, data, V, verbose = FALSE) {
	n <- ncol(V)
	invcholV <- t(backsolve(chol((1 - nugget) * V + nugget * diag(n)), diag(n)))
	z <- GLS.fit(formula, data = data, invcholV = invcholV)
	if (verbose == TRUE)
		show(c(z$logLik, nugget))
	return(z$logLik)
}


nugget.fit <- function(formula, data, V, nugget.tol = 1e-05, interval = c(0, 1), verbose = FALSE) {
	opt.nugget <- optimize(nugget.fit.funct, formula, data = data, V = V, interval = interval, maximum = T, tol = nugget.tol, verbose = verbose)
	# check at the zero boundary
	if (opt.nugget$maximum < nugget.tol) {
		nugget0.fit <- nugget.fit.funct(0, formula, data, V)
		if (nugget0.fit > opt.nugget$objective)
			opt.nugget$maximum <- 0
	}
	return(opt.nugget$maximum)
}

#################################################
GLS.fit <- function(formula, formula0 = NULL, data, V = NULL, invcholV = NULL, save.invcholV = F) {

	mf <- model.frame(formula = formula, data = data)
	x <- model.matrix(attr(mf, "terms"), data = mf)
	y <- model.response(mf)
	n <- length(y)

	if (is.null(invcholV)) {
		if (is.null(V)) {
			invcholV <- diag(n)
		} else {
			invcholV <- t(backsolve(chol(V), diag(n)))
		}
	}

	xx <- invcholV %*% x
	yy <- invcholV %*% y

	coef <- as.numeric(solve(crossprod(xx), crossprod(xx, yy)))
	names(coef) <- colnames(x)
	varX <- t(xx) %*% xx
	SSE <- as.numeric(crossprod(yy - xx %*% coef))
	MSE <- SSE/(n - ncol(xx))

	varcov <- MSE * solve(varX)
	se <- diag(varcov)^0.5
	t <- coef/se
	df.t <- n - ncol(xx)
	p.t <- 2 * pt(abs(t), df = df.t, lower.tail = F)

	logdetV <- -2 * sum(log(diag(invcholV)))
	logLik <- -0.5 * (n * log(2 * pi) + n * log((n - ncol(xx)) * MSE/n) + logdetV + n)

	if (is.null(formula0)) {
		x0 <- matrix(1, nrow = n, ncol = 1)
		xx0 <- invcholV %*% x0
		coef0 <- solve(crossprod(xx0), crossprod(xx0, yy))
		SSE0 <- as.numeric(crossprod(yy - xx0 %*% coef0))
		df0 <- ncol(coef0)
	} else {
		mf0 <- model.frame(formula = formula0, data = data)
		x0 <- model.matrix(attr(mf0, "terms"), data = mf0)
		xx0 <- invcholV %*% x0

		if (any(xx0 != 0)) {
			coef0 <- solve(crossprod(xx0), crossprod(xx0, yy))
			SSE0 <- as.numeric(crossprod(yy - xx0 %*% coef0))
			df0 <- ncol(coef0)
		} else {
			SSE0 <- as.numeric(crossprod(yy))
			df0 <- 1
			coef0 <- NA
		}
	}
	MSE0 <- SSE0/(n - ncol(xx))
	MSR <- (SSE0 - SSE)/(ncol(xx) - ncol(xx0))
	logLik0 <- -0.5 * (n * log(2 * pi) + n * log((n - df0) * MSE0/n) + logdetV + n)

	if (any(xx0 != 0)) {
		varX0 <- t(xx0) %*% xx0
		varcov0 <- MSE0 * solve(varX0)
		se0 <- diag(varcov0)^0.5
	} else {
		varcov0 <- NULL
		se0 <- NULL
	}

	if (ncol(xx) > 1) {
		FF <- (n - ncol(xx))/(ncol(xx) - ncol(xx0)) * (SSE0 - SSE)/SSE
		df1.F <- ncol(xx) - ncol(xx0)
		df2.F <- n - ncol(xx)
		p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
		df.F <- c(df1.F, df2.F)
	} else {
		FF <- (n - 1) * (SSE0 - SSE)/SSE
		df1.F <- 1
		df2.F <- n - 1
		p.F <- pf(FF, df1 = df1.F, df2 = df2.F, lower.tail = F)
		df.F <- c(df1.F, df2.F)
	}

	if (!save.invcholV)
		invcholV <- NULL

	return(list(coef = coef, se = se, t = t, df.t = df.t, p.t = p.t, F = FF, df1.F = df1.F, df2.F = df2.F, p.F = p.F, logLik = logLik, logLik0 = logLik0, MSE = MSE, MSE0 = MSE0, MSR = MSR,
		SSE = SSE, SSE0 = SSE0, SSR = SSE0 - SSE, coef0 = coef0, se0 = se0, varX = varX, varcov = varcov, varcov0 = varcov0, invcholV = invcholV, xx = xx, xx0 = xx0, yy = yy))
}


#################################################
GLS.partition.data <- function(formula, formula0 = NULL, data, spatial.autocor.FUN = "exponential-power", spatialcor = spatialcor, est.nugget = T, npart = 10, partition = NULL, nugget.interval = c(0,
	1), fixed.nugget = NULL, nugget.tol = 1e-05, min.num.cross.part = 5, verbose = F, rm.spatial.autocorrelation = F) {

	n <- nrow(data)
	if (!is.null(partition)) {
		npart <- nrow(partition)
		nn <- n - (n%%npart)
		n.p <- nn/npart
		pick <- partition
	} else {
		nn <- n - (n%%npart)
		n.p <- nn/npart
		pick <- matrix(sample(n)[1:nn], nrow = npart)
	}

	mf <- model.frame(formula = formula, data = data)
	df2 <- n.p - (ncol(model.matrix(attr(mf, "terms"), data = mf)) - 1)
	mf0 <- model.frame(formula = formula0, data = data)
	df0 <- n.p - (ncol(model.matrix(attr(mf0, "terms"), data = mf0)) - 1)
	df1 <- df0 - df2

	if (min.num.cross.part > npart)
		min.num.cross.part <- npart
	# NEED TO REMOVE
	# if(ncol(mf) <= 1){
# min.num.cross.part <- NA
# }

	SSR.part <- NULL
	SSE.part <- NULL
	SSE0.part <- NULL
	coef.part <- NULL
	coef0.part <- NULL
	se.part <- NULL
	se0.part <- NULL
	F.part <- NULL
	p.F.part <- NULL
	logLik.part <- NULL
	logLik0.part <- NULL
	nugget.part <- NULL
	invcholV.part <- list(NULL)
	xx.part <- list(NULL)
	xx0.part <- list(NULL)
	for (i in 1:npart) {

		data.part <- data[pick[i, ], ]

		if (rm.spatial.autocorrelation == F) {
			# create distance matrix in kilometers
			location <- data.part[, c("lng", "lat")]
			Dist.part <- geosphere::distm(location, fun = distGeo)/1000

			Vp <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)
			if (is.null(fixed.nugget) & est.nugget) {
				nugget <- nugget.fit(formula, data.part, Vp, interval = nugget.interval, verbose = verbose)
				nugget.interval <- c(0, max(1000 * nugget.tol, min(100 * nugget, 1)))
				Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
			} else {
				if (is.null(fixed.nugget)) {
					nugget <- 0
					Vp <- Vp
				} else {
					nugget <- fixed.nugget[i]
					Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
				}
			}
		} else {
			Vp <- diag(n.p)
		}
		invcholV <- t(backsolve(chol(Vp), diag(n.p)))
		z.part <- GLS.fit(formula, formula0, data = data.part, invcholV = invcholV, save.invcholV = F)

		SSR.part <- c(SSR.part, z.part$SSR)
		SSE.part <- c(SSE.part, z.part$SSE)
		SSE0.part <- c(SSE0.part, z.part$SSE0)
		coef.part <- cbind(coef.part, z.part$coef)
		coef0.part <- cbind(coef0.part, z.part$coef0)
		se.part <- cbind(se.part, z.part$se)
		se0.part <- cbind(se0.part, z.part$se0)
		F.part <- c(F.part, z.part$F)
		p.F.part <- c(p.F.part, z.part$p.F)
		logLik.part <- c(logLik.part, z.part$logLik)
		logLik0.part <- c(logLik0.part, z.part$logLik0)
		nugget.part <- c(nugget.part, nugget)

		if (!is.na(min.num.cross.part) && i <= min.num.cross.part) {
			invcholV.part[[i]] <- invcholV
			xx.part[[i]] <- z.part$xx
			xx0.part[[i]] <- z.part$xx0
		}
		if (verbose) {
			show(paste0("partition ", i, " of ", npart))
			show(z.part$coef)
			show(z.part$p.F)
		}
	}

	#	if(!is.na(min.num.cross.part)){  # NEED TO REMOVE
	########## start NEED TO ADD IN remoteSTAR ##########
	rSSE.part <- matrix(NA, nrow = npart, ncol = npart)
	rSSR.part <- matrix(NA, nrow = npart, ncol = npart)
	rcoef.part <- matrix(NA, nrow = min.num.cross.part * (min.num.cross.part - 1)/2, ncol = nrow(coef.part)) # NEED TO ADD
	counter <- 0
	for (i in 1:(min.num.cross.part - 1)) for (j in (i + 1):min.num.cross.part) {
		counter <- counter + 1
		# create distance matrix in kilometers
		location1 <- data[pick[i, ], c("lng", "lat")]
		location2 <- data[pick[j, ], c("lng", "lat")]
		Dist.part <- geosphere::distm(location1, location2, fun = distGeo)/1000
		V12 <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)
		V12 <- sqrt((1 - nugget.part[i]) * (1 - nugget.part[j])) * V12

		Rij <- crossprod(t(invcholV.part[[i]]), tcrossprod(V12, invcholV.part[[j]]))

		xx1 <- xx.part[[i]]
		xx2 <- xx.part[[j]]
		xx10 <- xx0.part[[i]]
		xx20 <- xx0.part[[j]]

		W1 <- solve(t(xx1) %*% xx1)
		W2 <- solve(t(xx2) %*% xx2)
		H1 <- xx1 %*% W1 %*% t(xx1)
		H2 <- xx2 %*% W2 %*% t(xx2)

		if (!is.na(xx10[1])) {
			H10 <- xx10 %*% solve(t(xx10) %*% xx10) %*% t(xx10)
			H20 <- xx20 %*% solve(t(xx20) %*% xx20) %*% t(xx20)
		} else {
			H10 <- 0
			H20 <- 0
		}

		S1R <- H1 - H10
		S2R <- H2 - H20

		S1E <- diag(n.p) - H1
		S2E <- diag(n.p) - H2

		rSSR.part[i, j] <- matrix(S1R, nrow = 1) %*% matrix(Rij %*% S2R %*% t(Rij), ncol = 1)/df1
		rSSE.part[i, j] <- matrix(S1E, nrow = 1) %*% matrix(Rij %*% S2E %*% t(Rij), ncol = 1)/df2

		Vcoefij <- W1 %*% (t(xx1) %*% Rij %*% xx2) %*% t(W2) # NEED TO ADD
		rcoef.part[counter, ] <- diag(Vcoefij) * (diag(W1) * diag(W2))^-0.5 # NEED TO ADD
		# NOTE: if we ever want to compute conditional se's for the coefficients, we'll need the full Vcoef12 matrix, which would have to be returned as a list of matricies. I don't know if it is worth doing this now just in case.
		}

	rSSR <- mean(rSSR.part, na.rm = T)
	rSSE <- mean(rSSE.part, na.rm = T)
	rcoef <- colMeans(rcoef.part, na.rm = T)
	Fmean <- mean(F.part)
	# } else{
	# rSSR.part <- NA
# rSSE.part <- NA
# rcoef.part <- NA

	# rSSR <- NA
	# rSSE <- NA
# rcoef <- NA

	# Fmean <- NA
	# }

	coef <- rowMeans(coef.part)
	coef0 <- rowMeans(coef0.part)

	return(list(coef = coef, Fmean = Fmean, df1 = df1, df2 = df2, SSR.part = SSR.part, SSE.part = SSE.part, SSE0.part = SSE0.part, logLik.part = logLik.part, logLik0.part = logLik0.part,
		nugget = mean(nugget.part), nugget.part = nugget.part, F.part = F.part, p.F.part = p.F.part, coef.part = coef.part, se.part = se.part, coef0.part = coef0.part, se0.part = se0.part,
		rSSR = rSSR, rSSE = rSSE, rcoef = rcoef, rSSR.part = rSSR.part, rSSE.part = rSSE.part, rcoef.part = rcoef.part, npart = npart, partition = pick, spatial.autocor.FUN = "exponential-power",
		spatialcor = spatialcor)) # NEED TO ADD rcoef and rcoef.part, and maybe a matrix for Vcoefij
	########## end NEED TO ADD IN remoteSTAR ##########
	}



#################################################
# If the nugget is the same for multiple formulae, then they can all be computed quickly
GLS.partition.data.multiformula <- function(formula, formula0 = NULL, data, spatial.autocor.FUN = "exponential-power", spatialcor = spatialcor, est.nugget = T, npart = 10, partition = NULL,
	nugget.interval = c(0, 1), fixed.nugget = NULL, nugget.tol = 1e-04, min.num.cross.part = 5, verbose = F, rm.spatial.autocorrelation = F) {

	n.formula <- length(formula)

	n <- nrow(data)
	if (!is.null(partition)) {
		npart <- nrow(partition)
		nn <- n - (n%%npart)
		n.p <- nn/npart
		pick <- partition
	} else {
		nn <- n - (n%%npart)
		n.p <- nn/npart
		pick <- matrix(sample(n)[1:nn], nrow = npart)
	}

	df0 <- list(numeric())
	df1 <- list(numeric())
	df2 <- list(numeric())

	if (min.num.cross.part > npart)
		min.num.cross.part <- npart

	dummy.part <- as.list(rep(NA, n.formula))
	dummy.part <- lapply(dummy.part, FUN = function(.) NULL)

	SSR.part <- dummy.part
	SSE.part <- dummy.part
	SSE0.part <- dummy.part
	SSE0.part <- dummy.part

	coef.part <- dummy.part
	coef0.part <- dummy.part
	se.part <- dummy.part
	se0.part <- dummy.part
	F.part <- dummy.part
	p.F.part <- dummy.part
	logLik.part <- dummy.part
	logLik0.part <- dummy.part
	xx.part <- dummy.part
	xx0.part <- dummy.part

	nugget.part <- array()
	invcholV.part <- list()

	for (i in 1:npart) {

		data.part <- data[pick[i, ], ]

		if (rm.spatial.autocorrelation == F) {
			# create distance matrix in kilometers
			location <- data.part[, c("lng", "lat")]
			Dist.part <- geosphere::distm(location, fun = distGeo)/1000

			Vp <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)

			if (is.null(fixed.nugget) & est.nugget) {
				nugget <- nugget.fit(formula[[1]], data.part, Vp, interval = nugget.interval, verbose = verbose)
				nugget.interval <- c(0, max(1000 * nugget.tol, min(100 * nugget, 1)))
				Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
			} else {
				if (is.null(fixed.nugget)) {
					nugget <- 0
					Vp <- Vp
				} else {
					nugget <- fixed.nugget[i]
					Vp <- (1 - nugget) * Vp + nugget * diag(n.p)
				}
			}
		} else {
			Vp <- diag(n.p)
		}

		invcholV <- t(backsolve(chol(Vp), diag(n.p)))
		for (i.formula in 1:n.formula) {

			mf <- model.frame(formula = formula[[i.formula]], data = data)
			df2[[i.formula]] <- n.p - (ncol(model.matrix(attr(mf, "terms"), data = mf)) - 1)
			mf0 <- model.frame(formula = formula0[[i.formula]], data = data)
			df0[[i.formula]] <- n.p - (ncol(model.matrix(attr(mf0, "terms"), data = mf0)) - 1)
			df1[[i.formula]] <- df0[[i.formula]] - df2[[i.formula]]

			z.part <- GLS.fit(formula[[i.formula]], formula0[[i.formula]], data = data.part, invcholV = invcholV, save.invcholV = F)

			SSR.part[[i.formula]] <- c(SSR.part[[i.formula]], z.part$SSR)
			SSE.part[[i.formula]] <- c(SSE.part[[i.formula]], z.part$SSE)
			SSE0.part[[i.formula]] <- c(SSE0.part[[i.formula]], z.part$SSE0)
			coef.part[[i.formula]] <- cbind(coef.part[[i.formula]], z.part$coef)
			se.part[[i.formula]] <- cbind(se.part[[i.formula]], z.part$se)
			coef0.part[[i.formula]] <- cbind(coef0.part[[i.formula]], z.part$coef0)
			se0.part[[i.formula]] <- cbind(se0.part[[i.formula]], z.part$se0)
			F.part[[i.formula]] <- c(F.part[[i.formula]], z.part$F)
			p.F.part[[i.formula]] <- c(p.F.part[[i.formula]], z.part$p.F)
			logLik.part[[i.formula]] <- c(logLik.part[[i.formula]], z.part$logLik)
			logLik0.part[[i.formula]] <- c(logLik0.part[[i.formula]], z.part$logLik0)
			if (i <= min.num.cross.part) {
				xx.part[[i.formula]][[i]] <- z.part$xx
				xx0.part[[i.formula]][[i]] <- z.part$xx0
			}
		}

		nugget.part[i] <- nugget
		if (i <= min.num.cross.part) {
			invcholV.part[[i]] <- invcholV
		}

		if (verbose) {
			print(paste0("partition ", i, " of ", npart))
		}
	}

	rSSE.part <- rep(list(matrix(NA, nrow = npart, ncol = npart)), n.formula)
	rSSR.part <- rep(list(matrix(NA, nrow = npart, ncol = npart)), n.formula)
	for (i in 1:(min.num.cross.part - 1)) for (j in (i + 1):min.num.cross.part) {

		# create distance matrix in kilometers
		location1 <- data[pick[i, ], c("lng", "lat")]
		location2 <- data[pick[j, ], c("lng", "lat")]
		Dist.part <- geosphere::distm(location1, location2, fun = distGeo)/1000
		V12 <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)
		# this is only approximate, because nugget.part is esstiamted for each partition, but here the overall mean nugget is used.
		V12 <- sqrt((1 - nugget.part[i]) * (1 - nugget.part[j])) * V12

		xx1 <- xx.part[[i]]
		xx2 <- xx.part[[j]]
		xx10 <- xx0.part[[i]]
		xx20 <- xx0.part[[j]]

		Rij <- crossprod(t(invcholV.part[[i]]), tcrossprod(V12, invcholV.part[[j]]))

		# create distance matrix in kilometers
		location1 <- data[pick[i, ], c("lng", "lat")]
		location2 <- data[pick[j, ], c("lng", "lat")]
		Dist.part <- geosphere::distm(location1, location2, fun = distGeo)/1000
		V12 <- V.fit(Dist.part, spatialcor = spatialcor, FUN = spatial.autocor.FUN)

		# this is only approximate, because nugget.part is estimated for each partition.
		V12 <- sqrt((1 - nugget.part[i]) * (1 - nugget.part[j])) * V12

		Rij <- crossprod(t(invcholV.part[[i]]), tcrossprod(V12, invcholV.part[[j]]))

		for (i.formula in 1:n.formula) {
			xx1 <- xx.part[[i.formula]][[i]]
			xx2 <- xx.part[[i.formula]][[j]]
			xx10 <- xx0.part[[i.formula]][[i]]
			xx20 <- xx0.part[[i.formula]][[j]]

			H1 <- xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1)
			H2 <- xx2 %*% solve(t(xx2) %*% xx2) %*% t(xx2)

			if (!is.na(xx10[1])) {
				H10 <- xx10 %*% solve(t(xx10) %*% xx10) %*% t(xx10)
				H20 <- xx20 %*% solve(t(xx20) %*% xx20) %*% t(xx20)
			} else {
				H10 <- 0
				H20 <- 0
			}

			S1R <- H1 - H10
			S2R <- H2 - H20

			S1E <- diag(n.p) - H1
			S2E <- diag(n.p) - H2

			rSSR.part[[i.formula]][i, j] <- matrix(S1R, nrow = 1) %*% matrix(Rij %*% S2R %*% t(Rij), ncol = 1)/df1[[i.formula]]
			rSSE.part[[i.formula]][i, j] <- matrix(S1E, nrow = 1) %*% matrix(Rij %*% S2E %*% t(Rij), ncol = 1)/df2[[i.formula]]
		}
	}
	coef <- dummy.part
	coef0 <- dummy.part
	rSSR <- dummy.part
	rSSE <- dummy.part
	Fmean <- dummy.part
	for (i.formula in 1:n.formula) {
		coef[[i.formula]] <- rowMeans(coef.part[[i.formula]])
		coef0[[i.formula]] <- rowMeans(coef0.part[[i.formula]])
		rSSR[[i.formula]] <- mean(rSSR.part[[i.formula]], na.rm = T)
		rSSE[[i.formula]] <- mean(rSSE.part[[i.formula]], na.rm = T)

		Fmean[[i.formula]] <- mean(F.part[[i.formula]])
	}

	return.list <- list()
	for (i.formula in 1:n.formula) {

		return.list[[i.formula]] <- list(coef = coef[[i.formula]], Fmean = Fmean[[i.formula]], df1 = df1[[i.formula]], df2 = df2[[i.formula]], SSR.part = SSR.part[[i.formula]], SSE.part = SSE.part[[i.formula]],
			SSE0.part = SSE0.part[[i.formula]], logLik.part = logLik.part[[i.formula]], logLik0.part = logLik0.part[[i.formula]], nugget = mean(nugget.part), nugget.part = nugget.part, F.part = F.part[[i.formula]],
			p.F.part = p.F.part[[i.formula]], coef.part = coef.part[[i.formula]], se.part = se.part[[i.formula]], coef0.part = coef0.part[[i.formula]], se0.part = se0.part[[i.formula]], rSSR = rSSR[[i.formula]],
			rSSE = rSSE[[i.formula]], rSSR.part = rSSR.part[[i.formula]], rSSE.part = rSSE.part[[i.formula]], npart = npart, partition = pick, spatial.autocor.FUN = "exponential-power", spatialcor = spatialcor)

	}
	return(return.list)
}

#################################################
# bootstrap F-test
correlated.F.bootstrap <- function(Fmean.obs, rSSR, rSSE, df1, df2, npart, nboot = 2000) {
	part <- rep(1:npart, each = df1)

	rZ <- rSSR^0.5/df1
	v.MSR <- diag(df1) - rZ
	v.MSR <- kronecker(diag(npart), v.MSR) + rZ
	D.MSR <- t(chol(v.MSR, pivot = T))
	if (attr(D.MSR, "rank") < npart) {
		rank.MSR <- attr(D.MSR, "rank")
		v.MSR <- diag(df1) - 0.99/df1
		v.MSR <- kronecker(diag(npart), v.MSR) + rZ
		D.MSR <- t(chol(v.MSR, pivot = T))
	} else {
		rank.MSR <- NA
	}

	v.MSE <- (1 - rSSE) * diag(npart) + rSSE
	D.MSE <- t(chol(v.MSE))

	count <- 0
	for (boot in 1:nboot) {
		Z1 <- D.MSR %*% rnorm(npart * df1)
		MSR.boot <- aggregate(Z1^2, by = list(part), FUN = sum)[, 2]/df1
		MSE.boot <- 1 + D.MSE %*% rnorm(npart, mean = 0, sd = (2 * df2)^0.5/df2)
		if (Fmean.obs < mean(MSR.boot/MSE.boot))
			count <- count + 1
	}
	return(list(pvalue = count/nboot, nboot = nboot, rank.MSR = rank.MSR))
}

#################################################
# analytical chisq-test
correlated.chisq <- function(Fmean.obs, rSSR, df1, npart) {

	require(CompQuadForm)

	rZ <- rSSR^0.5/df1
	v.MSR <- diag(df1) - rZ
	V.MSR <- kronecker(diag(npart), v.MSR) + rZ
	lambda <- eigen(V.MSR)$values
	pvalue <- suppressMessages(imhof(q = npart * df1 * Fmean.obs, lambda = lambda)$Qq)
	return(pvalue)
}

########## start NEED TO ADD IN remoteSTAR ##########

#################################################
# analytical t-test
correlated.t <- function(coef, se.part, rcoef, df2, npart) {

	secoef <- matrix(NA, length(coef), 1)
	for (i.coef in 1:length(coef)) {
		R <- (1 - rcoef[i.coef]) * diag(npart) + rcoef[i.coef] * matrix(1, npart, npart)
		secoef[i.coef, ] <- (se.part[i.coef, ] %*% R %*% se.part[i.coef, ])^0.5/npart
	}

	tscore <- coef/secoef
	pvalue <- 2 * pt(abs(tscore), df = df2 * npart, lower.tail = FALSE)

	ttest <- cbind(coef, secoef, tscore, pvalue)
	colnames(ttest) <- c("coef", "se", "tscore", "P")

	return(p.t = ttest)
}

#################################################
# wrapper for bootstrap test
GLS.partition.pvalue <- function(z, doFtest = F, nboot = 1e+05) {
	if (doFtest) {
		if (is.finite(z$rSSR)) {
			p.Fmean <- correlated.F.bootstrap(Fmean.obs = z$Fmean, rSSR = z$rSSR, rSSE = z$rSSE, df1 = z$df1, df2 = z$df2, npart = z$npart, nboot = nboot)
		} else {
			p.Fmean <- list(NA, NA, NA)
		}
	} else {
		p.Fmean <- list(NA, NA, NA)
	}

	p.chisq <- correlated.chisq(Fmean.obs = z$Fmean, rSSR = z$rSSR, df1 = z$df1, npart = z$npart)

	p.Fhochberg <- min(p.adjust(z$p.F.part, "hochberg"))
	p.Fhommel <- min(p.adjust(z$p.F.part, "hommel"))
	p.Ffdr <- min(p.adjust(z$p.F.part, "fdr"))

	p.t <- correlated.t(coef = z$coef, se.part = z$se.part, rcoef = z$rcoef, z$df2, npart = z$npart)

	return(list(p.Fmean = p.Fmean, p.chisq = p.chisq, p.Fhochberg = p.Fhochberg, p.Fhommel = p.Fhommel, p.Ffdr = p.Ffdr, p.t = p.t))
}
########## end NEED TO ADD IN remoteSTAR ##########
