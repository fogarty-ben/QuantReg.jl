# Fixing errors in quantreg w/ CI size and dyhat before QR factorization

library(quantreg)

rq.fit.br <- function (x, y, tau = 0.5, alpha = 0.1, ci = FALSE, iid = TRUE, 
                       interp = TRUE, tcrit = TRUE) 
{
  tol <- .Machine$double.eps^(2/3)
  eps <- tol
  big <- .Machine$double.xmax
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  ny <- NCOL(y)
  nsol <- 2
  ndsol <- 2
  if (qr(x)$rank < p) 
    stop("Singular design matrix")
  if (tau < 0 || tau > 1) {
    nsol <- 3 * n
    ndsol <- 3 * n
    lci1 <- FALSE
    qn <- rep(0, p)
    cutoff <- 0
    tau <- -1
  }
  else {
    if (p == 1) 
      ci <- FALSE
    if (ci) {
      lci1 <- TRUE
      if (tcrit) 
        cutoff <- qt(1 - alpha/2, n - p)
      else cutoff <- qnorm(1 - alpha/2)
      if (!iid) {
        h <- bandwidth.rq(tau, n, hs = T, alpha=alpha)
        if (tau + h > 1) 
          stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
          stop("tau - h < 0:  error in summary.rq")
        bhi <- rq.fit.br(x, y, tau + h, ci = FALSE)
        bhi <- coefficients(bhi)
        blo <- rq.fit.br(x, y, tau - h, ci = FALSE)
        blo <- coefficients(blo)
        dyhat <- x %*% (bhi - blo)
        if (any(dyhat <= 1e-7)) 
          warning(paste(sum(dyhat <= 1e-7), "non-positive fis"))
        dyhat[dyhat <= 1e-7] = 0
        f <- pmax(eps, (2 * h)/(dyhat - eps))
        qn <- rep(0, p)
        for (j in 1:p) {
          qnj <- lm(x[, j] ~ x[, -j] - 1, weights = f)$resid
          qn[j] <- sum(qnj * qnj)
        }
      }
      else qn <- 1/diag(solve(crossprod(x)))
    }
    else {
      lci1 <- FALSE
      qn <- rep(0, p)
      cutoff <- 0
    }
  }
  z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n + 
                                                                   5), as.integer(p + 3), as.integer(p + 4), as.double(x), 
                as.double(y), as.double(tau), as.double(tol), flag = as.integer(1), 
                coef = double(p), resid = double(n), integer(n), double((n + 
                                                                           5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol), 
                sol = double((p + 3) * nsol), dsol = double(n * ndsol), 
                lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn), 
                cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 * 
                                                                                 p), as.double(big), as.logical(lci1))
  if (z$flag != 0) 
    warning(switch(z$flag, "Solution may be nonunique", 
                   "Premature end - possible conditioning problem in x"))
  if (tau < 0 || tau > 1) {
    sol <- matrix(z$sol[1:((p + 3) * z$lsol)], p + 3)
    dsol <- matrix(z$dsol[1:(n * z$lsol)], n)
    vnames <- dimnames(x)[[2]]
    dimnames(sol) <- list(c("tau", "Qbar", "Obj.Fun", vnames), 
                          NULL)
    return(list(sol = sol, dsol = dsol))
  }
  if (!ci) {
    coef <- z$coef
    dual <- z$dsol[1:n]
    names(coef) <- dimnames(x)[[2]]
    return(list(coefficients = coef, x = x, y = y, residuals = y - 
                  x %*% z$coef, dual = dual))
  }
  if (interp) {
    Tn <- matrix(z$tnmat, nrow = 4)
    Tci <- matrix(z$ci, nrow = 4)
    Tci[3, ] <- Tci[3, ] + (abs(Tci[4, ] - Tci[3, ]) * (cutoff - 
                                                          abs(Tn[3, ])))/abs(Tn[4, ] - Tn[3, ])
    Tci[2, ] <- Tci[2, ] - (abs(Tci[1, ] - Tci[2, ]) * (cutoff - 
                                                          abs(Tn[2, ])))/abs(Tn[1, ] - Tn[2, ])
    Tci[2, ][is.na(Tci[2, ])] <- -big
    Tci[3, ][is.na(Tci[3, ])] <- big
    coefficients <- cbind(z$coef, t(Tci[2:3, ]))
    vnames <- dimnames(x)[[2]]
    cnames <- c("coefficients", "lower bd", "upper bd")
    dimnames(coefficients) <- list(vnames, cnames)
    residuals <- y - drop(x %*% z$coef)
    return(list(coefficients = coefficients, residuals = residuals))
  }
  else {
    Tci <- matrix(z$ci, nrow = 4)
    coefficients <- cbind(z$coef, t(Tci))
    residuals <- y - drop(x %*% z$coef)
    vnames <- dimnames(x)[[2]]
    cnames <- c("coefficients", "lower bound", "Lower Bound", 
                "upper bd", "Upper Bound")
    dimnames(coefficients) <- list(vnames, cnames)
    c.values <- t(matrix(z$tnmat, nrow = 4))
    c.values <- c.values[, 4:1]
    dimnames(c.values) <- list(vnames, cnames[-1])
    p.values <- if (tcrit) 
      matrix(pt(c.values, n - p), ncol = 4)
    else matrix(pnorm(c.values), ncol = 4)
    dimnames(p.values) <- list(vnames, cnames[-1])
    list(coefficients = coefficients, residuals = residuals, 
         c.values = c.values, p.values = p.values)
  }
}

summary.rq <- function (object, se = NULL, covariance = FALSE, hs = TRUE, 
                        U = NULL, gamma = 0.7, ...) 
{
  if (object$method == "lasso") 
    stop("no inference for lasso'd rq fitting: try rqss (if brave, or credulous)")
  mt <- terms(object)
  m <- model.frame(object)
  y <- model.response(m)
  dots <- list(...)
  if (object$method == "sfn") 
    x <- object$model$x
  else x <- model.matrix(mt, m, contrasts = object$contrasts)
  wt <- as.vector(model.weights(object$model))
  tau <- object$tau
  eps <- .Machine$double.eps^(1/2)
  coef <- coefficients(object)
  if (is.matrix(coef)) 
    coef <- coef[, 1]
  vnames <- dimnames(x)[[2]]
  resid <- object$residuals
  n <- length(resid)
  p <- length(coef)
  rdf <- n - p
  if (!is.null(wt)) {
    resid <- resid * wt
    x <- x * wt
    y <- y * wt
  }
  if (is.null(se)) {
    if (n < 1001 & covariance == FALSE) 
      se <- "rank"
    else se <- "nid"
  }
  if (se == "rank") {
    f <- rq.fit.br(x, y, tau = tau, ci = TRUE, ...)
  }
  if (se == "iid") {
    xxinv <- diag(p)
    xxinv <- backsolve(qr(x)$qr[1:p, 1:p, drop = FALSE], 
                       xxinv)
    xxinv <- xxinv %*% t(xxinv)
    pz <- sum(abs(resid) < eps)
    h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs, alpha = alpha)))
    ir <- (pz + 1):(h + pz + 1)
    ord.resid <- sort(resid[order(abs(resid))][ir])
    xt <- ir/(n - p)
    sparsity <- rq(ord.resid ~ xt)$coef[2]
    cov <- sparsity^2 * xxinv * tau * (1 - tau)
    scale <- 1/sparsity
    serr <- sqrt(diag(cov))
  }
  else if (se == "nid") {
    h <- bandwidth.rq(tau, n, hs = hs, alpha = alpha)
    if (tau + h > 1) 
      stop("tau + h > 1:  error in summary.rq")
    if (tau - h < 0) 
      stop("tau - h < 0:  error in summary.rq")
    bhi <- rq.fit.fnb(x, y, tau = tau + h)$coef
    blo <- rq.fit.fnb(x, y, tau = tau - h)$coef
    dyhat <- x %*% (bhi - blo)
    if (any(dyhat <= 1e-7)) 
      warning(paste(sum(dyhat <= 1e-7), "non-positive fis"))
    dyhat[dyhat <= 1e-7] = 0
    f <- pmax(0, (2 * h)/(dyhat - eps))
    fxxinv <- diag(p)
    fxxinv <- backsolve(qr.R(qr(sqrt(f) * x))[1:p, 1:p, drop = FALSE], 
                        fxxinv)
    fxxinv <- fxxinv %*% t(fxxinv)
    cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% 
      fxxinv
    scale <- mean(f)
    
    serr <- sqrt(diag(cov))
  }
  else if (se == "ker") {
    h <- bandwidth.rq(tau, n, hs = hs, alpha=alpha)
    if (tau + h > 1) 
      stop("tau + h > 1:  error in summary.rq")
    if (tau - h < 0) 
      stop("tau - h < 0:  error in summary.rq")
    uhat <- c(y - x %*% coef)
    h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
                                                 (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
    f <- dnorm(uhat/h)/h
    fxxinv <- diag(p)
    fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE], 
                        fxxinv)
    fxxinv <- fxxinv %*% t(fxxinv)
    cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*% 
      fxxinv
    scale <- mean(f)
    serr <- sqrt(diag(cov))
  }
  else if (se == "boot") {
    if ("cluster" %in% names(dots)) {
      bargs <- modifyList(list(x = x, y = y, tau = tau), 
                          dots)
      if (length(object$na.action)) {
        cluster <- dots$cluster[-object$na.action]
        bargs <- modifyList(bargs, list(cluster = cluster))
      }
      B <- do.call(boot.rq, bargs)
    }
    else B <- boot.rq(x, y, tau, ...)
    cov <- cov(B$B)
    serr <- sqrt(diag(cov))
  }
  else if (se == "BLB") {
    n <- length(y)
    b <- ceiling(n^gamma)
    S <- n%/%b
    V <- matrix(sample(1:n, b * S), b, S)
    Z <- matrix(0, NCOL(x), S)
    for (i in 1:S) {
      v <- V[, i]
      B <- boot.rq(x[v, ], y[v], tau, bsmethod = "BLB", 
                   blbn = n, ...)
      Z[, i] <- sqrt(diag(cov(B$B)))
    }
    cov <- cov(B$B)
    serr <- apply(Z, 1, mean)
  }
  if (se == "rank") {
    coef <- f$coef
  }
  else {
    coef <- array(coef, c(p, 4))
    dimnames(coef) <- list(vnames, c("Value", "Std. Error", 
                                     "t value", "Pr(>|t|)"))
    coef[, 2] <- serr
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- if (rdf > 0) 
      2 * (1 - pt(abs(coef[, 3]), rdf))
    else NA
  }
  object <- object[c("call", "terms")]
  if (covariance == TRUE) {
    object$cov <- cov
    if (se == "iid") 
      object$scale <- scale
    if (se %in% c("nid", "ker")) {
      object$Hinv <- fxxinv
      object$J <- crossprod(x)
      object$scale <- scale
    }
    else if (se == "boot") {
      object$B <- B$B
      object$U <- B$U
    }
  }
  object$coefficients <- coef
  object$residuals <- resid
  object$rdf <- rdf
  object$tau <- tau
  class(object) <- "summary.rq"
  object
}

# Generating results sets as individual files

data(mtcars)
for (tau in c(0.25, 0.5, 0.75)) {
  model <- rq(mpg ~ disp + hp + cyl, data=mtcars, tau=tau)
  for (invers in c(T, F)) {
    for (alpha in c(0.01, 0.05, 0.10)) {
      for (hs in c(T, F)) {
        for (iid in c(T, F)) {
          for (interp in c(T, F)) {
            se <- if(invers) "rank" else (if(iid) "iid" else("nid"))
            fn <- paste0(paste("tau", tau, "invers", invers, "alpha", alpha, "hs", hs, "iid", iid, "interp", interp, sep="_"), ".txt")
            tryCatch({
              model.summary <- summary(model, se=se, alpha=alpha, hs=hs, iid=iid, interp=interp)
              if (invers & interp){
                write.table(format(unname(model.summary[["coefficients"]][1:4, 2:3]), nsmall=10, scientific=F), file=fn, row.names=F, col.names=F, sep="\t")
              } else if (invers & !interp) {
                write.table(format(unname(model.summary[["coefficients"]][1:4, 2:5]), nsmall=10, scientific=F), file=fn, row.names=F, col.names=F, sep="\t")
              } else {
                write.table(format(unname(model.summary[["coefficients"]][1:4, 2:2]), nsmall=10, scientific=F), file=fn, row.names=F, col.names=F, sep="\t")
              }
            }, error = function(e) {print(fn); e})
          }
        }
      }
    }
  }
}