# cross validation using deviance approach

######################################
# Functions Produce Sums of Squares: #
######################################

# Using Kalbfleisch-Prentice's method (discrete estimation)
loss.deviance.KP <- function (y, eta) 
{
  ind <- order(y[, 1])
  d <- as.numeric(y[ind, 2])
  # d, eta, and w are all ordered
  # more than 2 lambda:
  if (is.matrix(eta)) {
    w <- exp(eta[ind, , drop = FALSE])
    val <- apply(w, 2, function(x){ 
      r <- rev(cumsum(rev(x)))
      a <- ifelse(d, (1 - x/r)^(1/x), 1)
      BH0 <- cumsum(1-a)
      M <- d - BH0*x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      #dev <- (- M - ifelse(d == 0, 0, d * log(d - M)))
      sum(dev^2)})
  }
  else {
    #construct baseline hazard
    w <- exp(eta[ind])
    r <- rev(cumsum(rev(w)))
    a <- ifelse(d, (1 - w/r)^(1/w), 1)
    BH0 <- cumsum(1-a)
    M <- d - BH0*w
    dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                             log(d - M))))
    #dev <- (- M - ifelse(d == 0, 0, d * log(d - M)))
    #dev <- - d - d*log(BH0) - d*eta[ind]
    val <- sum(dev^2)
    #val = list(BH0 = BH0, dev = dev)
  }
  return(val)
}


###################################################
# CV functions that takes the baseline estimation #
###################################################
cvf.surv <- function (i, XX, y, cv.ind, cv.args) 
{
  cv.args$X <- XX[cv.ind != i, , drop = FALSE]
  cv.args$y <- y[cv.ind != i, ]
  fit.i <- do.call("ncvsurv", cv.args)
  X2 <- XX[cv.ind == i, , drop = FALSE]
  y2 <- y[cv.ind == i, ]
  nl <- length(fit.i$lambda)
  yhat <- predict(fit.i, X2)
  list(nl = length(fit.i$lambda), yhat = yhat)
}

# KP
cv.ncvsurv.KP <- function (X, y, ..., cluster, nfolds = 10, seed, returnY = FALSE, 
                              trace = FALSE) 
{
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$returnX <- TRUE
  fit <- do.call("ncvsurv", fit.args)
  X <- fit$X
  y <- cbind(fit$time, fit$fail)
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) 
    fit$X <- NULL
  n <- nrow(X)
  if (!missing(seed)) 
    set.seed(seed)
  cv.ind <- ceiling(sample(1:n)/n * nfolds)
  Y <- matrix(NA, nrow = n, ncol = length(fit$lambda))
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$penalty.factor <- fit$penalty.factor
  if (!missing(cluster)) {
    if (!("cluster" %in% class(cluster))) 
      stop("cluster is not of class 'cluster'; see ?makeCluster")
    parallel::clusterExport(cluster, c("cv.ind", "fit", "X", 
                                       "y", "cv.args"), envir = environment())
    parallel::clusterCall(cluster, function() require(ncvreg))
    fold.results <- parallel::parLapply(cl = cluster, X = 1:nfolds, 
                                        fun = cvf.surv, XX = X, y = y, cv.ind = cv.ind, cv.args = cv.args)
  }
  for (i in 1:nfolds) {
    if (!missing(cluster)) {
      res <- fold.results[[i]]
    }
    else {
      if (trace) 
        cat("Starting CV fold #", i, sep = "", "\\n")
      res <- cvf.surv(i, X, y, cv.ind, cv.args)
    }
    Y[cv.ind == i, 1:res$nl] <- res$yhat
  }
  ind <- which(apply(is.finite(Y), 2, all))
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]
  cve <- as.numeric(loss.deviance.KP(y, Y)) #zzzzz
  min <- which.min(cve)
  val <- list(cve = cve, lambda = lambda, fit = fit, min = min, 
              lambda.min = lambda[min], null.dev = cve[1])
  if (returnY) 
    val$Y <- Y
  structure(val, class = c("cv.ncvsurv", "cv.ncvreg"))
}



