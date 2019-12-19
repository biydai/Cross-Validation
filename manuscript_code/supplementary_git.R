#--------------------------------------
# R code that generates the simluation 
# in the supplementary material:
# compares different baseline construction
# for CV deviance residuals
#--------------------------------------

library(ncvreg)
library(survival)

generatebeta <- function(p,nbeta,c = 1){
  pos <- rep(c, ceiling(nbeta/2))
  neg <- rep(-c, floor(nbeta/2))
  beta <- c(pos,neg,rep(0,p-nbeta))
}

genSurvData <- function(n = 100, beta = c(1,-1),cpct = 0.70,sd = 1,h = 1 # 70% means 70% of the data are censored; 30% are observed
){
  p <- length(beta)
  ## Generate column vector of coefficients
  beta <- matrix(beta, ncol = 1)
  x <- matrix(rnorm(n*p, sd = sd),n,p)
  #standardize x first
  x <- ncvreg::std(x)
  hazard <- h*exp(x %*% beta)
  ## Calculate survival times
  time <- rexp(n, hazard)
  ## observations are censored ("event = rbinom(n, 1, 0.30)")
  y <- Surv(time = time, event = rbinom(n, 1, 1-cpct))
  return(list(x=x,y=y,beta = beta))
}

#---------------------------------------------------------------------------------------------------------
# simulate data

n <- 150
p <- 1000
c <- 4:9/10
rep <- 200

result_lambda <- array(NA,dim = c(6,length(c),rep),dimnames = list(c("NA","Full","Full_breslow","CV","LP","Oracle"),
                                                                     as.character(c),
                                                                     NULL))
result_MSER <- array(NA,dim = c(6,length(c),rep),dimnames = list(c("NA","Full","Full_breslow","CV","LP","Oracle"),
                                                                 as.character(c),
                                                                 NULL))
set.seed(1)
for(R in 1:rep){
  for(j in 1:length(c)){
    beta <- generatebeta(p,nbeta = 10 ,c = c[j])
    data <- genSurvData(n=n , beta = beta, cpct = 0.1, sd = 1)

    #initial fit
    fit <- ncvsurv(X = data$x, y = data$y, penalty = "lasso")
    lambda <- fit$lambda
    l <- length(lambda)
    
    # oracle fit
    fit_oracle <- numeric(p)
    fit_oracle[1:10] <- coxph(data$y ~ data$x[,1:10])$coef
    
    MSE_oracle <- mean((fit_oracle - beta)^2)
    MSE_all <- apply(fit$beta,2,function(x){mean((x - beta)^2)})
    log_MSE_ratio <- log(MSE_all/MSE_oracle)
    
    # allocation for cross validation
    foldid <- numeric(n)
    eta.cv <- matrix(nrow = n, ncol = l)
    colnames(eta.cv) <- lambda
    foldid[data$y[,2] == 1] <- sample(1:10, sum(data$y[,2]), replace = TRUE)
    foldid[data$y[,2] == 0] <- sample(1:10, sum(1 - data$y[,2]), replace = TRUE)
    ind <- order(data$y[,1])
    d <- as.numeric(data$y[ind, 2])
    
    # get cross-validated linear predictors
    for(i in 1:10){
      ifold <- foldid == i
      cvfit <- ncvsurv(X = data$x[!ifold,], y = data$y[!ifold,], lambda = lambda, penalty = "lasso")
      cvlambda <- cvfit$lambda
      eta.cv[ifold,lambda%in%cvlambda] <- predict(cvfit,X = data$x[ifold,], type = "link")
    }
    
    eta.cv <- eta.cv[,apply(eta.cv,2,function(x){sum(is.na(x)) == 0})]
    ncol(eta.cv)
    w.cv <- exp(eta.cv[ind, , drop = FALSE])
    
    #-------------------------------------------------------
    # completely non-parametric (nelson aalen estimator)
    KM <- survfit(data$y~1)
    BH_NA <- cumsum(KM$n.event/KM$n.risk)
    
    val_NA <- apply(w.cv, 2, function(x){ 
      M <- d - BH_NA*x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      sum(dev^2)})
    
    #-------------------------------------------------------
    # with eta from full fit: KP
    ind <- order(data$y[,1])
    d <- as.numeric(data$y[ind, 2])  
    Eta <- predict(fit,data$x,type = "link")
    w <- exp(Eta[ind, , drop = FALSE])
    BH_full <- apply(w, 2, function(x){ 
      r <- rev(cumsum(rev(x)))
      a <- ifelse(d, (1 - x/r)^(1/x), 1)
      cumsum(1-a)
      }
    )
    BH_full <- BH_full[,1:ncol(w.cv)]
    
    val_full <- apply(w.cv*BH_full, 2, function(x){ 
      M <- d - x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      sum(dev^2)})
    
    #-------------------------------------------------------
    # with eta from full fit: KP
    BH_full_breslow <- apply(w, 2, function(x){ 
      r <- rev(cumsum(rev(x)))
      h <- d/r
      cumsum(h)
    }
    )
    BH_full_breslow <- BH_full_breslow[,1:ncol(w.cv)]
    
    val_full_breslow <- apply(w.cv*BH_full_breslow, 2, function(x){ 
      M <- d - x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      sum(dev^2)})
    
    
    #-------------------------------------------------------
    # with cross-validated eta (old KP)
    
    val_cv <- apply(w.cv, 2, function(x){ 
      r <- rev(cumsum(rev(x)))
      a <- ifelse(d, (1 - x/r)^(1/x), 1)
      BH0 <- cumsum(1-a)
      M <- d - BH0*x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      #dev <- (- M - ifelse(d == 0, 0, d * log(d - M)))
      sum(dev^2)})
    
    #-------------------------------------------------------
    # partial likelihood built from linear predictors
    eta.cv.ordered <- eta.cv[ind,,drop = FALSE]
    r <- apply(eta.cv.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
    pl <- -2 * (crossprod(d, eta.cv.ordered) - crossprod(d, log(r)))
    colnames(pl) <- colnames(eta.cv)

#-------------------------------------------------------
# extract results
  id <- c(which.min(val_NA), which.min(val_full), which.min(val_full_breslow),which.min(val_cv), which.min(pl),
          which.min(log_MSE_ratio))
  result_lambda[,j,R] <- lambda[id]
    
  result_MSER[,j,R] <- log_MSE_ratio[id]
    
  }
}

apply(result_lambda,c(1,2),mean)
apply(result_MSER,c(1,2),mean)

