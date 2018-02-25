################################
# File that include functions 
# that are used in simulations to: 
# 1. generate data
# 2. run penalized regression models 
# 3. extract simulation results
library(glmnet)
library(survival)
library(ncvreg)

##############################
# function that generates beta
# with p as the length
# nbeta as non-zero betas
##############################

generatebeta <- function(p,nbeta,c = 1){
  pos <- rep(c, ceiling(nbeta/2))
  neg <- rep(-c, floor(nbeta/2))
  beta <- c(pos,neg,rep(0,p-nbeta))
}

#########################################
# function that generates time-to-event
# data based on exponential distribution
#########################################
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

genSurvData_weibull <- function(n = 100, beta = c(1,-1),cpct = 0.70,sd = 1,h = 1,rho = 2 # 70% means 70% of the data are censored; 30% are observed
){
  p <- length(beta)
  ## Generate column vector of coefficients
  beta <- matrix(beta, ncol = 1)
  x <- matrix(rnorm(n*p, sd = sd),n,p)
  #standardize x first
  x <- ncvreg::std(x)
  v <- runif(n)
  time <- (-log(v)/(h*exp(x%*%beta)))^(1/rho)
  
  ## observations are censored ("event = rbinom(n, 1, 0.30)")
  y <- Surv(time = time, event = rbinom(n, 1, 1-cpct))
  return(list(x=x,y=y,beta = beta))
}

# for generating correlated datasets
genX <- function(n, J, K=1, rho=0, rho.g=rho, corr=corr) {
  a <- sqrt(rho/(1-rho.g))
  b <- sqrt((rho.g-rho)/(1-rho.g))
  Z <- rnorm(n)
  ZZ <- t(matrix(rep(rnorm(n*J), rep(K,n*J)), ncol=n))
  ZZZ <- matrix(rnorm(n*J*K),nrow=n)
  return(matrix(as.numeric(a*Z + b*ZZ + ZZZ),nrow=n)/sqrt(1+a^2+b^2))
}

genSurvDataCor <- function(n = 100, beta = c(1,-1),cpct = 0.70,sd = 0.5 # 70% means 70% of the data are censored; 30% are observed
                  # extra argument for genX:
                  ){
  p <- length(beta)
  ## Generate column vector of coefficients
  beta <- matrix(beta, ncol = 1)
  x <- genX(n=n, J=p, K=K, rho=rho, rho.g=rho.g)
  hazard <- exp(x %*% beta)
  ## Calculate survival times
  time <- rexp(n, hazard)
  ## observations are censored ("event = rbinom(n, 1, 0.30)")
  y <- Surv(time = time, event = rbinom(n, 1, 1-cpct))
  return(list(x=x,y=y,beta = beta))
}


#genSurvData <- function(n = 100, beta = c(1,-1),cpct = 0.70 # 70% means 70% of the data are censored; 30% are observed
#){
#  p <- length(beta)
#  ## Generate column vector of coefficients
#  beta <- matrix(beta, ncol = 1)
#  x <- matrix(rnorm(n*p),n,p)
#  hazard <- exp(x %*% beta)
#  ## Calculate survival times
#  #time <- rexp(n, hazard)
#  time <- - log(runif(n))/(1*(hazard))
#  ## observations are censored ("event = rbinom(n, 1, 0.30)")
#  y <- Surv(time = time, event = rbinom(n, 1, 1-cpct))
#  return(list(x=x,y=y,beta = beta))
#}

#data <- genSurvData (n=100 , beta = c(1,-1), cpct = 0.2)
#coxph(Surv(data$y[,1],data$y[,2])~ x,data)
#data <- genSurvData_0 (n=100 , beta = c(1,-1), cpct = 0.2)
#coxph(Surv(data$y[,1],data$y[,2])~ x,data)

#data <- genSurvData_weibull (n=100 , beta = c(1,-1), cpct = 0.2)
#coxph(Surv(data$y[,1],data$y[,2])~ x,data)

#beta <- generatebeta(p = 100, nbeta = 10, c= 15)
#data <- genSurvData (n=100 , beta = c(5,-5), cpct = 0.2)

##############
# Fit GLMNET #
##############
# function that fits survival data using glmnet package
survglmnet <- function(x,y,beta,...){
  # fit glmnet to get the "true" lambda
  fit <- glmnet(data$x, Surv(time = data$y[,1],event = data$y[,2]), family = "cox",...)
  lambda <- as.numeric(fit$lambda)
  beta_hat <- as.matrix(fit$beta)
  norm <- numeric(length(lambda))
  #for (i in 1:length(lambda)){
  #  norm[i] <- mean((beta_hat[,i]-beta)^2)
  #}
  norm <- apply(beta_hat,2,function(x){mean((x-beta)^2)})
  lambda_t <- lambda[which.min(norm)]
  # return:
  # seq of lambda that used by glmnet
  # the "true" lambda
  # distrance between beta_hat and beta
  # the "true" fitted model 
  list <- list(lambda = lambda,
               lambda_t = lambda_t,
               norm = norm,
               fit = fit)
  invisible(list)
}
#
#fit <- survglmnet(data$x,data$y,data$beta)
# got lambda from cv.glmnet:
#cv <- cv.glmnet(data$x,Surv(data$y[,1],data$y[,2]),nfold = 10,lambda = fit$lambda,family = "cox")
#plot(x = log(fit$lambda), y = scale(fit$norm),type = "l")
#lines(x = log(cv$lambda), y = scale(cv$cvm),type = "l",col ="red")

##############
# Fit ncvreg #
##############
survncvreg <- function(x,y,beta,...){
  # fit glmnet to get the "true" lambda
  fit <- ncvsurv(data$x, Surv(data$y[,1],data$y[,2]),penalty = "lasso",...)
  lambda <- as.numeric(fit$lambda)
  beta_hat <- as.matrix(fit$beta)
  norm <- numeric(length(lambda))
  #for (i in 1:length(lambda)){
  #  norm[i] <- mean((beta_hat[,i]-beta)^2)                                                                                                             
  #}
  norm <- apply(beta_hat,2,function(x){mean((x-beta)^2)})
  lambda_t <- lambda[which.min(norm)]
  
  #lambda_aic <- lambda[which.min(AIC(fit))]
  #lambda_bic <- lambda[which.min(BIC(fit))]
  
  aic <- AIC(fit)
  bic <- BIC(fit)
  
  list <- list(lambda = lambda,
               lambda_t = lambda_t,
               aic = aic,
               bic = bic,
               norm = norm,
               fit = fit)
  invisible(list)
}

#coef(fit$fit, s = fit$lambda_t)[1:50]
#cv_ncv$lambda.min
#fit2<- survncvreg(data$x, Surv(data$y[,1],data$y[,2]),beta)
#fit2$lambda_t


