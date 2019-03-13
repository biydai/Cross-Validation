#------------------------------
# Self-contained R code
# to reproduce simulation results
# in Figure 4: LOOCV
# scenario: increasing signal
# in n = 120, p = 1000
#------------------------------

library(survival)
library(ncvreg)

#-------------------------------------------------------------------------
# functions

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


#-------------------------------------------------------------------------
# scenarios: 
n <- 120
p <- 1000
N <- 200
c <- c(0.6,0.8,1,1.2)
nbeta <- 10
nfold <- 10
beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = nbeta ,c = c[i])
}

# initiation
gp <- array(NA,c(2,N,length(c)))
ug <- array(NA,c(2,N,length(c)))
ncv <- array(NA,c(2,N,length(c)))
dvr <- array(NA,c(2,N,length(c)))
gp.loo <- array(NA,c(2,N,length(c)))
ncv.loo <- array(NA,c(2,N,length(c)))
dvr.loo <- array(NA,c(2,N,length(c)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(c))

foldid <- numeric(n)


for(i in 1:N){
  set.seed(999 + i)
  for (j in 1:length(c)){
    # generate data
    data <- genSurvData(n=n, beta = beta[j,], cpct = 0.3, sd = 1, h = 1)
    data.test <- genSurvData(n=1000 , beta = beta[j,], cpct = 0, sd = 1, h = 1)
    nevent <- sum(data$y[,2])
    sde <- sqrt(.Machine$double.eps)
    foldid[data$y[,2] == 1] <-  ceiling(sample(1:nevent)/(nevent + sde) * nfold)
    foldid[data$y[,2] == 0] <-  ceiling(sample(1:(n-nevent))/(n - nevent + sde) * nfold)
    
    # initial fit
    fit <- ncvsurv(data$x, data$y,penalty = "lasso")
    lambda <- fit$lambda
    beta_hat <- fit$beta
    norm <- apply(beta_hat,2,function(x){mean((x-beta[j,])^2)})
    
    #oracle model
    data.oracle <- data.frame(data$y, data$x[,1:nbeta])
    fit.oracle <- coxph(data.y ~ ., data = data.oracle)
    oracle.mse[i,j] <- sum((fit.oracle$coef - beta[j,1:nbeta])^2)/p
     
    # cross validation fit
    eta.cv <- matrix(nrow = n, ncol = length(lambda))
    l_diff <- matrix(nrow = nfold , ncol = length(lambda))
    l_test <- matrix(nrow = nfold, ncol = length(lambda))
    
    ind <- order(data$y[,1])
    d <- data$y[ind,2]
    
    # 10 fold for ungrouped
    for(f in 1:nfold){
      ifold <- foldid == f
      cvfit <- ncvsurv(data$x[!ifold,], data$y[!ifold,],penalty = "lasso", lambda = lambda)
      cvlambda <- cvfit$lambda
      beta_train <- cvfit$beta
      
      # linear predictors
      eta.cv[ifold,lambda%in%cvlambda] <- predict(cvfit,X = data$x[ifold,], type = "link")
      
      # grouped
      eta.full.ordered <- apply(beta_train,2,function(b){
        data$x %*% b})[ind,]
      r.full.ordered <- apply(eta.full.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_full <- -2 * (crossprod(d, eta.full.ordered) - crossprod(d, log(r.full.ordered)))
      
      ind_train <- order(data$y[!ifold,1])
      d_train <- data$y[!ifold,2][ind_train]
      eta.train.ordered <- apply(beta_train,2,function(b){
        data$x[!ifold,] %*% b})[ind_train,]
      r.train.ordered <- apply(eta.train.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_train <- -2 * (crossprod(d_train, eta.train.ordered) - crossprod(d_train, log(r.train.ordered)))
      
      l_diff[f,lambda%in%cvlambda] <- l_full - l_train
      
      # ungrouped
      ind_test <- order(data$y[ifold,1])
      d_test <- data$y[ifold,2][ind_test]
      eta.test.ordered <- apply(beta_train,2,function(b){
        e <- data$x[ifold,] %*% b})[ind_test,]
      r.test.ordered <- apply(eta.test.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_test[f,lambda%in%cvlambda] <- -2 * (crossprod(d_test, eta.test.ordered) - crossprod(d_test, log(r.test.ordered)))
    }
    
    # grouped
    val_GP <- na.omit(apply(l_diff,2,mean))
    
    # ungrouped
    val_UG <- na.omit(apply(l_test,2,mean))
    
    # linear predictor
    eta.cv <- eta.cv[,apply(eta.cv,2,function(x){sum(is.na(x)) == 0})]
    
    eta.cv.ordered <- eta.cv[ind,,drop = FALSE]
    r <- apply(eta.cv.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
    val_LP <- -2 * (crossprod(d, eta.cv.ordered) - crossprod(d, log(r)))
    
    # deviance residual
    KM <- survfit(data$y~1)
    BH_NA <- cumsum(KM$n.event/KM$n.risk)
    w.cv <- exp(eta.cv[ind, , drop = FALSE])
    val_NA <- apply(w.cv, 2, function(x){ 
      M <- d - BH_NA*x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      sum(dev^2)})
    
#----------------------------------------------------------------------------------------------    
    
    eta.cv <- matrix(nrow = n, ncol = length(lambda))
    l_diff <- matrix(nrow = n , ncol = length(lambda))
    l_test <- matrix(nrow = n, ncol = length(lambda))
  
    # leave one out for all other three methods
    for(f in 1:n){
      cvfit <- ncvsurv(data$x[-f,], data$y[-f,],penalty = "lasso", lambda = lambda)
      cvlambda <- cvfit$lambda
      beta_train <- cvfit$beta
      
      # linear predictors
      eta.cv[f,lambda%in%cvlambda] <- predict(cvfit,X = data$x[f,], type = "link")
      
      # grouped
      eta.full.ordered <- apply(beta_train,2,function(b){
        data$x %*% b})[ind,]
      r.full.ordered <- apply(eta.full.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_full <- -2 * (crossprod(d, eta.full.ordered) - crossprod(d, log(r.full.ordered)))
      
      ind_train <- order(data$y[-f,1])
      d_train <- data$y[-f,2][ind_train]
      eta.train.ordered <- apply(beta_train,2,function(b){
        data$x[-f,] %*% b})[ind_train,]
      r.train.ordered <- apply(eta.train.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_train <- -2 * (crossprod(d_train, eta.train.ordered) - crossprod(d_train, log(r.train.ordered)))
      
      l_diff[f,lambda%in%cvlambda] <- l_full - l_train
    }
    
    # grouped
    val_GP_loo <- na.omit(apply(l_diff,2,mean))
    
    # linear predictor
    eta.cv <- eta.cv[,apply(eta.cv,2,function(x){sum(is.na(x)) == 0})]
    
    eta.cv.ordered <- eta.cv[ind,,drop = FALSE]
    r <- apply(eta.cv.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
    val_LP_loo <- -2 * (crossprod(d, eta.cv.ordered) - crossprod(d, log(r)))
    
    # deviance residual
    KM <- survfit(data$y~1)
    BH_NA <- cumsum(KM$n.event/KM$n.risk)
    w.cv <- exp(eta.cv[ind, , drop = FALSE])
    val_NA_loo <- apply(w.cv, 2, function(x){ 
      M <- d - BH_NA*x
      dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                               log(d - M))))
      sum(dev^2)})
    
    ###########
    # lambdas #
    ###########
    gp[1,i,j]<- lambda[which.min(val_GP)]
    ug[1,i,j]<- lambda[which.min(val_UG)]
    ncv[1,i,j]<- lambda[which.min(val_LP)]
    dvr[1,i,j]<- lambda[which.min(val_NA)]
    
    gp.loo[1,i,j]<- lambda[which.min(val_GP_loo)]
    ncv.loo[1,i,j]<- lambda[which.min(val_LP_loo)]
    dvr.loo[1,i,j]<- lambda[which.min(val_NA_loo)]
    
    ########
    # mses #
    ########
    gp[2,i,j] <- norm[gp[1,i,j] == lambda]
    ug[2,i,j] <- norm[ug[1,i,j] == lambda]
    ncv[2,i,j] <- norm[ncv[1,i,j] == lambda]
    dvr[2,i,j] <- norm[dvr[1,i,j] == lambda]
    
    gp.loo[2,i,j] <- norm[gp.loo[1,i,j] == lambda]
    ncv.loo[2,i,j] <- norm[ncv.loo[1,i,j] == lambda]
    dvr.loo[2,i,j] <- norm[dvr.loo[1,i,j] == lambda]
    
    print(c(i,j))
  }
}


#-----------------------------------------------------------------------------------
# load saved results and make plots
library(ggplot2)

coef <-  rep(c(0.6,0.8,1,1.2),9)

mse_ratio <- log(c(apply(gp[2,,]/oracle.mse,2,mean),
               apply(ncv[2,,]/oracle.mse,2,mean),
               apply(dvr[2,,]/oracle.mse,2,mean),
               apply(gp.loo[2,,]/oracle.mse,2,mean),
               apply(ncv.loo[2,,]/oracle.mse,2,mean),
               apply(dvr.loo[2,,]/oracle.mse,2,mean),
               apply(ug[2,,]/oracle.mse,2,mean),
               apply(ug[2,,]/oracle.mse,2,mean),
               apply(ug[2,,]/oracle.mse,2,mean)))

methods <- c(rep("V & VH",4),
             rep("Linear Predictor",4),
             rep("Deviance Residuals",4),
             rep("V & VH",4),
             rep("Linear Predictor",4),
             rep("Deviance Residuals",4),
             rep("V & VH",4),
             rep("Linear Predictor",4),
             rep("Deviance Residuals",4)
)

type <- c(rep("10 Fold CV",12),rep("LOOCV",12),rep("Standard",12))

df <- data.frame(mse_ratio,methods,type,coef)

p <- ggplot(df, aes(x=coef, y=mse_ratio, colour = methods)) + 
  geom_line(aes(color = type), stat="identity",size = 1.2)+
  geom_point(aes(color = type), stat="identity",size = 3)+
  facet_grid(.~methods)+
  labs(x = "s",
         #expression(paste("Magnitude of Non-Zero ", beta, "s", sep = ""))
       y = "log(MSE Ratio)")+
  theme_gray()+
  theme(legend.position = "top")

p

