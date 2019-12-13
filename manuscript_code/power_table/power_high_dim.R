#------------------------------
# Self-contained R code
# to reproduce part of the simulations
# in Table 1 of the Cross Validation 
# Manuscript:
# Explored 2 Scenarios here: 
# low dim low signal
# low dim high signal
#------------------------------

library(survival)
library(ncvreg)

args <- as.numeric(commandArgs(TRUE)[1])

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
N <- 10
p <- 10000
n <- 400
c <- c(0.3,0.6)
nbeta <- 20

nfold <- 10
beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = nbeta ,c = c[i])
}

# initiation
gp <- array(NA,c(4,N,length(c)))
ug <- array(NA,c(4,N,length(c)))
ncv <- array(NA,c(4,N,length(c)))
dvr <- array(NA,c(4,N,length(c)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(c))
foldid <- numeric(n)

set.seed(999 + args*100)
for(i in 1:N){
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
    true_positive <- apply(beta_hat,2,function(x){sum(beta[j,]*x != 0)})
    fdr <- apply(beta_hat,2,function(x){sum((beta[j,]== 0)*(x!=0))/max(sum(x!=0),1)})
    
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
      
      l_diff[f,] <- l_full - l_train
      
      # ungrouped
      ind_test <- order(data$y[ifold,1])
      d_test <- data$y[ifold,2][ind_test]
      eta.test.ordered <- apply(beta_train,2,function(b){
        e <- data$x[ifold,] %*% b})[ind_test,]
      r.test.ordered <- apply(eta.test.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_test[f,] <- -2 * (crossprod(d_test, eta.test.ordered) - crossprod(d_test, log(r.test.ordered)))
    }
  
    # grouped
    val_GP <- apply(l_diff,2,mean)
    
    # ungrouped
    val_UG <- apply(l_test,2,mean)
    
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
    
    ###########
    # lambdas #
    ###########
    gp[1,i,j]<- lambda[which.min(val_GP)]
    ug[1,i,j]<- lambda[which.min(val_UG)]
    ncv[1,i,j]<- lambda[which.min(val_LP)]
    dvr[1,i,j]<- lambda[which.min(val_NA)]
  
    ########
    # mses #
    ########
    gp[2,i,j] <- norm[gp[1,i,j] == lambda]
    ug[2,i,j] <- norm[ug[1,i,j] == lambda]
    ncv[2,i,j] <- norm[ncv[1,i,j] == lambda]
    dvr[2,i,j] <- norm[dvr[1,i,j] == lambda]
	
    ##################
    # true positives #
    ##################
    gp[3,i,j] <- true_positive[gp[1,i,j] == lambda]
    ug[3,i,j] <- true_positive[ug[1,i,j] == lambda]
    ncv[3,i,j] <- true_positive[ncv[1,i,j] == lambda]
    dvr[3,i,j] <- true_positive[dvr[1,i,j] == lambda]
    
    
    ##################
    # false discovers
    ##################
    gp[4,i,j] <- fdr[gp[1,i,j] == lambda]
    ug[4,i,j] <- fdr[ug[1,i,j] == lambda]
    ncv[4,i,j] <- fdr[ncv[1,i,j] == lambda]
    dvr[4,i,j] <- fdr[dvr[1,i,j] == lambda]
    
  }
}

results <- list()
for(i in 1:4){
  results[[i]] <- rbind(
    apply(gp[i,,],2,mean),
    apply(ncv[i,,],2,mean),
    apply(ug[i,,],2,mean),
    apply(dvr[i,,],2,mean)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- c
}
names(results) <- c("lambda",
                    "MSE",
                    "tp",
                    "fdr") 

results$lambda
results$MSE
results$tp
results$fdr

save(gp,ncv,dvr,ug,oracle.mse,
     file = paste("high_fdr_", args,".RData",sep = ""))

