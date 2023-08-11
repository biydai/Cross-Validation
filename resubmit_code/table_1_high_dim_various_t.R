#------------------------------
# Self-contained R code
# to reproduce part of the simulations
# in Table 1 of the Cross Validation 
# Manuscript:
# Explored 2 Scenarios here: 
# high dim low signal
# high dim high signal
#------------------------------

id <- as.numeric(commandArgs(TRUE)[1])
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


# Brier score for selected lambda
brierSurv_lambda <- function(data,time,fit,lambda) {
  predicted <- apply(data$x,1, 
                     function(x)
                     {predict(fit, x,type = "survival",
                              lambda = lambda)(time)})
  is_obs_after <- data$y[,1] > time
  mean((is_obs_after - predicted)^2)
}

# Kullback-Leibler for selected lambda
KL_lambda <- function(data,time,fit,lambda) {
  predicted <- apply(data$x,1, 
                     function(x)
                     {predict(fit, x,type = "survival",
                              lambda = lambda)(time)})
  is_obs_after <- data$y[,1] > time
  mean(-is_obs_after*log(predicted) + (is_obs_after-1)*log(1 - predicted),na.rm = TRUE)
}

# overall AUC score for selected lambda
AUC_lambda <- function(data,fit,lambda){
  lp <- predict(fit,X = data$x, type = "link",lambda = lambda)
  survConcordance(data$y ~ lp)$concordance
}

#-------------------------------------------------------------------------
# scenarios: 
N <- 5
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
gp <- array(NA,c(9,N,length(c)))
ug <- array(NA,c(9,N,length(c)))
ncv <- array(NA,c(9,N,length(c)))
dvr <- array(NA,c(9,N,length(c)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(c))
foldid <- numeric(n)


for(i in 1:N){
  set.seed(999*id + i)
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
	
    ####################################
    # Brier Score: 25% quantile #
    ####################################
    
    sfit <- survfit(data.test$y ~ 1)
    t1 <- quantile(sfit$time, 0.25)
    
    gp[3,i,j] <- brierSurv_lambda(data.test,t1,fit,gp[1,i,j])
    ug[3,i,j] <- brierSurv_lambda(data.test,t1,fit,ug[1,i,j])
    ncv[3,i,j] <- brierSurv_lambda(data.test,t1,fit, ncv[1,i,j])
    dvr[3,i,j] <- brierSurv_lambda(data.test,t1,fit,dvr[1,i,j])
    
    
    ####################################
    # Brier Score: median survival time#
    ####################################
    
    sfit <- survfit(data.test$y ~ 1)
    t2 <- summary(sfit)$table["median"]
    
    gp[4,i,j] <- brierSurv_lambda(data.test,t2,fit,gp[1,i,j])
    ug[4,i,j] <- brierSurv_lambda(data.test,t2,fit,ug[1,i,j])
    ncv[4,i,j] <- brierSurv_lambda(data.test,t2,fit, ncv[1,i,j])
    dvr[4,i,j] <- brierSurv_lambda(data.test,t2,fit,dvr[1,i,j])
    
    
    ####################################
    # Brier Score: 75% quantile #
    ####################################
    
    sfit <- survfit(data.test$y ~ 1)
    t3 <- quantile(sfit$time, 0.75)
    
    gp[5,i,j] <- brierSurv_lambda(data.test,t3,fit,gp[1,i,j])
    ug[5,i,j] <- brierSurv_lambda(data.test,t3,fit,ug[1,i,j])
    ncv[5,i,j] <- brierSurv_lambda(data.test,t3,fit, ncv[1,i,j])
    dvr[5,i,j] <- brierSurv_lambda(data.test,t3,fit,dvr[1,i,j])
    
    ###################
    # Kullbeck Leibler: 25% quantile time#
    ##################
    
    gp[6,i,j] <- KL_lambda(data.test,t1,fit,gp[1,i,j])
    ug[6,i,j] <- KL_lambda(data.test,t1,fit,ug[1,i,j])
    ncv[6,i,j] <- KL_lambda(data.test,t1,fit, ncv[1,i,j])
    dvr[6,i,j] <- KL_lambda(data.test,t1,fit,dvr[1,i,j])
    
    ###################
    # Kullbeck Leibler: median survival time#
    ##################
    
    gp[7,i,j] <- KL_lambda(data.test,t2,fit,gp[1,i,j])
    ug[7,i,j] <- KL_lambda(data.test,t2,fit,ug[1,i,j])
    ncv[7,i,j] <- KL_lambda(data.test,t2,fit, ncv[1,i,j])
    dvr[7,i,j] <- KL_lambda(data.test,t2,fit,dvr[1,i,j])
    
    ###################
    # Kullbeck Leibler: 75% survival time#
    ##################
    
    gp[8,i,j] <- KL_lambda(data.test,t3,fit,gp[1,i,j])
    ug[8,i,j] <- KL_lambda(data.test,t3,fit,ug[1,i,j])
    ncv[8,i,j] <- KL_lambda(data.test,t3,fit, ncv[1,i,j])
    dvr[8,i,j] <- KL_lambda(data.test,t3,fit,dvr[1,i,j])
    
    ##################
    # Overall C index#
    ##################
    gp[9,i,j] <- AUC_lambda(data.test,fit,gp[1,i,j])
    ug[9,i,j] <- AUC_lambda(data.test,fit,ug[1,i,j])
    ncv[9,i,j] <- AUC_lambda(data.test,fit, ncv[1,i,j])
    dvr[9,i,j] <- AUC_lambda(data.test,fit,dvr[1,i,j])
  
    print(c(i,j))
  }
}

results <- list()
for(i in 1:9){
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
                    "Brier (25% time)",
                    "Brier (median time)",
                    "Brier (75% time)",
                    "KL (25% time)",
                    "KL (median time)",
                    "KL (75% time",
                    "CIndex")
results

 save(gp,ncv,dvr,ug,oracle.mse,
      file = paste0("sim_high",id,".RData"))
