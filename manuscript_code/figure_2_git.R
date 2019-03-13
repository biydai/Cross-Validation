#------------------------------
# Self-contained R code
# to reproduce simulation results
# in Figure 2 of the Cross Validation 
# Manuscript
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
n <- 120
p <- 1000
N <- 200
c <- c(0.4,0.5,0.6,0.7,0.8,0.9)
nbeta <- 10
nfold <- 10
beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = nbeta ,c = c[i])
}

# initiation
gp <- array(NA,c(5,N,length(c)))
ug <- array(NA,c(5,N,length(c)))
ncv <- array(NA,c(5,N,length(c)))
dvr <- array(NA,c(5,N,length(c)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(c))
foldid <- numeric(n)


for(i in 1:N){
  set.seed(8888 + i)
  for (j in 1:length(c)){
    # generate data
    data <- genSurvData(n=n, beta = beta[j,], cpct = 0.1, sd = 1, h = 1)
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
    
    ##############
    # Brier Score#
    ##############
    
    sfit <- survfit(data.test$y ~ 1)
    t <- summary(sfit)$table["median"]
    
    gp[3,i,j] <- brierSurv_lambda(data.test,t,fit,gp[1,i,j])
    ug[3,i,j] <- brierSurv_lambda(data.test,t,fit,ug[1,i,j])
    ncv[3,i,j] <- brierSurv_lambda(data.test,t,fit, ncv[1,i,j])
    dvr[3,i,j] <- brierSurv_lambda(data.test,t,fit,dvr[1,i,j])
    
    ###################
    # Kullbeck Leibler#
    ##################
    
    gp[4,i,j] <- KL_lambda(data.test,t,fit,gp[1,i,j])
    ug[4,i,j] <- KL_lambda(data.test,t,fit,ug[1,i,j])
    ncv[4,i,j] <- KL_lambda(data.test,t,fit, ncv[1,i,j])
    dvr[4,i,j] <- KL_lambda(data.test,t,fit,dvr[1,i,j])
    
    ##################
    # Overall C index#
    ##################
    gp[5,i,j] <- AUC_lambda(data.test,fit,gp[1,i,j])
    ug[5,i,j] <- AUC_lambda(data.test,fit,ug[1,i,j])
    ncv[5,i,j] <- AUC_lambda(data.test,fit, ncv[1,i,j])
    dvr[5,i,j] <- AUC_lambda(data.test,fit,dvr[1,i,j])
  }
}

#-------------------------------------------------------------------------
# making Figure 2
library(ggplot2)
library(cowplot)

results <- list()
for(i in 1:5){
  results[[i]] <- c(
    apply(gp[i,,],2,mean),
    apply(ug[i,,],2,mean),
    apply(ncv[i,,],2,mean),
    apply(dvr[i,,],2,mean)
  )
}
names(results) <- c("lambda",
                    "MSE",
                    "Brier",
                    "KL",
                    "CIndex")
lambda <- results$lambda

mse <- log(c(apply(gp[2,,]/oracle.mse,2,mean),
             apply(ug[2,,]/oracle.mse,2,mean),
             apply(ncv[2,,]/oracle.mse,2,mean),
             apply(dvr[2,,]/oracle.mse,2,mean)))

Brier <- (results$Brier)
KL <- (results$KL)
CIndex <- results$CIndex
coef <- rep((4:9)/10,4)
method <- c(rep("V & VH",6),
            rep("Standard",6),rep("LinearPred",6),
            rep("DevResid",6))
df <- data.frame(lambda,mse,Brier,KL,CIndex,coef,method)
df$method <- factor(df$method,levels = c("V & VH","Standard","LinearPred","DevResid"))


lambda_p <- ggplot(data = df, aes(x = coef,y = lambda,group = method))+
  geom_line(aes(color=method),size=1.2,stat = 'identity')+
  geom_point(aes(color=method),size=3)+ theme_gray()+
  labs(x = "s", 
       y = expression(lambda))+
  theme(legend.position = "none")
mse_p <- ggplot(data = df, aes(x = coef,y = mse,group = method))+
  geom_line(aes(color=method),size=1.2)+
  geom_point(aes(color=method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$method))+
  labs(x = "s", 
       y = "log(MSE Ratio)") +
  theme(legend.position = "none")
Brier_p <- ggplot(data = df, aes(x = coef,y = Brier,group = method))+
  geom_line(aes(color=method),size=1.2)+
  geom_point(aes(color=method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$method))+
  labs(x = "s", 
       y = "Brier Score") +
  theme(legend.position = "none")
KL_p <- ggplot(data = df, aes(x = coef,y = KL,group = method))+
  geom_line(aes(color=method),size=1.2)+
  geom_point(aes(color=method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$method))+
  labs(x = "s", 
       y = "Kullback Leibler")+
  theme(legend.position = "none")
CIndex <- ggplot(data = df, aes(x = coef,y = CIndex,group = method))+
  geom_line(aes(color=method),size=1.2)+
  geom_point(aes(color=method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$method))+
  labs(x = "s" , 
       y = "Harrell's C Index") +
  theme(legend.position = "right")

plot_grid(mse_p,Brier_p,CIndex, ncol = 3,align = "h",axis = "b",rel_widths = c(1, 1, 1.5))



