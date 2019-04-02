#------------------------------
# Self-contained R code
# to reproduce part of the simulations
# in Figure 5 of the Cross Validation 
# Manuscript:
# Explored low dimension p =50
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


#-------------------------------------------------------------
N <- 200
n <- 150
p <- 50
c <- c(0.2,0.3,0.4,0.5,0.6,0.7)
nbeta <- 5
nfold <- 10
foldid <- numeric(n)

beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = nbeta ,c = c[i])
}

# initiation
ncv <- array(NA,c(3,N,length(c)))
auc <- array(NA,c(3,N,length(c)))

oracle.mse <- matrix(NA, nrow = N, ncol = length(c))

for(i in 1:N){
  set.seed(999 + i)
  
  for (j in 1:length(c)){
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
    
    cv.ncv <- try(cv.ncvsurv(X = data$x,
                             y = data$y,
                             fold = foldid,
                             lambda = lambda,
                             penalty = "lasso",
                             returnY = TRUE))
    
    cv.auc <- AUC(cv.ncv)
    
    #oracle model
    data.oracle <- data.frame(data$y, data$x[,1:nbeta])
    fit.oracle <- coxph(data.y ~ ., data = data.oracle)
    oracle.mse[i,j] <- sum((fit.oracle$coef - beta[j,1:nbeta])^2)/p
    
    ###########
    # lambdas #
    ###########
    ncv[1,i,j]<- cv.ncv$lambda.min
    auc[1,i,j]<- lambda[which.max(cv.auc)]
    
    ########
    # mses #
    ########
    ncv[2,i,j] <- norm[ncv[1,i,j] == lambda]
    auc[2,i,j] <- norm[auc[1,i,j] == lambda]
    
    ##############
    # Brier Score#
    ##############
    
    sfit <- survfit(data.test$y ~ 1)
    t <- summary(sfit)$table["median"]
    
    ncv[3,i,j] <- brierSurv_lambda(data.test,t,fit,ncv[1,i,j])
    auc[3,i,j] <- brierSurv_lambda(data.test,t,fit,auc[1,i,j])
    
  }
}

# store results
# lambda
ncv_lambda <- as.numeric(ncv[1,,])
auc_lambda <- as.numeric(auc[1,,])
# MSE ratio
ncv_mse <- log(as.numeric(ncv[2,,]/oracle.mse))
auc_mse <- log(as.numeric(auc[2,,]/oracle.mse))
# Brier
ncv_brier <- as.numeric(ncv[3,,])
auc_brier <- as.numeric(auc[3,,])

coef <- rep(c(rep("0.2",200),
              rep("0.3",200),
              rep("0.4",200),
              rep("0.5",200),
              rep("0.6",200),
              rep("0.7",200)),2)
method <- c(rep("CV-LP",200*6),rep("CV-AUC",200*6))

df <- data.frame(c(ncv_lambda,auc_lambda,
                   ncv_mse,auc_mse,
                   ncv_brier,auc_brier),
                 rep(coef,3),
                 rep(method,3))
colnames(df)<- c("Value","Coef","Method")
df$Type <- c(rep("lambda",2400),rep("log(MSE Ratio)",2400),rep("Brier Score",2400))


#------------------------------
# Explored high dimension p =1000
#------------------------------

#-------------------------------------------------------------
N <- 200
n <- 120
p <- 1000
c <- c(0.3,0.4,0.5,0.6,0.7,0.8)
nbeta <- 10
nfold <- 10
foldid <- numeric(n)

beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = nbeta ,c = c[i])
}

# initiation
ncv <- array(NA,c(3,N,length(c)))
auc <- array(NA,c(3,N,length(c)))

oracle.mse <- matrix(NA, nrow = N, ncol = length(c))

for(i in 1:N){
  set.seed(999 + i)
  
  for (j in 1:length(c)){
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
    
    cv.ncv <- try(cv.ncvsurv(X = data$x,
                             y = data$y,
                             fold = foldid,
                             lambda = lambda,
                             penalty = "lasso",
                             returnY = TRUE))
    
    cv.auc <- AUC(cv.ncv)
    
    #oracle model
    data.oracle <- data.frame(data$y, data$x[,1:nbeta])
    fit.oracle <- coxph(data.y ~ ., data = data.oracle)
    oracle.mse[i,j] <- sum((fit.oracle$coef - beta[j,1:nbeta])^2)/p
    
    ###########
    # lambdas #
    ###########
    ncv[1,i,j]<- cv.ncv$lambda.min
    auc[1,i,j]<- lambda[which.max(cv.auc)]
    
    ########
    # mses #
    ########
    ncv[2,i,j] <- norm[ncv[1,i,j] == lambda]
    auc[2,i,j] <- norm[auc[1,i,j] == lambda]
    
    ##############
    # Brier Score#
    ##############
    
    sfit <- survfit(data.test$y ~ 1)
    t <- summary(sfit)$table["median"]
    
    ncv[3,i,j] <- brierSurv_lambda(data.test,t,fit,ncv[1,i,j])
    auc[3,i,j] <- brierSurv_lambda(data.test,t,fit,auc[1,i,j])
    
  }
}


# store results

# lambda
ncv_lambda <- as.numeric(ncv[1,,])
auc_lambda <- as.numeric(auc[1,,])
# MSE ratio
ncv_mse <- log(as.numeric(ncv[2,,]/oracle.mse))
auc_mse <- log(as.numeric(auc[2,,]/oracle.mse))
# Brier
ncv_brier <- as.numeric(ncv[3,,])
auc_brier <- as.numeric(auc[3,,])

coef <- rep(c(rep("0.3",200),
              rep("0.4",200),
              rep("0.5",200),
              rep("0.6",200),
              rep("0.7",200),
              rep("0.8",200)),2)

df2 <- data.frame(c(ncv_lambda,auc_lambda,
                    ncv_mse,auc_mse,
                    ncv_brier,auc_brier),
                  rep(coef,3),
                  rep(method,3))
colnames(df2)<- c("Value","Coef","Method")
df2$Type <- c(rep("lambda",2400),rep("log(MSE Ratio)",2400),rep("Brier Score",2400))

#-----------------------------------------------------------
# make plots
library(ggplot2)

df_combined <- rbind(df,df2)
df_combined$dim <- c(rep("p = 50",7200),rep("p = 1000",7200))

df_combined$dim <- factor(df_combined$dim,levels = c("p = 50","p = 1000"))
df_combined$Type <- factor(df_combined$Type,levels = c("lambda","log(MSE Ratio)","Brier Score"))

auc_plot <- ggplot(df_combined, aes(Coef, Value,fill=Method)) + 
  geom_boxplot() + 
  facet_grid(Type ~ dim,scales = "free") +
  labs(x = "signal strength (s)") +
  theme_grey()

auc_plot
