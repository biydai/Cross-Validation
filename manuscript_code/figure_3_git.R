#------------------------------
# Self-contained R code
# to reproduce simulations
# in Figure 3 of the Cross Validation 
# Manuscript:
# Explore Stability of the Standard
# CV Method
#------------------------------

# functions
library(survival)
library(ncvreg)

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
  
  ## observations are censored
  id <- sample(n,round(cpct * n))
  event <- rep(1,n)
  event[id] <- 0
  
  y <- Surv(time = time, event = event)
  return(list(x=x,y=y,beta = beta))
}



# increase number of folds (plot on the right)
#----------------------------------------------------------------------------------------

n <- 120
p <- 1000
cpct <- 0.5
nbeta <- 10
beta <- generatebeta(p,nbeta = nbeta ,c = 0.5)
N <- 200
Nfold <- c(10,12,15,20,30,60)

# initiation
gp <- array(NA,c(2,N,length(Nfold)))
ug <- array(NA,c(2,N,length(Nfold)))
ncv <- array(NA,c(2,N,length(Nfold)))
dvr <- array(NA,c(2,N,length(Nfold)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(Nfold))
foldid <- numeric(n)

for(i in 1:N){
  set.seed(666 + i)
  
  for (j in 1:length(Nfold)){
    nfold <- Nfold[j]
    # generate data
    data <- genSurvData(n=n, beta = beta, cpct = cpct, sd = 1, h = 1)
    data.test <- genSurvData(n=1000 , beta = beta, cpct = 0, sd = 1, h = 1)
    nevent <- sum(data$y[,2])
    sde <- sqrt(.Machine$double.eps)
    foldid[data$y[,2] == 1] <-  ceiling(sample(1:nevent)/(nevent + sde) * nfold)
    foldid[data$y[,2] == 0] <-  ceiling(sample(1:(n-nevent))/(n - nevent + sde) * nfold)
    
    # initial fit
    fit <- ncvsurv(data$x, data$y,penalty = "lasso")
    lambda <- fit$lambda
    beta_hat <- fit$beta
    norm <- apply(beta_hat,2,function(x){mean((x-beta)^2)})
    
    #oracle model
    data.oracle <- data.frame(data$y, data$x[,1:nbeta])
    fit.oracle <- coxph(data.y ~ ., data = data.oracle)
    oracle.mse[i,j] <- sum((fit.oracle$coef - beta[1:nbeta])^2)/p
    
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
  
    print(c(i,j))
  }
}

results <- list()
for(i in 1:2){
  results[[i]] <- rbind(
    apply(gp[i,,],2,mean),
    apply(ncv[i,,],2,mean),
    apply(ug[i,,],2,mean),
    apply(dvr[i,,],2,mean)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- Nfold
}

for(i in 3:4){
  results[[i]] <- rbind(
    apply(gp[i-2,,],2,sd),
    apply(ncv[i-2,,],2,sd),
    apply(ug[i-2,,],2,sd),
    apply(dvr[i-2,,],2,sd)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- Nfold
}

names(results) <- c("lambda",
                    "MSE",
                    "lambda_sd",
                    "MSE_sd")
results$lambda
results$MSE
results$lambda_sd
results$MSE_sd

# making plots

lambda_sd <- c(apply(gp[1,,],2,sd),apply(ug[1,,],2,sd),apply(ncv[1,,],2,sd),apply(dvr[1,,],2,sd))
MSE_sd <- c(apply(gp[2,,]*120,2,sd),apply(ug[2,,]*120,2,sd),apply(ncv[2,,]*120,2,sd),apply(dvr[2,,]*120,2,sd))
method <- c(rep("V & VH",6),rep("Standard",6),rep("LinearPred",6),rep("DevResid",6))
Nfold <- rep(c(6,5,4,3,2,1),4)

df <- data.frame(lambda_sd,MSE_sd,method,Nfold)
df$method <- factor(df$method,levels = c("V & VH","Standard","LinearPred","DevResid"))

MSE_p <- ggplot(data = df, aes(x = Nfold,y = MSE_sd,group = method))+
  geom_line(aes(color=method),size=1.2,stat = 'identity')+
  geom_point(aes(color=method),size=3)+ theme_gray()+
  scale_x_reverse(breaks=c(6,5,4,3,2,1))+
  
  labs(x = "Number of Events Per Fold", 
       y = "SD(Squared Error Loss)")+
  theme(legend.position = "right")
MSE_p

# simulation on imbalance (plot on the left)
#--------------------------------------------------------------------------------------------------------
n <- 100
Cpct <- c(0.4,0.6,0.5,0.7,0.8)

N <- 1000
result <- matrix(NA,nrow = N, ncol = length(Cpct))
set.seed(777)
for(i in 1:length(Cpct)){
  nfold <- 10
  cpct <- Cpct[i]
  for(r in 1:N){
    id <- sample(n,round(cpct * n))
    event <- rep(1,n)
    event[id] <- 0
    sde <- sqrt(.Machine$double.eps)
    foldid <- ceiling(sample(1:n)/(n+sde) * nfold)
    result[r,i]<- prod((1:nfold) %in% foldid[event == 1])>=1
  }
}

balance <- (1 - apply(result,2,mean))*100
df <- data.frame(balance,Cpct)

balance_p <- ggplot(data = df, aes(x = Cpct,y = balance))+
  geom_line(color="#7CAE00",size=1.2,stat = 'identity')+
  geom_point(color="#7CAE00",size=3)+ theme_gray()+
  scale_x_continuous(breaks=c(0.4,0.6,0.5,0.7,0.8))+
  labs(x = "Censoring", 
       y = "% of Undefined PL")+
  theme(legend.position = "right")
balance_p

# combine two plots and output
#----------------------------------------------------------------------------------------------------------
png(filename="figure_3.png",
    units="px",
    width=1000*1.3,
    height=350*1.5,
    pointsize=12,
    res = 120*1.3)
plot_grid(balance_p,MSE_p, ncol = 2,align = "h",axis = "b",rel_widths = c(1, 1.5))
dev.off()
