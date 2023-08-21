#------------------------------
# Self-contained R code
# to reproduce simulation results
# for the correlated case in Figure 2 
# of the Cross Validation Manuscript
# scenario: increasing signal
# in n = 120, p = 1000
#------------------------------

rerun <- FALSE

if (rerun == FALSE){
  load("cor_exchange.RData")
}else{


library(survival)
library(ncvreg)

#-------------------------------------------------------------------------
# functions
genS <- function(p, rho, corr) {
  if (corr=='exchangeable') {
    S <- matrix(rho, p, p) + (1-rho)*diag(p)
  } else if (corr=='autoregressive') {
    RHO <- matrix(rho^(0:(p-1)), p, p, byrow=TRUE)
    S <- Matrix::bandSparse(p, k=0:(p-1), diagonals=RHO, symmetric=TRUE)
  }
  S
}

genX <- function(n, p, S) {
  R <- chol(S)
  as.matrix(matrix(rnorm(n*p), n, p) %*% R)
}

genSurvData <- function(n=100, p=60, a=6, b=2, rho=0.5, noise=c('exchangeable', 'autoregressive'),
                        rho.noise=0, beta = 0.5,
                        h = 1,
                        cpct = 0.70,# 70% means 70% of the data are censored; 30% are observed
                        corr = c("exchangeable","autoregressive")
){
  noise <- match.arg(noise)
  corr <- match.arg(corr)
  K <- b + 1
  
  # Gen X, S
  sigmaList <- vector('list', a+1)
  if(corr == "exchangeable"){
    for (i in 1:a) {
      sigmaList[[i]] <- matrix(rho, K, K) + (1-rho)*diag(K)
    }
  }
  if(corr == "autoregressive"){
    for (i in 1:a) {
      RHO <- matrix(rho^(0:(K-1)), K, K, byrow=TRUE)
      sigmaList[[i]] <- Matrix::bandSparse(K, k=0:(K-1), diagonals=RHO, symmetric=TRUE)
    }
  }
  sigmaList[[a+1]] <- genS(p-K*a, rho.noise, noise)
  S <- Matrix::.bdiag(sigmaList)
  X <- genX(n, p, S)
  
  # Gen beta
  bb <- beta*(c(-1,1)[(1:a)%%2+1])
  beta <- numeric(p)
  beta[((1:a)-1)*K+1] <- bb
  
  #standardize X first
  x <- ncvreg::std(X)
  hazard <- h*exp(x %*% beta)
  ## Calculate survival times
  time <- rexp(n, hazard)
  ## observations are censored ("event = rbinom(n, 1, 0.30)")
  y <- Surv(time = time, event = rbinom(n, 1, 1-cpct))
  InA <- ((1:a)-1)*K+1
  return(list(x=x,y=y,beta = beta,InA = InA))
}
# test genSurvData
# data <- genSurvData(n=100, p=10, a=2, b=3, rho=0.5, noise='exchangeable',
#                         rho.noise=0, beta = 0.5,
#                         h = 1,
#                         cpct = 0.70,# 70% means 70% of the data are censored; 30% are observed
#                         corr = "exchangeable")
# coxph(data$y~data$x)


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
N <- 500
rho_all <- c(0,0.2,0.4,0.6,0.8)
nbeta <- 10
signal_s <- 0.7
nfold <- 10
b <- 10

# initiation
gp <- array(NA,c(5,N,length(rho_all)))
ug <- array(NA,c(5,N,length(rho_all)))
ncv <- array(NA,c(5,N,length(rho_all)))
dvr <- array(NA,c(5,N,length(rho_all)))
oracle.mse <- matrix(NA, nrow = N, ncol = length(rho_all))
foldid <- numeric(n)


for(i in 1:N){
  set.seed(123456 + 3*i)
  for (j in 1:length(rho_all)){
    # generate data
    data <- genSurvData(n=n, p=p, a=nbeta, b=b, rho=rho_all[j], noise= 'exchangeable', beta = signal_s,
                        h = 1, cpct = 0.10,
                        corr = "exchangeable")
    data.test <- genSurvData(n=n, p=p, a=nbeta, b=b, rho=rho_all[j], noise= 'exchangeable', beta = signal_s,
                        h = 1, cpct = 0.10,
                        corr = "exchangeable")
    nevent <- sum(data$y[,2])
    sde <- sqrt(.Machine$double.eps)
    foldid[data$y[,2] == 1] <-  ceiling(sample(1:nevent)/(nevent + sde) * nfold)
    foldid[data$y[,2] == 0] <-  ceiling(sample(1:(n-nevent))/(n - nevent + sde) * nfold)
    
    # initial fit
    fit <- ncvsurv(data$x, data$y,penalty = "lasso")
    lambda <- fit$lambda
    beta_hat <- fit$beta
    norm <- apply(beta_hat,2,function(x){mean((x-data$beta)^2)})

    #oracle model
    data.oracle <- data.frame(data$y, data$x[,data$InA])
    fit.oracle <- coxph(data.y ~ ., data = data.oracle)
    oracle.mse[i,j] <- sum((fit.oracle$coef -data$beta[data$InA])^2)/p
    
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
    
    #print(c(i,j))
  }
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
coef <- c(0,0.2,0.4,0.6,0.8)
Method <- c(rep("V & VH",5),
            rep("Basic",5),rep("LinearPred",5),
            rep("DevResid",5))
df <- data.frame(lambda,mse,Brier,KL,CIndex,coef,Method)
df$Method <- factor(df$Method,levels = c("V & VH","Basic","LinearPred","DevResid"))


lambda_p <- ggplot(data = df, aes(x = coef,y = lambda,group = Method))+
  geom_line(aes(color=Method),size=1.2,stat = 'identity')+
  geom_point(aes(color=Method),size=3)+ theme_gray()+
  labs(x = "Correlation", 
       y = expression(lambda))+
  theme(legend.position = "none")
mse_p <- ggplot(data = df, aes(x = coef,y = mse,group = Method))+
  geom_line(aes(color=Method),size=1.2)+
  geom_point(aes(color=Method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$Method))+
  labs(x = "Correlation", 
       y = "log(MSE Ratio)") +
  theme(legend.position = "none")
Brier_p <- ggplot(data = df, aes(x = coef,y = Brier,group = Method))+
  geom_line(aes(color=Method),size=1.2)+
  geom_point(aes(color=Method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$Method))+
  labs(x = "Correlation", 
       y = "Brier Score") +
  theme(legend.position = "none")
KL_p <- ggplot(data = df, aes(x = coef,y = KL,group = Method))+
  geom_line(aes(color=Method),size=1.2)+
  geom_point(aes(color=Method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$Method))+
  labs(x = "Correlation", 
       y = "Kullback Leibler")+
  theme(legend.position = "none")
CIndex <- ggplot(data = df, aes(x = coef,y = CIndex,group = Method))+
  geom_line(aes(color=Method),size=1.2)+
  geom_point(aes(color=Method),size=3)+
  theme_gray()+
  scale_colour_discrete(drop=TRUE,limits = levels(df$Method))+
  labs(x = "Correlation" , 
       y = "Harrell's C Index") +
  theme(legend.position = "right")


# png(filename="figure_2_exchangeable.png",
#    units="px",
#    width=1000*1.5,
#    height=350*1.5,
#    pointsize=12,
#    res = 120*1.3)
plot_grid(mse_p,Brier_p,CIndex, ncol = 3,align = "h",axis = "b",rel_widths = c(1, 1, 1.5))
# dev.off()

# tiff(filename="figure_2_new.tiff",
#     units="px",
#     width=1000*1.5,
#     height=350*1.5,
#     pointsize=12,
#     res = 120*1.3)
# plot_grid(mse_p,Brier_p,CIndex, ncol = 3,align = "h",axis = "b",rel_widths = c(1, 1, 1.5))
# dev.off()


