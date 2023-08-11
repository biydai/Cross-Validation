#------------------------
#
# Self contained R code
# to analyze the Shedden Lung Cancer data set
#
#------------------------

#------------------------
# load/define functions that will be used in the analysis
library(survival)
library(ncvreg)
library(glmnet)


# Brier score for selected lambda

brierSurv_lambda <- function(data,time,fit,lambda) {
  predicted <- apply(data$x,1, 
                     function(x)
                     {predict(fit, x,type = "survival",
                              lambda = lambda)(time)})
  is_obs_after <- data$y[,1] > time
  mean((is_obs_after - predicted)^2)
}

#------------------------
# load data

# Shedden Lung Cancer Data
data <- readRDS(url("https://s3.amazonaws.com/pbreheny-data-sets/Shedden2008.rds"))
X <- data$X
S <- data$S
Z <- data$Z

# censor S by 72 month
plot(survfit(S ~ 1), xlab = "month", ylab = "Survival Probability")
abline(h = 0.5, col = "grey", lty = 2)

# clean data
ncolX <- ncol(X)
n <- nrow(X)
Z <- as.data.frame(Z)

data_low <- data.frame(S,Z)
data_low$Age <- as.numeric(data_low$Age)
data_low$Sex[data_low$Sex == "Female"] <- 0
data_low$Sex[data_low$Sex == "Male"] <- 1
data_low$Sex <- as.numeric(data_low$Sex)
data_low$AdjChemo <- as.numeric(data_low$AdjChemo == "Yes")

x <- as.matrix(data.frame(data_low[,c(2,3,5)],X))
x <- ncvreg::std(x)
data <- list(y=S,x=x)

# prespecify fold id
set.seed(123456)
nfold <- 10
nevent <- sum(data$y[,2])
sde <- sqrt(.Machine$double.eps)
foldid <- numeric(0)
foldid[data$y[,2] == 1] <-  ceiling(sample(1:nevent)/(nevent + sde) * nfold)
foldid[data$y[,2] == 0] <-  ceiling(sample(1:(n-nevent))/(n - nevent + sde) * nfold)


Alpha <- c(1,0.9,0.7,0.5,0.3,0.2,0.1,0.05)
A <- length(Alpha)

#------------------------------------------------------------------
for(a in 1:A){
  # initial fit
    alpha <- Alpha[a]
    fit <-  glmnet(data$x, data$y,family = "cox",alpha = alpha,
                 penalty.factor = c(rep(0,3),rep(1,ncolX)))
        #ncvsurv(data$x, data$y,penalty = "lasso", alpha = 0.7,
        #         penalty.factor = c(rep(0,3),rep(1,ncolX)))

    lambda <- fit$lambda
    beta_hat <- fit$beta
    
    eta.cv <- matrix(nrow = n, ncol = length(lambda))
    l_diff <- matrix(nrow = nfold , ncol = length(lambda))
    l_test <- matrix(nrow = nfold, ncol = length(lambda))
    
    ind <- order(data$y[,1])
    d <- data$y[ind,2]

# cross validation

    for(f in 1:nfold){
      ifold <- foldid == f
      cvfit <- glmnet(data$x[!ifold,], data$y[!ifold,],family = "cox",alpha = alpha,
                      penalty.factor = c(rep(0,3),rep(1,ncolX)),
                      lambda = lambda,intercept = FALSE)
      cvlambda <- cvfit$lambda
      beta_train <- cvfit$beta
      
      # linear predictors
      eta.cv[ifold,lambda%in%cvlambda] <- as.matrix(data$x[ifold,] %*% beta_train)
      
      # grouped
      eta.full.ordered <- apply(beta_train,2,function(b){
        data$x %*% b})[ind,]
      r.full.ordered <- apply(eta.full.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
      l_full <- -2 * (crossprod(d, eta.full.ordered) - crossprod(d, log(r.full.ordered)))
      
      ind_train <- order(data$y[!ifold,1])
      d_train <- data$y[!ifold,2][ind_train]
      eta.train.ordered <- apply(beta_train,2,function(b){
        data$x[!ifold,] %*% b})[ind_train,]
      r.train.ordered <- apply(eta.train.ordered, 2, 
                               function(x) rev(cumsum(rev(exp(x)))))
      l_train <- -2 * (crossprod(d_train, eta.train.ordered) - crossprod(d_train, log(r.train.ordered)))
      
      l_diff[f,lambda%in%cvlambda] <- l_full - l_train
      
      # ungrouped
      ind_test <- order(data$y[ifold,1])
      d_test <- data$y[ifold,2][ind_test]
      eta.test.ordered <- apply(beta_train,2,function(b){
        e <- data$x[ifold,] %*% b})[ind_test,]
      r.test.ordered <- apply(eta.test.ordered, 2, 
                              function(x) rev(cumsum(rev(exp(x)))))
      l_test[f,lambda%in%cvlambda] <- -2 * (crossprod(d_test, eta.test.ordered) - crossprod(d_test, log(r.test.ordered)))
    }

  # pool cross-validation results together
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
    # this one changed because of ties
    KM <- survfit(data$y~1)
    BH_tied <- cumsum(KM$n.event/KM$n.risk)
    BH_tied_time <- data.frame(KM$time,BH_tied)
    time <- data.frame(data$y[,1])
    BH_NA <- rep(NA,n)
    BH_NA[!duplicated(time)] <- BH_tied
    for(i in 1:n){
      if(is.na(BH_NA[i])){
        BH_NA[i] <- BH_NA[i-1]
      }
    }

    w.cv <- exp(eta.cv[ind, , drop = FALSE])
    val_NA <- apply(w.cv, 2, function(x){ 
    M <- d - BH_NA*x
    dev <- sign(M) * sqrt(-2 * (M + ifelse(d == 0, 0, d * 
                                             log(d - M))))
    sum(dev^2)})
  
  # put results into tables:
    lambda_min <- c(lambda[which.min(val_GP)],
                  lambda[which.min(val_UG)],
                  lambda[which.min(val_LP)],
                  lambda[which.min(val_NA)])
    log_lambda <- log(lambda_min)
  
  # compute AUC
    AUC_all <- apply(eta.cv,2, function(lp){survConcordance(data$y ~ lp)$concordance})
    AUC <- c(AUC_all[which.min(val_GP)],
           AUC_all[which.min(val_UG)],
           AUC_all[which.min(val_LP)],
           AUC_all[which.min(val_NA)])
  
  # use the actual survival prediction to compute the brier score
    ncvfit <- ncvsurv(X=data$x,y = data$y, penalty = "lasso",alpha = alpha,lambda= lambda,
                   penalty.factor = c(rep(0,3),rep(1,ncolX)))
    
    brier <- c(brierSurv_lambda(data = data, 
                              time =70.5, fit = ncvfit,
                              lambda = lambda_min[1]),
               brierSurv_lambda(data = data, 
                              time =70.5, fit = ncvfit,
                              lambda = lambda_min[2]),
               brierSurv_lambda(data = data, 
                              time =70.5, fit = ncvfit,
                              lambda = lambda_min[3]),
               brierSurv_lambda(data = data, 
                              time =70.5, fit = ncvfit,
                              lambda = lambda_min[4]))

  result <- data.frame(lambda_min,log_lambda,brier,AUC)
  result
  save(result,ncvfit,fit,lambda,val_GP,val_UG,val_LP,val_NA, file = paste0("elasticnet_result_",alpha,".RData"))

  range <- function(x){(x-min(x))/(max(x)-min(x))}

}

# make tables 


Alpha <- c(0.05,0.1,0.2,0.3,0.5,0.7,0.9,1)

C_Index <- matrix(NA, nrow = 4, ncol = length(Alpha))
Brier <- matrix(NA, nrow = 4, ncol =length(Alpha))

colnames(C_Index) <- Alpha
rownames(C_Index) <- c("GP", "UG","LP","DR")

colnames(Brier) <- Alpha
rownames(Brier) <- c("GP", "UG","LP","DR")


load("elasticnet_result_0.05.RData")
C_Index[,1] <- result[,4]
Brier[,1] <- result[,3]


load("elasticnet_result_0.1.RData")
C_Index[,2] <- result[,4]
Brier[,2] <- result[,3]

load("elasticnet_result_0.2.RData")
C_Index[,3] <- result[,4]
Brier[,3] <- result[,3]


load("elasticnet_result_0.3.RData")
C_Index[,4] <- result[,4]
Brier[,4] <- result[,3]

load("elasticnet_result_0.5.RData")
C_Index[,5] <- result[,4]
Brier[,5] <- result[,3]

load("elasticnet_result_0.7.RData")
C_Index[,6] <- result[,4]
Brier[,6] <- result[,3]

load("elasticnet_result_0.9.RData")
C_Index[,7] <- result[,4]
Brier[,7] <- result[,3]

load("elasticnet_result_1.RData")
C_Index[,8] <- result[,4]
Brier[,8] <- result[,3]

C_Index
Brier

# make figure
load("elasticnet_result_0.3.RData")
range <- function(x){(x-min(x))/(max(x)-min(x))}

ID <- 1:25
library(ggplot2)
library(dplyr)
df <- rbind(
  data.frame(lambda = log(lambda[ID]), CVE = range(val_GP[ID]), Method = rep("V VH",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_UG[ID]), Method = rep("Standard",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_LP[ID]), Method = rep("Linear Pred",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_NA[ID]), Method = rep("Dev Resid", length(ID)))
)
min <- group_by(df,Method)%>%summarise(lambdamin = lambda[which.min(CVE)])

 png(filename="shedden.png", 
     units="px", 
     width=1500, 
     height=1000, 
     pointsize=15,
     res = 120*2)
 ggplot(data = df, mapping = aes(x = lambda,y = CVE, group = Method)) + 
   geom_line(aes(color = Method),size = 1.25)+
   scale_x_reverse()+
   labs(x = paste0("log(", '\u03BB',")"))+
   geom_vline(data = min, aes(xintercept = lambdamin,colour = Method),linetype = 2,size = 0.75)
dev.off()

# numbder of genes selected by the model
GP <- ncvfit$beta[,which.min(val_GP)]!=0
LP <- ncvfit$beta[,which.min(val_LP)]!=0
UG <- ncvfit$beta[,which.min(val_UG)]!=0
DV <- ncvfit$beta[,which.min(val_NA)]!=0

sum(GP);sum(LP);sum(UG);sum(DV)

# make figure
load("elasticnet_result_0.05.RData")
range <- function(x){(x-min(x))/(max(x)-min(x))}

ID <- 1:25
library(ggplot2)
library(dplyr)
df <- rbind(
  data.frame(lambda = log(lambda[ID]), CVE = range(val_GP[ID]), Method = rep("V VH",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_UG[ID]), Method = rep("Standard",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_LP[ID]), Method = rep("Linear Pred",length(ID))),
  data.frame(lambda = log(lambda[ID]), CVE = range(val_NA[ID]), Method = rep("Dev Resid", length(ID)))
)
min <- group_by(df,Method)%>%summarise(lambdamin = lambda[which.min(CVE)])

png(filename="shedden.png", 
    units="px", 
    width=1500, 
    height=1000, 
    pointsize=15,
    res = 120*2)
ggplot(data = df, mapping = aes(x = lambda,y = CVE, group = Method)) + 
  geom_line(aes(color = Method),size = 1.25)+
  scale_x_reverse()+
  labs(x = paste0("log(", '\u03BB',")"))+
  geom_vline(data = min, aes(xintercept = lambdamin,colour = Method),linetype = 2,size = 0.75)
dev.off()

# numbder of genes selected by the model
GP <- ncvfit$beta[,which.min(val_GP)]!=0
LP <- ncvfit$beta[,which.min(val_LP)]!=0
UG <- ncvfit$beta[,which.min(val_UG)]!=0
DV <- ncvfit$beta[,which.min(val_NA)]!=0

sum(GP);sum(LP);sum(UG);sum(DV)
