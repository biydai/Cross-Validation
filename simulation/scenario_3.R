##############################################
# compare the stability of different methods #
##############################################

# ungrouped
# grouped
# linear predictor
# deviance residuals

# 3.1 vary censoring %
# 3.2 vary fold numbers


# load source files
source("survival_functions.R")

# load modified ungrouped method
source("modify_fold.R")

# load deviance residuals
source("deviance.R")



#######################
### vary censoring ####
#######################

n <- 100
p <- 100
cpct <- seq(0,0.5,by = 0.1)
beta <- generatebeta(p,nbeta = 10 ,c = 1)
N <- 200
# initiation
true.lambda.cen<- matrix(NA,nrow = N ,ncol = length(cpct))
gp.lambda.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.lambda.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.lambda.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.lambda.cen <- matrix(NA,nrow = N ,ncol = length(cpct))

true.mse.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.mse.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.mse.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.mse.cen <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.mse.cen <- matrix(NA,nrow = N ,ncol = length(cpct))


for(i in 1:N){
  set.seed(666 + i)
  foldid <- sample(1:10,size=n,replace=TRUE)
  for (j in 1:length(cpct)){
    data <- genSurvData(n=100 , beta = beta, cpct = cpct[j])
    f <- try(survncvreg(data$x,data$y,data$beta))
    cv.gp <- try(cv.glmnet(x = data$x,y = data$y,nfolds = 10,foldid=foldid,lambda = f$lambda,family = "cox"))
    cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,foldid=foldid,lambda = f$lambda,grouped = FALSE,family = "cox"))
    cv.ncvreg <- try(cv.ncvsurv(X = data$x,y = data$y,nfolds = 10,foldid=foldid,lambda = f$lambda,penalty = "lasso"))
    cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,foldid = foldid,lambda = f$lambda,penalty = "lasso"))
    # true
    true.lambda.cen[i,j] <- f$lambda_t
    true.mse.cen[i,j] <- min(f$norm)
    
    # gp
    if (typeof(cv.gp) == "character"){
      gp.lambda.cen[i,j] <- NA
      gp.mse.cen[i,j] <- NA
    }else{
      gp.lambda.cen[i,j]<- cv.gp$lambda.min
      gp.mse.cen[i,j] <- f$norm[cv.gp$lambda.min == f$lambda]
    }
    
    # ug
    if (typeof(cv.ug) == "character"){
      ug.lambda.cen[i,j] <- NA
      ug.mse.cen[i,j] <- NA
    }else{
      ug.lambda.cen[i,j]<- cv.ug$lambda.min
      ug.mse.cen[i,j] <- f$norm[cv.ug$lambda.min == f$lambda]
    }
    
    # ncv
    if (typeof(cv.ncvreg) == "character"){
      ncv.lambda.cen[i,j] <- NA
      ncv.mse.cen[i,j] <- NA
    }else{
      ncv.lambda.cen[i,j] <- cv.ncvreg$lambda.min
      ncv.mse.cen[i,j] <- f$norm[cv.ncvreg$lambda.min == f$lambda]
    }   
    
    # deviance  #ncv
    if (typeof(cv.dvr) == "character"){
      dvr.lambda.cen[i,j] <- NA
      dvr.mse.cen[i,j] <- NA
    }else if(length(cv.dvr$lambda.min) == 0){
      dvr.lambda.cen[i,j] <- NA
      dvr.mse.cen[i,j] <- NA
    }else{
      dvr.lambda.cen[i,j] <- try(cv.dvr$lambda.min)
      dvr.mse.cen[i,j] <- f$norm[cv.dvr$lambda.min == f$lambda]
    }  
  }
}

save(true.lambda.cen, true.mse.cen,
     gp.lambda.cen, gp.mse.cen,
     ug.lambda.cen, ug.mse.cen,
     ncv.lambda.cen, ncv.mse.cen,
     dvr.lambda.cen, dvr.mse.cen, file = "scenario_3_1.RData")


### vary number of folds
# need to run modify_fold to change the automatic correct in ungrouped
n <- 100
p <- 100
sd <- 1
nbeta <- 10
N <- 200
nfold <- c(5,10,15,20)


cv.gp <- NULL
cv.ug <- NULL
cv.ncvreg <- NULL
cv.dvr <- NULL

# generate beta:
library(dplyr)
set.seed(888)
beta <- generatebeta(p,nbeta)

true.lambda.fold<- matrix(NA,nrow = N ,ncol = length(nfold))
gp.lambda.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
ug.lambda.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
ncv.lambda.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
dvr.lambda.fold <- matrix(NA,nrow = N ,ncol = length(nfold))

true.mse.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
gp.mse.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
ug.mse.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
ncv.mse.fold <- matrix(NA,nrow = N ,ncol = length(nfold))
dvr.mse.fold <- matrix(NA,nrow = N ,ncol = length(nfold))

for(j in 1:length(nfold)){
  #foldid <- sample(1:nfold[j],size=n,replace=TRUE)
  for (i in 1:N){
    set.seed(666 + i)
    data <- genSurvData(n=100 , beta = beta, cpct = 0.2)
    f <- try(survncvreg(data$x,data$y,data$beta))
    cv.gp <- try(cv.glmnet(x = data$x,y = data$y,nfolds = nfold[j],lambda = f$lambda,family = "cox"))#foldid=foldid,
    cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = nfold[j],lambda = f$lambda,grouped = FALSE,family = "cox"))#foldid=foldid
    cv.ncvreg <- try(cv.ncvsurv(X = data$x,y = data$y,nfolds = nfold[j],lambda = f$lambda,penalty = "lasso"))#foldid=foldid
    cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = nfold[j],lambda = f$lambda,penalty = "lasso"))#foldid = foldid
    # true
    true.lambda.fold[i,j] <- f$lambda_t
    true.mse.cen[i,j] <- min(f$norm)
    
    # gp
    if (typeof(cv.gp) == "character"){
      gp.lambda.fold[i,j] <- NA
      gp.mse.fold[i,j] <- NA
    }else{
      gp.lambda.fold[i,j]<- cv.gp$lambda.min
      gp.mse.fold[i,j] <- f$norm[cv.gp$lambda.min == f$lambda]
    }
    
    # ug
    if (typeof(cv.ug) == "character"){
      ug.lambda.fold[i,j] <- NA
      ug.mse.fold[i,j] <- NA
    }else{
      ug.lambda.fold[i,j]<- cv.ug$lambda.min
      ug.mse.fold[i,j] <- f$norm[cv.ug$lambda.min == f$lambda]
    }
    
    # ncv
    if (typeof(cv.ncvreg) == "character"){
      ncv.lambda.fold[i,j] <- NA
      ncv.mse.fold[i,j] <- NA
    }else{
      ncv.lambda.fold[i,j] <- cv.ncvreg$lambda.min
      ncv.mse.fold[i,j] <- f$norm[cv.ncvreg$lambda.min == f$lambda]
    }   
    
    # deviance  #ncv
    if (typeof(cv.dvr) == "character"){
      dvr.lambda.fold[i,j] <- NA
      dvr.mse.fold[i,j] <- NA
    }else if(length(cv.dvr$lambda.min) == 0){
      dvr.lambda.fold[i,j] <- NA
      dvr.mse.fold[i,j] <- NA
    }else{
      dvr.lambda.fold[i,j] <- try(cv.dvr$lambda.min)
      dvr.mse.fold[i,j] <- f$norm[cv.dvr$lambda.min == f$lambda]
    }  
  
  }
}

save(true.lambda.fold, true.mse.fold,
     gp.lambda.fold, gp.mse.fold,
     ug.lambda.fold, ug.mse.fold,
     ncv.lambda.fold, ncv.mse.fold,
     dvr.lambda.fold, dvr.mse.fold, file = "scenario_3_2.RData")

################
# plot results #
################

cen <- rbind(
  apply(ncv.lambda.cen,2,function(x){sum(is.na(x))}),
  apply(gp.lambda.cen,2,function(x){sum(is.na(x))}),
  apply(ug.lambda.cen,2,function(x){sum(is.na(x))})
  #,apply(dvr.lambda.cen,2,function(x){sum(is.na(x))})
)
cen <- (200-cen)/200

fold <- rbind(
  apply(ncv.lambda.fold,2,function(x){sum(is.na(x))}),
  apply(gp.lambda.fold,2,function(x){sum(is.na(x))}),
  apply(ug.lambda.fold,2,function(x){sum(is.na(x))})
  #,apply(dvr.lambda.fold,2,function(x){sum(is.na(x))})
)
fold <- (200-fold)/200

#par(mfrow = c(1,2))
png(filename= "sim_3_1.png")
cpct <-seq(0,0.5,by = 0.1)
plot(x = cpct, y = cen[1,], type = "l", ylim = c(0.5,1), lwd = 2,
     ylab = "Probability of Defined CVE", xlab = "Censored Portion")
lines(x = cpct, y = cen[2,], col = 1, lwd = 2)
lines(x = cpct, y = cen[3,], col = 2, lwd = 2)
legend("bottomleft",legend = c("Linear Predictor","Grouped","Ungrouped"), col = c(1,1,2), lwd = 2)
dev.off()

png(filename= "sim_3_2.png")
nfold <- c(5,10,15,20)
plot(x = nfold, y = fold[1,], type = "l", ylim = c(0.5,1), lwd = 2,
     ylab = "Convergence %", xlab = "Number of Folds")
lines(x = nfold, y = fold[2,], col = 1, lwd = 2)
lines(x = nfold, y = fold[3,], col = 2, lwd = 2)
dev.off()

