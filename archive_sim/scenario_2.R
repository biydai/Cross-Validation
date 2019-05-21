# Scenario 1 (plots the lambda selection and MSE ratio)
# compare four cross validation methods: 
#  ungrouped, grouped, linear predictor, deviance residuals
# varying dimension (n = 100, p = 100, nonzero = 10)

# load source files
source("survival_functions.R")

# load modified ungrouped method
source("modify_fold.R")

# load deviance residuals
source("deviance.R")

#############################################################################
# scenario 2_1 : n = 100, p = 20, nbeta = 10, censoring = 10% , nfold = 10
#############################################################################
# change the magnitude of betas
n <- 100
p <- 20
#cpct <- seq(0,0.3,by = 0.1)
cpct <- 0.1
beta <- generatebeta(p,nbeta = 10 ,c = 1)
N <- 200

# initiation
true.lambda.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.lambda.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.lambda.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.lambda.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.lambda.1 <- matrix(NA,nrow = N ,ncol = length(cpct))

true.mse.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.mse.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.mse.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.mse.1 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.mse.1 <- matrix(NA,nrow = N ,ncol = length(cpct))

for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(1:10,size=n,replace=TRUE)
  data <- genSurvData(n=n , beta = beta, cpct = cpct,sd = 1)
  data$x <- unclass(ncvreg::std(data$x))
  
  f <- try(survncvreg(data$x,data$y,data$beta))
  cv.gp <- try(cv.glmnet(x = data$x,y = data$y,foldid=foldid,lambda = f$lambda,family = "cox"))
  cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,lambda = f$lambda,grouped = FALSE,family = "cox"))
  cv.ncv <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
  cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,lambda = f$lambda,penalty = "lasso")) 
  
  # lambdas
  true.lambda.1[i] <- f$lambda_t
  gp.lambda.1[i]<-cv.gp$lambda.min
  ncv.lambda.1[i]<- cv.ncv$lambda.min
  dvr.lambda.1[i]<- cv.dvr$lambda.min
  
  # mses
  true.mse.1[i] <- min(f$norm)
  gp.mse.1[i] <- f$norm[cv.gp$lambda.min == f$lambda]
  ncv.mse.1[i] <- f$norm[cv.ncv$lambda.min == f$lambda]
  dvr.mse.1[i] <- f$norm[cv.dvr$lambda.min == f$lambda]
  
  if (typeof(cv.ug) == "character"){
    ug.lambda.1[i] <- NA
    ug.mse.1[i] <- NA
  }else{
    ug.lambda.1[i]<- cv.ug$lambda.min
    ug.mse.1[i] <- f$norm[cv.ug$lambda.min == f$lambda]
  }  
  
}

mean(true.lambda.1)
mean(gp.lambda.1)
mean(ug.lambda.1)
mean(ncv.lambda.1)
mean(dvr.lambda.1)

#save(true.lambda.1,gp.lambda.1,ug.lambda.1,ncv.lambda.1,dvr.lambda.1,
#     true.mse.1,gp.mse.1,ug.mse.1,ncv.mse.1,dvr.mse.1,
#     file = "sim_2_1.RData")

#############################################################################
# scenario 2_2 : n = 100, p = 100, nbeta = 10, censoring = 10% , nfold = 10
#############################################################################
# change the magnitude of betas
n <- 100
p <- 100
#cpct <- seq(0,0.3,by = 0.1)
cpct <- 0.1
beta <- generatebeta(p,nbeta = 10 ,c = 1)
N <- 200

# initiation
true.lambda.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.lambda.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.lambda.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.lambda.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.lambda.2 <- matrix(NA,nrow = N ,ncol = length(cpct))

true.mse.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.mse.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.mse.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.mse.2 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.mse.2 <- matrix(NA,nrow = N ,ncol = length(cpct))

for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(1:10,size=n,replace=TRUE)
  data <- genSurvData(n=n , beta = beta, cpct = cpct,sd = 1)
  data$x <- unclass(ncvreg::std(data$x))
  
  f <- try(survncvreg(data$x,data$y,data$beta))
  cv.gp <- try(cv.glmnet(x = data$x,y = data$y,foldid=foldid,lambda = f$lambda,family = "cox"))
  cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,lambda = f$lambda,grouped = FALSE,family = "cox"))
  cv.ncv <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
  cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,lambda = f$lambda,penalty = "lasso")) 
  
  # lambdas
  true.lambda.2[i] <- f$lambda_t
  gp.lambda.2[i]<-cv.gp$lambda.min
  ncv.lambda.2[i]<- cv.ncv$lambda.min
  dvr.lambda.2[i]<- cv.dvr$lambda.min
  
  # mses
  true.mse.2[i] <- min(f$norm)
  gp.mse.2[i] <- f$norm[cv.gp$lambda.min == f$lambda]
  ncv.mse.2[i] <- f$norm[cv.ncv$lambda.min == f$lambda]
  dvr.mse.2[i] <- f$norm[cv.dvr$lambda.min == f$lambda]
  
  if (typeof(cv.ug) == "character"){
    ug.lambda.2[i] <- NA
    ug.mse.2[i] <- NA
  }else{
    ug.lambda.2[i]<- cv.ug$lambda.min
    ug.mse.2[i] <- f$norm[cv.ug$lambda.min == f$lambda]
  }  
  
}

mean(true.lambda.2)
mean(gp.lambda.2)
mean(ug.lambda.2)
mean(ncv.lambda.2)
mean(dvr.lambda.2)

#save(true.lambda.2,gp.lambda.2,ug.lambda.2,ncv.lambda.2,dvr.lambda.2,
#     true.mse.2,gp.mse.2,ug.mse.2,ncv.mse.2,dvr.mse.2,
#     file = "sim_2_2.RData")

#############################################################################
# scenario 2_3 : n = 100, p = 1000, nbeta = 10, censoring = 10% , nfold = 10
#############################################################################
# change the magnitude of betas
n <- 100
p <- 1000
#cpct <- seq(0,0.3,by = 0.1)
cpct <- 0.1
beta <- generatebeta(p,nbeta = 10 ,c = 1)
N <- 200

# initiation
true.lambda.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.lambda.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.lambda.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.lambda.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.lambda.3 <- matrix(NA,nrow = N ,ncol = length(cpct))

true.mse.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.mse.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.mse.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.mse.3 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.mse.3 <- matrix(NA,nrow = N ,ncol = length(cpct))

for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(1:10,size=n,replace=TRUE)
  data <- genSurvData(n=n , beta = beta, cpct = cpct,sd = 1)
  data$x <- unclass(ncvreg::std(data$x))
  
  f <- try(survncvreg(data$x,data$y,data$beta))
  cv.gp <- try(cv.glmnet(x = data$x,y = data$y,foldid=foldid,lambda = f$lambda,family = "cox"))
  cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,lambda = f$lambda,grouped = FALSE,family = "cox"))
  cv.ncv <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
  cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,lambda = f$lambda,penalty = "lasso")) 
  
  # lambdas
  true.lambda.3[i] <- f$lambda_t
  gp.lambda.3[i]<-cv.gp$lambda.min
  ncv.lambda.3[i]<- cv.ncv$lambda.min
  dvr.lambda.3[i]<- cv.dvr$lambda.min
  
  # mses
  true.mse.3[i] <- min(f$norm)
  gp.mse.3[i] <- f$norm[cv.gp$lambda.min == f$lambda]
  ncv.mse.3[i] <- f$norm[cv.ncv$lambda.min == f$lambda]
  dvr.mse.3[i] <- f$norm[cv.dvr$lambda.min == f$lambda]
  
  if (typeof(cv.ug) == "character"){
    ug.lambda.3[i] <- NA
    ug.mse.3[i] <- NA
  }else{
    ug.lambda.3[i]<- cv.ug$lambda.min
    ug.mse.3[i] <- f$norm[cv.ug$lambda.min == f$lambda]
  }  
  
}

mean(true.lambda.3)
mean(gp.lambda.3)
mean(ug.lambda.3)
mean(ncv.lambda.3)
mean(dvr.lambda.3)

#save(true.lambda.3,gp.lambda.3,ug.lambda.3,ncv.lambda.3,dvr.lambda.3,
#     true.mse.3,gp.mse.3,ug.mse.3,ncv.mse.3,dvr.mse.3,
#     file = "sim_2_3.RData")

#############################################################################
# scenario 2_4 : n = 500, p = 10000, nbeta = 20, censoring = 10% , nfold = 10
#############################################################################
# change the magnitude of betas
n <- 500
p <- 10000
#cpct <- seq(0,0.3,by = 0.1)
cpct <- 0.1
beta <- generatebeta(p,nbeta = 20 ,c = 1)
N <- 200

# initiation
true.lambda.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.lambda.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.lambda.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.lambda.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.lambda.4 <- matrix(NA,nrow = N ,ncol = length(cpct))

true.mse.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
gp.mse.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
ug.mse.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
ncv.mse.4 <- matrix(NA,nrow = N ,ncol = length(cpct))
dvr.mse.4 <- matrix(NA,nrow = N ,ncol = length(cpct))

for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(1:10,size=n,replace=TRUE)
  data <- genSurvData(n=n , beta = beta, cpct = cpct,sd = 1)
  data$x <- unclass(ncvreg::std(data$x))
  
  f <- try(survncvreg(data$x,data$y,data$beta))
  cv.gp <- try(cv.glmnet(x = data$x,y = data$y,foldid=foldid,lambda = f$lambda,family = "cox"))
  cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,lambda = f$lambda,grouped = FALSE,family = "cox"))
  cv.ncv <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
  cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,lambda = f$lambda,penalty = "lasso")) 
  
  # lambdas
  true.lambda.4[i] <- f$lambda_t
  gp.lambda.4[i]<-cv.gp$lambda.min
  ncv.lambda.4[i]<- cv.ncv$lambda.min
  dvr.lambda.4[i]<- cv.dvr$lambda.min
  
  # mses
  true.mse.4[i] <- min(f$norm)
  gp.mse.4[i] <- f$norm[cv.gp$lambda.min == f$lambda]
  ncv.mse.4[i] <- f$norm[cv.ncv$lambda.min == f$lambda]
  dvr.mse.4[i] <- f$norm[cv.dvr$lambda.min == f$lambda]
  
  if (typeof(cv.ug) == "character"){
    ug.lambda.4[i] <- NA
    ug.mse.4[i] <- NA
  }else{
    ug.lambda.4[i]<- cv.ug$lambda.min
    ug.mse.4[i] <- f$norm[cv.ug$lambda.min == f$lambda]
  }  
  
}

mean(true.lambda.4)
mean(gp.lambda.4)
mean(ug.lambda.4)
mean(ncv.lambda.4)
mean(dvr.lambda.4)

#save(true.lambda.4,gp.lambda.4,ug.lambda.4,ncv.lambda.4,dvr.lambda.4,
#     true.mse.4,gp.mse.4,ug.mse.4,ncv.mse.4,dvr.mse.4,
#     file = "sim_2_4.RData")

#making tables that combines the results here
load("sim_2_1.RData")
load("sim_2_2.RData")
load("sim_2_3.RData")
load("sim_2_4.RData")
col <- c(expression(lambda),"SE","MSR Ratio")
row <- c("optimal","grouped","ungrouped","linear predictor","deviance residuals")
  
tb1<- cbind(
  c(mean(true.lambda.1),
  mean(gp.lambda.1),
  mean(ug.lambda.1),
  mean(ncv.lambda.1),
  mean(dvr.lambda.1)),
  c(sd(true.lambda.1),
    sd(gp.lambda.1),
    sd(ug.lambda.1),
    sd(ncv.lambda.1),
    sd(dvr.lambda.1)),
  c(NA,
    mean(gp.mse.1/true.mse.1),
    mean(ug.mse.1/true.mse.1),
    mean(ncv.mse.1/true.mse.1),
    mean(dvr.mse.1/true.mse.1))
)

tb2<- cbind(
  c(mean(true.lambda.2),
    mean(gp.lambda.2),
    mean(ug.lambda.2),
    mean(ncv.lambda.2),
    mean(dvr.lambda.2)),
  c(sd(true.lambda.2),
    sd(gp.lambda.2),
    sd(ug.lambda.2),
    sd(ncv.lambda.2),
    sd(dvr.lambda.2)),
  c(NA,
    mean(gp.mse.2/true.mse.2),
    mean(ug.mse.2/true.mse.2),
    mean(ncv.mse.2/true.mse.2),
    mean(dvr.mse.2/true.mse.2))
)

tb3 <- cbind(
  c(mean(true.lambda.3),
    mean(gp.lambda.3),
    mean(ug.lambda.3),
    mean(ncv.lambda.3),
    mean(dvr.lambda.3)),
  c(sd(true.lambda.3),
    sd(gp.lambda.3),
    sd(ug.lambda.3),
    sd(ncv.lambda.3),
    sd(dvr.lambda.3)),
  c(NA,
    mean(gp.mse.3/true.mse.3),
    mean(ug.mse.3/true.mse.3),
    mean(ncv.mse.3/true.mse.3),
    mean(dvr.mse.3/true.mse.3))
)

tb4 <- cbind(
  c(mean(true.lambda.4),
    mean(gp.lambda.4),
    mean(ug.lambda.4),
    mean(ncv.lambda.4),
    mean(dvr.lambda.4)),
  c(sd(true.lambda.4),
    sd(gp.lambda.4),
    sd(ug.lambda.4),
    sd(ncv.lambda.4),
    sd(dvr.lambda.4)),
  c(NA,
    mean(gp.mse.4/true.mse.4),
    mean(ug.mse.4/true.mse.4),
    mean(ncv.mse.4/true.mse.4),
    mean(dvr.mse.4/true.mse.4))
)

colnames(tb1) <- col
rownames(tb1) <- row

colnames(tb2) <- col
rownames(tb2) <- row

colnames(tb3) <- col
rownames(tb3) <- row

colnames(tb4) <- col
rownames(tb4) <- row

knitr::kable(tb1, caption = "n = 100, p = 20")
knitr::kable(tb2, caption = "n = 100, p = 100")
knitr::kable(tb3, caption = "n = 100, p = 1000")
knitr::kable(tb4, caption = "n = 100, p = 10000")
