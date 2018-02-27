# Scenario 1 (plots the lambda selection and MSE ratio)
# compare four cross validation methods: 
#  ungrouped, grouped, linear predictor, deviance residuals
# varying magnitude of betas at higher dimension (n = 100, p = 1000, nonzero = 10)


# load source files
source("survival_functions.R")

# load modified ungrouped method
source("modify_fold.R")

# load deviance residuals
source("deviance.R")

#############################################################################
# comparison 1 : n = 200, p = 1000, nbeta = 10, censoring = 10% , nfold = 10
#############################################################################
# change the magnitude of betas
n <- 100
p <- 1000
c <-  c(0.5,0.75,1,1.25,1.5)
beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = 10 ,c = c[i])
}
N <- 200
# initiation
true.lambda <- matrix(NA,nrow = N ,ncol = length(c))
gp.lambda <- matrix(NA,nrow = N ,ncol = length(c))
ug.lambda <- matrix(NA,nrow = N ,ncol = length(c))
ncv.lambda <- matrix(NA,nrow = N ,ncol = length(c))
dvr.lambda <- matrix(NA,nrow = N ,ncol = length(c))


true.mse <- matrix(NA,nrow = N ,ncol = length(c))
gp.mse <- matrix(NA,nrow = N ,ncol = length(c))
ug.mse <- matrix(NA,nrow = N ,ncol = length(c))
ncv.mse <- matrix(NA,nrow = N ,ncol = length(c))
pd.mse <- matrix(NA,nrow = N ,ncol = length(c))
dvr.mse <- matrix(NA,nrow = N ,ncol = length(c))

#i <- 1
#j <- 1
for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(rep(seq(10),length = n))
  for (j in 1:length(c)){
    data <- genSurvData(n=n , beta = beta[j,], cpct = 0.1, sd = 1)
    f <- try(survncvreg(data$x,data$y,data$beta))
    cv.gp <- try(cv.glmnet(x = data$x,y = data$y,foldid=foldid,lambda = f$lambda,family = "cox"))
    cv.ug <- try(cv.glmnet(data$x,data$y,nfolds = 10,lambda = f$lambda,grouped = FALSE,family = "cox"))
    cv.ncv <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
    cv.dvr <- try(cv.ncvsurv.KP(X = data$x,y = data$y,nfolds = 10,lambda = f$lambda,penalty = "lasso")) 

    # lambdas
    true.lambda[i,j] <- f$lambda_t
    gp.lambda[i,j]<-cv.gp$lambda.min
    ncv.lambda[i,j]<- cv.ncv$lambda.min
    dvr.lambda[i,j]<- cv.dvr$lambda.min
    
    # mses
    true.mse[i,j] <- min(f$norm)
    gp.mse[i,j] <- f$norm[cv.gp$lambda.min == f$lambda]
    ncv.mse[i,j] <- f$norm[cv.ncv$lambda.min == f$lambda]
    dvr.mse[i,j] <- f$norm[cv.dvr$lambda.min == f$lambda]
    
    if (typeof(cv.ug) == "character"){
      ug.lambda[i,j] <- NA
      ug.mse[i,j] <- NA
    }else{
      ug.lambda[i,j]<- cv.ug$lambda.min
      ug.mse[i,j] <- f$norm[cv.ug$lambda.min == f$lambda]
    }  
  }
}

#apply(na.omit(true.lambda),2,mean)
#apply(na.omit(gp.lambda),2,mean)
#apply(na.omit(ug.lambda),2,mean)
#apply(na.omit(ncv.lambda),2,mean)
#apply(na.omit(dvr.lambda),2,mean)

save(true.lambda,gp.lambda,ug.lambda,ncv.lambda,dvr.lambda,
     gp.mse,ug.mse,ncv.mse,dvr.mse,true.mse,
     file = "sim_1_4.RData")

#load(file = "sim_1_4.RData")
# lambda
com_result1 <- matrix(NA,nrow = 5,ncol = length(c))
com_result1[1,] <- apply(true.lambda,2,mean)
com_result1[2,] <- apply(na.omit(gp.lambda),2,mean)
com_result1[3,] <- apply(na.omit(ug.lambda),2,mean)
com_result1[4,] <- apply(na.omit(ncv.lambda),2,mean)
com_result1[5,] <- apply(na.omit(dvr.lambda),2,mean)

colnames(com_result1) <- c
rownames(com_result1) <- c("true","grouped","ungrouped","ncvreg","deviance")

# lambda sd
sd_result1 <- matrix(NA,nrow = 5,ncol = length(c))
sd_result1[1,] <- apply(true.lambda,2,sd)
sd_result1[2,] <- apply(na.omit(gp.lambda),2,sd)
sd_result1[3,] <- apply(na.omit(ug.lambda),2,sd)
sd_result1[4,] <- apply(na.omit(ncv.lambda),2,sd)
sd_result1[5,] <- apply(na.omit(dvr.lambda),2,sd)


# ratio of lambda:
com_ratio1 <- matrix(NA,nrow = 4,ncol = length(c))
com_ratio1[1,] <- apply(gp.lambda/true.lambda,2,mean)
com_ratio1[2,] <- apply(na.omit(ug.lambda/true.lambda),2,mean)
com_ratio1[3,] <- apply(ncv.lambda/true.lambda,2,mean)
com_ratio1[4,] <- apply(dvr.lambda/true.lambda,2,mean)
colnames(com_ratio1) <- c
rownames(com_ratio1) <- c("grouped/true","ungrouped/true","ncvreg/true","dvr/true")

# MSE
com_mse1 <- matrix(NA,nrow = 5,ncol = length(c))
com_mse1[1,] <- apply(na.omit(true.mse),2,mean)
com_mse1[2,] <- apply(na.omit(gp.mse),2,mean)
com_mse1[3,] <- apply(na.omit(ug.mse),2,mean)
com_mse1[4,] <- apply(na.omit(ncv.mse),2,mean)
com_mse1[5,] <- apply(na.omit(dvr.mse),2,mean)

colnames(com_mse1) <- c
rownames(com_mse1) <- c("true","grouped","ungrouped","ncvreg","deviance")

# ratio of MSE
mse_ratio <- matrix(NA,nrow = 4,ncol = length(c))
mse_ratio[1,] <- apply(na.omit(gp.mse/true.mse),2,mean)
mse_ratio[2,] <- apply(na.omit(ug.mse/true.mse),2,mean)
mse_ratio[3,] <- apply(na.omit(ncv.mse/true.mse),2,mean)
mse_ratio[4,] <- apply(na.omit(dvr.mse/true.mse),2,mean)

colnames(mse_ratio) <- c
rownames(mse_ratio) <- c("grouped","ungrouped","ncvreg","deviance")

#plot lambda
png(filename= "sim_1_4_lambda.png")
maintitle <- expression(paste(lambda," vs the magnitude of non-zero ", beta, "s", sep = ""))
xlab <- expression(paste("the magnitude of non-zero ", beta, "s", sep = ""))
plot(x = c, y = com_result1[1,],ylim = c(0, 0.40), type = "l",lwd = 2,main = maintitle,
     ylab = expression(lambda),xlab = xlab, axes = FALSE,cex = 2)
axis(side = 1, at = c,labels = c)
axis(side = 2,las = 1)
lines(x = c, y = com_result1[2,],lty = 1,col = 2,lwd = 2)
lines(x = c, y = com_result1[3,],lty = 1,col = 3,lwd = 2)
lines(x = c, y = com_result1[4,],lty = 1,col = 4,lwd = 2)
lines(x = c, y = com_result1[5,],lty = 1,col = 6,lwd = 2)
legend("topright",legend = c("Optimal","Grouped","Ungrouped","Linear Pred","Deviance"),lty =1,col = c("black",2:4,6))
dev.off()

#plot MSE ratio
png(filename= "sim_1_4_mse.png")
maintitle <- expression(paste("MSE Ratio"," vs the magnitude of non-zero ", beta, "s", sep = ""))
xlab <- expression(paste("the magnitude of non-zero ", beta, "s", sep = ""))
plot(x = c, y = mse_ratio[1,],ylim = c(0, 6), type = "l",lwd = 2,main = maintitle,
     ylab = "MSE Ratio",xlab = xlab, axes = FALSE,col = 2,cex = 2)
axis(side = 1, at = c,labels = c)
axis(side = 2,las = 1)
lines(x = c, y = mse_ratio[2,],lty = 1,col = 3,lwd = 2)
lines(x = c, y = mse_ratio[3,],lty = 1,col = 4,lwd = 2)
lines(x = c, y = mse_ratio[4,],lty = 1,col = 6,lwd = 2)
legend("topleft",legend = c("Grouped","Ungrouped","Linear Pred","Deviance"),lty =1,col = c(2:4,6))
dev.off()


