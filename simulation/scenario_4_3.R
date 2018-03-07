# comparison with other model selection criteria
# AIC
# BIC
# GCV


#n <- 100
#p <- 100
#c <-  1
#beta <- matrix(NA,nrow = length(c), ncol = p)
#for (i in 1:length(c)){
#  beta[i,]<- generatebeta(p,nbeta = 10 ,c = c[i])
#}
#data <- genSurvData(n=n , beta = beta, cpct = 0.1)

# try out AIC computations
# with splines
# fit <- ncvsurv(data$x, Surv(data$y[,1],data$y[,2]),penalty = "lasso")
# which.min(AIC(fit))
# which.min(BIC(fit))

source("survival_functions.R")
#############################################
# n = 100, p = 50, # of none-zero beta = 10 #
#############################################

n <- 100
p <- 10
c <-  c(0.5,0.75,1,1.25,1.5)
beta <- matrix(NA,nrow = length(c), ncol = p)
for (i in 1:length(c)){
  beta[i,]<- generatebeta(p,nbeta = 2 ,c = c[i])
}
N <- 200
# initiation
true.lambda <- matrix(NA,nrow = N ,ncol = length(c))
ncv.lambda <- matrix(NA,nrow = N ,ncol = length(c))
aic.lambda <- matrix(NA,nrow = N ,ncol = length(c))
bic.lambda <- matrix(NA,nrow = N ,ncol = length(c))
#gcv.lambda <- matrix(NA,nrow = N ,ncol = length(c))

true.mse <- matrix(NA,nrow = N ,ncol = length(c))
ncv.mse <- matrix(NA,nrow = N ,ncol = length(c))
aic.mse <- matrix(NA,nrow = N ,ncol = length(c))
bic.mse <- matrix(NA,nrow = N ,ncol = length(c))
#gcv.mse <- matrix(NA,nrow = N ,ncol = length(c))

#i <- 1
#j <- 2
for(i in 1:N){
  set.seed(999 + i)
  foldid <- sample(rep(seq(10),length = n))
  for (j in 1:length(c)){
    data <- genSurvData(n=n , beta = beta[j,], cpct = 0.1, sd = 1)
    f <- try(survncvreg(data$x,data$y,data$beta))
    cv.ncvreg <- try(cv.ncvsurv(X = data$x,y = data$y,nfold = 10,lambda = f$lambda,penalty = "lasso"))
    
    # lambdas
    true.lambda[i,j] <- f$lambda_t
    ncv.lambda[i,j]<- cv.ncvreg$lambda.min
    
    #########################
    # select lambda for aic #
    #########################
    
    n.l <- length(f$lambda)
    aic.ss <- smooth.spline(log(f$lambda),f$aic,df = 6)
    d <- diff(rev(aic.ss$y))
    if (all(d<0)){
      ind <- n.l
    }else {
      ind <- min(which(d>0))-1
    }
    if (ind ==0) {ind <- 1}
    aic.lambda[i,j] <- f$lambda[ind]
    
    #########################
    # select lambda for bic #
    #########################
    
    bic.ss <- smooth.spline(log(f$lambda),f$bic,df = 6)
    d <- diff(rev(bic.ss$y))
    if (all(d<0)){
      ind <- n.l
    }else {
      ind <- min(which(d>0))-1
    }
    if (ind ==0) {ind <- 1}
    bic.lambda[i,j] <- f$lambda[ind]
    
    ###########################
    # select lambda using GCV #
    ###########################
    #e <- apply(f$fit$beta,2,function(x){sum(x!=0)})
    #gcv <-  f$fit$loss/((1-e/n)^2*n)
    
    #gcv.ss <- smooth.spline(log(f$lambda),gcv,df = 8)
    #d <- diff(rev(gcv.ss$y))
    #if (all(d<0)){
    #  ind <- n.l
    #}else {
    #  ind <- min(which(d>0))-1
    #}
    #if (ind ==0) {ind <- 1}
    #gcv.lambda[i,j] <- f$lambda[ind]
    
    
    # mses
    true.mse[i,j] <- min(f$norm)
    aic.mse[i,j] <- f$norm[f$lambda == aic.lambda[i,j]]
    bic.mse[i,j] <- f$norm[f$lambda == bic.lambda[i,j]]
    #gcv.mse[i,j] <- f$norm[f$lambda == gcv.lambda[i,j]]
    ncv.mse[i,j] <- f$norm[cv.ncvreg$lambda.min == f$lambda]
    
  }
}

apply(na.omit(true.lambda),2,mean)
apply(na.omit(ncv.lambda),2,mean)
apply(na.omit(aic.lambda),2,mean)
apply(na.omit(bic.lambda),2,mean)
#apply(na.omit(gcv.lambda),2,mean)

save(true.lambda,ncv.lambda,aic.lambda,bic.lambda,#gcv.lambda,
     true.mse,ncv.mse,aic.mse,bic.mse,gcv.mse,
     file = "scenario_4_3.RData")

############
#plotting #
############

#plot lambda
com_result1 <- matrix(NA,nrow = 4,ncol = length(c))
com_result1[1,] <- apply(na.omit(true.lambda),2,mean)
com_result1[2,] <- apply(na.omit(ncv.lambda),2,mean)
com_result1[3,] <- apply(na.omit(aic.lambda),2,mean)
com_result1[4,] <- apply(na.omit(bic.lambda),2,mean)
#com_result1[5,] <- apply(na.omit(gcv.lambda),2,mean)

colnames(com_result1) <- c
rownames(com_result1) <- c("true","ncvreg","aic","bic")

png(filename= "sim_4_3_lambda.png")
maintitle <- expression(paste(lambda," vs the magnitude of non-zero ", beta, "s", sep = ""))
xlab <- expression(paste("the magnitude of non-zero ", beta, "s", sep = ""))
plot(x = c, y = com_result1[1,], type = "l",lwd = 2,main = maintitle,ylim = c(0,0.2),
     ylab = expression(lambda),xlab = xlab, axes = FALSE)
axis(side = 1, at = c,labels = c)
axis(side = 2,las = 1)
lines(x = c, y = com_result1[2,],lty = 1,col = 2,lwd = 2)
lines(x = c, y = com_result1[3,],lty = 1,col = 3,lwd = 2)
lines(x = c, y = com_result1[4,],lty = 1,col = 4,lwd = 2)
#lines(x = c, y = com_result1[5,],lty = 1,col = 5,lwd = 2)
legend("topright",legend = c("Optimal",
                             "Linear Pred","AIC","BIC"),lty = 1, col = 1:4)
dev.off()

#plot MSE ratio
# MSE
mse_ratio <- matrix(NA,nrow = 3,ncol = length(c))
mse_ratio[1,] <- apply(na.omit(ncv.mse/true.mse),2,mean)
mse_ratio[2,] <- apply(na.omit(aic.mse/true.mse),2,mean)
mse_ratio[3,] <- apply(na.omit(bic.mse/true.mse),2,mean)
#mse_ratio[4,] <- apply(na.omit(gcv.mse/true.mse),2,mean)

colnames(mse_ratio) <- c
rownames(mse_ratio) <- c("ncvreg","aic","bic")

png(filename= "sim_4_3_mse.png")
maintitle <- expression(paste("MSE Ratio"," vs the magnitude of non-zero ", beta, "s", sep = ""))
xlab <- expression(paste("the magnitude of non-zero ", beta, "s", sep = ""))
plot(x = c, y = mse_ratio[1,],ylim = c(0, 5), type = "l",lwd = 2,main = maintitle,
     ylab = "MSE Ratio",xlab = xlab, axes = FALSE,col = 2)
axis(side = 1, at = c,labels = c)
axis(side = 2,las = 1)
lines(x = c, y = mse_ratio[2,],lty = 1,col = 3,lwd = 2)
lines(x = c, y = mse_ratio[3,],lty = 1,col = 4,lwd = 2)
#lines(x = c, y = mse_ratio[4,],lty = 1,col = 5,lwd = 2)
legend("topleft",legend = c("Linear Pred","AIC","BIC"),lty =1,col = 2:4)
dev.off()
