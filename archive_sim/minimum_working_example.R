#------------------------------  
# a minimal working example
# to validate cv.ncvreg 
#------------------------------
library(survival)
library(ncvreg)

# generate data
beta <- c(rep(1,5),rep(0,20))
n <- 100
p <- length(beta)
x <- matrix(rnorm(n*p, sd = 1),n,p)
x <- ncvreg::std(x)
hazard <- exp(x %*% beta)
time <- rexp(n, hazard)
y <- Surv(time = time, event = rbinom(n, 1, 0.7))

data <- list(y=y,x=x)

# inital fit

fit <- ncvsurv(data$x,data$y,penalty = "lasso")
lambda <- fit$lambda

# cross-validation
nfold <- 10
nevent <- sum(data$y[,2])
sde <- sqrt(.Machine$double.eps)
foldid <- numeric(n)
foldid[data$y[,2] == 1] <-  ceiling(sample(1:nevent)/(nevent + sde) * nfold)
foldid[data$y[,2] == 0] <-  ceiling(sample(1:(n-nevent))/(n - nevent + sde) * nfold)

# use cv.ncvsurv:
cv_ncv <- cv.ncvsurv(data$x,data$y,penalty = "lasso",
                     fold = foldid, lambda = lambda,
                     returnY = TRUE,returnX=TRUE)
cv_ncv$fit$X[1:10,1:10]
data$x[1:10,1:10]

sort(cv_ncv$fit$X[,1])
sort(data$x[,1])

# # use Biyue's code:
# eta.cv <- matrix(nrow = n, ncol = length(lambda))
# l_diff <- matrix(nrow = nfold , ncol = length(lambda))
# l_test <- matrix(nrow = nfold, ncol = length(lambda))
# 
# ind <- order(data$y[,1])
# d <- data$y[ind,2]
# 
# for(f in 1:nfold){
#   ifold <- foldid == f
#   cvfit <- ncvsurv(cv_ncv$fit$X[!ifold,], data$y[!ifold,],
#                    penalty = "lasso", lambda = lambda)
#   cvlambda <- cvfit$lambda
#   beta_train <- cvfit$beta
#   
#   # linear predictors
#   eta.cv[ifold,lambda%in%cvlambda] <- predict(cvfit,X = data$x[ifold,], type = "link")
# }
# 
# eta.cv <- eta.cv[,apply(eta.cv,2,function(x){sum(is.na(x)) == 0})]
# eta.cv.ordered <- eta.cv[ind,,drop = FALSE]
# r <- apply(eta.cv.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
# val_LP <- -2 * (crossprod(d, eta.cv.ordered) - crossprod(d, log(r)))
# 
# 
# #------------------------------------------------------------------
# # compare results:
# 
# which.min(cv_ncv$cve)
# which.min(val_LP)
# 
# # replace eta.cv.ordered with Y,
# # then they would be pretty much the same
# eta.cv.ordered <- cv_ncv$Y
# r <- apply(eta.cv.ordered, 2, function(x) rev(cumsum(rev(exp(x)))))
# val_LP <- -2 * (crossprod(d, eta.cv.ordered) - crossprod(d, log(r)))
# which.min(val_LP)
# as.numeric(val_LP/nevent)
# cv_ncv$cve
# 
# # X matrices are different:
# cv_ncv$fit$
