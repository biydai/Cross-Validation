setwd("A:\\cv_manuscript_code\\fdr_power")

library(abind)

gp_combined <- NULL
dvr_combined <- NULL
ncv_combined <- NULL
ug_combined <- NULL

#--------------------------------------------------------------
# low dimension
#--------------------------------------------------------------

load("low_fdr.RData")

results <- list()
for(i in 1:4){
  results[[i]] <- rbind(
    apply(gp[i,,],2,mean),
    apply(ncv[i,,],2,mean),
    apply(ug[i,,],2,mean),
    apply(dvr[i,,],2,mean)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- c(0.3,0.6)
}
names(results) <- c("lambda",
                    "MSE",
                    "tp",
                    "fdr")

results$lambda
results$MSE
results$tp
results$fdr

results_low <- results

#--------------------------------------------------------------
# middle dimension
#--------------------------------------------------------------

gp_combined <- NULL
dvr_combined <- NULL
ncv_combined <- NULL
ug_combined <- NULL

for(i in 1:20){
  
  load(paste("mid_fdr_",i,".RData",sep=""))
  gp_1 <- gp
  gp_combined <- abind(gp_combined,gp_1,along = 2)

  dvr_1 <- dvr
  dvr_combined <- abind(dvr_combined,dvr_1,along = 2)
  
  ncv_1 <- ncv
  ncv_combined <- abind(ncv_combined,ncv_1,along = 2)
  
  ug_1 <- ug
  ug_combined <- abind(ug_combined,ug_1,along = 2)
  
}

results <- list()
for(i in 1:4){
  results[[i]] <- rbind(
    apply(gp_combined[i,,],2,mean),
    apply(ncv_combined[i,,],2,mean),
    apply(ug_combined[i,,],2,mean),
    apply(dvr_combined[i,,],2,mean)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- c(0.2,0.6)
}
names(results) <- c("lambda",
                    "MSE",
                    "tp",
                    "fdr")

results$lambda
results$MSE
results$tp
results$fdr

results_mid <- results

#--------------------------------------------------------------
# high dimension
#--------------------------------------------------------------

gp_combined <- NULL
dvr_combined <- NULL
ncv_combined <- NULL
ug_combined <- NULL

for(i in 1:20){
  
  load(paste("high_fdr_",i,".RData",sep=""))
  gp_1 <- gp
  gp_combined <- abind(gp_combined,gp_1,along = 2)
  
  dvr_1 <- dvr
  dvr_combined <- abind(dvr_combined,dvr_1,along = 2)
  
  ncv_1 <- ncv
  ncv_combined <- abind(ncv_combined,ncv_1,along = 2)
  
  ug_1 <- ug
  ug_combined <- abind(ug_combined,ug_1,along = 2)
  
}

results <- list()
for(i in 1:4){
  results[[i]] <- rbind(
    apply(gp_combined[i,,],2,mean),
    apply(ncv_combined[i,,],2,mean),
    apply(ug_combined[i,,],2,mean),
    apply(dvr_combined[i,,],2,mean)
  )
  rownames(results[[i]]) <- c("VVH","LP","ST","DVR")
  colnames(results[[i]]) <- c(0.3,0.6)
}
names(results) <- c("lambda",
                    "MSE",
                    "tp",
                    "fdr")

results$lambda
results$MSE
results$tp
results$fdr

results_high <- results