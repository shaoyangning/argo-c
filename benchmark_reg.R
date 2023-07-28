####20230412 
### benchmarks: ARGO2, VAR, GFT

###ARGOX
##GTstate API raw
## 20211001 

options(echo=TRUE)
library(abind)
library(argo)
library(boot)
library(cluster)
library(factoextra)
library(fossil)
library(glmnet)
library(parallel)
library(SGL)
library(xts)
## parameters

seed <- 1000

##work space
##need to change to local directory
path_ws <- "/Users/sn9/Dropbox/Williams/Research/SGL_ARGO/ARGOX_shared2022"
path_data <- file.path(path_ws, "data")

## data version
GT_all <- "api_raw_results-2023-01-29"
ili_version <- "ili20230207/"


## ili
ili.folder <- file.path(path_data, ili_version)
## GT
GT_version <- paste0(GT_all, "/")
## output
out.folder <- paste0(ili.folder, GT_version)
if(!exists(out.folder)){
  dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)
}  

path_out_reg <- file.path(out.folder, "Regions")
if(!exists(path_out_reg)){
  dir.create(path_out_reg, showWarnings = FALSE, recursive = TRUE)
}  


### data input
source("argo_data_input_script.R")


ili_national <- state_data$ili_national
GT_national <- state_data$GT_national

ili_regional <- state_data$ili_regional
GT_regional <- state_data$GT_regional


#### first-step argo prediction ####
transY <- function(y){
  logit((y+1e-1) / 100)
}

inv_transY <- function(y){
  100*logit_inv(y)-1e-1
}

transYnat <- function(y){
  logit((y+1e-1) / 100)
}

inv_transYnat <- function(y){
  100*logit_inv(y)-1e-1
}

transYstate <- function(y){
  logit((y+1e-1) / 100)
}

inv_transYstate <- function(y){
  100*logit_inv(y)-1e-1
}

######################################################################################################
##############################################################################
###regional ARGO
### Only need to run once, can change to False later
if(T){
  get_argo1_regional <- function(terms, period, reglag = NULL, seed=1000){
    common_idx <- period
    set.seed(seed)
    
    # national 
    idx.nat1 <- c(index(GT_national)[1]-(52:1)*7, index(GT_national))
    GT_national1 <- xts(matrix(nrow=length(idx.nat1), ncol=ncol(GT_national)), order.by = as.Date(idx.nat1))
    colnames(GT_national1) <- colnames(GT_national)
    GT_national1[index(GT_national), ] <- GT_national
    common_idx_nat <- c(common_idx[1]-(52:1)*7, common_idx)
    argo_national <- argo(data = transY(ili_national[common_idx_nat]),
                          exogen = log(GT_national1[common_idx_nat, terms]+1),
                          mc.cores = 8)
    
    # regional
    argo_result <- list()
    all.pred <- list()
    for(region.id in 1:10){
      j <- paste0("Region.", region.id)
      set.seed(seed)
      argo_result[[j]] <- argo(transY(ili_regional[common_idx, j]),
                               log(GT_regional[[region.id]][common_idx, terms]+1),
                               mc.cores = 8,
                               N_lag = reglag)
      
      pred_xts_blend <- inv_transY(argo_result[[j]]$pred)
      pred_xts_blend <- merge(ili_regional[,j], pred_xts_blend, all=FALSE)
      pred_xts_blend$naive <- c(NA, as.numeric(pred_xts_blend[1:(nrow(pred_xts_blend)-1), j]))
      names(pred_xts_blend)[1] <- "CDC.data"
      all.pred[[j]] <- pred_xts_blend
      print(j)
    }
    list(argo_national=argo_national,
         argo_regional_result=argo_result,
         all.pred=all.pred)
  }
  
  common_idx <- index(merge(ili_national, GT_national, all=FALSE))
  common_idx09 <- common_idx[common_idx < as.Date("2010-05-22")]
  terms <- colnames(GT_national)
  terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
  
  # first-step argo for two periods
  argo1.09 <- NULL
  if(length(common_idx09)>104){
    argo1.09 <- get_argo1_regional(terms09, common_idx09, reglag = NULL, seed=1000)
  }
  argo1.post10 <- get_argo1_regional(terms, common_idx, reglag = NULL, seed=1000)
  
  #### blend 09 data to all data ####
  argo1.reg.coef <- list()
  all.pred.reg <- argo1.post10$all.pred
  
  for(region.id in 1:10){
    j <- paste0("Region.", region.id)
    all.pred.reg[[j]][index(argo1.09$all.pred[[j]]),"predict"] <- argo1.09$all.pred[[j]][,"predict"]
    argo1.reg.coef[[j]] <- argo1.post10$argo_regional_result[[j]]$coef
    argo1.reg.coef[[j]][,colnames(argo1.09$argo_regional_result[[j]]$coef)] <- NA
    argo1.reg.coef[[j]][rownames(argo1.09$argo_regional_result[[j]]$coef),colnames(argo1.09$argo_regional_result[[j]]$coef)] <- 
      argo1.09$argo_regional_result[[j]]$coef
  }
  argo.nat.p <- inv_transY(argo1.post10$argo_national$pred)
  if(length(index(argo1.09$argo_national$pred))>0){
    argo.nat.p[index(argo1.09$argo_national$pred)] <- inv_transY(argo1.09$argo_national$pred)
  }
  argo.nat.coef <- argo1.post10$argo_national$coef
  argo.nat.coef[,colnames(argo1.09$argo_national$coef)] <- NA
  argo.nat.coef[rownames(colnames(argo1.09$argo_national$coef)),colnames(argo1.09$argo_national$coef)] <- 
    argo1.09$argo_national$coef
  
  
  argo1.09.reg <- argo1.09
  argo1.post10.reg <- argo1.post10
  
  #### argo second step ####
  argo.reg.pred <- lapply(all.pred.reg, function(x) x[,"predict"])
  argo.reg.pred <- do.call(merge, argo.reg.pred)
  colnames(argo.reg.pred) <- paste0(colnames(argo.reg.pred)[1], ".", 1:10)
  
  ili_regional2003post <- ili_regional["2003/"]
  argo2_reg_result <- argo2(ili_regional2003post, argo.reg.pred, argo.nat.p)
  
  
  
  
  save(ili_national, ili_regional, argo2_reg_result,
       argo.nat.p, argo.reg.pred,
       argo1.post10.reg, argo1.09.reg,
       file=file.path(path_out_reg, "argo2_reg_all_logit01_lag0a.Rdata"))
  
  #print(paste0("Nat/Reg: ", GTpub_version, " Done."))
}




list_argo_var_reg_pred <- list()
for(reg in 1:10){
  path_out <- file.path(path_out_reg, paste0("Region", reg))
  prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl,  "lag", length(lag_cur), sep = "_")
  out_name <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata"))
  load(out_name)
  list_argo_var_reg_pred[[method]][[reg]] <- var.pred
}  

var_sgl_reg_pred <- do.call(merge, list_argo_var_reg_pred[[method]])
colnames(var_sgl_reg_pred) <- paste0("Region", ".", 1:10)


save(var_sgl_reg_pred,
     file=file.path(path_out_reg, "var_reg_all_lag1.Rdata"))
