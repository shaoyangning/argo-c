### 20230207
### ARGO_SGL regional pipeline

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
if(!file.exists(out.folder)){
  dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)
}  

path_out_reg <- file.path(out.folder, "Regions")
if(!file.exists(path_out_reg)){
  dir.create(path_out_reg, showWarnings = FALSE, recursive = TRUE)
}  

source("argo_data_input_script.R")

##############
####regional
ili_national <- state_data$ili_national
GT_national <- state_data$GT_national

ili_regional <- state_data$ili_regional
GT_regional <- state_data$GT_regional

#save(ili_national, GT_national, ili_regional, GT_regional,
#     seed, GT_all, ili_version, file = file.path(path_out_reg, "data_ili_GT.Rdata"))

load(file.path(path_out_reg, "data_ili_GT.Rdata"))

## clustering
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")
# covid clustering period
covid_train <- seq(as.Date("2015-03-14"), as.Date("2020-03-14"), by = "week")


#############clustering
# transposing and scaling the matrix
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
periods_all <- c("pre09", "post10", "covid_train")
method <- "hc_corr_ave"
k_stats_all <- c("elbow", "gap", "silhouette")

idx_optim_k <- NULL ## 
#idx_optim_k <-"kopt_regional_subopt" ## manually screened for larger k
####idx_optim_k <- "kopt_regional" ## silhouette mostly

# source("GT_clust_region.R")

###########ARGO-SGL on regional: <-> step 1 of ARGO2
lambdas_all <- c(10^seq(-4.0, -2.8, length.out = 20), 10^seq(-2.8, -2, length.out = 10)) #lambda2
alpha_sgl <- 0.95

lag_cur <-  NULL # consistent with argo2
#lag_cur <- 1:52
seed <- 1000
idx_precovid_clust <- T

## step 1 + step 2
#source("argo_sgl_script_reg.R")

## summary
source("argo_sgl_summary_reg.R")
## rmse is saved as MSE here!
prefix_out1 <- paste("argo_sgl_reg", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
prefix_out1 <- paste0(paste(prefix_out1, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
       ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata")
out_name_result <- paste0("summary", prefix_out1)

save(tab.allregion, interval.coverage, file = file.path(path_out_reg, out_name_result))
