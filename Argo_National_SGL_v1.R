## 20221021
## ARGO + SGL, national pipeline
## S. Ning
## based on codes by A. Hussain

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
population_file <- file.path(path_data, "Population.csv") 
gft_file <- file.path(path_data, "GFT.txt")

## ili
ili.folder <- file.path(path_data, ili_version)
## GT
GT_version <- paste0(GT_all, "/")
## output
out.folder <- paste0(ili.folder, GT_all)
if(!exists(out.folder)){
  dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)
}  

source("argo_data_input_script.R")

##############
####national

## national ili data
ili_national <- state_data$ili_national
GT_national <- state_data$GT_national


#####!!! need to work on this 
## 1. automated/yearly update cluster

## clustering
## periods
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")

covid_2020 <- seq(as.Date("2020-03-14"), as.Date("2021-04-03"), by = "week")
covid_train <- seq(as.Date("2015-03-14"), as.Date("2020-03-14"), by = "week")


# transposing and scaling the matrix
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
pre09 <- t(scale(GT_national[pre_2009, terms09]))
post10 <- t(scale(GT_national[post_2010, ]))


## Covid period is from 2020-03-14 to 2021-04-03
covid_GT <- GT_national[covid_train,]
covid <- t(scale(covid_GT))

## loading cluster results

load("GT_clust_optim_k.Rdata")

method_all <- c("hc_corr_ave"  #"pam", "kmeans", "hc_corr_comp", "hc_corr_single"
)
periods_all <- c("pre09", "post10", #"covid",
                 "covid_train")
k_stats_all <- c("elbow", "gap", "silhouette")

### change method!
method <- "hc_corr_ave"

cluster_pre09 <- list_cluster_results[[method]]$pre09
cluster_post10 <- list_cluster_results[[method]]$post10
cluster_covid <- list_cluster_results[[method]]$covid_train


## lambdas: need length 20
#lambdas_all <- 10^seq(-4.0, -2.8, length.out = 20) # lambda1

lambdas_all <- c(10^seq(-4.0, -2.8, length.out = 20), 10^seq(-2.8, -2, length.out = 10)) #lambda2
alpha_sgl <- 0.95
idx_precovid_clust <- T

source("argo_sgl_script_nat.R")
## need to streamline this
### 1. customize training + prediction range for each period
### 2. make it to weekly predict?

prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, sep = "_")
prefix_out <- paste0(prefix_out, ifelse(idx_precovid_clust, "_precovidclust", ""))
path_out <- file.path(out.folder, "national")
out_name <- file.path(path_out, paste0(prefix_out, ".Rdata"))
if(!exists(path_out)){
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
}
if(F){
  save(seed, lambdas_all, method, cluster_pre09, cluster_post10, cluster_covid, 
       argo_SGL_09, argo_SGL_10, argo_SGL_covid,
       alpha_sgl, var.pred,
       argo.nat.p, argo_nat_sgl, file = out_name)
}


## check lambda selection index

####################################################
####################################################
## organize
list_argo_sgl_pred <- list()

#for(method in method_all){
prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, sep = "_")
prefix_out <- paste0(prefix_out, ifelse(idx_precovid_clust, "_precovidclust", ""))
out_name_curr <- paste0(prefix_out, ".Rdata")
load(file.path(path_out, out_name_curr))
list_argo_sgl_pred[[method]] <- argo_nat_sgl
#}
method_print <- sapply(names(list_argo_sgl_pred), function(x){paste("ARGO_SGL", "covid", x, sep="_")})
names(list_argo_sgl_pred) <- method_print

## GFT
source("ARGOX_shared2022/argo2_states_function.R")
list_gft <- load_gft_data1(population.file=population_file, gft.file=gft_file)
list_argo_sgl_pred[["GFT"]] <- list_gft$GFT_nat

## VAR1
list_argo_sgl_pred[["VAR1"]] <- var.pred

colnames(var.pred) <- paste0("var1")
argo_sgl_pred_all <- do.call(merge, list_argo_sgl_pred)
method_print <- names(list_argo_sgl_pred)
colnames(argo_sgl_pred_all) <- method_print

source("argo_sgl_summary.R")

out_name_result <- paste0(paste("summary", prefix_out), ".Rdata")
save(tab, file = file.path(path_out, out_name_result))

## print out 


##### interval estimate
prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, sep = "_")
prefix_out <- paste0(prefix_out, ifelse(idx_precovid_clust, "_precovidclust", ""))
path_out <- file.path(out.folder, "national")
out_name_curr <- paste0(prefix_out, ".Rdata")
load(file.path(path_out, out_name_curr))


