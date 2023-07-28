## ARGO-SGL: cluster step -- regional
## pipeline + K determine
## 20230207
## S. Ning
## Input: idx_optim_k, terms09, method

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

source("argo_data_input_script.R")

##############
####regional
ili_national <- state_data$ili_national
GT_national <- state_data$GT_national

ili_regional <- state_data$ili_regional
GT_regional <- state_data$GT_regional

## clustering
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")
# covid clustering period
covid_train <- seq(as.Date("2015-03-07"), as.Date("2020-02-29"), by = "week")


#############clustering
# transposing and scaling the matrix
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
path_out_reg <- file.path(out.folder, "Regions")
if(!exists(path_out_reg)){
  dir.create(path_out_reg, showWarnings = FALSE, recursive = TRUE)
}  

## organize optim k (from manually input table) 
## + save cluster results 
## + fill back group cluster id for all zero entries
# idx == NULL: use national level optim k
# idx: take in table of manually identified optim k
#idx_optim_k <-"kopt_regional_subopt"
#idx_optim_k <- NULL

if(T){
  for(reg in 1:10){
    path_out <- file.path(path_out_reg, paste0("Region", reg))
    if(!exists(path_out)){
      dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
    }  
    ## filter out all zero entries
    pre09 <- t(scale(GT_regional[[reg]][pre_2009, terms09]))
    idx_pre09 <- apply(is.na(pre09), 1, sum) == 0
    post10 <- t(scale(GT_regional[[reg]][post_2010, ]))
    idx_post10 <- apply(is.na(post10), 1, sum) == 0
    covid_GT <- GT_regional[[reg]][covid_train,]
    covid <- t(scale(covid_GT))
    idx_covid <- apply(is.na(covid), 1, sum) == 0
    ## input to GT_cluster.R
    list_data <- list(pre09 = pre09[idx_pre09, ], post10 = post10[idx_post10, ], covid = covid, covid_train = covid[idx_covid, ])
    list_idx_na <- list(idx_pre09, idx_post10, idx_covid)
    #print(sum(idx_pre09))
    source("GT_cluster.R")
  }
}



