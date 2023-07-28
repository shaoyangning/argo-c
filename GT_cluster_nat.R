## ARGO-SGL: cluster step
## 20221101 covid period clustering
## S. Ning
## 20221025 organizing cluster results by Ahmed
##

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


path_ws <- "/Users/sn9/Dropbox/Williams/Research/SGL_ARGO/ARGOX_shared2022"

## data version
GT_all <- "api_raw_results-2023-01-29"
ili_version <- "ili20230207/"

GT_version <- paste0(GT_all, "/")
## output
out.folder <- paste0(ili.folder, GT_version)
if(!exists(out.folder)){
  dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)
}  


## data input
source("argo_data_input_script.R")
## GT data

ili_national <- state_data$ili_national
GT_national <- state_data$GT_national


seed <- 1000

## clustering
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")
# covid clustering period
covid_train <- seq(as.Date("2015-03-07"), as.Date("2020-02-29"), by = "week")


# transposing and scaling the matrix
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))


path_out <- file.path(out.folder, "national")
if(!exists(path_out)){
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
}  


pre09 <- t(scale(GT_national[pre_2009, terms09]))
post10 <- t(scale(GT_national[post_2010, ]))


## Covid period 
covid <- t(scale(GT_national[covid_train, ]))

list_data <- list(pre09 = pre09, post10 = post10, covid_train = covid)


source("GT_cluster.R")

##results.pdf

##
#list_optim_k_results <- list()
#load("GT_clust_optim_k.Rdata")


