## 20230516
### ARGO_SGL state pipeline

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

path_out_state <- file.path(out.folder, "States")
if(!file.exists(path_out_state)){
  dir.create(path_out_state, showWarnings = FALSE, recursive = TRUE)
}  

source("argo_data_input_script.R")

##################
####state
ili_national <- state_data$ili_national
GT_national <- state_data$GT_national

ili_regional <- state_data$ili_regional
GT_regional <- state_data$GT_regional

ili_state <- state_data$ili_state
GT_state <- state_data$GT_state

## state ili organize
states.abb.all <- colnames(ili_state)
states.rm <- c("US.FL", "US.Virgin Islands",
               "US.Puerto Rico", "US.Commonwealth of the Northern Mariana Islands")
states.abb.all <- states.abb.all[-match(c(states.rm),states.abb.all)]
ili_state <- ili_state[, states.abb.all]
## impute NA
idx_na1 <- which(is.na(ili_state$US.DC))
ili_state$US.DC[idx_na1] <- ili_state$US.DC[idx_na1-1]

## state GT organize
GT_state[["US"]] <- NULL
GT_state[["US.FL"]] <- NULL
names(GT_state)[which(names(GT_state)=="501")] <-  "US.NYC"
colnames(ili_state)[which(colnames(ili_state)=="US.New York City")] <-  "US.NYC"

states.abb.all <- colnames(ili_state)

## state info
tab.states <- read.csv(file.path(path_data, "Population.csv"))
states.info <- data.frame(names=states.abb.all, 
                          abbre=unlist(lapply(strsplit(states.abb.all, split = "[.]"), function(x){x[2]})))

states.info$region <- tab.states$Region[match(states.info$abbre, tab.states$Abbre)]
states.info$region[states.info$abbre=="NYC"] <- 2
states.info$states <- tab.states$State[match(states.info$abbre, tab.states$Abbre)]

states.info$states[states.info$abbre=="NYC"] <- "New York NY"
#save(ili_national, GT_national, ili_regional, GT_regional,
#     ili_state, GT_state, tab.states, states.abb.all,states.info,
#    seed, GT_all, ili_version, file = file.path(path_out_state, "data_ili_GT.Rdata"))
 

#load(file.path(path_out_state, "data_ili_GT.Rdata"))

###################
## clustering    
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")
# covid clustering period
covid_train <- seq(as.Date("2015-03-14"), as.Date("2020-03-14"), by = "week")

# transposing and scaling the matrix
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
periods_all <- c("pre09", "post10", "covid_train")
method <- "hc_corr_ave"
k_stats_all <- c("elbow", "gap", "silhouette")

idx_optim_k <- NULL ## 
## use national idx optim

## GT with regmix 

idx.GTregmix <- F
# source("GT_clust_state.R")
idx.GTregmix <- T
reg.weight.mix <- 1/3
# source("GT_clust_state.R")

###########ARGO-SGL on state: <-> step 1 of ARGOX
lambdas_all <- c(10^seq(-4.0, -2.8, length.out = 20), 10^seq(-2.8, -2, length.out = 10)) #lambda2
alpha_sgl <- 0.95

lag_cur <-  NULL # consistent with argo2
seed <- 1000
idx_precovid_clust <- T

#source("argo_sgl_script_state.R")

## summary
## benchmarks
#source("argo_sgl_summary_state.R")

prefix_out1 <- paste("argo_sgl_state", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
out_name_result <- paste0("summary_", paste(prefix_out1, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
                      ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata")
#save(tab.allregion, file = file.path(path_out_state, out_name_result))

