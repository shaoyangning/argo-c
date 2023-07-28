## 20221031
## S. Ning
## argo + sgl, national
## input
# ili_national: ili
# GT_national: GT
# cluster_post10, cluster_pre09
# ready to run

source("argo3.R")

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

##############################################################################
## shared dates between GT and ili data
common_idx <- index(merge(ili_national, GT_national, all=FALSE))
terms <- colnames(GT_national)

## fix random seed
seed <- 1000
set.seed(seed)

## national 
# post 2010
idx.nat1 <- c(index(GT_national)[1]-(52:1)*7, index(GT_national))
GT_national1 <- xts(matrix(nrow=length(idx.nat1), ncol=ncol(GT_national)), 
                    order.by = as.Date(idx.nat1))
colnames(GT_national1) <- colnames(GT_national)
GT_national1[index(GT_national), ] <- GT_national
common_idx_nat <- c(common_idx[1]-(52:1)*7, common_idx)

# pre-2009 set of predictors
common_idx09 <- common_idx[common_idx < as.Date("2010-05-22")]
terms09 <- intersect(colnames(GT_national), colnames(load_data()$GC09))
common_idx09_nat <- c(common_idx09[1]-(52:1)*7, common_idx09)

## covid
covid_period <- "2020-03-14/2021-04-03"

common_covid <- common_idx[common_idx >= as.Date("2020-03-14")]
common_covid_nat <- c(common_covid[1]-((2+52*3):1)*7, common_covid)

# post2010
set.seed(seed)
argo_SGL_10 <- argo3(data = transY(ili_national[common_idx_nat,]),
                     exogen = log(GT_national1[common_idx_nat, terms]+1), 
                     clusterGroups = cluster_post10$cluster, lambdas_sgl = lambdas_all,
                     alpha = alpha_sgl,
                     mc.cores = 8)

print("post10 done")

# pre-2009
set.seed(seed)
argo_SGL_09 <- argo3(data = transY(ili_national[common_idx09_nat,]),
                     exogen = log(GT_national1[common_idx09_nat, terms09]+1), 
                     clusterGroups = cluster_pre09$cluster, lambdas_sgl = lambdas_all,
                     alpha = alpha_sgl,
                     mc.cores = 8)

print("pre09 done")

# covid_train
argo_SGL_covid <- NULL
if(idx_precovid_clust){
  set.seed(seed)
  argo_SGL_covid <- argo3(data = transY(ili_national[common_covid_nat,]),
                          exogen = log(GT_national1[common_covid_nat, terms]+1), 
                          clusterGroups = cluster_covid$cluster, lambdas_sgl = lambdas_all,
                          alpha = alpha_sgl,
                          mc.cores = 8)
  
  print("covid done")
}


## concatenate
argo_nat_sgl <- inv_transY(argo_SGL_10$pred)
if(length(index(argo_SGL_09$pred))>0){
  argo_nat_sgl[index(argo_SGL_09$pred)] <- inv_transY(argo_SGL_09$pred)
}

if(idx_precovid_clust){
  argo_nat_sgl[covid_period] <- inv_transY(argo_SGL_covid$pred[covid_period])
}

## vanilla argo
set.seed(seed)
argo_national <- argo(data = transY(ili_national[common_idx_nat]),
                      exogen = log(GT_national1[common_idx_nat, terms]+1),
                      mc.cores = 8)

argo_national09 <- argo(data = transY(ili_national[common_idx09_nat]),
                        exogen = log(GT_national1[common_idx09_nat, terms09] + 1),
                        mc.cores = 8)

argo.nat.p <- inv_transY(argo_national$pred)
if(length(index(argo_national09$pred))>0){
  argo.nat.p[index(argo_national09$pred)] <- inv_transY(argo_national09$pred)
}


common_idx_var <- index(ili_national)

var.each <- argo(data = transY(ili_national)[common_idx_var],
                 exogen = transY(stats::lag(ili_national,1))[common_idx_var],
                 N_lag = NULL, N_training = 104, alpha = NA, mc.cores = 8)
var.pred <- inv_transY(var.each$pred)
