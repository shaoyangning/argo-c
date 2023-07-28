## 20221031
## S. Ning
## argo + sgl, general
## input
# ili_curr: ili
# GT_curr: GT
# lag_cur: lag
# cluster_post10, cluster_pre09
# ready to run

source("argo3.R")

##############################################################################
## shared dates between GT and ili data
common_idx <- index(merge(ili_curr, GT_curr, all=FALSE))
terms <- colnames(GT_curr)

## fix random seed
seed <- 1000
set.seed(seed)

# post 2010
idx.nat1 <- c(index(GT_curr)[1]-(52:1)*7, index(GT_curr))
GT_curr1 <- xts(matrix(nrow=length(idx.nat1), ncol=ncol(GT_curr)), 
                    order.by = as.Date(idx.nat1))
colnames(GT_curr1) <- colnames(GT_curr)
GT_curr1[index(GT_curr), ] <- GT_curr
common_idx_nat <- c(common_idx[1]-(52:1)*7, common_idx)
if(length(lag_cur)==52){
  common_idx_curr <- common_idx_nat
} else{
  common_idx_curr <- common_idx
}

# pre-2009 set of predictors
common_idx09 <- common_idx[common_idx < as.Date("2010-05-22")]
terms09 <- intersect(colnames(GT_curr), colnames(load_data()$GC09))
common_idx09_nat <- c(common_idx09[1]-(52:1)*7, common_idx09)

if(length(lag_cur)==52){
  common_idx09_curr <- common_idx09_nat 
} else{
  common_idx09_curr <- common_idx09
}

## covid
covid_period <- "2020-03-14/2021-04-03"

common_covid <- common_idx[common_idx >= as.Date("2020-03-14")]
common_covid_nat <- c(common_covid[1]-((2+52*3):1)*7, common_covid)

if(length(lag_cur)==52){
  common_covid_curr <- common_covid_nat 
} else{
  common_covid_curr <- c(common_covid[1]-((2+52*2):1)*7, common_covid)
}

# post2010
set.seed(seed)
argo_SGL_10 <- argo3(data = transY(ili_curr[common_idx_curr,]),
                     exogen = log(GT_curr1[common_idx_curr, terms]+1), 
                     clusterGroups = cluster_post10$cluster, lambdas_sgl = lambdas_all,
                     alpha = alpha_sgl,
                     N_lag = lag_cur,
                     mc.cores = 8)

print("post10 done")

# pre-2009
argo_SGL_09 <- NULL
if(length(common_idx09) > 0){
  set.seed(seed)
  argo_SGL_09 <- argo3(data = transY(ili_curr[common_idx09_curr,]),
                       exogen = log(GT_curr1[common_idx09_curr, terms09]+1), 
                       clusterGroups = cluster_pre09$cluster, lambdas_sgl = lambdas_all,
                       alpha = alpha_sgl,
                       N_lag = lag_cur,
                       mc.cores = 8)
  
  print("pre09 done")
}


# covid_train
set.seed(seed)

argo_SGL_covid <- argo3(data = transY(ili_curr[common_covid_curr,]),
                        exogen = log(GT_curr1[common_covid_curr, terms]+1), 
                        clusterGroups = cluster_covid$cluster, lambdas_sgl = lambdas_all,
                        alpha = alpha_sgl,
                        N_lag = lag_cur,
                        mc.cores = 8)

print("covid done")

## concatenate
argo_sgl_final <- inv_transY(argo_SGL_10$pred)
if(length(index(argo_SGL_09$pred))>0){
  argo_sgl_final[index(argo_SGL_09$pred)] <- inv_transY(argo_SGL_09$pred)
}

if(idx_precovid_clust){
  argo_sgl_final[covid_period] <- inv_transY(argo_SGL_covid$pred[covid_period])
}

## vanilla argo
set.seed(seed)
argo_vanilla <- argo(data = transY(ili_curr[common_idx_curr]),
                      exogen = log(GT_curr1[common_idx_curr, terms]+1),
                      N_lag = lag_cur,
                      mc.cores = 8)

argo_vanilla09 <- NULL
if(length(common_idx09) > 0){
argo_vanilla09 <- argo(data = transY(ili_curr[common_idx09_curr]),
                        exogen = log(GT_curr1[common_idx09_curr, terms09] + 1),
                        N_lag = lag_cur,
                        mc.cores = 8)
}

argo.vanilla.p <- inv_transY(argo_vanilla$pred)
if(length(index(argo_vanilla09$pred))>0){
  argo.vanilla.p[index(argo_vanilla09$pred)] <- inv_transY(argo_vanilla09$pred)
}


common_idx_var <- index(ili_curr)

var.each <- argo(data = transY(ili_curr)[common_idx_var],
                 exogen = transY(stats::lag(ili_curr,1))[common_idx_var],
                 N_lag = NULL, N_training = 104, alpha = NA, mc.cores = 8)
var.pred <- inv_transY(var.each$pred)
