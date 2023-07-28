## ARGO-SGL: cluster step -- state
## pipeline + K determine
## 20230516
## S. Ning
## Input: idx_optim_k, terms09, method

## state only
common_idx <- index(merge(ili_national, GT_national, ili_state, all=FALSE))
#terms <- colnames(GT_national)

if(!idx.GTregmix){
  GT_cur <- GT_state
} else{
  GT_cur <- GT_state
  for(j in 1:length(GT_state)){
    state_curr <- names(GT_state)[j]
    region_curr <- states.info$region[states.info$names==state_curr]
    GT.state.curr <- GT_state[[state_curr]]
    GT.reg.curr <- GT_regional[[paste0("Region.", region_curr)]]
    GT_cur[[state_curr]] <-  (1-reg.weight.mix)*GT.state.curr + reg.weight.mix* GT.reg.curr
  }
}


for(reg in states.abb.all){
  path_out <- file.path(path_out_state, reg)
  if(!exists(path_out)){
    dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  }  
  ## filter out all zero entries
  pre09 <- t(scale(GT_cur[[reg]][pre_2009, terms09]))
  idx_pre09 <- apply(is.na(pre09), 1, sum) == 0
  post10 <- t(scale(GT_cur[[reg]][post_2010, ]))
  idx_post10 <- apply(is.na(post10), 1, sum) == 0
  covid_GT <- GT_cur[[reg]][covid_train,]
  covid <- t(scale(covid_GT))
  idx_covid <- apply(is.na(covid), 1, sum) == 0
  ## input to GT_cluster.R
  list_data <- list(pre09 = pre09[idx_pre09, ], post10 = post10[idx_post10, ], covid = covid, covid_train = covid[idx_covid, ])
  list_idx_na <- list(idx_pre09, idx_post10, idx_covid)
  #print(sum(idx_pre09))
  source("GT_cluster.R")
}


