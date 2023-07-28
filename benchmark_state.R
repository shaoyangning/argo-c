####20230529 
### benchmarks: ARGOX, VAR, GFT

###ARGOX

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
path_ws <- "/Users/sn9/Dropbox/Williams/Research/SGL_ARGO/ARGOX_shared2022"
path_data <- file.path(path_ws, "data")
GT_all <- "api_raw_results-2023-01-29"
ili_version <- "ili20230207/"
ili.folder <- file.path(path_data, ili_version)
GT_version <- paste0(GT_all, "/")
out.folder <- paste0(ili.folder, GT_version)
path_out_state <- file.path(out.folder, "States")
load(file.path(path_out_state, "data_ili_GT.Rdata"))

method <- "hc_corr_ave"
alpha_sgl <- 0.95
lag_cur <-  NULL # consistent with argo2
seed <- 1000
idx_precovid_clust <- T
reg.weight.mix <- 1/3
idx_optim_k <- NULL

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
### vanilla ARGOX 2nd step
### Only need to run once, can change to False later
for(idx.GTregmix in c(F, T)){
  # organize step 1 results
  list_argox_state_pred <- list()
  for(reg in 1:length(GT_state)){
    state_curr <- names(GT_state)[reg]
    path_out <- file.path(path_out_state, state_curr)
    prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
    out_name <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
                                           ifelse(idx_precovid_clust, "_precovidclust", ""), 
                                           ifelse(idx.GTregmix, "_regmix", ""),
                                           ".Rdata"))
    load(out_name)
    list_argox_state_pred[[method]][[state_curr]] <- argo.vanilla.p
  }  
  
  argox_state_pred <- do.call(merge, list_argox_state_pred[[method]])
  colnames(argox_state_pred) <- names(list_argox_state_pred[[method]])
  argo.state.pred <- argox_state_pred[,names(ili_state)]
  
  # national & regional
  path_out_reg <- file.path(out.folder, "Regions")
  out_name1 <- file.path(path_out_reg, "argo2_reg_all_logit01_lag0a.Rdata")
  load(out_name1)
  
  ### step 2
  idx.lagstate <- NULL
  truth <-  ili_state
  naive.p <- stats::lag(truth, 1)
  argo1.p <- argo.state.pred
  common_idx <- zoo::index(na.omit(merge(truth, naive.p, argo1.p, 
                                         argo.nat.p,argo.reg.pred)))
  if (ncol(argo.nat.p) == 1) {
    argo.nat.p <- do.call(cbind, lapply(1:ncol(argo1.p), function(i) argo.nat.p))
  }
  colnames(argo.nat.p) <- colnames(truth)
  
  argo.reg.p <- do.call(cbind, lapply(states.abb.all, function(x){argo.reg.pred[,states.info$region[states.info$names==x]]}))
  
  {
    ##Z W
    argo2.pred <- naive.p
    argo2.pred[] <- NA
    var.est.block44 <-  naive.p
    var.est.block44[] <- NA
    
    idx.state.select <- c("US.HI", "US.AK", "US.MT", "US.VT", "US.SD", "US.ND", "US.ME")
    idx.states <- which(!states.abb.all%in%idx.state.select)
    
    #idx.states <- 1:dim(truth)[2]
    print(length(idx.states))
    Z <- (truth - naive.p)[,idx.states]
    W <- merge(stats::lag(Z, 1), (argo1.p - naive.p)[,idx.states], (argo.nat.p - naive.p)[,idx.states], (argo.reg.p - naive.p)[,idx.states])
    ZW <- as.matrix(na.omit(merge(Z, W)))
    
    Z.hat <- Z
    Z.hat[] <- NA
    var.est <- Z
    var.est[] <- NA
    w <- 0.5
    
    #w=0.8
    n.train <- 52*2
    n.all <- nrow(ZW) - n.train
    
    for(it in (n.train+1):nrow(ZW)){
      idx.train <- it-(1:n.train)
      t.now <- rownames(as.matrix(ZW))[it]
      
      mu.hat.z <- colMeans(ZW[idx.train,(1:ncol(Z))])
      mu.hat.w <- colMeans(ZW[idx.train,-(1:ncol(Z))])
      
      Sigma.zz <- var(ZW[idx.train,(1:ncol(Z))]) 
      Sigma.nat <- var(ZW[idx.train, (1:ncol(Z))+ncol(Z)*3]-ZW[idx.train, (1:ncol(Z))+ncol(Z)*0])
      Sigma.step1 <- var(ZW[idx.train, (1:ncol(Z))+ncol(Z)*2]-ZW[idx.train, (1:ncol(Z))+ncol(Z)*0])
      
      D.argo1 <- diag(diag(Sigma.step1))
      D.argo1 <- Sigma.step1
      
      Sigma.reg <- var(ZW[idx.train, (1:ncol(Z))+ncol(Z)*4]-ZW[idx.train, (1:ncol(Z))+ncol(Z)*0])
      
      ##rho: same
      rho.hat <- sum(cor(ZW[idx.train, (1:ncol(Z))],ZW[idx.train, (1:ncol(Z))+ncol(Z)*1])*
                       cor(ZW[idx.train, (1:ncol(Z))]))/sum(cor(ZW[idx.train, (1:ncol(Z))])^2)
      rho.hat <- diag(rho.hat, ncol(Z))
      
      Sigma.zz1 <- rho.hat %*% Sigma.zz  ##cov(t,t-1)
      Sigma.step1.nat <-  matrix(0, ncol=ncol(Z), nrow=ncol(Z))
      Sigma.step1.z1 <- matrix(0, ncol=ncol(Z), nrow=ncol(Z))
      Sigma.step1.z <- matrix(0, ncol=ncol(Z), nrow=ncol(Z))
      Sigma.nat.z1 <- matrix(0, ncol=ncol(Z), nrow=ncol(Z))
      Sigma.nat.z <- matrix(0, ncol=ncol(Z), nrow=ncol(Z))
      Sigma.GT.z <- Sigma.zz + Sigma.step1.z
      Sigma.GT.z1 <- Sigma.zz1 + Sigma.step1.z1
      
      Sigma.ww.str <-  kronecker(matrix(1, 4, 4), Sigma.zz)
      Sigma.ww.str[(1:ncol(Z))+ncol(Z)*1, (1:ncol(Z))+ncol(Z)*1] <- Sigma.zz + D.argo1 + Sigma.step1.z + t(Sigma.step1.z)
      Sigma.ww.str[(1:ncol(Z))+ncol(Z)*2, (1:ncol(Z))+ncol(Z)*2] <-  Sigma.zz + Sigma.nat + Sigma.nat.z + t(Sigma.nat.z)
      Sigma.ww.str[(1:ncol(Z))+ncol(Z)*3, (1:ncol(Z))+ncol(Z)*3] <-  Sigma.zz + Sigma.reg
      Sigma.ww.str[(1:ncol(Z)), (1:ncol(Z)) + ncol(Z)*1] <- t(Sigma.GT.z1) #t(Sigma.zz1) + t(Sigma.step1.z1)
      Sigma.ww.str[(1:ncol(Z)), (1:ncol(Z)) + ncol(Z)*2] <- t(Sigma.zz1) + t(Sigma.nat.z1)
      Sigma.ww.str[(1:ncol(Z)), (1:ncol(Z)) + ncol(Z)*3] <- t(Sigma.zz1) 
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*1, (1:ncol(Z))] <- Sigma.GT.z1 #Sigma.zz1 + Sigma.step1.z1
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*2, (1:ncol(Z))] <- Sigma.zz1 + Sigma.nat.z1
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*3, (1:ncol(Z))] <- Sigma.zz1 
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*2, (1:ncol(Z)) + ncol(Z)*1] <- Sigma.zz + t(Sigma.step1.z) + (Sigma.nat.z) + t(Sigma.step1.nat)
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*1, (1:ncol(Z)) + ncol(Z)*2] <- Sigma.zz + Sigma.step1.z + t(Sigma.nat.z) + (Sigma.step1.nat)
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*3, (1:ncol(Z)) + ncol(Z)*1] <- t(Sigma.GT.z)#Sigma.zz + t(Sigma.step1.z)  
      Sigma.ww.str[(1:ncol(Z)) + ncol(Z)*1, (1:ncol(Z)) + ncol(Z)*3] <- Sigma.GT.z#Sigma.zz + Sigma.step1.z
      
      Sigma.zw.str <- cbind(Sigma.zz1, t(Sigma.GT.z), Sigma.zz + t(Sigma.nat.z), Sigma.zz)
      
      D.ww <- diag(diag(var(ZW[idx.train,-(1:ncol(Z))])))
      
      Z.hat[t.now, ] <- (w*Sigma.zw.str) %*% solve(w*Sigma.ww.str + (1-w)*D.ww) %*% (ZW[it,-(1:ncol(Z))] - mu.hat.w) + mu.hat.z
      var.est[t.now, ] <- diag(Sigma.zz - (w*Sigma.zw.str) %*% solve(w*Sigma.ww.str + (1-w)*D.ww) %*% t((w*Sigma.zw.str)))
    }
    argo2.pred[, idx.states] <- Z.hat + naive.p[, idx.states]
    
    argo2.pred.block44 <- argo2.pred
    idx.state.select44 <- idx.state.select
    var.est.block44[, idx.states] <- var.est
  } #44block
  {
    argo2.pred[] <- NA
    var.est.ind <- naive.p
    var.est.ind[] <- NA
    
    for(idx.state.select in states.abb.all){
      idx.states <- which(states.abb.all==idx.state.select)
      
      #print(length(idx.states))
      Z <- (truth - naive.p)[, idx.states, drop=F]
      W <- merge(stats::lag(Z, 1), (argo1.p - naive.p)[, idx.states, drop=F], (argo.nat.p - naive.p)[, idx.states, drop=F], (argo.reg.p - naive.p)[, idx.states, drop=F])
      
      #No regional
      W <- merge(stats::lag(Z, 1), (argo1.p - naive.p)[, idx.states, drop=F], (argo.nat.p - naive.p)[, idx.states, drop=F])
      ZW <- as.matrix(na.omit(merge(Z, W)))
      
      Z.hat <- Z
      Z.hat[] <- NA
      var.est <- Z
      var.est[] <- NA
      
      w <- 0.5
      
      #w=0.8
      n.train <- 52*2
      n.all <- nrow(ZW) - n.train
      
      
      for(it in (n.train+1):nrow(ZW)){
        idx.train <- it-(1:n.train)
        t.now <- rownames(as.matrix(ZW))[it]
        Sigma.zz <- cov(ZW[idx.train, (1:ncol(Z)), drop=F], ZW[idx.train, (1:ncol(Z)), drop=F])
        mu.hat.z <- colMeans(ZW[idx.train, (1:ncol(Z)), drop=F])
        mu.hat.w <- colMeans(ZW[idx.train, -(1:ncol(Z)), drop=F])
        Sigma.ww.str <- cov(ZW[idx.train, -(1:ncol(Z)), drop=F], ZW[idx.train, -(1:ncol(Z)), drop=F]) 
        Sigma.zw.str <- cov(ZW[idx.train, (1:ncol(Z)), drop=F], ZW[idx.train, -(1:ncol(Z)), drop=F])
        
        D.ww <- diag(diag(var(ZW[idx.train, -(1:ncol(Z)), drop=F])))
        
        Z.hat[t.now, ] <- (w * Sigma.zw.str) %*% solve(w*Sigma.ww.str + (1-w)*D.ww) %*% t(ZW[it,-(1:ncol(Z)), drop=F] - mu.hat.w) + mu.hat.z
        var.est[t.now, ] <- (Sigma.zz - (w*Sigma.zw.str) %*% solve(w*Sigma.ww.str + (1-w)*D.ww) %*% t((w*Sigma.zw.str)))
        
      }
      argo2.pred[, idx.states] <- Z.hat + naive.p[, idx.states]
      var.est.ind[, idx.states] <- var.est
    }
    argo2.pred.ind <- argo2.pred
  } #individual
  
  if(idx.GTregmix){
    gt.out <- paste0("GTregmix", round(reg.weight.mix, 2), "lag", length(idx.lagstate))
  }
  if(!idx.GTregmix){
    gt.out <- paste0("GTstate", "lag", length(idx.lagstate))
  }
  
  file.name <- paste0("argox_states_all_logit01_", gt.out, "a_step2_44b.Rdata")
  
  #print(index(argo2.pred.block44[!is.na(argo2.pred.block44[,1]), ]))
  save(argo2.pred.ind, var.est.ind,
    argo2.pred.block44, idx.state.select44, var.est.block44,
    file=file.path(path_out_state, file.name))
}


idx.state.select <- idx.state.select44
idx.GTregmix <- T
gt.out <- paste0("GTregmix", round(reg.weight.mix, 2), "lag", length(idx.lagstate))
file.name <- paste0("argox_states_all_logit01_", gt.out, "a_step2_44b.Rdata")
load(file.path(path_out_state, file.name))

argox_state_pred_final <- argo2.pred.block44
argox_state_var_est_final <- var.est.block44

idx.GTregmix <- F
gt.out <- paste0("GTstate", "lag", length(idx.lagstate))
file.name <- paste0("argox_states_all_logit01_", gt.out, "a_step2_44b.Rdata")
load(file.path(path_out_state, file.name))

argox_state_pred_final[, idx.state.select] <- argo2.pred.ind[, idx.state.select]
argox_state_var_est_final[, idx.state.select] <- var.est.ind[, idx.state.select]

file.name <- paste0("argox_states_all_logit01_final.Rdata")
save(argox_state_pred_final, argox_state_var_est_final, ili_state, file=file.path(path_out_state, file.name))

###### var
idx.GTregmix <- F
list_argo_var_state_pred <- list()
for(reg in 1:length(GT_state)){
  state_curr <- names(GT_state)[reg]
  ili_curr <- ili_state[, state_curr]
  GT_curr <- GT_state[[state_curr]]
  path_out <- file.path(path_out_state, state_curr)
  prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
  out_name <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
                                         ifelse(idx_precovid_clust, "_precovidclust", ""), 
                                         ifelse(idx.GTregmix, "_regmix", ""),
                                         ".Rdata"))
  load(out_name)
  list_argo_var_state_pred[[method]][[reg]] <- var.pred
}  

var_sgl_state_pred <- do.call(merge, list_argo_var_state_pred[[method]])
colnames(var_sgl_state_pred) <- names(GT_state)
save(var_sgl_state_pred,
     file=file.path(path_out_state, "var_state_all_lag1.Rdata"))



