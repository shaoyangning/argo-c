## 20230328
## S. Ning
## argo + sgl, regional
## input
# seed
# alpha_sgl
# lag_cur
# lambdas_all

source("argo3.R")

## shared dates between GT and ili data

## fix random seed

# reglag = NULL # consistent with Argo2
for(reg in 1:10){
  ili_curr <- ili_regional[, reg]
  GT_curr <- GT_regional[[reg]]
  ## read in cluster info
  path_out <- file.path(path_out_reg, paste0("Region", reg))
  if(!file.exists(path_out)){
    dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  }  
  out_name_clust <- paste0(paste("GT_clust_optim_k", idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), ".Rdata")
  load(file.path(path_out, out_name_clust))
  cluster_pre09 <- list_cluster_results[[method]]$pre09
  cluster_post10 <- list_cluster_results[[method]]$post10
  cluster_covid <- list_cluster_results[[method]]$covid_train
  
  ## argo-sgl
  set.seed(seed)
  prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
  out_name_precovid <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), "_precovidclust", ".Rdata"))
  if(!file.exists(out_name_precovid)){
    source("argo_sgl_script.R")
  } else{
    load(out_name_precovid)
    if(!idx_precovid_clust){
      argo_sgl_final <- inv_transY(argo_SGL_10$pred)
      if(length(index(argo_SGL_09$pred))>0){
        argo_sgl_final[index(argo_SGL_09$pred)] <- inv_transY(argo_SGL_09$pred)
      }
    }
  }
  #out_name <- file.path(path_out, paste0(prefix_out, ".Rdata"))
  out_name <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata"))
  save(seed, lambdas_all, method, cluster_pre09, cluster_post10, cluster_covid, 
       argo_SGL_09, argo_SGL_10, argo_SGL_covid,
       alpha_sgl, var.pred,
       argo.vanilla.p, argo_sgl_final, file = out_name)
  gc()
}

## argo regional step 2
## organize step 1 argo-sgl results

## national
path_out_nat <- file.path(out.folder, "national")
prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl, sep = "_")
#prefix_out <- paste0(prefix_out, ifelse(idx_precovid_clust, "_precovidclust", ""))
out_name <- paste0(prefix_out, ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata")
load(file.path(path_out_nat, out_name))


## regional
list_argo_sgl_reg_pred <- list()

for(reg in 1:10){
  path_out <- file.path(path_out_reg, paste0("Region", reg))
  prefix_out <- paste("argo_sgl", method, "lambda2", "alpha", alpha_sgl,  "lag", length(lag_cur), sep = "_")
  out_name <- file.path(path_out, paste0(paste(prefix_out, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata"))
  load(out_name)
  list_argo_sgl_reg_pred[[method]][[reg]] <- argo_sgl_final
}  


argo_sgl_reg_pred <- do.call(merge, list_argo_sgl_reg_pred[[method]])
colnames(argo_sgl_reg_pred) <- paste0("Region", ".", 1:10)

ili_regional2003post <- ili_regional["2003/"]
argo_sgl_reg_result <- argo2(ili_regional2003post, argo_sgl_reg_pred, argo_nat_sgl)

#prefix_out1 <- paste0(paste("argo_sgl_reg", method, "lambda2", "alpha", alpha_sgl, ".Rdata")
 
prefix_out1 <- paste("argo_sgl_reg", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
out_name1 <- file.path(path_out_reg, paste0(paste(prefix_out1, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
                                       ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata"))
                                           
#out_name1 <- file.path(path_out_reg, prefix_out1)
save(argo_sgl_reg_pred, argo_nat_sgl, argo_sgl_reg_result, file = out_name1)
