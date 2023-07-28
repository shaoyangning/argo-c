## ARGO-SGL: cluster step
## 20230207 generalize codes to all loc
### input: 
#### path_out: output folder name
#### list_data: list of datasets (already transposed, unfiltered to zeros)
#### idx_optim_k
### output: plots and clustering results

## 20221101 covid period clustering
## S. Ning
## 20221025 organizing cluster results by Ahmed
##

options(echo=TRUE)
library(cluster)
library(factoextra)
library(fossil)
library(xts)

## 
method_all <- c("hc_corr_ave" #, "hc_corr_comp", "hc_corr_single", "kmeans" , "pam",
)
periods_all <- c("pre09", "post10", #"covid",
                 "covid_train")
k_stats_all <- c("elbow", "gap", "silhouette")


seed <- 1000

## clustering periods
# pre 2009 search terms clustered on data from Jan 10 2010 to March 28 2009
pre_2009 <- seq(as.Date("2004-01-10"), as.Date("2009-03-28"), by = "week")
# post-2010 search terms clustered on data from Jan 10 2004 to May 22 2010
post_2010 <- seq(as.Date("2004-01-10"), as.Date("2010-05-22"), by = "week")
# covid clustering period
covid_train <- seq(as.Date("2015-03-14"), as.Date("2020-03-14"), by = "week")


## clustering
hc_average_pearson <- function (x, k){ hcut(x,k,hc_method = "average", hc_metric= "pearson")}
hc_average_euclidian<- function (x, k){ hcut(x,k,hc_method = "average", hc_metric= "euclidian")}
list_clust_fn <- list(kmeans = kmeans, hc_corr_ave = hc_average_pearson, 
                      hc_euc_ave = hc_average_pearson, pam = pam)



list_stats <- list(elbow = "wss", gap = "gap", silhouette = "silhouette")


##
list_optim_k_results <- list()
#load("GT_clust_optim_k.Rdata")

## optim K
if(F){
  for(method in method_all){
    for(stat in k_stats_all){
      for(period in periods_all){
        if(period == "pre09"){
          kmax_cur <- 65
        } else {
          kmax_cur <- 140
        }
        set.seed(seed)
        list_optim_k_results[[method]][[period]][[stat]] <- 
          fviz_nbclust(list_data[[period]], 
                       list_clust_fn[[method]], method = list_stats[[stat]], k.max = kmax_cur) 
        print(paste(method, stat, period, "done"))
      }
    }
  }
}


## plot
if(F){
  pdf(file.path(path_out, "optim_k.pdf"), width = 20, height =10)
  for(period in periods_all){
    for(method in method_all){
      for(stat in k_stats_all){
        fig_title <- paste(period, method, stat)
        print(list_optim_k_results[[method]][[period]][[stat]] 
              + ggtitle(fig_title))
      }
    }
  }
  dev.off()
}

#load(file.path(path_out, "GT_clust_optim_k.Rdata"))


list_optim_k <- list()
## saving up optim K
out_name_clust <- paste0("GT_clust_optim_k", ".Rdata")
if(T){
  ##hc_corr_ave
  list_optim_k[["hc_corr_ave"]] <- list()
  list_optim_k[["hc_corr_ave"]][["post10"]] <- list(kmax = 140, kopt = 45, sil = 2, elbow = 19, gap = 88)
  list_optim_k[["hc_corr_ave"]][["pre09"]] <- list(kmax = 65, kopt=53, gap = 2, sil = 53)
  ## in## post 10: elbow 19 both 45 gap 88
  list_optim_k[["hc_corr_ave"]][["covid"]] <- list(kmax = 140, kopt = 30, gap = 1, sil = 2, elbow = 30)
  ## elbow 11 = 18 > 30,
  list_optim_k[["hc_corr_ave"]][["covid_train"]] <- list(kmax = 140, kopt = 15, gap = 15, sil = 2, elbow = 15)
  ## elbow 14 >>,
  
  ##hc_euc_ave
  list_optim_k[["hc_euc_ave"]] <- list()
  list_optim_k[["hc_euc_ave"]][["post10"]] <- list(kmax = 140, kopt = 45, gap = 88, sil = 2, elbow= 35)
  list_optim_k[["hc_euc_ave"]][["pre09"]] <- list(kmax = 65, kopt=53, gap = 1, sil = 53)
  ###ea3 46 53
  ### same as corr_ave?
  
  ##pam
  list_optim_k[["pam"]] <- list()
  list_optim_k[["pam"]][["post10"]] <- list(kmax = 140, kopt = 70, gap = 140, sil = 70)
  list_optim_k[["pam"]][["pre09"]] <- list(kmax = 65, kopt = 48, gap = 1, sil = 48)
  list_optim_k[["pam"]][["covid"]] <- list(kmax = 140, kopt = 2, gap = NA, sil = 2, elbow = 1)
  list_optim_k[["pam"]][["covid_train"]] <- list(kmax = 140, kopt = 13, gap = 13, sil = 2, elbow = 13)
  
}

if(length(idx_optim_k) > 0){
  out_name_clust <- paste0(paste("GT_clust_optim_k", idx_optim_k, sep = "_"), ".Rdata")
  print(out_name_clust)
  list_optim_k <- list()
  file_optim_k <- file.path(path_out_cur, paste0(idx_optim_k, ".csv"))
  tab_optim_k <- read.csv(file_optim_k)
  rownames(tab_optim_k) <- tab_optim_k[,1]
  for(period in periods_all){
    list_optim_k[["hc_corr_ave"]][[period]] <- list(kopt = tab_optim_k[reg, period]) 
  }
} 

## recluster and save cluster results
## cluster_harm: harmonized clusters, zero/NA removed
## cluster_back: fill-back search terms, singleton clusters
list_cluster_results <- list()
if(T){  
  for(period in periods_all){
    for(method in method_all){
      set.seed(seed)
      list_cluster_results[[method]][[period]] <- 
        list_clust_fn[[method]](list_data[[period]],  
                                min(dim(list_data[[period]])[1]-1, list_optim_k[[method]][[period]]$kopt))
      if(period == "pre09"){
        ## update
        cluster_back <- numeric(length(terms09))
        names(cluster_back) <- terms09
      } else{
        cluster_back <- numeric(length(colnames(GT_cur[[reg]])))
        names(cluster_back) <- colnames(GT_cur[[reg]])
      }
      ## save clusters with harmonized terms
      cluster_harm <- (list_cluster_results[[method]][[period]]$cluster)
      list_cluster_results[[method]][[period]]$cluster_harm <- cluster_harm
      ## fill back group id (as singleton) for all zero entries
      if(length(cluster_back)!=length(cluster_harm)){
          cluster_back[names(cluster_harm)] <- list_cluster_results[[method]][[period]]$cluster
          cluster_back[cluster_back==0] <- (length(unique(cluster_harm))+1)+(1:sum(cluster_back==0))
          list_cluster_results[[method]][[period]]$cluster <- cluster_back
      }
      print(paste(method, period, "done"))
      }
    }
}

if(exists("idx.GTregmix")){
  out_name_clust <- paste0("GT_clust_optim_k", 
                           ifelse(length(idx_optim_k)==0, "", paste0("_", idx_optim_k)),
                           ifelse(idx.GTregmix, "_regmix", ""), 
                           ".Rdata")
  
}

save(list_optim_k, list_optim_k_results, list_cluster_results, file = file.path(path_out, out_name_clust))
