## 20230425 summary argo-sgl regional
## S. Ning
## input
# ili_national: ili
# GT_national: GT
# ready to run

## read in argo sgl regional data
prefix_out1 <- paste("argo_sgl_reg", method, "lambda2", "alpha", alpha_sgl, "lag", length(lag_cur), sep = "_")
out_name1 <- file.path(path_out_reg, paste0(paste(prefix_out1, idx_optim_k, sep = ifelse(length(idx_optim_k)==0, "", "_")), 
                                        ifelse(idx_precovid_clust, "_precovidclust", ""), ".Rdata"))
load(out_name1)

## argo2 
load(file.path(path_out_reg, "argo2_reg_all_logit01_lag0a.Rdata"))
load(file.path(path_out_reg, "var_reg_all_lag1.Rdata"))

## GFT
population_file <- file.path(path_data, "Population.csv") 
gft_file <- file.path(path_data, "GFT.txt")
source("ARGOX_shared2022/argo2_states_function.R")
list_gft <- load_gft_data1(population.file=population_file, gft.file=gft_file)


## Zoom periods
zoom_periods <- c("2009-03-29/2020-02-29", # overall pre-covid
                  "2009-03-29/2023-01-28", # overall including covid
                  "2020-03-07/2023-01-28", # Covid-period
                  "2014-10-11/2020-03-21", # post10
                  "2009-10-10/2010-05-22", # year-wise
                  "2010-10-03/2011-05-22",
                  "2011-10-02/2012-05-20",
                  "2012-09-30/2013-05-19",
                  "2013-10-05/2014-05-17",
                  "2014-10-04/2015-05-23",
                  "2015-10-10/2016-05-22",
                  "2016-10-08/2017-05-21",
                  "2017-10-07/2018-05-19",
                  "2018-10-06/2019-05-18",
                  "2019-10-06/2020-05-16",
                  "2020-10-05/2021-05-22",
                  "2021-10-09/2022-05-21",
                  "2022-10-08/2023-01-28",
                  "2019-10-06/2020-03-21",
                  "2020-03-21/2020-05-16")

eval.period <- "2009-03-29/2015-07-11"


reg.method.p1 <- sapply(1:10, function(region.id){
  pred.region <- 
    merge(ili_regional[,region.id], argo_sgl_reg_result$twostep[,region.id],
          argo2_reg_result$twostep[, region.id], var_sgl_reg_pred[, region.id],
          list_gft$GFT[, region.id],
          stats::lag(ili_regional[,region.id], 1))
  colnames(pred.region) <- c("CDC.data", "ARGO-SGL", "Argo2", "VAR1", "GFT", "Naive")
  pred.region <- pred.region[is.finite(pred.region$Argo2),]
  data.matrix(pred.region)
}, simplify = "array")

#print(dim(reg.method.p1))

tab.allregion <- mclapply(1:10, function(region.id){
  tab <- summary_argo(xts(reg.method.p1[,,region.id], as.Date(dimnames(reg.method.p1)[[1]])),
                      dimnames(reg.method.p1)[[2]], dimnames(reg.method.p1)[[2]],
                      zoom_periods, eval.period)
}, mc.cores = 2)

tab.allregion <- sapply(c("rmse","abse","rmspe","mape","corr", "corr_diff"), function(type){
  sapply(1:10, function(region.id){
    if(type=="rmse"){
      tab.allregion[[region.id]][[type]]^2
    }else{
      tab.allregion[[region.id]][[type]]  
    }
  }, simplify = "array")
}, simplify = "array")

dimnames(tab.allregion)[[3]] <- paste0("Region", 1:10, sep=".")

eval.period1 <- "2009-03-29/2023-01-28"
tt <- sqrt(rbind(t(tab.allregion[,eval.period1,,"rmse"]), 
            colMeans(t(tab.allregion[,eval.period1,,"rmse"]),na.rm=T)))
#View(sqrt(tt))
eval.period1 <- "2009-03-29/2020-02-29"
tt1_rmse <- sqrt(rbind(t(tab.allregion[,eval.period1,,"rmse"]), 
            colMeans(t(tab.allregion[,eval.period1,,"rmse"]),na.rm=T)))


tt1_mae <- (rbind(t(tab.allregion[,eval.period1,,"abse"]), 
                       colMeans(t(tab.allregion[,eval.period1,,"abse"]),na.rm=T)))

tt1_corr  <- (rbind(t(tab.allregion[,eval.period1,,"corr"]), 
                        colMeans(t(tab.allregion[,eval.period1,,"corr"]),na.rm=T)))

##View((tt1_rmse))



eval.period1 <- "2009-03-29/2015-07-11"
tt2_rmse <- sqrt(rbind(t(tab.allregion[,eval.period1,,"rmse"]), 
                       colMeans(t(tab.allregion[,eval.period1,,"rmse"]),na.rm=T)))

######################################################
### CI
var.est <- sapply(1:dim(argo_sgl_reg_result$projection.mat)[3], 
                  function(it) {
                    stopifnot(argo_sgl_reg_result$projection.mat[, , it] == 
                                argo_sgl_reg_result$sigma_zwzw.structured[1:10, -(1:10), 
                                                                   it] %*% solve(argo_sgl_reg_result$sigma_zwzw.structured[-(1:10), 
                                                                                                                    -(1:10), it] + diag(diag(argo_sgl_reg_result$sigma_zwzw.empirical[-(1:10), 
                                                                                                                                                                               -(1:10), it]))))
                    argo_sgl_reg_result$sigma_zwzw.structured[1:10, 1:10, it] - 
                      0.5 * (argo_sgl_reg_result$sigma_zwzw.structured[1:10, 
                                                                -(1:10), it] %*% solve(argo_sgl_reg_result$sigma_zwzw.structured[-(1:10), 
                                                                                                                          -(1:10), it] + diag(diag(argo_sgl_reg_result$sigma_zwzw.empirical[-(1:10), 
                                                                                                                                                                                     -(1:10), it]))) %*% argo_sgl_reg_result$sigma_zwzw.structured[-(1:10), 
                                                                                                                                                                                                                                            1:10, it])
                  }, simplify = "array")
dimnames(var.est)[[3]] <- as.character(index(na.omit(argo_sgl_reg_result$twostep)))
err.realized <- na.omit(argo_sgl_reg_result$truth - argo_sgl_reg_result$twostep)

err.realized <- err.realized["2009-03-29/2020-02-29"]
var.est <- var.est[,,dimnames(var.est)[[3]] %in% as.character(index(err.realized))]

interval.coverage <- list()
for (region.id in 1:10) {
  interval.coverage[[region.id]] <- 
    c(mean( abs( data.matrix(err.realized)[, region.id]/
                   sqrt(var.est[region.id, region.id, ])) <= 1.96), 
      var(data.matrix(err.realized)[, region.id]/sqrt(var.est[region.id, region.id, ])))
}                                                                                                            

interval.coverage <- do.call(rbind, interval.coverage)
rownames(interval.coverage) <- colnames(err.realized)


