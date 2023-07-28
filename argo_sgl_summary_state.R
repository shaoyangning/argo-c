## 20230529 summary argo-sgl state
## S. Ning
## input
# ili_state: ili
# states.abb.all: states names abv
# ready to run
## ARGO-C
file.name <- paste0("argo_sgl_states_all_logit01_final.Rdata")
load(file.path(path_out_state, file.name))
## VAR
load(file.path(path_out_state, "var_state_all_lag1.Rdata"))
## ARGOX
load(file.path(path_out_state, "argox_states_all_logit01_final.Rdata"))

## GFT
population_file <- file.path(path_data, "Population.csv") 
gft_file <- file.path(path_data, "GFT.txt")
source("ARGOX_shared2022/argo2_states_function.R")
list_gft <- load_gft_data1(population.file=population_file, gft.file=gft_file)
gft_state <- list_gft$GFT_state
colnames(gft_state) <- states.info$names[match(colnames(gft_state), states.info$states)]
gft_state <- gft_state[,colnames(ili_state)]


reg.method.p <- sapply(1:length(states.abb.all), function(region.id){
  pred.region <- 
    merge(ili_state[,region.id],
          na.omit(argo_sgl_state_pred_final[,region.id]), 
          argox_state_pred_final[,region.id], 
          gft_state[,region.id],
          var_sgl_state_pred[,region.id],
          stats::lag(ili_state,1)[,region.id])
  colnames(pred.region) <- c("CDC.data", "ARGO-C", "ARGOX", "GFT", 
                             "var1","naive")
  pred.region <- pred.region[is.finite(pred.region$ARGOX),]
  data.matrix(pred.region)
}, simplify = "array")
dimnames(reg.method.p)[[3]] <- colnames(ili_state)

zoom_periods <- c("2014-10-11/2020-02-29", # overall pre-covid
                  "2014-10-11/2023-01-28", # overall including covid
                  "2020-03-07/2023-01-28", # Covid-period
                  "2014-10-11/2020-03-21", # post10
                  "2014-10-04/2015-05-23", # year wise
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

eval.period <- "2014-10-11/2015-07-11"

tab.allregion <- mclapply(1:length(states.abb.all), function(region.id){
  tab <- summary_argo(xts(reg.method.p[,,region.id], as.Date(dimnames(reg.method.p)[[1]])),
                      dimnames(reg.method.p)[[2]], dimnames(reg.method.p)[[2]],
                      zoom_periods, eval.period)
}, mc.cores = 2)

tab.allregion <- sapply(c("rmse","abse","rmspe","mape","corr", "corr_diff"), function(type){
  sapply(1:length(states.abb.all), function(region.id){
    if(type=="rmse"){
      tab.allregion[[region.id]][[type]]^2
    }else{
      tab.allregion[[region.id]][[type]]  
    }
  }, simplify = "array")
}, simplify = "array")

dimnames(tab.allregion)[[3]] <- colnames(ili_state)
## rmse = mse


eval.period1 <- "2014-10-11/2020-02-29"
tt1_mse <- sqrt(rbind(t(tab.allregion[,eval.period1,,"rmse"]), 
                      colMeans(t(tab.allregion[,eval.period1,,"rmse"]),na.rm=T)))
sum(tt1_mse[,2] <= tt1_mse[,3])
eval.period1 <- "2014-10-11/2023-01-28"
tt1_mse1 <- sqrt(rbind(t(tab.allregion[,eval.period1,,"rmse"]), 
                       colMeans(t(tab.allregion[,eval.period1,,"rmse"]),na.rm=T)))

sum(tt1_mse1[,2] <= tt1_mse1[,3])


