## 20221021 summary argo-sgl
## S. Ning
## argo + sgl, national
## input
# ili_national: ili
# GT_national: GT
# ready to run


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

pred <- merge(ili_national,
              argo_sgl_pred_all,
              argo.nat.p,
              ## included other concatenated results as needed
              stats::lag(ili_national,1))
colnames(pred) <- c("CDC.data",  method_print, "Argo", "naive")

tab <- summary_argo(pred, colnames(pred), colnames(pred),
                    zoom_periods, eval.period)

### CI by bootstrap
pred1 <- na.omit(pred[,1:2])
err1 <- pred1[,1] - pred1[,2]
nlag <- 104
nboot <- 500
err1.bound1 <- err1[,1]
err1.bound1[,] <- NA
for(i in 2:dim(err1)[1])  {
  err1.bound1[i,] <- sqrt(mean(sample(err1[max(1,i-nlag):(i-1),], size = nboot, replace = T)^2))
}
  #sapply(2:dim(err1)[1], function(i){sample(err1[max(1,i-nlag):(i-1),], size = nboot, replace = T)})

err1.bound1 <- lapply(err1.boot, function(x){sqrt(mean(x^2))})
err1.bound1 <- unlist(err1.bound1)
err1.bound1 <- xts(err1.bound1, order.by = index(pred1)[-1])

ci1 <- cbind(pred1[-1, 2] + 1.96 * err1.bound1, pred1[-1, 2] - 1.96 * err1.bound1)
idx_ci1 <- pred1[-1, 1] > ci1[,2] & pred1[-1, 1] < ci1[,1]
mean(idx_ci1["2009-03-29/2020-02-29"])
