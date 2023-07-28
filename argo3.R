## 20220717
## ARGO + SGL
## v1.5
## Edited by S. Ning, based on codes by A. Hussain
## update v1.5: lambda tracking
## update: bug-free, taking care fixed all zero entry 
## alpha must be finite
## 20221208: alpha, nlam change -> fix previous bug on lambda

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

argo3 <- function (data, exogen = xts::xts(NULL), clusterGroups, N_lag = 1:52, N_training = 104, 
          alpha = 1, use_all_previous = FALSE, mc.cores = 1, schedule = list(), lambdas_sgl = NULL) 
{
  if (is.null(schedule$y_gap)) {
    schedule$y_gap <- 1
  }
  if (is.null(schedule$forecast)) {
    schedule$forecast <- 0
  }
  parm <- list(N_lag = N_lag, N_training = N_training, alpha = alpha, 
               use_all_previous = use_all_previous, schedule = schedule)
  if (ncol(data) == 1) {
    data_mat <- matrix(rep(data, nrow(data)), nrow = nrow(data))
    colnames(data_mat) <- as.character(index(data))
    for (i in schedule$y_gap:min(ncol(data_mat), ncol(data_mat) - 
                                 1 + schedule$y_gap)) {
      data_mat[(i - schedule$y_gap + 1):nrow(data_mat), 
               i] <- NA
    }
    data_mat <- xts::xts(data_mat, zoo::index(data))
    data <- data_mat
  }
  lasso.pred <- rep()
  lasso.coef <- list()
  lasso.lambda_idx <- c()
  lasso.lambda_list <- list()
  if (length(exogen) > 0) 
    if (!all(zoo::index(data) == zoo::index(exogen))) 
      stop("error in data and exogen: their time steps must match")
  starttime <- N_training + max(c(N_lag, 0)) + 2 * schedule$y_gap + 
    schedule$forecast
  endtime <- nrow(data)
  clusterGroups_all <- c(N_lag+max(clusterGroups), clusterGroups)
  each_iteration <- function(i) {
    if (use_all_previous) {
      training_idx <- (schedule$y_gap + schedule$forecast + 
                         max(c(N_lag, 0))):(i - schedule$y_gap)
    } else {
      training_idx <- (i - N_training + 1):i - schedule$y_gap
    }
    lagged_y <- sapply(N_lag, function(l) as.numeric(diag(data.matrix(data[training_idx - 
                                                                             l + 1 - schedule$y_gap - schedule$forecast, training_idx - 
                                                                             schedule$forecast]))))
    if (length(lagged_y) == 0) {
      lagged_y <- NULL
    } else {
      colnames(lagged_y) <- paste0("lag_", N_lag + schedule$y_gap - 
                                     1)
    }
    if (length(exogen) > 0 && schedule$y_gap > 0) {
      xmat <- lapply(1:schedule$y_gap, function(l) as.matrix(exogen[training_idx - 
                                                                      schedule$forecast - l + 1, ]))
      xmat <- do.call(cbind, xmat)
      design_matrix <- cbind(lagged_y, xmat)
    } else {
      design_matrix <- cbind(lagged_y)
    }
    
    ### can not take column with all zeros
    idx_nonzero <- which(apply(design_matrix^2, 2, sum) != 0)
    
    y.response <- data[training_idx, i]
    data_t <- list(x = design_matrix[, idx_nonzero], y = as.matrix(y.response))
    ## ignore infinite alpha
    lasso_reg <- cvSGL(data = data_t, index = clusterGroups_all[idx_nonzero], type = "linear", lambdas = lambdas_sgl, alpha = alpha, nlam = length(lambdas_sgl))
    idx_lam <- which.min(lasso_reg$lldiff)
    lasso_reg_refit <- SGL(data = data_t, index = clusterGroups_all[idx_nonzero], type = "linear", lambdas = lasso_reg$lambdas, alpha = alpha, nlam = length(lambdas_sgl))
    coef_cur <- numeric(ncol(design_matrix))
    coef_cur[idx_nonzero] <- lasso_reg_refit$beta[,idx_lam]
    coef_cur <- c(lasso_reg_refit$intercept, coef_cur)
    lasso.coef[[i]] <- as.matrix((coef_cur))
    
    lasso.lambda_idx[i] <- idx_lam
    lasso.lambda_list[[i]] <- lasso_reg$lambdas
    
    lagged_y_next <- matrix(sapply(N_lag, function(l) as.numeric(data[i - 
                                                                        schedule$y_gap + 1 - l, i])), nrow = 1)
    if (length(lagged_y_next) == 0) 
      lagged_y_next <- NULL
    if (length(exogen) > 0 && schedule$y_gap > 0) {
      xmat.new <- lapply(1:(schedule$y_gap), function(l) data.matrix(exogen[i - 
                                                                              l + 1, ]))
      xmat.new <- do.call(cbind, xmat.new)
      newx <- cbind(lagged_y_next, xmat.new)
    } else {
      newx <- lagged_y_next
    }
    #if (is.finite(alpha)) {
    pred1 <- predictSGL(lasso_reg_refit, newX = as.matrix(newx[, idx_nonzero, drop=F]), lam = idx_lam)
    lasso.pred[i] <- ifelse(length(pred1)==0, NA, pred1)
    #}
    #else {
    #  colnames(newx) <- names(design_matrix)
    #  newx <- as.data.frame(newx)
    #  lasso.pred[i] <- predictSGL(lasso_reg_refit, newdata = newx, lam = which.min(lasso_reg$lldiff))
    #}
    
    result_i <- list()
    result_i$pred <- lasso.pred[i]
    result_i$coef <- lasso.coef[[i]]
    result_i$lambda_idx <- lasso.lambda_idx[i]
    result_i$lambda_list<-  lasso.lambda_list[[i]]
    print(index(data)[i])
    #print(result_i)
    return(result_i)
  }
  
  result_all <- parallel::mclapply(starttime:endtime, each_iteration, 
                                   mc.cores = mc.cores, mc.set.seed = FALSE)
  #print(length(result_all))
  lasso.pred[starttime:endtime] <- sapply(result_all, function(x) x$pred)
  lasso.coef <- lapply(result_all, function(x) x$coef)
  data$predict <- NA
  data$predict[1:endtime] <- lasso.pred
  lasso.coef <- do.call("cbind", lasso.coef)
  colnames(lasso.coef) <- as.character(zoo::index(data))[starttime:endtime]
## add lambda idx, list
  lasso.lambda_idx[starttime:endtime] <- sapply(result_all, function(x) x$lambda_idx)
  lasso.lambda_list <- lapply(result_all, function(x) x$lambda_list)
  lasso.lambda_list <- do.call("cbind", lasso.lambda_list)
  colnames(lasso.lambda_list) <- as.character(zoo::index(data))[starttime:endtime]
  argo <- list(pred = data$predict, coef = lasso.coef, parm = parm, 
               lambda.idx = lasso.lambda_idx, lambda.all = lasso.lambda_list)
  class(argo) <- "argo"
  argo
}
