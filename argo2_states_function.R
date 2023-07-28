load_state_data1 <- function (gt.folder, ili.folder, population.file, gft.file) 
{
  gtfiles <- list.files(gt.folder)
  gtfiles <- gtfiles[which(sapply(strsplit(gtfiles, "[.]"), 
                                  function(x) x[2]) == "csv")]
  GTdata <- list()
  GTdata_week <- list()
  GTdata_month <- list()
  state.info <- c()
  for (f in gtfiles) {
    fscaned <- scan(file.path(gt.folder, f), what = character(), 
                    sep = "\n", quiet = TRUE)
    if (length(fscaned) < 4 || all(fscaned != "Interest over time")) 
      next
    state.info[f] <- substr(fscaned[2], 1, gregexpr(" 2004", 
                                                    fscaned[2])[[1]] - 1)
    endline <- which(!grepl(",", fscaned))[4] - 1
    if (is.na(endline)) 
      endline <- length(fscaned)
    con <- textConnection(paste(fscaned[4:endline], collapse = "\n"))
    GTdata[[f]] <- read.csv(con)
    if (substr(fscaned[4], 1, 4) == "Week") {
      GTdata_week[[f]] <- xts(GTdata[[f]][, 2, drop = FALSE], 
                              as.Date(GTdata[[f]]$Week) + 6)
    }
    else {
      GTdata_month[[f]] <- xts(GTdata[[f]][, 2, drop = FALSE], 
                               as.yearmon(GTdata[[f]]$Month))
    }
    if (!grepl(colnames(GTdata[[f]])[2], gsub(" ", ".", 
                                              f))) {
      cat(f, colnames(GTdata[[f]])[2], "\n")
    }
    #print(f)
    close(con)
  }
  GTdata_week.state <- tapply(GTdata_week, state.info, function(gt.eachstate) {
    tab <- do.call(merge, gt.eachstate)
    tab[, !grepl("1", colnames(tab))]
  })
  names(GTdata_week.state) <- sapply(names(GTdata_week.state), 
                                     function(x) {
                                       if (x == "United States") 
                                         return("US")
                                       paste0("US-", state.abb[gsub("\\s+$", "", strsplit(x, 
                                                                                          "\\(")[[1]][1]) == state.name])
                                     })
  GTdata_week.state <- GTdata_week.state[order(names(GTdata_week.state))]
  GT_national <- GTdata_week.state$US
  GT_national <- na.omit(GT_national)
  GT_data <- list()
  for (j in state.abb) {
    empty_GT <- GT_national
    empty_GT[] <- 0
    GT_data[[j]] <- GTdata_week.state[[paste0("US-", j)]]
    empty_GT[, colnames(GT_data[[j]])] <- GT_data[[j]][index(GT_national)]
    GT_data[[j]] <- empty_GT
    #print(j)
  }
  
  ###ili
  ili <- read.csv(file.path(ili.folder, "ILINet_regional.csv"), 
                  skip = 1)
  ili_idx <- subset(ili, REGION == "Region 1")[, c("YEAR", 
                                                   "WEEK")]
  ili_idx$DATE <- as.Date("1997-10-04") + (1:nrow(ili_idx) - 
                                             1) * 7
  ili_idx$SEASON <- ifelse(ili_idx$WEEK < 40, paste0(ili_idx$YEAR - 
                                                       1, "-", ili_idx$YEAR), paste0(ili_idx$YEAR, "-", ili_idx$YEAR + 
                                                                                       1))
  ili_regional <- list()
  for (j in 1:10) {
    ili_regional[[paste0("Region.", j)]] <- as.numeric(as.character(ili$X..WEIGHTED.ILI))[ili$REGION == 
                                                                                            paste("Region", j)]
  }
  ili_regional <- do.call(cbind, ili_regional)
  ili_regional <- xts(ili_regional, ili_idx$DATE)
  
  ili_national <- read.csv(file.path(ili.folder, "ILINet_nat.csv"), 
                           skip = 1)
  ili_national <- xts(as.numeric(as.character(ili_national$X..WEIGHTED.ILI)), 
                      as.Date("1997-10-04") + (1:nrow(ili_national) - 1) * 
                        7)
  states <- read.csv(population.file, stringsAsFactor = F)
  states$Population <- as.numeric(gsub(",", "", states$Population))
  states <- subset(states, Abbre != "DC")
  
  
  ##state
  ili1 <- read.csv(file.path(ili.folder, "ILINet_state.csv"), 
                  skip = 1)
  
  ili1_idx <- subset(ili1, REGION == "Alabama")[, c("YEAR", 
                                                   "WEEK")]
  # week ending 20101009
  ili1_idx$DATE <- as.Date("2010-10-09") + (1:nrow(ili1_idx) - 
                                             1) * 7
  ili1_idx$SEASON <- ifelse(ili1_idx$WEEK < 40, paste0(ili1_idx$YEAR - 
                                                       1, "-", ili1_idx$YEAR), paste0(ili1_idx$YEAR, "-", ili1_idx$YEAR + 
                                                                                       1))
  
  ili_state <- list()
  for (j in 1:length(states$State)) {
    ili_state[[states$Abbre[j]]] <- as.numeric(as.character(ili1$X.UNWEIGHTED.ILI))[ili1$REGION ==states$State[j]]
  }
  ili_state <- do.call(cbind, ili_state)
  ili_state <- xts(ili_state, ili1_idx$DATE)
  
  ili_state <- ili_state[,-match("FL", colnames(ili_state))]
  
  ##GT regional
  GT_regional <- lapply(1:10, function(region_number) {
    states_id <- which(states$Region == region_number)
    GT_states_wgted <- lapply(states_id, function(j) {
      GT_data[[states$Abbre[j]]] * states$Population[j]/sum(states$Population[states_id])
    })
    Reduce("+", GT_states_wgted)
  })
  names(GT_regional) <- paste0("Region.", 1:10)
  
  ##GT state
  GT_state <- list()
  for(j in 1:length(names(GT_data))){
    if(states$Abbre[[j]]!="FL"){
      GT_state[[states$Abbre[[j]]]] <- GT_data[[states$Abbre[[j]]]]
    }
  }
  
  ##GFT state
  gft.pred <- read.csv(gft.file, skip = 10)
  gft.pred <- xts(gft.pred[, -1], as.Date(gft.pred[, 1]) + 
                    6)
  colnames(gft.pred) <- lapply(strsplit(colnames(gft.pred), "[.][.]"), 
                               function(x) paste(x, collapse = " "))
  colnames(gft.pred) <- lapply(strsplit(colnames(gft.pred), "[.]"), 
                               function(x) paste(x, collapse = " "))
  
  gft.pred <- gft.pred[, match(c(states$State, "New York NY"), colnames(gft.pred))]
  gft.pred <- gft.pred/1000
  #names(gft.pred) <- sapply(strsplit(names(gft.pred), "[.][.]"), 
  #                          function(x) x[1])
  
  list(ili_national = ili_national, ili_regional = ili_regional, 
       GT_national = GT_national, GT_regional = GT_regional, GT_regional = GT_regional, 
       GT_state =GT_data, ili_state =ili_state,
       GFT = gft.pred)
}

load_gft_data1  <- function (population.file, gft.file) 
{
  states <- read.csv(population.file, stringsAsFactor = F)
  states$Population <- as.numeric(gsub(",", "", states$Population))
  state.abb <- states$Abbre
  state.name <- states$State
  
  gft.pred <- read.csv(gft.file, skip = 10)
  gft.pred <- xts(gft.pred[, -1], as.Date(gft.pred[, 1]) + 
                    6)
  gft.pred.nat <-  gft.pred[, 1]
  gft.pred.nat <- gft.pred.nat/1000
  
  gft.pred <- gft.pred[, grep("HHS.Region", colnames(gft.pred))]
  gft.pred <- gft.pred/1000
  names(gft.pred) <- sapply(strsplit(names(gft.pred), "[.][.]"), 
                            function(x) x[1])
  

  
  gft.pred.state <- read.csv(gft.file, skip = 10)
  gft.pred.state <- xts(gft.pred.state[, -1], as.Date(gft.pred.state[, 1]) + 
                    6)
  colnames(gft.pred.state) <- lapply(strsplit(colnames(gft.pred.state), "[.][.]"), 
                               function(x) paste(x, collapse = " "))
  colnames(gft.pred.state) <- lapply(strsplit(colnames(gft.pred.state), "[.]"), 
                               function(x) paste(x, collapse = " "))
  
  gft.pred.state <- gft.pred.state[, match(c(states$State, "New York NY"), colnames(gft.pred.state))]
  gft.pred.state <- gft.pred.state/1000
  
  list(GFT = gft.pred, GFT_state = gft.pred.state, GFT_nat = gft.pred.nat)
}


argo2State <- function (truth, argo1.p, argo.nat.p, states.abb.all1=states.abb.all) 
{
  naive.p <- lag(truth, 1)
  common_idx <- zoo::index(na.omit(merge(truth, naive.p, argo1.p, 
                                         argo.nat.p)))
  if (ncol(argo.nat.p) == 1) {
    argo.nat.p <- do.call(cbind, lapply(1:ncol(argo1.p), function(i) argo.nat.p))
  }
  colnames(argo.nat.p) <- colnames(truth)
  Z <- truth - naive.p
  W <- merge(lag(Z, 1), argo1.p - naive.p, argo.nat.p - naive.p)
  arg2zw_result <- argo2zwState(as.matrix(na.omit(merge(Z, W))))
  Z.hat <- arg2zw_result$Z.hat
  Z.hat <- xts(Z.hat, as.Date(rownames(Z.hat)))
  argo2.p <- Z.hat + naive.p
  heat.vec <- na.omit(merge(Z, lag(Z), argo1.p - truth, argo.nat.p - 
                              truth))
  colnames(heat.vec) <- paste0(rep(c("detla.ili.", "err.argo.", 
                                     "err.nat.", "lag.detla.ili."), each = ncol(Z)), rep(states.abb.all1, 
                                                                                    3))
  result_wrapper <- list(onestep = argo1.p, twostep = argo2.p, 
                         naive = naive.p, truth = truth, heat.vec = heat.vec)
  c(result_wrapper, arg2zw_result)
}


argo2zwState <- function (zw.mat) 
{
  if (is.null(rownames(zw.mat))) 
    stop("row name must exist as index")
  projection.mat <- list()
  mean.mat <- list()
  zw_used <- list()
  sigma_ww.structured <- sigma_ww.empirical <- sigma_zw.structured <- sigma_zw.empirical <- heat.vec.structured <- sigma_zwzw.structured <- sigma_zwzw.empirical <- list()
  epsilon <- list()
  Z <- zw.mat[, 1:(ncol(zw.mat)/4)]
  W <- zw.mat[, -(1:(ncol(zw.mat)/4))]
  W1 <- W[, 1:(ncol(zw.mat)/4)]
  W2 <- W[, 1:(ncol(zw.mat)/4)+ (ncol(zw.mat)/4)]
  W3 <- W[, 1:(ncol(zw.mat)/4)+ (ncol(zw.mat)/4)*2]
  Z.hat <- Z
  Z.hat[] <- NA
  for (it in 105:nrow(zw.mat)) {
    training_idx <- (it - 104):(it - 1)
    t.now <- rownames(zw.mat)[it]
    if (is.null(t.now)) 
      t.now <- it
    epsilon[[as.character(t.now)]] <- list()
    sigma_zz <- var(Z[training_idx, ])
    sigma_zz_chol <- chol(sigma_zz)
    epsilon[[as.character(t.now)]]$z <- solve(t(sigma_zz_chol), 
                                              t(data.matrix(Z[training_idx, ])))
    zw_used[[as.character(t.now)]] <- cbind(Z[training_idx, 
                                              ], W[training_idx, ])
    m1 <- cor(Z[training_idx, ], W[training_idx, 1:(ncol(zw.mat)/4)])
    m2 <- cor(Z[training_idx, ])
    rho <- sum(m1 * m2)/sum(m2^2)
    d.gt <- diag(diag(var(W2[training_idx, ] - Z[training_idx, 
                                                 ])))
    epsilon[[as.character(t.now)]]$argoreg <- t(data.matrix((W2 - 
                                                               Z)[training_idx, ]))
    epsilon[[as.character(t.now)]]$argoreg <- epsilon[[as.character(t.now)]]$argoreg/sqrt(diag(d.gt))
    sigma.nat <- var((W3 - Z)[training_idx, ])
    epsilon[[as.character(t.now)]]$argonat <- t(data.matrix((W3 - 
                                                               Z)[training_idx, ]))
    sigma.nat_chol <- chol(sigma.nat)
    epsilon[[as.character(t.now)]]$argonat <- solve(t(sigma.nat_chol), 
                                                    epsilon[[as.character(t.now)]]$argonat)
    sigma_ww <- rbind(cbind(sigma_zz, rho * sigma_zz, rho * 
                              sigma_zz), cbind(rho * sigma_zz, sigma_zz + d.gt, 
                                               sigma_zz), cbind(rho * sigma_zz, sigma_zz, sigma_zz + 
                                                                  sigma.nat))
    sigma_zw <- cbind(rho * sigma_zz, sigma_zz, sigma_zz)
    mu_w <- colMeans(W[training_idx, ])
    mu_z <- colMeans(Z[training_idx, ])
    d_ww <- diag(diag(var(W[training_idx, ])))
    pred.blp <- mu_z + sigma_zw %*% solve(sigma_ww, W[t.now, 
                                                      ] - mu_w)
    Kzz <- solve((1 - rho^2) * sigma_zz)
    Kgt <- diag(1/diag(var((W2 - Z)[training_idx, ])))
    Knat <- solve(var((W3 - Z)[training_idx, ]))
    pred.bayes <- mu_z + solve(Kzz + Knat + Kgt, Knat %*% 
                                 (W[t.now, 1:(ncol(zw.mat)/4)+(ncol(zw.mat)/4)*2] - 
                                    mu_w[1:(ncol(zw.mat)/4)+(ncol(zw.mat)/4)*2]) + Kgt %*% 
                                 (W[t.now, 1:(ncol(zw.mat)/4)+(ncol(zw.mat)/4)] - 
                                    mu_w[1:(ncol(zw.mat)/4)+(ncol(zw.mat)/4)]) + 
                                 Kzz %*% (rho * (W[t.now,1:(ncol(zw.mat)/4)] - 
                                                   mu_w[1:(ncol(zw.mat)/4)])))
    if (all(is.finite(pred.blp))) {
      stopifnot(all(abs(pred.blp - pred.bayes) < 1e-08))
    }
    z.hat <- mu_z + sigma_zw %*% solve(sigma_ww + d_ww, 
                                       (W[t.now, ] - mu_w))
    Z.hat[t.now, ] <- t(z.hat)
    projection.mat[[as.character(t.now)]] <- sigma_zw %*% 
      solve(sigma_ww + d_ww)
    mean.mat[[as.character(t.now)]] <- c(mu_z, mu_w)
    sigma_ww.structured[[as.character(t.now)]] <- sigma_ww
    sigma_ww.empirical[[as.character(t.now)]] <- var(W[training_idx, 
                                                       ])
    sigma_zw.structured[[as.character(t.now)]] <- sigma_zw
    sigma_zw.empirical[[as.character(t.now)]] <- cov(Z[training_idx, 
                                                       ], W[training_idx, ])
    sigma_zwzw.structured[[as.character(t.now)]] <- rbind(cbind(sigma_zz, 
                                                                sigma_zw), cbind(t(sigma_zw), sigma_ww))
    sigma_zwzw.empirical[[as.character(t.now)]] <- var(cbind(Z, 
                                                             W)[training_idx, ])
    heat.vec.struc <- rbind(cbind(sigma_zz, rho * sigma_zz), 
                            cbind(rho * sigma_zz, sigma_zz))
    heat.vec.struc <- Matrix::bdiag(heat.vec.struc, d.gt, 
                                    sigma.nat)
    heat.vec.structured[[as.character(t.now)]] <- as.matrix(heat.vec.struc)
  }
  projection.mat <- sapply(projection.mat, identity, simplify = "array")
  mean.mat <- sapply(mean.mat, identity, simplify = "array")
  sigma_ww.structured <- sapply(sigma_ww.structured, identity, 
                                simplify = "array")
  sigma_ww.empirical <- sapply(sigma_ww.empirical, identity, 
                               simplify = "array")
  sigma_zw.structured <- sapply(sigma_zw.structured, identity, 
                                simplify = "array")
  sigma_zw.empirical <- sapply(sigma_zw.empirical, identity, 
                               simplify = "array")
  sigma_zwzw.structured <- sapply(sigma_zwzw.structured, identity, 
                                  simplify = "array")
  sigma_zwzw.empirical <- sapply(sigma_zwzw.empirical, identity, 
                                 simplify = "array")
  zw_used <- sapply(zw_used, identity, simplify = "array")
  heat.vec.structured <- sapply(heat.vec.structured, identity, 
                                simplify = "array")
  return(list(Z.hat = Z.hat, heat.vec.structured = heat.vec.structured, 
              projection.mat = projection.mat, mean.mat = mean.mat, 
              sigma_ww.structured = sigma_ww.structured, sigma_ww.empirical = sigma_ww.empirical, 
              sigma_zw.structured = sigma_zw.structured, sigma_zw.empirical = sigma_zw.empirical, 
              sigma_zwzw.structured = sigma_zwzw.structured, sigma_zwzw.empirical = sigma_zwzw.empirical, 
              zw_used = zw_used, epsilon = epsilon))
}
