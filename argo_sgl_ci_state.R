## 20230627 summary argo-sgl state CI
## S. Ning
## input
## path_out_state
# ready to run

load(file.path(path_out_state, "argo_sgl_states_all_logit01_final.Rdata"))

truth <-  ili_state
ind.CI <- (truth > argo_sgl_state_pred_final - 1.96 * sqrt(argo_sgl_state_var_est_final)) &
  (truth < argo_sgl_state_pred_final + 1.96 * sqrt(argo_sgl_state_var_est_final))

coverage_CI <- apply(na.omit(ind.CI), 2, mean)
names(coverage_CI) <- states.info$abbre

library(xtable)

print(xtable(t(as.matrix(coverage_CI[1:15])), digits = 3), include.rownames = F)
print(xtable(t(as.matrix(coverage_CI[1:15+15])), digits = 3), include.rownames = F)
print(xtable(t(as.matrix(coverage_CI[1:15+15*2])), digits = 3), include.rownames = F)
print(xtable(t(as.matrix(coverage_CI[46:51])), digits = 3), include.rownames = F)

mean(coverage_CI)
