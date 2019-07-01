#################################################
## Permutation test simulations
## Author : Adam Elder
#################################################

source("cv_test_neo.R")
source("one_fold_est.R")
source("other_functions.R")
source("est_cv.R")
source("ZL_test.R")
library(MASS) 

sim_mat <- expand.grid("samp_size" = c(100), 
                       "marg_indp" = c(TRUE),
                       "xx_cor" = c(0.5),
                       "prop_x_cor" = c(0, 0.0001, 0.6),
                       "xy_cor" = c(0.1, 0.5), "dim" = c(10, 50),
                       "dist_acc" = c(250, 500, 750),
                       "nrm_typ" = c("l2", "lp", "max", "ZL", "ssq"),
                       "proc_type" = c(0))

make_cov_mat <- function(dim, prop_cor, x_cor, xy_cor){
  mat <- matrix(x_cor, nrow = dim + 1, ncol = dim + 1)
  diag(mat) <- 1
  num_cor <- ceiling(prop_cor * dim)
  mat[1, ] <- mat[, 1] <- c(1, rep(xy_cor, num_cor), rep(0, dim - num_cor))
  return(mat)
}

c_r <- rep(NA, nrow(sim_mat))

for(i in 1:nrow(sim_mat)){
  c_r[i] <- det(make_cov_mat(dim = sim_mat[i, "dim"], prop_cor = sim_mat[i, "prop_x_cor"],
                             x_cor = sim_mat[i, "xx_cor"], xy_cor = sim_mat[i, "xy_cor"]))
}

sim_mat <- sim_mat[c_r > 0, ]

syst.envr <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if(!is.na(syst.envr)){
  j_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
}else{
  cat("No Job ID found, choosing the answer to
        life the universe and everything \n")
  j_id <- 32
}

this_sim <- sim_mat[j_id, ]
#this_sim <- sim_mat[1, ];rownames(sim_mat) <- 1:nrow(sim_mat); subset(sim_mat, xx_cor == 0 & prop_x_cor == 0.1 & xy_cor == 0.3 & dim == 10)
#this_sim <- sim_mat[146, ]

cor_mat <- make_cov_mat(dim = this_sim[1, "dim"], prop_cor = this_sim[1, "prop_x_cor"], 
                        x_cor = this_sim[1, "xx_cor"], xy_cor = this_sim[1, "xy_cor"])

if(this_sim[1, "nrm_typ"] == "lp"){
  norm_type <- "lp"
  norms_chosen <- c(1, 2, 4, 6, "max")
}else if (this_sim[1, "nrm_typ"] == "ssq"){
  norm_type <- "ordl2"
  norms_chosen <- unique(round(seq(from = 1, to = this_sim[1, "dim"], length.out = 6)))
}else if (this_sim[1, "nrm_typ"] == "l2"){
  norm_type <- "lp"
  norms_chosen <- c(2)
}else if (this_sim[1, "nrm_typ"] == "max"){
  norm_type <- "lp"
  norms_chosen <- "max"
}else if (this_sim[1, "nrm_typ"] == "ZL"){
  norm_type <- "ZL"
  norms_chosen <- "tstat"
}

if(this_sim[1, "proc_type"] == 0){proc_type <- "mag"}else{proc_type <- "est_pow"}
mi <- this_sim[1, "marg_indp"]
dist_acc <- this_sim[1, "dist_acc"]

num_correlated <- ceiling(this_sim[1, "prop_x_cor"] * this_sim[1, "dim"])
cat("Correlation Matrix : \n" )
samp_data <- gen_pois_samp(10000, 4, xycor = this_sim[1, "xy_cor"],
                           xxcor = this_sim[1, "xx_cor"],
                           num_cor = num_correlated,
                           num_uncor = this_sim[1, "dim"] - num_correlated,
                           marg_indep = this_sim[1, "marg_indp"])

print(cor(samp_data))

cat("norm type is ", norm_type, "\n")
cat("Each test uses ", 1, "fold \n")
cat("The possible norms are ", norms_chosen, "\n")
cat("sample size is ", this_sim[1, 1], "\n")
cat("data is generated", if(!mi){" not "}else{" "}, "with marginal independence \n")
cat("The measure of performance is ", proc_type, " \n")

this_res <- matrix(NA, nrow = 400, ncol = length(norms_chosen) + 1)
bonf_res <- rep(NA, 400)

for(bb in 1:400){
  if(bb %% 4 == 0){cat(bb, " ")}
  num_correlated <- ceiling(this_sim[1, "prop_x_cor"] * this_sim[1, "dim"])
  vld_data <- FALSE
  while(!vld_data){
    gen_data <- gen_pois_samp(this_sim[1, "samp_size"], 4, xycor = this_sim[1, "xy_cor"],
                              xxcor = this_sim[1, "xx_cor"],
                              num_cor = num_correlated,
                              num_uncor = this_sim[1, "dim"] - num_correlated,
                              marg_indep = mi)
    vld_data <- (sum(gen_data[, 1]) != 0)
  }
  if(norm_type != "ZL"){
    this_res[bb, ] <- cv_test_neo(obs_data = gen_data, pos_lp_norms = norms_chosen,
                                  ts_ld_bs_samp = dist_acc,
                                  num_folds = 1, n_bs_smp = dist_acc, 
                                  nrm_type = norm_type, f_estimate = est_pearson,
                                  f_cv_summary = mean, test_stat_func = proc_type,
                                  incl_chsn_norm = TRUE, test_type = "perm",
                                  perf_meas = proc_type) 
  }else{
    this_res[bb, ] <- ZL(observed_data = gen_data, dist_acc * 4, dist_acc * 4)
  }

  bonf_res[bb] <- as.integer(min(bonf_test(obs_data = gen_data)) * this_sim[1, "dim"] <= 0.05)
}

this_res <- cbind(this_res, bonf_res)

cat("\n")
if (!dir.exists("res")){dir.create("res")}
write.csv(this_res,
          file = paste0("res/sims_n_samp_", this_sim[1, "samp_size"],
                        "_mar_indep_", this_sim[1, "marg_indp"], 
                        "_norm_type_", this_sim[1, "nrm_typ"],
                        "_prop_cor_", this_sim[1, "prop_x_cor"], 
                        "_dist_acc_", this_sim[1, "dist_acc"],
                        "_meas_", proc_type, 
                        "_xy_cor_", this_sim[1, "xy_cor"], 
                        "_dim_", this_sim[1, "dim"],"_tsq_", ".csv"))

sub_idx <- function(data, sub_chsn_var){
  subs_needed <- length(sub_chsn_var)
  subset_idx <- rep(1, nrow(data))
  for(subs in 1:subs_needed){
    sub_info <- sub_chsn_var[[subs]]
    sub_mlt <- as.numeric(data[, sub_info$var_name] %in% sub_info$var_vals)
    subset_idx <- subset_idx * sub_mlt
  }
  return(subset_idx)
}
