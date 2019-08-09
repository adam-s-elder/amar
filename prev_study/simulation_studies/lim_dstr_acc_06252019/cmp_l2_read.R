#############################
#### read sims
#### december 6th 2018
#### author : adam elder
#############################

library(tidyverse)
library(kableExtra)

sim_mat <- expand.grid("samp_size" = c(100),
                       "marg_indp" = c(TRUE),
                       "prop_x_cor" = c(0, 0.0001, 0.6),
                       "xy_cor" = c(0.1, 0.5), "dim" = c(10, 50),
                       "dist_acc" = c(250, 500, 750),
                       "nrm_typ" = c("l2", "lp", "max", "ZL", "ssq"),
                       "proc_type" = c(0))


res_name <- function(n, marg_indp, norm_type, prop_cor, dist_acc,
                     proc_type, xy_cor, dim, res_fold = "res"){
  if(prop_cor == 0.0001){prop_cor <- "1e-04"}
  paste0("res/sims_n_samp_", n, "_mar_indep_", marg_indp,
         "_norm_type_", norm_type, "_prop_cor_", prop_cor,
         "_dist_acc_", dist_acc, "_meas_", proc_type,
         "_xy_cor_", xy_cor, "_dim_", dim,"_tsq_", ".csv")
}

one_sim <- function(n, marg_indp, norm_type, prop_cor, dist_acc, proc_type, xy_cor, dim, res_fold = "res"){
  res_name <- res_name(n, marg_indp, norm_type, prop_cor, dist_acc,
                       proc_type, xy_cor, dim, res_fold = res_fold)
  if(file.exists(res_name)){
    result <- read.csv(res_name)
    num_cols <- ncol(result)
    return(round(c(as.numeric(mean(result[, 2] <= 0.05)),
                   mean(result[, num_cols])), 3))
  }else{
    return(NA)
  }
}

one_sim_df <- function(df_row, res_fold = "res"){
  proc_type <- ifelse(df_row[1, "proc_type"] == 0, "mag", "est_pow")
  one_sim(df_row[1, "samp_size"], df_row[1, "marg_indp"], df_row[1, "nrm_typ"],
          df_row[1, "prop_x_cor"], df_row[1, "dist_acc"], proc_type,
          df_row[1, "xy_cor"], df_row[1, "dim"], res_fold = res_fold)
}

get_res <- function(df, res_fold = "res"){
  nrows <- nrow(df)
  inf_df <- one_sim_df(df[1, ])
  if(length(inf_df) > 1){
    for(bb in 2:nrows){
      inf_df <- rbind(inf_df, one_sim_df(df[bb, ], res_fold = res_fold))
    }
    }else{
      for(bb in 2:nrows){
        inf_df <- c(inf_df, one_sim_df(df[bb, ], res_fold = res_fold)[1])
      }
    }
  return(cbind(df, inf_df))
}

good_sort <- function(df, list_of_cols){
  loc <- rev(list_of_cols)
  num_cols <- length(loc)
  fin_df <- df
  col_idx <- c()
  for(bb in 1:num_cols){
    col_idx <- c(which(colnames(df) == loc[bb]), col_idx)
    df_ord <- order(fin_df[, loc[bb]])
    fin_df <- fin_df[df_ord, ]
  }
  fin_df <- cbind(fin_df[, col_idx], fin_df[, -col_idx])
  return(fin_df)
}

rename_col_vals <- function(df, col_names, pre_val_str){
  num_repl <- length(col_names)
  fin_df <- df
  for(col_idx in 1:num_repl){
    new_vals <- paste0(pre_val_str[col_idx], fin_df[, col_names[col_idx]])
    fin_df[, col_names[col_idx]] <- new_vals
  }
  return(fin_df)
}

## functions for making exciting tables

make_kable_table <- function(data_table, form_table){
  if(all(dim(data_table) != dim(form_table))){
    stop("The data_table and form_table must have matching dimensions")
  }
  fin_df <- data.frame(matrix(NA, nrow = nrow(data_table), ncol = ncol(data_table)))
  n_cols <- ncol(data_table)
  for(col_idx in 1:n_cols){
    temp_col <- kableExtra::cell_spec(data_table[, col_idx], "latex",bold = TRUE,
                                      color = factor(form_table[, col_idx],
                                                     c(0, 1, 2),
                                                     c("black", "white", "red")),
                                      background = ifelse(form_table[, col_idx] == 0,
                                                          "white", "black"))
    fin_df[, col_idx] <- temp_col
  }
  return(fin_df)
}

get_form_table <- function(dat_table, tol_lev){
  n_col <- ncol(dat_table)
  data_table <- matrix(as.numeric(dat_table), ncol = n_col)
  form_table <- matrix(NA, nrow = nrow(data_table), ncol = n_col)
  for(row_idx in 1:nrow(data_table)){
    row_vals <- rep(0, n_col)
    sub_x <- data_table[row_idx, ]
    x_ord <- order(sub_x, decreasing = TRUE)
    if(!is.na(sub_x[x_ord[2]])){
      if(sub_x[x_ord[1]] - sub_x[x_ord[2]] > tol_lev){
        row_vals[x_ord[1]] <- 2
      }else{
        row_vals[x_ord[1]] <- 1
      }
    }
    form_table[row_idx, ] <- row_vals
  }
  return(form_table)
}

color_tbl <- function(table, cutoff){
  make_kable_table(table, get_form_table(table, cutoff))
}

make_table <- function(res_tbl, grp_vars, grp_var_add, splt1, splt2, val, sg_df = 0.1, nd = 2){
  n_res_tbl <- res_tbl[, c(grp_vars, splt1, splt2, val)]
  # colnames(n_res_tbl)[ncol(n_res_tbl)] <- "val"
  # colnames(n_res_tbl)[ncol(n_res_tbl) - 2] <- "splt1"
  splt_1_vals <- unique(n_res_tbl[, splt1])
  splt_2_vals <- unique(n_res_tbl[, splt2])
  n_grp_vars <- length(grp_vars)
  n_s1_vars <- length(splt_1_vals)
  n_s2_vars <- length(splt_2_vals)
  split_dfs <- list()
  mis_idx <- matrix(NA, nrow = nrow(n_res_tbl)/(n_s1_vars * n_s2_vars), ncol = n_s2_vars)
  for(s2_idx in 1:n_s2_vars){
    sub_df <- n_res_tbl[n_res_tbl[, splt2] == splt_2_vals[s2_idx], c(grp_vars, splt1, val)]
    s_sub_df <- spread(sub_df, key = splt1, value = val, drop = FALSE)
    exl_cols <- which(colnames(s_sub_df) %in% grp_vars)
    splt_1_vals <- colnames(s_sub_df[, -exl_cols])
    s_sub_df <- good_sort(s_sub_df, grp_vars)
    s_sub_df <- rename_col_vals(s_sub_df, grp_vars, grp_var_add)
    s_mis_idx <- apply(s_sub_df[, -exl_cols], 1, function(x) all(is.na(x)))
    if(nrow(mis_idx) != length(s_mis_idx)){
      mis_idx <- matrix(NA, nrow = length(s_mis_idx), ncol = n_s2_vars)
    }
    mis_idx[, s2_idx] <- s_mis_idx
    mat_s <- round(as.matrix(s_sub_df[, -exl_cols]), nd)
    s_sub_df[, -exl_cols] <- color_tbl(mat_s, sg_df)
    split_dfs[[as.character(splt_2_vals[s2_idx])]] <- s_sub_df
  }
  fin_fin_df <- split_dfs[[1]]
  for(bb in 2:length(split_dfs)){
    fin_fin_df <- cbind(fin_fin_df, split_dfs[[bb]][, -exl_cols])
  }
  cut_idx <- which(apply(mis_idx, 1, function(x) all(x)))
  fin_fin_df <- fin_fin_df[-cut_idx, ]
  spl_1 <- rep(1, n_s1_vars)
  names(spl_1) <- splt_1_vals
  spl_2 <- rep(n_s1_vars, n_s2_vars)
  names(spl_2) <- paste0(splt2, " = ", splt_2_vals)
  colnames(fin_fin_df) <- NULL
  knitr::kable(fin_fin_df, format = "latex", align = "c", #longtable = TRUE,
               booktabs = TRUE, escape = FALSE, linesep = "", row.names = FALSE) %>%
    kable_styling(latex_options = "scale_down") %>%
    # row_group_label_position = "stack"
    collapse_rows(columns = c(1:n_grp_vars)) %>%
    add_header_above(c(" "  = n_grp_vars, rep(spl_1, n_s2_vars))) %>%
    add_header_above(c(" " = n_grp_vars, spl_2))
}

initial_table <- get_res(sim_mat, res_fold = "res")


bonf_table <- expand.grid("samp_size" = c(100),  "marg_indp" = c(TRUE),
                          "prop_x_cor" = c(0, 0.0001, 0.6),
                          "xy_cor" = c(0.1, 0.5), "dim" = c(10, 50),
                          "dist_acc" = c(250, 500, 750),
                          "nrm_typ" = c("bonf"), "proc_type" = c(0), "1" = 1)

  expand.grid("samp_size" = c(100),  "marg_indp" = c(TRUE, FALSE),
                          "prop_x_cor" = c(0, 0.0001, 0.1, 0.2),
                          "xy_cor" = c(0.1, 0.3, 0.5), "dim" = c(10, 50, 100),
                          "dist_acc" = c(250, 500, 750),
                          "nrm_typ" = c("bonf"), "proc_type" = c(0), '1' = 1)

for(bon_idx in 1:nrow(bonf_table)){
  sub_idx <- which(apply(initial_table[, 1:6], 1,
                         function(x) all(x == bonf_table[bon_idx, 1:6])))
  sub_res <- initial_table[sub_idx,]
  bonf_table[bon_idx, "1"] <- mean(sub_res[, "2"], na.rm = TRUE)
}

final_table <- rbind(initial_table[, 1:9], bonf_table)
colnames(final_table)[9] <- "val"

mar_ind_tb <- subset(final_table, marg_indp == TRUE)

mar_dep_tb <- subset(final_table, marg_indp == FALSE)

make_table(res_tbl = mar_ind_tb, grp_vars = c("prop_x_cor", "xy_cor", "dim"),
           splt1 = "dist_acc", splt2 = "nrm_typ",
           grp_var_add = c("Prop of X cor = ", "X Y Correlation = ", ""),
           val = "val")

make_table(res_tbl = mar_dep_tb, grp_vars = c("prop_x_cor", "xy_cor", "dim"),
           splt1 = "nrm_typ", splt2 = "xx_cor",
           grp_var_add = c("Prop of X cor = ", "X Y Correlation = ", ""),
           val = "val")







