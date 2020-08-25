######################################################################
#### Estimating Influence Curve and parameter for various parameters.
#### Adam Elder
#### This script contains the various proceedures
#### used for estimating both the parameter, and the influence
#### curves used for estimating the projected risk ratio.
#### These functions will be given an entire dataset, and will
#### return the estimate for the parameter, and for each
#### estimated influence curve at each observation.
######################################################################

# The following function will take a set of observations, and return the
# estimated covariance matrix for the estimates of the pearson correlation.
# Estimates of the covariance are generated using the empirical influence
# function.  The first column of your data should correspond to the
# variable of interest (the variable for which pearson correlation is
# calculated).
#

#' Obtain an estimator of the probability delta = 1 given w
#' @param delt_vec Vector of missingess indicators
#' @param obs_ws The observed w's
#' @export

las_dw_nl <- function(delt_vec, obs_ws){
  sl_fit <- SuperLearner::SuperLearner(
    Y = delt_vec, X = obs_ws,
    ## Maybe try SL.polymars instead??
    SL.library = c("SL.glm", "SL.loess", "SL.ranger"),
    family = "gaussian")

  new_funct <- function(ws){
    pmax(stats::predict(sl_fit, newdata = ws,
                        onlySL = TRUE)$pred[, , drop = TRUE],
         0.01)
  }
  return(new_funct)
}

las_ydw_nl <- function(y_vec, delt_vec, obs_ws){
  ys <- y_vec[delt_vec == 1]
  ws <- obs_ws[delt_vec == 1, ]
  sl_fit <- SuperLearner::SuperLearner(
    Y = ys, X = ws,
    ## Maybe try SL.polymars instead??
    SL.library = c("SL.glm", "SL.loess", "SL.ranger"),
    family = "gaussian")
  l_mod_ydw <- function(ws){
    pmax(stats::predict(sl_fit, newdata = ws,
                   onlySL = TRUE)$pred[, , drop = TRUE],
         0.01)
  }
  return(l_mod_ydw)
}

sl_me_nl <- function(m_e_v, w_j){
  sl_mod <- SuperLearner::SuperLearner(Y = m_e_v,
                                       X = data.frame("WJ" = w_j),
                                       ## Maybe try SL.polymars instead??
                                       SL.library = c("SL.glm", "SL.loess",
                                                      "SL.ranger"),
                                       family = "gaussian")

  sl_exp_wj <- function(w_j){
    new_data <- data.frame("WJ" = w_j)
    pmax(0.01, stats::predict(sl_mod, newdata = new_data,
            onlySL = TRUE)$pred[, , drop = TRUE])
  }
  return(sl_exp_wj)
}

get_Dmat_nl <- function(obs_data){
  num_cov <- ncol(obs_data) - 2
  num_obs <- nrow(obs_data)

  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))

  las_del <- las_dw_nl(obs_data$delta, obs_data[, w_cols])
  las_y  <- las_ydw_nl(obs_data$y, obs_data$delta,
                    obs_data[, w_cols])
  Dmat <- array(NA,  dim = c(num_obs, 4, num_cov))
  psi_mat <- matrix(nrow = num_cov, ncol = 4)

  for(cov_idx in 1:num_cov){

    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj_nl <- suppressWarnings(sl_me_nl(las_y(obs_data[, w_cols]), cov_col))
    mar_exp_vec <- log(marg_exp_wj_nl(cov_col)) # This subset may be needed [obs_data$delta == 1]
    mme_wj <- mean(mar_exp_vec * cov_col)
      # This subset may be needed [obs_data$delta == 1]
    mme_one <- mean(mar_exp_vec)

    psi_1 <- mme_wj
    psi_2 <- mme_one
    psi_3 <- mean(cov_col)
    psi_4 <- mean(cov_col ** 2)
    psi_mat[cov_idx, ] <- c(psi_1, psi_2, psi_3, psi_4)

    dmat12 <- psi_12(y = obs_data$y, d = obs_data$delta,
                     w = obs_data[, w_cols],
                     wj = obs_data[, w_cols[cov_idx]],
                     pr_dw = las_del, epr_ywj = marg_exp_wj_nl,
                     pr_yw = las_y, exp_wj = mme_wj,
                     exp_one = mme_one)
    dmat3 <- cov_col - mean(cov_col)
    cov_col_sq <- cov_col ** 2
    dmat4 <- cov_col_sq - mean(cov_col_sq)
    Dmat[, ,cov_idx] <- cbind(dmat12[, 1], dmat12[, 2], dmat3, dmat4)
  }
  return(list("Dmat" = Dmat, "psi_mat" = psi_mat))
}

est_psi_nl <- function(obs_data){
  num_cov <- ncol(obs_data) - 2

  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))

  las_y  <- las_ydw_nl(obs_data$y, obs_data$delta,
                    obs_data[, w_cols])

  psi <- rep(NA, num_cov)

  for(cov_idx in 1:num_cov){
    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj_nl <- sl_me_nl(las_y(obs_data[, w_cols]), cov_col)
    mar_exp_vec <- log(marg_exp_wj_nl(cov_col)) #[obs_data$delta == 1]
    mme_wj <- mean(mar_exp_vec * cov_col) # [obs_data$delta == 1]
    mme_one <- mean(mar_exp_vec)

    psi_1 <- mme_wj
    psi_2 <- mme_one
    psi_3 <- mean(cov_col)
    psi_4 <- mean(cov_col ** 2)

    psi[cov_idx] <- (psi_1 - psi_2 * psi_3)/(psi_4 - psi_3 ** 2)
  }
  return(psi)
}

# exp_wj <- function(c_idx, ws){
#   betas <- c(2, 0.5, -3, rep(0, 7))
#   oth_vals <- betas[-c_idx]
#   inside <- as.matrix(ws[, - c_idx]) %*% oth_vals
#   marg_exp <- function(wj){
#     f_vec <- rep(NA, length(wj))
#     for(j_idx in 1:length(wj)){
#       f_vec[j_idx] <- mean(expit(inside + betas[c_idx] * wj[j_idx]))
#     }
#     return(f_vec)
#   }
#   return(marg_exp)
# }

psi_12_nl <- function(y, d, w, wj, pr_dw, epr_ywj, pr_yw,
                   exp_wj, exp_one){
  piece_one <- (y * d / (pr_dw(w)) + pr_yw(w)) / epr_ywj(wj) +
    log(epr_ywj(wj)) - 1 -
    d * pr_yw(w) / (pr_dw(w) * epr_ywj(wj))
  piece_one_w <- wj * piece_one - exp_wj
  piece_one <- piece_one - exp_one
  return(cbind(piece_one_w, piece_one))
}

# For reference:
# y = obs_data$y, d = obs_data$delta,
# w = obs_data[, w_cols],
# wj = obs_data[, w_cols[cov_idx]],
# pr_dw = las_del, epr_ywj = marg_exp_wj,
# pr_yw = las_y, exp_wj = mme_wj,
# exp_one = mme_one

delt_methd <- function(xes){
  x_1 <- xes[1] ; x_2 <- xes[2]
  x_3 <- xes[3] ; x_4 <- xes[4]
  grad <- c(
    1/(x_4 - x_3 ** 2),
    -x_3/((x_4 - (x_3) ** 2)),
    (2 * x_1 * x_3 - x_2 * x_4 - x_2 * x_3 ** 2)/((x_4 - x_3 ** 2) ** 2),
    -(x_1 - x_2 * x_3)/((x_4 - x_3 ** 2) ** 2)
  )
  return(grad)
}

calc_IC_mat_nl <- function(obs_data){
  dmat_and_psi <- get_Dmat_nl(obs_data = obs_data)
  D_array <- dmat_and_psi$Dmat
  Psi_mat <- dmat_and_psi$psi_mat
  grad_mat <- matrix(NA, nrow = nrow(Psi_mat), ncol = 4)
  for(c_idx in 1:nrow(Psi_mat)){
    grad_mat[c_idx, ] <- delt_methd(Psi_mat[c_idx, ])
  }
  IC <- matrix(NA, nrow = nrow(obs_data), ncol = nrow(grad_mat))
  for(p_idx in 1:nrow(grad_mat)){
    IC[, p_idx] <- D_array[, , p_idx] %*% t(grad_mat[c_idx, , drop = FALSE])
  }
  return(IC)
}

de_2_nl <- function(est_or_IC){
  if(est_or_IC == "est"){return(est_psi_nl)}
  if(est_or_IC == "IC"){return(calc_IC_mat_nl)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}




