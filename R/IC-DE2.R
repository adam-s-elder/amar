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
#' @param simp If the lasso function should simply return 0.5
#' @export
#'

las_dw <- function(delt_vec, obs_ws, simp = FALSE){
  if(simp){
    new_funct <- function(ws){
      return(rep(0.5, nrow(ws)))
    }
    return(new_funct)
  }else{
    lasso_fit <- glmnet::cv.glmnet(x = as.matrix(obs_ws), y= delt_vec,
                                   family = "binomial")
    new_funct <- function(ws){
      stats::predict(lasso_fit, newx = as.matrix(ws),
                     type = "response", s = "lambda.1se")[, ,drop = TRUE]
    }
    return(new_funct)
  }
}

smth_func <- function(x){
  ifelse(x > 0.9928119 | x < 0.007188064,
         exp(10 * (x - 0.5))/(1 + exp(10 * (x - 0.5))),
         x)
}

las_ydw <- function(y_vec, delt_vec, obs_ws, simp = FALSE){
  if(simp){
    l_mod_ydw <- function(ws){
      return(rep(0.5, nrow(ws)))
    }
    return(l_mod_ydw)
  }else{
    ys <- y_vec[delt_vec == 1]
    ws <- obs_ws[delt_vec == 1, ]
    lasso_fit <- glmnet::cv.glmnet(x = as.matrix(ws), y = ys,
                                   family = "binomial")
    l_mod_ydw <- function(ws){
      stats::predict(lasso_fit, newx = as.matrix(ws),
                     type = "response",
                     s = "lambda.1se")[, , drop = TRUE]
    }
    return(l_mod_ydw)
  }
}

sl_me <- function(m_e_v, w_j, simp = FALSE){
  if (simp) {
    sl_exp_wj <- function(w_j){
      return(rep(0.5, length(w_j)))
    }
    return(sl_exp_wj)
  }else{
    sl_mod <- stats::loess(m_e_v ~ w_j)

    sl_exp_wj <- function(w_j){
      new_data <- data.frame("wj" = w_j)
      stats::predict(sl_mod, newdata = new_data)
    }
    return(sl_exp_wj)
  }
}

sl_me_new <- function(m_e_v, w_j, simp = FALSE){
  if (simp){
    sl_exp_wj <- function(w_j){
      return(rep(0.5, length(w_j)))
    }
    return(sl_exp_wj)
  }else{
    sl_mod <- SuperLearner::SuperLearner(Y = m_e_v,
                                         X = data.frame("WJ" = w_j),
                                         ## Maybe try SL.polymars instead??
                                         SL.library = c("SL.glm", "SL.loess",
                                                        "SL.ranger"),
                                         family = "gaussian")

    sl_exp_wj <- function(w_j){
      new_data <- data.frame("WJ" = w_j)
      stats::predict(sl_mod, newdata = new_data,
                     onlySL = TRUE)$pred[, , drop = TRUE]
    }
    return(sl_exp_wj)
  }
}

get_Dmat <- function(obs_data, fit_delt = FALSE, fit_jnt = FALSE,
                     fit_marg = FALSE, hav_mis_prob = NULL){
  if (!is.null(hav_mis_prob)){
    true_mis_prob <- obs_data$pr_d_mis
    obs_data$pr_d_mis <- NULL
  }else{
    true_mis_prob <- NULL
  }
  num_cov <- ncol(obs_data) - 2
  num_obs <- nrow(obs_data)


  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))

  las_del <- las_dw(obs_data$delta, obs_data[, w_cols], simp = fit_delt)
  las_y  <- las_ydw(obs_data$y, obs_data$delta,
                    obs_data[, w_cols], simp = fit_jnt)
  Dmat <- array(NA,  dim = c(num_obs, 4, num_cov))
  psi_mat <- matrix(nrow = num_cov, ncol = 4)

  for(cov_idx in 1:num_cov){

    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj <- suppressWarnings(sl_me(las_y(obs_data[, w_cols]),
                                          cov_col, simp = fit_marg))
    mar_exp_vec <- log(smth_func(marg_exp_wj(cov_col))) # This subset may be needed [obs_data$delta == 1]
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
                     pr_dw = las_del, epr_ywj = marg_exp_wj,
                     pr_yw = las_y, exp_wj = mme_wj,
                     exp_one = mme_one,
                     true_mis_prob = true_mis_prob)
    dmat3 <- cov_col - mean(cov_col)
    cov_col_sq <- cov_col ** 2
    dmat4 <- cov_col_sq - mean(cov_col_sq)
    Dmat[, ,cov_idx] <- cbind(dmat12[, 1], dmat12[, 2],
                              dmat3, dmat4)
  }
  return(list("Dmat" = Dmat, "psi_mat" = psi_mat))
}

est_psi <- function(obs_data, comps = FALSE,
                    fit_marg = FALSE,
                    fit_jnt = FALSE,
                    p_m_k = FALSE){
  if (p_m_k){
    obs_data$pr_d_mis <- NULL
  }

  num_cov <- ncol(obs_data) - 2

  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))

  las_y  <- las_ydw(obs_data$y, obs_data$delta,
                    obs_data[, w_cols], simp = fit_jnt)

  psi <- rep(NA, num_cov)
  psi_mat <- matrix(NA, nrow = num_cov, ncol = 4)

  for (cov_idx in 1:num_cov) {
    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj <- sl_me(las_y(obs_data[, w_cols]),
                         cov_col, simp = fit_marg)
    mar_exp_vec <- log(smth_func(marg_exp_wj(cov_col))) #[obs_data$delta == 1]
    mme_wj <- mean(mar_exp_vec * cov_col) # [obs_data$delta == 1]
    mme_one <- mean(mar_exp_vec)

    psi_1 <- mme_wj
    psi_2 <- mme_one
    psi_3 <- mean(cov_col)
    psi_4 <- mean(cov_col ** 2)

    psi_mat[cov_idx, ] <- c(psi_1, psi_2, psi_3, psi_4)

    psi[cov_idx] <- (psi_1 - psi_2 * psi_3)/(psi_4 - psi_3 ** 2)
  }
  if(comps){
    return(list("psi" = psi, "psi_mat" = psi_mat))
  }else{
    return(psi)
  }
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

psi_12 <- function(y, d, w, wj, pr_dw, epr_ywj, pr_yw,
                   exp_wj, exp_one, true_mis_prob = NULL){
  if (!is.null(true_mis_prob)){
    prob_dw <- true_mis_prob
  }else{
    prob_dw <- pmax(0.25, pr_dw(w))
  }
  piece_one <- (y * d / prob_dw  + pr_yw(w)) /
    epr_ywj(wj) +
    log(smth_func(epr_ywj(wj))) - 1 -
    d * pr_yw(w) / (prob_dw * epr_ywj(wj))
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
    -x_3/(x_4 - x_3 ** 2),
    (2 * x_1 * x_3 - x_2 * x_4 - x_2 * x_3 ** 2)/((x_4 - x_3 ** 2) ** 2),
    -(x_1 - x_2 * x_3)/((x_4 - x_3 ** 2) ** 2)
  )
  return(grad)
}

calc_IC_mat <- function(obs_data, comps = FALSE, fit_delt = FALSE,
                        fit_jnt = FALSE, fit_marg = FALSE, hv_ms_pb = NULL){
  dmat_and_psi <- get_Dmat(obs_data = obs_data, fit_delt = fit_delt,
                           fit_jnt = fit_jnt, fit_marg = fit_marg,
                           hav_mis_prob = hv_ms_pb)
  D_array <- dmat_and_psi$Dmat
  Psi_mat <- dmat_and_psi$psi_mat
  grad_mat <- matrix(NA, nrow = nrow(Psi_mat), ncol = 4)
  for(c_idx in 1:nrow(Psi_mat)){
    grad_mat[c_idx, ] <- delt_methd(Psi_mat[c_idx, ])# +
                                      # apply(D_array[, , c_idx], 2, mean))
  }
  IC <- matrix(NA, nrow = nrow(obs_data), ncol = nrow(grad_mat))
  for(p_idx in 1:nrow(grad_mat)){
    IC[, p_idx] <- D_array[, , p_idx] %*% t(grad_mat[p_idx, , drop = FALSE])
  }
  if(comps){
    return(list("IC" = IC, "d_and_psi" = dmat_and_psi))
  }else{
    return(IC)
  }
}

de_2 <- function(est_or_IC){
  if(est_or_IC == "est"){return(est_psi)}
  if(est_or_IC == "IC"){return(calc_IC_mat)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}

oth_calc_IC <- function(obs_data){
  calc_IC_mat(obs_data = obs_data,
              comps = FALSE, fit_delt = FALSE,
              fit_jnt = FALSE, fit_marg = FALSE,
              hv_ms_pb = TRUE)
}

oth_est_psi <- function(obs_data){
  est_psi(obs_data = obs_data,
              comps = FALSE,
              fit_jnt = FALSE, fit_marg = FALSE,
              p_m_k = TRUE)
}

de_2_know_mis <- function(est_or_IC){
  if(est_or_IC == "est"){return(oth_est_psi)}
  if(est_or_IC == "IC"){return(oth_calc_IC)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}

# check_psi <- function(delt, y, wj, f = TRUE){
#   if (f){
#     return((delt * (2 * y - 1) * wj / 0.7) - log(2) * wj)
#     }
#   else{
#     return(delt * (2 * y - 1) / 0.7)
#   }
# }



