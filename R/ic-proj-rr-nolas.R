#' Estimate both the parameter, and or the influence
#' curves used for estimating the projected risk ratio (not using lasso).
#' The first column of your data should correspond to the variable of interest.
#'
#' @param obs_data the observed data.  The first column should be the outcome.
#' @param what the desired return value. Should be one of \code{"ic"}
#' (influence curve), \code{"est"} (estimate), or \code{"both"}.
#' @param control any other control parameters to be passed to the estimator.
#' @return If \code{what} is
#'
#' - \code{"est"}, then return the estimated parameter.
#'
#' - \code{"ic"}, then return the estimated IC of the parameter estimate.
#'
#' - \code{"both"}, then return both the parameter estimate and
#' corresponding estimated IC.
#'
#' @examples
#' # not run (make sure to load in SuperLearner if running)
#' # set.seed(1010)
#' # fake_dat <- data.frame(y = rbinom(100, size = 1, prob = 0.5),
#' #                       delta = rbinom(100, size =  1, prob = 0.5),
#' #                      w = matrix(rnorm(500), ncol = 5))
#' # ic.proj.rr.nolas(fake_dat)
#'
#' @export

ic.proj.rr.nolas <-  function(obs_data, what = "both", control = NULL){
  if (!(what %in% c("ic", "est", "both"))) {
    stop("what must be one of ic (influence curve), est (estimate), or both")
  }
  ret <- list()
  num_cov <- ncol(obs_data) - 2
  num_obs <- nrow(obs_data)
  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))

  las_del <- las_dw_nl(obs_data$delta, obs_data[, w_cols])
  las_y  <- las_ydw_nl(obs_data$y, obs_data$delta, obs_data[, w_cols])
  Dmat <- array(NA,  dim = c(num_obs, 4, num_cov))
  psi_mat <- matrix(NA, nrow = num_cov, ncol = 4)

  est <- rep(NA, num_cov)

  for(cov_idx in 1:num_cov){

    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj_nl <- suppressWarnings(
      sl_me_nl(las_y(obs_data[, w_cols]), cov_col)
    )
    mar_exp_vec <- log(marg_exp_wj_nl(cov_col)) # This subset may be needed [obs_data$delta == 1]
    mme_wj <- mean(mar_exp_vec * cov_col)
    # This subset may be needed [obs_data$delta == 1]
    mme_one <- mean(mar_exp_vec)

    psi_1 <- mme_wj
    psi_2 <- mme_one
    psi_3 <- mean(cov_col)
    psi_4 <- mean(cov_col ** 2)
    if (what %in% c("est", "both")) {
      est[cov_idx] <- (psi_1 - psi_2 * psi_3)/(psi_4 - psi_3 ** 2)
    }

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
  grad_mat <- matrix(NA, nrow = nrow(psi_mat), ncol = 4)
  for(c_idx in seq_len(nrow(psi_mat))) {
    grad_mat[c_idx, ] <- ic.delt_methd(psi_mat[c_idx, ])
  }
  IC <- matrix(NA, nrow = nrow(obs_data), ncol = nrow(grad_mat))
  for(p_idx in seq_len(nrow(grad_mat))) {
    IC[, p_idx] <- Dmat[, , p_idx] %*% t(grad_mat[c_idx, , drop = FALSE])
  }
  if (what %in% c("est", "both")) {
    os_cor <- as.vector(apply(IC, 2, mean))
    ret$est <- os_cor + as.vector(est)
    ret$onestep.correction <- os_cor
  }
  if (what %in% c("ic", "both")) {
    ret$ic <- IC
  }
  return(ret)
}


# For reference:
# y = obs_data$y, d = obs_data$delta,
# w = obs_data[, w_cols],
# wj = obs_data[, w_cols[cov_idx]],
# pr_dw = las_del, epr_ywj = marg_exp_wj,
# pr_yw = las_y, exp_wj = mme_wj,
# exp_one = mme_one

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

ic.delt_methd <- function(xes){
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

las_dw_nl <- function(delt_vec, obs_ws) {
  sl_fit <- SuperLearner::SuperLearner(
    Y = delt_vec, X = obs_ws,
    ## Maybe try SL.polymars instead??
    SL.library = c("SL.glm", "SL.loess", "SL.ranger"),
    family = "gaussian")

  new_funct <- function(ws) {
    pmax(stats::predict(sl_fit, newdata = ws,
                        onlySL = TRUE)$pred[, , drop = TRUE],
         0.01)
  }
  return(new_funct)
}

las_ydw_nl <- function(y_vec, delt_vec, obs_ws) {
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

