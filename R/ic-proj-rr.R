#' Estimate both the parameter, and the influence
#' curves used for estimating the projected risk ratio. The first column
#' of your data should correspond to the variable of interest.
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
#'
#' @export

ic.proj.rr <-  function(obs_data, what = "both", control = NULL) {
  if (!(what %in% c("ic", "est", "both"))) {
    stop("what must be one of ic (influence curve), est (estimate), or both")
  }

  ret <- list()
  fit_delt <- control$ic_fit_delt
  fit_jnt <- control$ic_fit_jnt
  fit_marg <- control$ic_marg
  fit_delt <- control$ic_fit_delt
  hav_mis_prob <- control$ic_hav_mis_prop
  prob_mis_known <- control$ic_prob_mis_known
  if (is.null(fit_marg)) fit_marg <- FALSE
  if (is.null(fit_jnt)) fit_jnt <- FALSE
  if (is.null(fit_delt)) fit_delt <- FALSE
  if (is.null(prob_mis_known)) prob_mis_known <- FALSE
  if (prob_mis_known) {obs_data$pr_d_mis <- NULL}
  if (!is.null(hav_mis_prob)){
    true_mis_prob <- obs_data$pr_d_mis
    obs_data$pr_d_mis <- NULL
  }else{
    true_mis_prob <- NULL
  }

  num_obs <- nrow(obs_data)
  num_cov <- ncol(obs_data) - 2
  w_cols <- which(!colnames(obs_data) %in% c("delta", "y"))
  las_y  <- las_ydw(obs_data$y, obs_data$delta,
                    obs_data[, w_cols], simp = fit_jnt)
  las_del <- las_dw(obs_data$delta, obs_data[, w_cols], simp = fit_delt)
  las_y  <- las_ydw(obs_data$y, obs_data$delta,
                    obs_data[, w_cols], simp = fit_jnt)
  est <- rep(NA, num_cov)
  Dmat <- array(NA,  dim = c(num_obs, 4, num_cov))
  psi_mat <- matrix(NA, nrow = num_cov, ncol = 4)

  for (cov_idx in 1:num_cov) {
    cov_col <- obs_data[, w_cols[cov_idx]]
    marg_exp_wj <- suppressWarnings(sl_me(las_y(obs_data[, w_cols]),
                                          cov_col, simp = fit_marg))
    mar_exp_vec <- log(smth_func(marg_exp_wj(cov_col)))
    mme_wj <- mean(mar_exp_vec * cov_col)
    mme_one <- mean(mar_exp_vec)

    psi_1 <- mme_wj
    psi_2 <- mme_one
    psi_3 <- mean(cov_col)
    psi_4 <- mean(cov_col ** 2)
    psi_mat[cov_idx, ] <- c(psi_1, psi_2, psi_3, psi_4)
    est[cov_idx] <- (psi_1 - psi_2 * psi_3)/(psi_4 - psi_3 ** 2)
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
  grad_mat <- matrix(NA, nrow = nrow(psi_mat), ncol = 4)
  for (c_idx in seq_len(nrow(psi_mat))) {
    grad_mat[c_idx, ] <- ic_delt_methd(psi_mat[c_idx, ])
  }
  IC <- matrix(NA, nrow = nrow(obs_data), ncol = nrow(grad_mat))
  for (p_idx in seq_len(nrow(grad_mat))) {
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

ic_delt_methd <- function(xes){
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

las_dw <- function(delt_vec, obs_ws, simp = FALSE){
  if (simp) {
    new_funct <- function(ws){
      return(rep(0.5, nrow(ws)))
    }
    return(new_funct)
  }else{
    lasso_fit <- glmnet::cv.glmnet(
      x = as.matrix(obs_ws), y = delt_vec, family = "binomial"
      )
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
  if (simp) {
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
