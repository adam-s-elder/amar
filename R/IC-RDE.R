expit <- function(x) exp(x) / (1 + exp(x))

rde_ic <- function(obs_data){
  y_vc <- obs_data[, "y"]
  wts <- obs_data[, "wt"]
  ws <- obs_data[, - which(colnames(obs_data) %in% c("y", "wt"))]
  fin_IC <- matrix(NA, nrow = length(y_vc), ncol = ncol(ws))

  for(cov_idx in 1:ncol(ws)){
    w_j <- ws[, cov_idx]
    glm_fit <- suppressWarnings(glm(y_vc ~ w_j, weights = wts,
                   family = binomial))
    beta_0 <- glm_fit$coefficients[1]
    beta_1 <- glm_fit$coefficients[2]
    jst_expit <- expit(beta_0 + beta_1 * w_j)
    expit_fracs <-
      jst_expit / (
        1 + exp(beta_0 + beta_1 * w_j)
      )

    exp_w2 <- mean(expit_fracs * (w_j ** 2) * wts)
    exp_w <- mean(expit_fracs * w_j * wts)
    exp_1 <- mean(expit_fracs * wts)

    c_b <- (exp_w2 * exp_1 - (exp_w ** 2)) ** (-1)
    ic_vals <- c_b * wts * (y_vc - jst_expit) * (w_j * exp_1 - exp_w)
    fin_IC[, cov_idx] <- ic_vals
  }
  return(fin_IC)
}

rde_est <- function(obs_data){
  y_vc <- obs_data[, "y"]
  wts <- obs_data[, "wt"]
  ws <- obs_data[, - which(colnames(obs_data) %in% c("y", "wt"))]
  returnVal <- rep(NA, ncol(ws))

  for(cov_idx in 1:ncol(ws)){
    w_j <- ws[, cov_idx]

    glm_fit <- suppressWarnings(glm(y_vc ~ w_j, weights = wts,
                   family = binomial))
    returnVal[cov_idx] <- glm_fit$coefficients[2]
  }
  return(returnVal)
}


rde <- function(est_or_IC){
  if(est_or_IC == "est"){
    return(rde_est)
  }else if(est_or_IC == "IC"){
        return(rde_ic)
  }else{
    stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}

## Check

# ws <- matrix(rnorm(30000), ncol = 3)
# probs <- expit(ws  %*% c(-1, 0, 2))
# y <- rbinom(n = nrow(probs), size = 1, prob = probs[, 1])
# wts <-   abs(rnorm(length(y))) + 1
# wts <- length(wts) * wts / sum(wts)
# obs_dat <- cbind(y, "wt" = wts, ws)
# my_est <- rde_est(obs_data = obs_dat)
# my_est
# my_ic <- rde_ic(obs_dat) / nrow(obs_dat)
# var_mat <- t(my_ic) %*% my_ic
# sqrt(diag(var_mat))
#
# for(cov_idx in 1:ncol(ws)){
#   print(summary(glm(y ~ ws[, cov_idx], weights = obs_dat[, "wt"],
#                     family = binomial))$coefficients[2, 1:2])
# }


