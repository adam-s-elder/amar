################################################################
#### Estimating Influence Curve and parameter
#### for pearson corrlation
#### Adam Elder
#### This script will contain all of the various proceedures
#### used for estimating both the parameter, and the influence
#### curves used for estimating spearman correlation.
#### These functions will be given an entire dataset, and will
#### return the estimate for the parameter, and for each
#### estimated influence curve at each observation.
################################################################

#### The following function will take a set of observations, and return the
#### estimated covariance matrix for the estimates of the pearson correlation.
#### Estimates of the covariance are generated using the empirical influence
#### function.  The first column of your data should correspond to the
#### variable of interest (the variable for which pearson correlation is
#### calculated.

est_influence_pearson <- function(observ, trans = "none"){
  n <- nrow(observ)
  num_cov <- ncol(observ) - 1
  means <- colMeans(observ)
  y_mean <- means[1]
  cent_obs <- observ - matrix(rep(colMeans(observ), each = n), nrow = n)
  cent_obs_sqrd <- cent_obs^2
  sigmas <- colSums(cent_obs_sqrd)/(n - 1)
  covs <- as.numeric(crossprod(cent_obs[, 1], cent_obs[, -1]))/(n - 1)
  y_var <- sigmas[1]
  cent_cov <- cent_obs[, -1] * cent_obs[, 1] - rep(covs, each = n)
  cent_var <- cent_obs_sqrd - rep(sigmas, each = n)
  psi_1 <- matrix(rep(sigmas[-1], each = n), nrow = n) *
    cent_var[, 1] + y_var * cent_var[, -1]
  ic <- (cent_cov - psi_1 * rep(covs/(y_var * sigmas[-1]), each = n)) *
    rep((1/sqrt(y_var * sigmas[-1])), each = n)
  if(trans == "none"){
    return(ic)
  }else if(trans == "tsqd"){
    num_var <- ncol(observ)
    rho <- cor(observ[, 1], observ[, -1], method = "pearson")[1, ]
    mat_mult <- diag(2 * rho / (1 - rho ** 2) ** 2)
    return(ic  %*% t(mat_mult))
  }
}
pearson <- list("est_IC" = est_influence_pearson, "est_param" = est_pearson)

est_pearson <- function(observ, trans = NULL){
  num_var <- ncol(observ)
  rho <- cor(observ[, 1], observ[, -1], method = "pearson")
  if(trans == "tsqd"){return(rho ** 2 / (1 - rho ** 2))}
  else(return(rho))
}









