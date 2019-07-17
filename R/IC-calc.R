######################################################################
#### Estimating Influence Curve and parameter for various parameters. 
#### Adam Elder
#### This script will contain all of the various proceedures
#### used for estimating both the parameter, and the influence
#### curves used for estimating spearman correlation, 
#### pearson correlation and mean.
#### These functions will be given an entire dataset, and will
#### return the estimate for the parameter, and for each
#### estimated influence curve at each observation.
######################################################################

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

est_pearson <- function(observ, trans = NULL){
    num_var <- ncol(observ)
    rho <- cor(observ[, 1], observ[, -1], method = "pearson")
    if(trans == "tsqd"){
      return(rho ** 2 / (1 - rho ** 2))
    }else{
      return(rho)
    }
}

pearson <- list("est_IC" = est_influence_pearson, "est_param" = est_pearson)


#### The following function will take a set of observations, and return the
#### estimated covariance matrix for the estimates of the spearman correlation.
#### Estimates of the covariance are generated using the empirical influence
#### function.  The first column of your data should correspond to the
#### variable of interest (the variable for which spearman correlation is
#### calculated.


est_influence_spearman <- function(observ){
  n <- nrow(observ)
  num_cov <- ncol(observ) - 1
  est_IC <- matrix(NA, nrow = n, ncol = num_cov)
  y_vals <- observ[, 1]
  Y_ecdf <- ecdf(y_vals)
  y_ecdf_vals <- Y_ecdf(y_vals)

  for(i in 2:(num_cov + 1)){
    x_vals <- observ[, i]
    X_ecdf <- ecdf(x_vals)
    x_ecdf_vals <- X_ecdf(x_vals)
    xy_constant <- -3 * mean(x_ecdf_vals * y_ecdf_vals)
    ic_fun <- function(x, y){xy_constant + X_ecdf(x) * Y_ecdf(y) +
        mean(x_ecdf_vals * ifelse(y_vals >= y, 1, 0)) +
        mean(y_ecdf_vals * ifelse(x_vals >= x, 1, 0))}
    est_IC[, i - 1] <- mapply(ic_fun, x_vals, y_vals)
  }
  return( 12 * est_IC)
}

est_spearman <- function(observ){
  return(cor(observ[, 1], observ[, -1], method = "spearman"))
}

spearman <- list("est_IC" = est_influence_spearman, "est_param" = est_spearman)

## Influence function and parameter esitmation for the mean.

est_influence_mean <- function(observ){
  col_means <- colMeans(x = observ)
  infl <- sweep(x = observ, MARGIN = 2,
                STATS = col_means, FUN = "-")
  return(infl)
}

est_mean <- function(observ){
  col_means <- colMeans(x = observ)
  return(col_means)
}

mean_f <- list("est_IC" = est_influence_mean, "est_param" = est_mean)



#### Placing all of the different influence function and parameter estimation
#### inside of a single function

est_infl_funcs <- list("spearman" = spearman, "pearson" = pearson,
                       "mean" = mean_f)


