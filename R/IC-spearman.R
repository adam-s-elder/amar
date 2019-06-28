################################################################
#### Estimating Influence Curve and parameter for spearman corrlation
#### Adam Elder
#### This script will contain all of the various proceedures
#### used for estimating both the parameter, and the influence
#### curves used for estimating spearman correlation.
#### These functions will be given an entire dataset, and will
#### return the estimate for the parameter, and for each
#### estimated influence curve at each observation.
################################################################

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
  num_var <- ncol(observ)
  return(cor(observ[, 1], observ[, -1], method = "spearman"))
}

spearman <- list("est_IC" = est_influence_spearman, "est_param" = est_spearman)

