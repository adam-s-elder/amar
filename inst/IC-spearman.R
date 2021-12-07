#### The following function will take a set of observations, and return the
#### estimated covariance matrix for the estimates of the spearman correlation.
#### Estimates of the covariance are generated using the empirical influence
#### function.  The first column of your data should correspond to the
#### variable of interest (the variable for which spearman correlation is
#### calculated.


est_influence_spearman <- function(observ) {
  n <- nrow(observ)
  num_cov <- ncol(observ) - 1
  est_IC <- matrix(NA, nrow = n, ncol = num_cov)
  y_vals <- observ[, 1]
  Y_ecdf <- stats::ecdf(y_vals)
  y_ecdf_vals <- Y_ecdf(y_vals)
  for(i in 2:(num_cov + 1)){
    x_vals <- observ[, i]
    X_ecdf <- stats::ecdf(x_vals)
    x_ecdf_vals <- X_ecdf(x_vals)
    xy_constant <- -3 * mean(x_ecdf_vals * y_ecdf_vals)
    ic_fun <- function(x, y){xy_constant + X_ecdf(x) * Y_ecdf(y) +
        mean(x_ecdf_vals * ifelse(y_vals >= y, 1, 0)) +
        mean(y_ecdf_vals * ifelse(x_vals >= x, 1, 0))}
    est_IC[, i - 1] <- mapply(ic_fun, x_vals, y_vals)
  }
  return(12 * est_IC)
}

est_spearman <- function(observ){
  return(stats::cor(observ[, 1], observ[, -1], method = "spearman"))
}

spearman <- function(est_or_IC){
  if(est_or_IC == "est"){return(est_spearman)}
  if(est_or_IC == "IC"){return(est_influence_spearman)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}
