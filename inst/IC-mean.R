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

mean_f <- function(est_or_IC){
  if(est_or_IC == "est"){return(est_influence_mean)}
  if(est_or_IC == "IC"){return(est_mean)}
  else{stop("You must specify if you want the estimate of the parameter (est),
            or the influence curve (IC)")}
}
