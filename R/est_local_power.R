#' Generate an estimate of the limiting distribution of
#' the vector covariate specific parameter estimates
#'
#' @param boot_data bootstrap estimates of the (centered) limiting distribution. Data where
#'  columns correspond to different covariates, and rows are independent observations.
#' @param lp the lp norm to be used
#' @param method either "average" for an average local power over a sphere of local
#' alternatives, or "min" to get the minimum power over a sphere of local alternatives
#' @param tries the number of places to check local alternatives either for averages or minima
#' @param magnitude geometric distance away from zero of the local alternative
#'
#' @return Either the average, or minimum power for the specified
#' lp norm for all local alternatives.
#'
#' @export

est_loc_pow <- function(boot_data, lp = "max", method = "average",
                        tries, magnitude = 3){
  n_covs <- ncol(boot_data)
  n_obs  <- nrow(boot_data)
  crits <- rep(0, length(lp))
  for(p_indx in 1:length(lp)){
    ##Solving for the estimated 95% quantile of the l_p norm of
    ##the centered limiting distribution of parameter estiamtes.
    norm_boot <- apply(boot_data, 1, l_p_norm, p = lp[p_indx])
    crits[p_indx] <- stats::quantile(norm_boot, 0.95)
  }
  if (method == "average"){
   pows <- matrix(-1, nrow = tries, ncol = length(lp))
   for(t_indx in 1:tries){
     local_alt_t <- (magnitude) * r_unif_sphere(n_covs)
     shift_distr <- sweep(boot_data, 2, local_alt_t, "+")
     for(p_indx in 1:length(lp)){
       norm_shift_distr <- apply(shift_distr, 1, l_p_norm, p = lp[p_indx])
       pows[t_indx, p_indx] <- mean(as.numeric(norm_shift_distr > crits[p_indx]))
     }
   }
   return(apply(pows, 2, mean))

  }else if (method == "minimum"){
    minpows <- rep(0, length(lp))
    for(p_indx in 1:length(lp)){
      obj_fun <- function(dir){
        n_dir <- dir/(sqrt(sum(dir ** 2)))
        shift <- n_dir * magnitude
        shift_distr <- sweep(boot_data, 2, shift, "+")
        norm_shift_distr <- apply(shift_distr, 1, l_p_norm, p = lp[p_indx])
        return( - mean(as.numeric(norm_shift_distr > crits[p_indx])))
      }
      opt_solv <- GA::ga(type = "real-valued", fitness = obj_fun, lower = rep(-1, n_covs),
         upper = rep(1, n_covs), popSize = tries, run = 10, monitor = FALSE)
      minpows[p_indx] <- opt_solv@fitnessValue
    }
    return(- minpows)
  }
}

