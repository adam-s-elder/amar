#' Control function for the adaptive norm test
#'
#' @param n_peld_mc_samples Number of samples to be used in approximating the
#' estimated limiting distribution of the parameter estimate under the null.
#' Increasing this value reduces the approximation error of the test statistic.
#' @param nrm_type The type of norm to be used for the test.
#' Generally the l_p norm
#' @param perf_meas the prefered measure used to generate the test statistic.
#' @param pos_lp_norms The index of the norms to be considered.  For example if
#' we use the l_p norm, norms_indx specifies the different p's to try.
#' @param ld_est_meth String indicating method for estimating the limiting
#' distribution of the test statistic parametric bootstrap or permutation.
#' @param ts_ld_bs_samp The number of test statistic limiting distribution
#' bootstrap samples to be drawn.
#' @param more_info Boolean indicating if the function should return more
#' information that just the p-value.  When true, the chosen norm index,
#' the bonferroni based p-value, the test statistic,
#' and the estimated distribution of the test statistic will be returned.
#' @param ret_cov_mat Boolean indicating if the function should also return
#' the estimated covariance matrix.  This can somewhat slow computation and
#' result in large file sizes at higher dimensions. The argument `more_info`
#' must be set to `TRUE` for this option to be possible.
#' @param ... Other arguments needed in other places.
#'
#' @export

test.control <- function(
  n_peld_mc_samples = 300,
  nrm_type = "lp",
  perf_meas = "est_acc",
  pos_lp_norms = c(1, 2, 3, "max"),
  ld_est_meth = "par_boot",
  ts_ld_bs_samp = 250,
  more_info = TRUE,
  ret_cov_mat = FALSE,
  ... ## RENAME
) {
  formal_args <- formals(sys.function())
  dot_args <- list(...)
  p <- .get.args(formal_args, dot_args)
  if (is.character(p$perf_meas)) {
    if (!(p$perf_meas %in% c("pval", "est_acc", "mag"))) {
      stop(paste0("The control argument perf_meas must be either ",
                  "'pval', 'mag', or 'est_acc' if it is a string."))
    }
    p$perf_meas <- list(
      "est_acc" = amp::accept_rate,
      "pval" = amp::pval_for_mag,
      "mag" = amp::mag_for_pow
    )[[p$perf_meas]]
  }
  return(p)
}

.get.args <- function(formal.args, dot.args) {
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in seq_len(length(dot.args))) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}
