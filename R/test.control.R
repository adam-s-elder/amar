#' Control function for the adaptive norm test
#'
#' @param more_info Boolean indicating if the test should return more
#' information that just the p-value.  When true, the chosen norm index,
#' the bonferroni based p-value, the test statistic,
#' and the estimated distribution of the test statistic will be returned.
#' @param n_mc_samples Number of samples to be used in estimating the limiting
#' distribution of the test statistic under the null.
#' @param nrm_type The type of norm to be used for the test.
#' Generally the l_p norm
#' @param perf_meas the prefered measure used to generate the test statistic.
#' @param pos_lp_norms The index of the norms to be considered.  For example if
#' we use the l_p norm, norms_indx specifies the different p's to try.
#' @param show_hist set to True to see a histogram of the test statistic's
#' limiting distribution MC draws compared to the estimated parameter.
#' @param ts_ld_bs_samp The number of test statistic limiting distribution
#' bootstrap samples to be drawn.
#' @param ld_est_meth String indicating method for estimating the limiting
#' distribution of the test statistic parametric bootstrap or permutation.
#' @param ... Other arguments needed in other places.
#'
#' @export

test.control <- function(
  more_info = TRUE,
  n_mc_samples = 300, ## RENAME
  nrm_type = "lp",
  perf_meas = "est_acc",
  pos_lp_norms = c(1, 2, 3, "max"),
  show_hist = FALSE,
  ld_est_meth = "par_boot",
  ts_ld_bs_samp = 250, ... ## RENAME
) {
  formal_args <- formals(sys.function())
  dot_args <- list(...)
  p <- .get.args(formal_args, dot_args)
  if (is.character(p$perf_meas)) {
    p$perf_meas <- list(
      "est_acc" = amar::accept_rate,
      "pval" = amar::pval_for_mag,
      "mag" = amar::mag_for_pow
    )[[p$perf_meas]]
  }
  return(p)
}


## Pulled from epiModel
.get.args <- function (formal.args, dot.args) {
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}
