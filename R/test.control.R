#' Control function for the adaptive norm test
#'
#' @param n_peld_mc_samples Number of samples to be used in approximating the
#' estimated limiting distribution of the parameter estimate under the null.
#' Increasing this value reduces the approximation error of the test statistic.
#' @param nrm_type The type of norm to be used for the test.
#' Generally the l_p norm
#' @param perf_meas the preferred measure used to generate the test statistic.
#' @param pos_lp_norms The index of the norms to be considered.  For example if
#' we use the l_p norm, norms_indx specifies the different p's to try.
#' @param ld_est_meth String indicating method for estimating the limiting
#' distribution of the test statistic parametric bootstrap or permutation.
#' @param ts_ld_bs_samp The number of test statistic limiting distribution
#' bootstrap samples to be drawn.
#' @param other_output A vector indicating additional data that should be
#' returned. Currently only \code{"var_est"} and \code{data} is supported.
#' @param ... Other arguments needed in other places.
#' @return A list that provide controls for \code{mv_pn_test} (specified by the
#' arguments passed to \code{test.control}).
#'
#' @examples
#' test.control()
#'
#' @export

test.control <- function(
  n_peld_mc_samples = 300,
  nrm_type = "lp",
  perf_meas = "est_acc",
  pos_lp_norms = c(1, 2, 3, "max"),
  ld_est_meth = "par_boot",
  ts_ld_bs_samp = 250,
  other_output = c(),
  ... ## RENAME
) {
  if (length(other_output) > 0) {
    w_n_m <-
      which(is.na(match(other_output, c("var_est", "obs_data"))))
    if (length(w_n_m) > 0) {
      warning(paste0(
        "The following output options are not supported",
        " by other_output: \n ",
        paste0(other_output[w_n_m], ", ", collapse = ""),
        "\n These arguments will be ignored."))
    }
  }
  formal_args <- formals(sys.function())
  dot_args <- list(...)
  p <- .get.args(formal_args, dot_args)
  if (is.character(p$perf_meas)) {
    if (!(p$perf_meas %in% c("pval", "est_acc", "mag"))) {
      stop(paste0("The control argument perf_meas must be either ",
                  "'pval', 'mag', or 'est_acc' if it is a string."))
    }
    p$perf_meas <- list(
      "est_acc" = accept_rate,
      "pval" = pval_for_mag,
      "mag" = mag_for_pow
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
