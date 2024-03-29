% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rr-msm-ic.R
\name{rr.msm.ic}
\alias{rr.msm.ic}
\title{Estimate both the parameter, and the influence
curves used for estimating the projected risk ratio from a working
marginal structural mean model.}
\usage{
rr.msm.ic(obs_data, what = "both", control = NULL)
}
\arguments{
\item{obs_data}{the observed data.  The first column should be the outcome.}

\item{what}{the desired return value. Should be one of `"ic"`
(infludence curve), `"est"` (estimate), or `"both"`.}

\item{control}{any other control parameters to be passed to the estimator.}
}
\description{
The function takes a dataset, and will returns the estimate for the
parameter, and or the estimated influence curve at each observation.
The first column of obs_data should be the binary outcome of interest.
}
\details{
Most of this code has been copied and slightly modified from the
ltmle package and is copyrighted by Joshua Schwab under the terms of
the GPL-2 license.

Obtain an estimator of the probability delta = 1 given w
}
