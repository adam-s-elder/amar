% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_test.R
\name{cv_test}
\alias{cv_test}
\title{CV-test runs a cross-validation, parametric bootstrap test. This function returns
an approximate p-value for the specified test statistic.}
\usage{
cv_test(obs_data, param_est = NULL, control = test.control())
}
\arguments{
\item{obs_data}{The observed data to be used for finding the optimal
norm (training), and finding the test statistic (testing).  Similar to
above, each row is an observation and each column corresponds to either
the outcome (first column) or a covariate.}

\item{param_est}{Function used to estimate the parameter and corresponding
influence curve.}

\item{control}{List used to define controls for test.}
}
\value{
learned test statistic for a single fold of data
}
\description{
CV-test runs a cross-validation, parametric bootstrap test. This function returns
an approximate p-value for the specified test statistic.
}
