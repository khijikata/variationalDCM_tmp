% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmdcm_vb_diff_q.R
\name{hmdcm_diff_q}
\alias{hmdcm_diff_q}
\title{estimates attributes mastery patterns for hidden Markov DCM.}
\usage{
hmdcm_diff_q(
  X,
  Q,
  A_0 = NULL,
  B_0 = NULL,
  delta_0 = NULL,
  ommega_0 = NULL,
  max_it = 500,
  epsilon = 0.001,
  Test_versions,
  test_order,
  model = "General",
  random_start = FALSE
)
}
\arguments{
\item{X}{N by J matrix containing 0 or 1}

\item{Q}{J by k matrix containing 0 or 1}

\item{max_it}{Maximum number of iterations}

\item{epsilon}{convergence tolerance for iterations}

\item{random_start}{}
}
\description{
\code{hmdcm_diff_q()} returns variational Bayesian estimates for hidden
Markov DCM.
}
\references{
Yamaguchi, K., & Martinez, A. J. (2021). Variational Bayesian
  Inference Posterior Approximation Algorithm for Hidden Markov Diagnostic
  Classification Models. \url{https://doi.org/10.31234/osf.io/28jkf}.
}
