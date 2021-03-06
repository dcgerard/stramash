% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mix_like.R
\name{quantile_to_mix}
\alias{quantile_to_mix}
\title{Approximate the density function of any distribution given its
quantile function.}
\usage{
quantile_to_mix(q_func, gridsize = 20, eps = 10^-6)
}
\arguments{
\item{q_func}{A quantile function.}

\item{gridsize}{The number of points to approximate.}

\item{eps}{A positive numeric. The tolerance from 0 or 1 to start
the grid.}
}
\description{
Approximate the density function of any distribution given its
quantile function.
}

