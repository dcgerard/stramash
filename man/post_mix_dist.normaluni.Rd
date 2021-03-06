% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mix_like.R
\name{post_mix_dist.normaluni}
\alias{post_mix_dist.normaluni}
\title{Normal-mixture (prior) - uniform-mixture (error) posterior calculation.}
\usage{
post_mix_dist.normaluni(g, betahat, errordist)
}
\arguments{
\item{g}{The prior distribution. Either of class \code{"normalmix"}
or \code{"unimix"}.}

\item{betahat}{The observations, a vector of numerics.}

\item{errordist}{A list of mixture distributions.}
}
\description{
Normal-mixture (prior) - uniform-mixture (error) posterior calculation.
}
\author{
David Gerard
}

