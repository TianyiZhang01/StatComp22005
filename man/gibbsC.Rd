% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbsC}
\alias{gibbsC}
\title{a Gibbs sampler using Rcpp}
\usage{
gibbsC(inix, iniy, mu1, mu2, sigma1, sigma2, rho, n = 10000L, burn = 1000L)
}
\arguments{
\item{inix}{the initial value of the markov chain}

\item{iniy}{the initial value of the markov chain}

\item{mu1}{the mean of X of the target distribution}

\item{mu2}{the mean of Y of the target distribution}

\item{sigma1}{the standard deviation of X of the target distribution}

\item{sigma2}{the standard deviation of Y of the target distribution}

\item{rho}{the correlation between X and Y of the target distribution}

\item{n}{the quantity of samples would like to generate, with initial 10000}

\item{burn}{the quantity of samples dropped before the chain go to steady}
}
\value{
a random sample of size \code{n+burn}
}
\description{
a Gibbs sampler using Rcpp, which comes from my homework 10.
}
\examples{
\dontrun{
temp = gibbsC(100, 10, 0, 0, 1, 1, 0.9)[-(1:1000), ]
par(mfrow=c(2,1));
plot(temp[-(1:1000),1],type='l')
plot(temp[-(1:1000),2],type='l')
}
}
