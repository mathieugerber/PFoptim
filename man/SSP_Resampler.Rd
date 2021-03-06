% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SSP_Resampler}
\alias{SSP_Resampler}
\title{SSP resampling}
\usage{
SSP_Resampler(U,W)
}
\arguments{
\item{W}{A vector of normalized weights.}

\item{U}{A vector of points in (0,1) such that \code{length(U)=length(W)}.}
}
\value{
A vector of length N with elements in the set \eqn{\{1,...,N\}}, with \code{N=length(U)=length(W)}.
}
\description{
This function implements the SSP resampling algorithm \insertCite{gerber2019negative}{PFoptim}.
}
\details{
For efficiency reasons, \code{SSP_Resampler} does not perform checks on the supplied arguments.
}
\examples{
N<-100
W<-rbeta(N,0.5,2)
W<-W/sum(W)
J<-SSP_Resampler(runif(N),W)
}
\references{
\insertAllCited{}
}
\keyword{algorithms}
\keyword{resampling}
