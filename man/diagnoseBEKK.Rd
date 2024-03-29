% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEKK.R
\name{diagnoseBEKK}
\alias{diagnoseBEKK}
\title{Diagnose BEKK process estimation}
\usage{
diagnoseBEKK(estimation)
}
\arguments{
\item{estimation}{The return value of the \code{mvBEKK.est} function}
}
\value{
Nothing special
}
\description{
Provides diagnostics for a BEKK process estimation.
}
\details{
This procedure provides console output and browsable plots for a
given BEKK process estimation. Therefore, it is meant to be
interactive as the user needs to proceed by pressing \code{c} on
the keyboard to see each plot one-by-one.
}
\examples{
## Simulate series:
simulated = simulateBEKK(2, 1000, c(1,1))

## Prepare the matrix:
simulated = do.call(cbind, simulated$eps)

## Estimate with default arguments:
estimated = BEKK(simulated)

\dontrun{
## Show diagnostics:
diagnoseBEKK(estimated)
}

}
