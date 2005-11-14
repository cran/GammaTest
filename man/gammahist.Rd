\name{gammahist}
\alias{gammahist}

\title{The Gamma Histogram}

\description{
Plots the Gamma statistics returned by a full embedding search in the form of a histogram.
}

\usage{
gammahist(fe.results, ...)
}

\arguments{
  \item{fe.results}{The returned results from a full embedding search.}
  \item{...}{Graphical parameters (e.g. setting the colour)}
}
\details{
The Gamma histogram is a useful way to represent a full embedding search. Here,
the range of Gamma statistics is divided into bins and the frequency of the Gamma statistics
per-bin is plotted vertically. The left-most area on the histogram is known as the low Gamma region, 
conversly, the right-most area on the histogram is known as the high Gamma region.
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{fesearch}} \code{\link{hist}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
gfts.sim <- dvec(ts.sim, 3)
my.fe	<- fesearch(gfts.sim) # Run full embedding on 10 previous lags
gammahist(my.fe, col="red")
}
\keyword{hplot}
