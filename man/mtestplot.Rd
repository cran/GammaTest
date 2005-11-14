\name{mtest.plot}
\alias{mtest.plot}

\title{Plotting the results of the M-test}

\description{
Plots the M-test graph.
}

\usage{
mtest.plot(mt.results, stat="Gamma", ...)
}

\arguments{
  \item{mt.results}{The returned results from an M-test.}
  \item{stat}{Either \code{"Gamma"} to plot the Gamma statistics, or \code{"Vratio"} to plot the V ratios.} 
  \item{...}{Graphical parameters (e.g. setting the colour)}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}



\examples{
# Noisy Sine wave example
x <- seq(length=500, from=-2*pi, to=2*pi)
y <- sin(x) + rnorm(500, sd=sqrt(0.075)) # Set the variance of the noise to 0.075.
xy <- data.frame(x,y)			 # Create an input/output dataset.
mt <- mtest(xy)
mtest.plot(mt)
}
\keyword{hplot}

