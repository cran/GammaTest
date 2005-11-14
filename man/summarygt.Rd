\name{gtsummary}
\alias{gtsummary}

\title{The Gamma Test: Summary of Results}

\description{
Shows the key statistics calculated by a Gamma test including the Mask used, Gamma statistic, gradient
of the regression line and the V ratio.
}

\usage{
gtsummary(gt.list)
}

\arguments{
  \item{gt.list}{The returned results from a Gamma test.}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
gfts.sim <- dvec(ts.sim, 1)
g <- gammatest(gfts.sim)
gtsummary(g)

# Noisy Sine wave example
x <- seq(length=500, from=-2*pi, to=2*pi)
y <- sin(x) + rnorm(500, sd=sqrt(0.075)) # Set the variance of the noise to 0.075.
xy <- data.frame(x,y)			 # Create an input/output dataset.
noise.est <- gammatest(xy)		 # Run Gamma test.
gtsummary(noise.est)			 # Produce a summary of the results
}
\keyword{ts}
