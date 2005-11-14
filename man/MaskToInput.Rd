\name{mask2input}
\alias{mask2input}

\title{Create an Input/Output file for time series.}

\description{
To run a Gamma test on a time series it needs to be converted into an input/output
\code{data.frame}, where the output is in the last column.
}

\usage{
mask2input(mask, timeseries, multiple=FALSE)
}

\arguments{
  \item{mask}{A vector of 1's and 0's representing the input inclusion/exclusion.} 
  \item{timeseries}{The time series we wish to make into an input/output dataset.}
  \item{multiple}{Logical indicating whether we are dealing with muliple time series. If we are then 
  the output is not included in the last column}
}

\value{
  \item{newIO}{The input/output dataset.}
}

\references{For papers and Gamma test related material visit \url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{fesearch}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
mask <- c(1)
ts.sim.io <- mask2input(mask, as.matrix(ts.sim))	

# Example on an AR(2) process
ts.sim 	<- arima.sim(500, model=list(ar=c(0.9, -0.8)), sd=sqrt(1))
mask <- c(1,1)
ts.sim.io <- mask2input(mask, as.matrix(ts.sim))

# Example on an AR(3) process
X 	<- arima.sim(500, model=list(ar=c(-0.2, 0.3, -0.2)), sd=sqrt(1))
mask <- c(0,1,1)	# Leave out X_{t-1}
X.io <- mask2input(mask, as.matrix(X)) # Write a data.frame with {(X_{t-2}, X_{t-3}), X_t}
}
\keyword{ts}
