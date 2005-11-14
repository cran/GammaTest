\name{dvec}
\alias{dvec}

\title{Creating delay vectors from time series}

\description{
Converts a time series dataset into an input/output format, where the output is in the last column.
This is required to run a gamma test analysis on time series.
}

\usage{
dvec(time.series, lag)
}

\arguments{
  \item{time.series}{Array, matrix or ts data to be converted.}
  \item{lag}{The maximum lag with which to create the inputs e.g. if lag=2 then the inputs will
  be lags 1 and 2.}
}

\value{
  \item{io}{The input/output dataset.}
}

\references{For papers and Gamma test related material visit \url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{fesearch}}}

\examples{
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
dvar <- dvec(ts.sim, 1)
}
\keyword{ts}
