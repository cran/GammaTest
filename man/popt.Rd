\name{popt}
\alias{popt}

\title{Selecting p}

\description{
Calculates the MSE of the Gamma regression fit over increasing p.
}

\usage{
popt(data, pmax)
}

\arguments{
  \item{data}{The input/output dataset.}
  \item{pmax}{The Maximum p one wishes to consider}
}
\details{
Determines if the local density of points is sufficient for a Gamma test, consult the three graphs shown in Figure 3.8 on Page 65 of Kemp's thesis. 
}

\value{
  \item{gammaMSE}{The returned Gamma statistics with the MSE of the associated Gamma regression line.}
}

\references{
For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}}}

\examples{
# Example on an AR(1) process
tssim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
dvtssim <- dvec(tssim, 1)
ts.popt	<- popt(dvtssim, pmax=100)

}
\keyword{ts}
