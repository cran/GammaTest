\name{durrantsmethod}
\alias{durrantsmethod}

\title{Durrant's Feature Selection Method}

\description{
Durrant's Gamma test based feature selection method can be used to identify the relevent inputs to a smooth model based on the results of a full embedding search.
}

\usage{
durrantsmethod(fe.results, percentage=10)
}

\arguments{
  \item{fe.results}{The returned results from a full embedding search.}
  \item{percentage}{The proportion of masks from the high and low Gamma regions required to conduct
  a frequency analysis on the input inclusion/exclusion. The default is 12.5 percent.}
  }
\details{
 	Consult Durrant's thesis from the web link below for more details. }

\value{
  \item{results}{The normalised frequency counts.}
}

\references{

Durrant P. J (2002), winGamma: a non-linear data analysis and modelling tool with
applications to flood prediction, PhD Thesis, Department of Computer Science, Cardiff University, Wales,
UK.

For papers, theses and other Gamma test related material please visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{fesearch}} \code{\link{dvec}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
gfts.sim <- dvec(ts.sim, 8)
my.fe	<- fesearch(gfts.sim)
durrantsmethod(my.fe)	
}
\keyword{ts}

