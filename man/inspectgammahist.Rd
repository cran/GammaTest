	\name{inspectgammahist}
\alias{inspectgammahist}

\title{Model Identification: Inspecting The Peaks of a Gamma Histogram}

\description{
Identifies the relevent inputs to a smooth model based on the analysing the peaks in the Gamma histogram
from a full embedding search.
}

\usage{
inspectgammahist(fe.results, gamma.range, ones=TRUE)
}

\arguments{
  \item{fe.results}{The returned results from a full embedding search.}
  \item{gamma.range}{The range of the Gamma statistics we wish to inspect.}
  \item{ones}{Logical indicating whether to count the occurances of input inclusion/exclusion. \code{TRUE}
  indicates to count input exclusion, \code{FALSE} means to count input exclusion.}
}
\details{
Counts the number of input inclusions or exclusions within a particular Gamma range. For more information
see Durrant (2002).
}

\value{
  \item{results}{The normalised frequency counts.}
}

\references{
Kemp S. E. (2004), Gamma Test Based Data Analysis Applied to
Modelling and Forecasting Crime Rates, MPhil to PhD Transfer Report, School of Computing,
University of Glamorgan, Wales, UK.

Durrant P. J (2002), winGamma: a non-linear data analysis and modelling tool with
applications to flood prediction, PhD Thesis, Department of Computer Science, Cardiff University, Wales,
UK.

Jones A. J (2003), New tools in non-linear modelling and prediction.
Computational Management Science, 1(1):xx.

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
inspectgammahist(my.fe, gamma.range=c(0,1), ones=TRUE)	
}
\keyword{ts}

