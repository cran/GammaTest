\name{FCM}
\alias{FCM}

\title{FCM: Frequency Combination Method}

\description{
Identifies the relevent inputs to a smooth model based on the results of an embedding search.
}

\usage{
FCM(mask.array, percentage=12.5)
}

\arguments{
  \item{mask.array}{The returned masks from an embedding search.}
  \item{percentage}{The proportion of masks from the high and low Gamma regions required to conduct
  a frequency analysis on the input inclusion/exclusion. The default is 12.5 percent.}
}
\details{
	FCM is a simple function which takes the mask array from the full embedding and given the proportion
 	of low and high Gammas counts the number of times each input was included in the proportion of low
 	Gammas AND the number of times each input was excluded in the proportion of high Gammas. The frequency 
 	of each input is normalised, so that it lies in the range [0,1].
 	
 	Non-parametric box-plot confidence limits are also pplied, thus the light blue dashed line represents the inner-fence (mild outliers). In the context of model identification the general rule of thumb is: an input should be included if it is greater than the inner fence.
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
FCM(my.fe$mask.array)	
}
\keyword{ts}

