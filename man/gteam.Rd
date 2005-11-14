\name{gteam}
\alias{gteam}

\title{GTEAM: Gamma Test Embedding Analysis Method}

\description{
Identifies the relevent inputs to a smooth model based on the results of a full embedding search.
}

\usage{
gteam(fe.results, percentage=12.5, cl="p", alpha=1.96, multiple=FALSE, ...)
}

\arguments{
  \item{fe.results}{The returned results from a full embedding search.}
  \item{percentage}{The proportion of masks from the high and low Gamma regions required to conduct
  a frequency analysis on the input inclusion/exclusion. The default is 12.5 percent.}
  \item{cl}{Whether to use non-parametric or parametric confidence limits. The default is parametric.To change to parametric use "p".}
  \item{alpha}{Confidence limit.}
  \item{multiple}{logical indicating whether there are multiple time series inputs.} 
  \item{...}{If dealing with multiple time series then one can send the names, so that a bar plot
  produces the correct time series names for each lag.}
}
\details{
In the case of time series model identification, GTEAM (at present) can only be applied to de-trended and de-seasonalised data.

 	GTEAM is a simple function which takes the mask array from the full embedding and given the proportion
 	of low and high Gammas counts the number of times each input was included in the proportion of low
 	Gammas AND the number of times each input was excluded in the proportion of high Gammas. The frequency 
 	of each input is normalised, so that it lies in the range [0,1].
 	
 	This function also plots the normalised frequencies and provides parametric confidence limits.  Non-parametric box-plot confidence limits can also be applied, thus the light blue dashed line represents the inner-fence (mild outliers) and the
 	navy dashed line represents the outer-fence (extreme outliers). In the context of model identification the general rule of thumb is: an input should be included if it is greater than the outer fence and an input may be
 	included if it lies between the inner and outer fences.
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
gteam(my.fe)	
}
\keyword{ts}

