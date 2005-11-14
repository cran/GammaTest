\name{fesearch}
\alias{fesearch}

\title{The Full Embedding Search}

\description{
Calculates the Gamma statistic for each possible subset of inputs and plots the results in the form
of a Gamma histogram.
}

\usage{
fesearch(data, plot=TRUE, ...)
}

\arguments{
  \item{data}{Either an input-output dataset, where the outputs are in the last column.}
  \item{plot}{Logical indicating whether the Gamma histogram is plotted. Defaults to TRUE.}
  \item{...}{Additional \code{Gamma.Test} parameters. Use if you wish to change the default values 
  of \code{p} and \code{eps} when computing the Gamma Tests.}
}

\details{
 Given an input dimension d, there are \deqn{2^d - 1} possible input subsets. A full embedding search calculates 
 the Gamma statistic for each subset. The input subsets are then ordered in increasing Gamma value. The 
 full embedding search is used to identify relevent variables for a smooth model.
}

\value{
  \item{Gammas}{An array of Gamma statistics ordered in increasing value.}
  \item{mask.array}{An array of ordered masks by Gamma value.}
  \item{input.names}{Character array specify the input names.}
}

\references{
Durrant P. J (2002), winGamma: a non-linear data analysis and modelling tool with
applications to flood prediction, PhD Thesis, Department of Computer Science, Cardiff University, Wales,
UK.

Jones A. J (2003), New tools in non-linear modelling and prediction.
Computational Management Science, 1(1):xx.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{gammahist}} \code{\link{gteam}}}

\examples{
# Example on an AR(1) process
ar1 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
dvar1 <- dvec(ar1, 5)
my.fe	<- fesearch(dvar1)
}
\keyword{ts}
