\name{research}
\alias{research}

\title{The Random Embedding Search}

\description{
Calculates the Gamma statistic for a random selection of possible subset of inputs and plots the results in the form
of a Gamma histogram.
}

\usage{
research(data, percentage=10, plot=TRUE, ...)
}

\arguments{
  \item{data}{Either an input-output dataset, where the outputs are in the last column.}
  \item{percentage}{The random sample to take as a percentage of the total number of embeddings. Defaults to 10.}
  \item{plot}{Logical indicating whether the Gamma histogram is plotted. Defaults to TRUE.}
  \item{...}{Additional \code{Gamma.Test} parameters. Use if you wish to change the default values 
  of \code{p} and \code{eps} when computing the Gamma Tests.}
}

\details{
 Given an input dimension d, there are \deqn{2^d - 1} possible input subsets. A full embedding search calculates 
 the Gamma statistic for each subset. The input subsets are then ordered in increasing Gamma value. The 
 random embedding search is used to identify relevent variables for a smooth model when the number of candidate inputs is large i.e. around 20.
}

\value{
  \item{Gammas}{An array of Gamma statistics ordered in increasing value.}
  \item{mask.array}{An array of ordered masks by Gamma value.}
}

\references{
For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{gammahist}} \code{\link{FCM}}}

\examples{
# Example on an AR(1) process
ar1 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
dvar1 <- dvec(ar1, 10)
my.fe	<- research(dvar1)
}
\keyword{ts}
