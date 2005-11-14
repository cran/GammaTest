\name{iesearch}
\alias{iesearch}

\title{Calculating the Embedding Dimension: The Increasing Embedding Search}

\description{
Calculates and plots the Gamma statistic over an increasing embedding of a time series.
}

\usage{
iesearch(data, ...)
}

\arguments{
  \item{data}{The input/output dataset.}
  \item{...}{Additional \code{Gamma.Test} parameters. Use if you wish to change the default values of 
  \code{p} and \code{eps} when computing the Gamma tests.}
}
\details{
 The increasing embedding technique has been used as an alternative to the false near neighbours (FNN)
 method proposed by Kennel. The increasing embedding technique has been successfully applied to determine
 the embedding dimension of many chaotic time series - see references below.
}

\value{
  \item{Gammas}{The returned Gamma statistics.}
}

\references{
Otani M and Jones A. J (1997). Guiding chaotic orbits.
Technical report, Department of Computing, Imperial College of Science, Technology and Medicine, 
and Department of Computer Science, Cardiff University.

Durrant P. J (2002), winGamma: a non-linear data analysis and modelling tool with
applications to flood prediction, PhD Thesis, Department of Computer Science, Cardiff University, Wales,
UK.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}}}

\examples{
# Example on an AR(1) process
tssim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
dvtssim <- dvec(tssim, 15)
my.ie	<- iesearch(dvtssim) # Run inc. embedding on 20 previous lags

# Henon Map example
data(HenonMap)
dvhm <- dvec(HenonMap[,1], 15)
a <- iesearch(dvhm)
}
\keyword{ts}
