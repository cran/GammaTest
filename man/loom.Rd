\name{loom}
\alias{loom}

\title{Calculating the Embedding Dimension: The Leave One Out Method}

\description{
Calculates and plots the Gamma statistic, where successive inputs are left out of the embedding. For
time series inputs are left out in reverse time order. For example, a 5 input mask iterates as follows:
1 1 1 1 0, 1 1 1 0 0 , 1 1 0 0 0, 1 0 0 0 0, 0 0 0 0 1. Notice for the last iteration where we exclude
lag one we include the last lag because otherwise we have no data on which to calculate the effect of leaving
out the lag!
}

\usage{
loom(data,...)
}

\arguments{
  \item{data}{An input/output dataset.}
  \item{...}{Additional \code{gammatest} parameters. Use if you wish to change the default values of 
  \code{p} and \code{eps} when computing the Gamma tests.}
}
\details{
 The leave one out technique is a similar procedure to that of the increasing embedding search.
}

\value{
  \item{results}{The returned Gamma statistics.}
}

\references{
Otani M and Jones A. J (1997). Guiding chaotic orbits.
Technical report, Department of Computing, Imperial College of Science, Technology and Medicine, 
and Department of Computer Science, Cardiff University.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{iesearch}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
gfts.sim <- dvec(ts.sim, 15)
my.loom	<- loom(gfts.sim)

# Henon Map example
data(HenonMap)
gfHenonMap <- dvec(HenonMap[,1], 20)
hm.loom <- loom(gfHenonMap)
}
\keyword{ts}
