\name{gammatest}
\alias{gammatest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating the Gamma Statistic}
\description{
  Estimates the variance of the noise in an input-output dataset, provided the underlying function
  is smooth (i.e. bound by partial derivatives).
}
\usage{
gammatest(data, mask=seq(from=1, to=1, length=(length(data[1,])-1)), p=10, eps=0.00, plot=TRUE,
summary=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An input-output dataset, where the output is in the last column.}
  \item{mask}{A vector of 1's and 0's representing input inclusion/exclusion. The default mask is all
  1's (i.e. include all inputs in the test).}
  \item{p}{The maximum number of near neighbours to calculate the Gamma statistic. The default value is set
  to 10.}
  \item{eps}{The error bound of calculating the near neighbour distances. The defualt is 0 meaning the
  exact near neighbours are calculated.}
  \item{plot}{A logical indicating whether a plot is drawn.}
  \item{summary}{Logical indicating whether to output the returned Gamma test statistics. Default
  is set to \code{TRUE}.}
  \item{...}{Graphical options for plotting the Gamma regression line.}
}
\details{
 The Gamma test is an algorithm capable of estimating the lowest
 attainable mean squared error (or noise variance) in an input/output
 dataset. The algorithm does not assume anything regarding the
 parametric equations governing the system, to this end the only
requirement is that the underlying model is smooth.

The algorithm itself works by calculating the nearest neighbour distances in input space. This calculation
is achieved in O(M log M) time, where M is the number of data points using Bentley's kd-tree. The
\code{GammaTest} package utilizes the Approximate Near Neighbor (ANN) C++ library, which can give the exact near
neighbours or (as the name suggests) approximate near neighbours to within a specified error bound. For more
information on the ANN library please visit \url{http://www.cs.umd.edu/~mount/ANN/}.
}
\value{
  \item{Mask}{The mask used to calculate the Gamma statistic}
  \item{deltas.gammas}{The values of the deltas and gammas that are 
  regressed (least squares) on in order to calculate the Gamma statistic.}
  \item{Gamma}{The Gamma statistic}
  \item{Gradient}{The Gradient of the regression line}
  \item{Vratio}{The Gamma statistic divided by the variance of the output.}
}
\references{
Stefansson A., Koncar N. and Jones A. J. (1997), A note on the
gamma test. Neural Computing Applications, 5:131-133.

Evans D., Jones A. J. and Schmidt W. (2002). Asymptotic moments of near neighbour distance 
distributions. Proc. Roy. Soc. Lond. Series A,
458(2028):2839-2849.

Evans D. and Jones A. J. (2002), A proof of the gamma
test. Proc. Roy. Soc. Lond. Series A, 458(2027):2759-2799.

Bentley J. L. (1975), Multidimensional binary search trees used for associative
search. Communication ACM, 18:309-517.

Arya S. and Mount D. M. (1993),
Approximate nearest neighbor searching,
Proc. 4th Ann. ACM-SIAM Symposium on Discrete Algorithms (SODA'93), 271-280.

Arya S., Mount D. M., Netanyahu N. S., Silverman R. and Wu A. Y (1998),
An optimal algorithm for approximate nearest neighbor searching,
Journal of the ACM, 45, 891-923.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{fesearch}} \code{\link{Mtest}} \code{\link{iesearch}} \code{\link{dvec}}}
\examples{
# A noisy sine wave example
x <- seq(length=500, from=-2*pi, to=2*pi)
y <- sin(x) + rnorm(500, sd=sqrt(0.075)) # Set the variance of the noise to 0.075.
xy <- data.frame(x,y)			 # Create an input/output dataset.
nest <- gammatest(xy)		 # Run Gamma test.

# Zero noise time series data: The Henon Map
data(HenonMap)
dvhm <- dvec(HenonMap[,1], 2)
a <- gammatest(dvhm)
}
\keyword{ts}
