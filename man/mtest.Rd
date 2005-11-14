\name{mtest}
\alias{mtest}

\title{The M-test: Judging the Reliability of the Gamma Statistic}

\description{
Computes and plots successive Gamma statistics and V ratios over an increasing number of 
data points \emph{M}.
}

\usage{
mtest(data, start=20, mask=seq(from=1, to=1, length=(length(data[1,])-1)), step=1, p=10, eps=0.00, ...)
}

\arguments{
  \item{data}{Input/Output \code{data.frame}, where the outputs are in the last column.}
  \item{start}{Burn-in period, must be greater than the number of near neighbours \code{p}.}
  \item{mask}{The input inclusion/exclusion on which to calculate the Gamma statistics. The default
  is all the 1s i.e. include all inputs.}
  \item{step}{The number of data points to increment.}
  \item{p}{The sample of near neighbours}
  \item{eps}{The error bound for calculating approximate near neighbours. The default is 0 meaning
  the exact near neighbours are calculated.}
  \item{...}{Graphical device parameters (\emph{e.g.} colour, labelling, titles, etc)}
}

\details{
The proof of the Gamma test assures us that as
\emph{M} tends to infinity the noise estimate given by the Gamma test will converge on the
true noise variance with probability one. However,
this does not give an indication of how many data points are
required to give an accurate estimate of the noise. In order
to achieve this, the Gamma statistic is computed for increasing \emph{M}. By
plotting the Gamma statistics over \emph{M} it can be seen whether the graph
appears to be approaching a stable asymptote.

If \emph{M}-test graph stabilises then this indicates that there is enough data to
construct a smooth model, which is of considerable benefit neural network users because the number of
data points required for the training set can be determined prior to modelling.
}

\value{
{results}{The returned Gamma statistics}
{vratios}{The returned V ratios}
}

\references{
Stefansson A., Koncar N. and Jones A. J. (1997), A note on the
gamma test. Neural Computing Applications, 5:131-133.

Evans D., Jones A. J. and Schmidt W. (2002). Asymptotic moments of near neighbour distance 
distributions. Proc. Roy. Soc. Lond. Series A,
458(2028):2839-2849.

Evans D. and Jones A. J. (2002), A proof of the gamma
test. Proc. Roy. Soc. Lond. Series A, 458(2027):2759-2799.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{mask2input}}}

\examples{
# Example on an AR(1) process
ts.sim 	<- arima.sim(500, model=list(ar=0.9), sd=sqrt(1))
ts.io	<- mask2input(1, as.matrix(ts.sim)) # Create input/output dataset
mt		<- mtest(ts.io)

# Noisy Sine wave example
x <- seq(length=500, from=-2*pi, to=2*pi)
y <- sin(x) + rnorm(500, sd=sqrt(0.075)) # Set the variance of the noise to 0.075.
xy <- data.frame(x,y)			 # Create an input/output dataset.
mt <- mtest(xy)
}
\keyword{ts}
