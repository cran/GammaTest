\name{Mtest}
\alias{Mtest}

\title{The M-test: Judging the Reliability of the Gamma Statistic}

\description{
Computes and plots successive Gamma statistics and V ratios over an increasing number of data points \emph{M}. Now includes heuristic confidence limits for the Gamma statistic.}

\usage{
Mtest(data, start=20, step=10, cl=.90, p=10, eps=0.00, ...)
}

\arguments{
  \item{data}{Input/Output \code{data.frame}, where the outputs are in the last column.}
  \item{start}{Burn-in period, must be greater than the number of near neighbours \code{p}. Defaults as 20.}
  \item{step}{The number of data points to increment. Default is 10.}
  \item{cl}{The confidence level for the heuristic confidence intervals. Defaults to .90 i.e. 90\% level.}
  \item{p}{The number of near neighbours to use in the Gamma test analysis.}
  \item{eps}{The error bound for calculating approximate near neighbours. The default is 0 meaning the exact near neighbours are calculated.}
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
{Mlist}{A data.frame containing M and its associated Gamma statistic.}
}

\references{
Stefansson A., Koncar N. and Jones A. J. (1997), A note on the
gamma test. Neural Computing Applications, 5:131-133.

A. J. Jones and S. E. Kemp (2006), Heuristic confidence limits for the Gamma test. Submitted paper available to download from url below.

For papers and Gamma test related material visit 
\url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}
}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{mask2input}}}

\examples{
# Noisy Sine wave example
x <- runif(500, 0, 2*pi)
y <- sin(x) + rnorm(500, sd=sqrt(0.075)) 
xy <- data.frame(x,y)
plot(xy)
mt <- Mtest(xy)
}
\keyword{ts}
