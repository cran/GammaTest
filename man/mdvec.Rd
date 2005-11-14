\name{mdvec}
\alias{mdvec}

\title{Converting Multiple Time Series into Input/Output Datasets}

\description{
Converts multiple time series into an input/output format, where the output is in the last column and the
input time series are in the other columns.
This is required to run a gamma test analysis on multiple time series.
}

\usage{
mdvec(data, lag)
}

\arguments{
  \item{data}{A data.frame where the last column is the output and the other columns are input time
  series. It is recommended that the names of the variables in the data.frame are meaningful. For
  example, if one of the input time series is the FTSE then use FTSE and not X as a variable names.}
  \item{lag}{The maximum lag with which to create the inputs FOR EACH TIME SERIES.}
}

\value{
  \item{d}{The input/output dataset.}
}

\references{For papers and Gamma test related material visit \url{http://users.cs.cf.ac.uk:81/Antonia.J.Jones/GammaArchive/IndexPage.htm}}

\author{Samuel E. Kemp. To report any bugs or suggestions please email: \email{sekemp@glam.ac.uk}}

\seealso{\code{\link{gammatest}} \code{\link{fesearch}}}

\examples{
# One input time series example
X <- arima.sim(500, model=list(ar=0.7, ma=0.6), sd=sqrt(1))
Y <- NULL
Noise <- rnorm(500, sd=sqrt(1))
for(i in 2 : 500){
Y[i] <- 0.6*X[i-1] + Noise[i]}
input.output <- data.frame(X=X[2:500],Y=Y[2:500])
multi.io     <- mdvec(input.output, lag=10)

# Three input time series example
X <- arima.sim(500, model=list(ar=0.7, ma=0.6), sd=sqrt(1))
R <- arima.sim(500, model=list(ar=0.7), sd=sqrt(1))
Q <- arima.sim(500, model=list(ma=0.6), sd=sqrt(1))
Y <- NULL
Noise <- rnorm(500, sd=sqrt(1))
for(i in 3 : 500){
Y[i] <- (0.6*X[i-1] + 0.7*R[i-2] + 0.8*Q[i-1]) + Noise[i]}
input.output <- data.frame(X=X[3:495],R=R[3:495],Q=Q[3:495],Y=Y[3:495])
multi.io     <- mdvec(input.output, lag=3)
}
}
\keyword{ts}
