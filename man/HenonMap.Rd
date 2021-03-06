\name{HenonMap}
\docType{data}
\alias{HenonMap}

\title{The Henon Map: A Zero Noise Chaotic Time Series}

\description{
This is a zero noise chaotic time series dataset. Generated by 
\deqn{X_t = -X^2_{t-1} + bX_{t-2} + a}{X_t = -X^2_{t-1} + bX_{t-2} + a}
where \deqn{X_0 = 0}, \deqn{X_1=0}, deqn{a=1.4} and deqn{b=0.3}. 
}

\usage{
data(HenonMap)
}

\format{A \code{data.frame} containing 1000 observations}

\references{
Henon M. (1976), A two-dimensional mapping with a strange attractor, Commun. Math. Phys., 50, 69-77.
}

\keyword{datasets}
