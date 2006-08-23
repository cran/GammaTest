# HeuristicIntervals.R
#==========================================================================
# 				The M-test with Confidence Intervals
#==========================================================================
# 
# We shall assume that the Gamma distribution is approximately normal 
# and the SD scales like c/sqrt(M) as M varies. With this
# assumption we 'normalise' the estimates for different M values (as in
# the M-test) and apply the usual Student t-test approach to estimate the
# confidence intervals. The only barrier is to estimate c in c/sqrt(M).
#
# This file contains all the neccessary functions to compute the M-test
# with heuristic confidence intervals.

# Run Mtest and get M and the associated Gammas


"Mtest" <- function(data, start=20, step=10, cl=.90, p=10, eps=0.00, ...)
{
	# SUPPORT FUNCTIONS
	"getMlist" <- function(data, start, step)
	{
		M <- length(data[,1])
		gammalist <- NULL
		f <- (M-start)/step
		l<- NULL
		len <- start
		for(i in 1 : (f+1))
		{
			l[i] <- len
			newdata <- data[1:len,]
			gammalist[i] <- gammatest(newdata, p=p, eps=eps, summary=FALSE, plot=FALSE)$Gamma
			len <- len + step
		}
	
		return(data.frame(M=l, Gamma=gammalist))
	}

	"GammaCI" <- function(Mlist, cl)
	{
		L         <- length(Mlist[,1])
		sqrtMLs   <- sqrt(Mlist[,1])
		ML        <- Mlist[L,1]
		Gammas    <- Mlist[,2]
		mtmean    <- sum(sqrtMLs*Gammas)/(sqrt(ML)*L)
		temp      <- ((sqrtMLs*Gammas)-(sqrtMLs*mtmean))^2
	
		mtsd      <- sqrt(sum(temp)/(ML*(L-1)))
		stderr    <- mtsd/sqrt(L)
		alpha     <- 1-cl
		
		limits <- Mlist[L,2]+c(-1,1)*stderr*qt(1-(alpha/2), df=(L-1))
		return(limits)
	}

	Mlist <- getMlist(data, start, step)
	L     <- length(Mlist[,1])
	ebars <- array(dim=c((L-2),2)) # store the upper and lower CI
	for(i in 2 : (L-1))
	{
		temp <- GammaCI(Mlist[1:(i+1),], cl)
		ebars[i-1,1]<-temp[1]
		ebars[i-1,2]<-temp[2]
	}
	
	ylimits <- c(min(ebars[,1]), max(ebars[,2])) 
	plot(x=Mlist[3:L,1], Mlist[3:L,2], type="l", ylim=ylimits, ylab=expression(Gamma), xlab="M", ...)
	segments(Mlist[3:L,1], ebars[,1],Mlist[3:L,1], ebars[,2], col="red")
	
	cat("\n")
	cat("Gamma test (Heuristic) confidence intervals:\n")
	cat("	Lower Limit = ", ebars[(L-2),1], "\n")
	cat("	Upper Limit = ", ebars[(L-2),2], "\n")
	cat("\n")
	
	return(Mlist)	
}

