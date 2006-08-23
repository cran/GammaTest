#=================================================================
#=================================================================
# res.R
#
# This file provides the function for a Random Embedding search
# Written by Samuel E. Kemp, February 2006.
#=================================================================
#=================================================================


"research" <- function(data, percentage=10, plot=TRUE, ...)
{
	create.binary <- function(num, num.of.bits)
	{
    .C("createBinary", as.integer(num), as.integer(num.of.bits), 
       	mask=integer(num.of.bits), PACKAGE="GammaTest")$mask
	}

	q           <- (length(data[1,]) - 1)
	pop         <- c(1:(2^q -1))
	sampleSize  <- as.integer(percentage*((2^q - 1)/100))
	samp        <- sample(x=pop, size=sampleSize)   
    gtarray     <- NULL
    temp        <- array(dim=c(sampleSize, q))
    
    
    start.time <- Sys.time()
    for(i in 1 : sampleSize)
    {
		m <- create.binary(samp[i], q)
		temp[i,] <- m 
		gtarray[i] <- gammatest(data, mask=m, plot=FALSE, summary=FALSE, ...)$Gamma
	}
    
    y <- sort(abs(gtarray), index.return=TRUE)
    r <- y$ix
    bins <- array(dim=c(sampleSize, q))
    for(i in 1 : sampleSize)
    {
    	bins[i,] <- temp[r[i],]
    }
    
    if(plot==TRUE)
    	gammahist(list(Gammas=y, mask.array=bins))
    
    finish.time <- Sys.time()
    	
    cat("\n")
	cat("================================================\n")
	cat("            Random Embedding Search\n")
	cat("================================================\n")
	cat("Number of data points:     ", length(data[,1]), "\n")
	cat("Number of inputs:          ", q, "\n")
	cat("Search time:               ", finish.time - start.time, 
		attr(finish.time - start.time,"units"),"\n")
	cat("Population size:           ", (2^q - 1), "\n")
	cat("Sample size:               ", sampleSize)
	cat("\n")
 
 	return(list(Gammas=y, mask.array=bins))
}