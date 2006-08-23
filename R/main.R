#===========================================================
# GammaTest package file
#===========================================================

#===========================================================
gammatest <- function(data, mask=seq(from=1, to=1, length=(length(data[1,])-1)), p=10, eps=0.00, plot=TRUE, summary=TRUE, ...)
#===========================================================
{
    # Coerce to a data.frame
    if(is.data.frame(data) == FALSE)
    	data <- data.frame(data)
	
    # Check that this is an input/output dataset
    if(length(data[1,]) <= 1)
	stop("Please make this an input/output dataset.")	
		
	
    if(sum(mask) == 0) # Check that we are going to have some inputs!!
	stop("Invalid mask: You must have at least one input!")
	
    if(length(mask) != (length(data[1,])-1))	# Check the mask length.
	stop("Invalid mask length")
	
	
    for(i in 1 : length(mask)) # Check that the mask only contains 1s and 0s
    {
	if(mask[i] > 1 | mask[i] < 0)
	    stop("Mask vector can only contain zeros or ones, refer to help pages")
    }
		
    # Check that p is not too large.
    if(p > (length(data[,1]) - 1))			
    	stop("There is not enough data to calculate the p near neighbours you are 
              selecting.")
	
    num.inputs 	<- sum(mask)
    dimension	<- length(data[1,])
    M			<- length(data[,1])
	
    # Run the C++ code to get the deltas and gammas
    results <- .C("Gamma_Test_Main",
		as.matrix(data),
		as.integer(mask),
		as.integer(num.inputs),
		as.integer(p),
		as.integer(dimension),
		as.integer(M),
		as.double(eps),
		d = double(p),
		g = double(p), PACKAGE="GammaTest")
		
    dg 	<- data.frame(delta=results$d, gamma=results$g)
    z 	<- lm(dg$gamma ~ dg$delta, dg)
    gs 	<- z$coefficients[[1]]
    grad <- z$coefficients[[2]]
    z.sum <- summary(z)
    vr <- gs/var(data[,dimension])
    mse <- var(z$residuals)
    
    if(plot == TRUE)
    {
		plot(dg$delta, dg$gamma, xlab=expression(delta[M](k)), ylab=expression(gamma[M](k)), ...)
		abline(z, lwd=2)
    }
    
    gtlist <- list(mask=mask, deltas.gammas=dg, Gamma=gs, Gradient=grad, Vratio=vr, MSE=mse)
	
    if(summary == TRUE)
		gtsummary(gtlist)
	
    return(gtlist)
}

"popt" <- function(data, pmax)
{
	criterion <- NULL
	gammas    <- NULL
	for(i in 3 : pmax)
	{
		temp <- gammatest(data, p=i, summary=FALSE, plot=FALSE)
		criterion[i-2] <- temp$MSE
		gammas[i-2]    <- temp$Gamma
	}
	
	return(gammaMSE<-data.frame(Gammas=gammas, MSE=criterion))
}


#===========================================================
gtsummary <- function(gt.list)
#===========================================================
{
    cat("==========================================\n")
    cat("     Gamma Test: Summary of Results\n")
    cat("==========================================\n")
    cat("\n")
    cat("Mask:                  ", gt.list$mask, "\n")
    cat("Gamma Statistic:       ", gt.list$Gamma, "\n")
    cat("Gradient:              ", gt.list$Gradient, "\n")
    cat("V Ratio:               ", gt.list$Vratio, "\n")
    cat("MS Error of Reg. fit:  ", gt.list$MSE, "\n")
    cat("\n")
}
#===========================================================
fesearch <- function(data, plot=TRUE, annerror=0.0, ...)
#===========================================================
{
	#===========================================================
	create.binary <- function(num, num.of.bits)
	#===========================================================
	{
    	.C("createBinary", as.integer(num), 
       	as.integer(num.of.bits), 
       	mask=integer(num.of.bits),
       	, PACKAGE="GammaTest")$mask
	}

    maxlag <- length(data[1,]) - 1  # Don't include output!!
    n <- names(data)[1:maxlag]
    
    data <- as.matrix(data)
	
    if(length(data[1,]) == 1)
    {
    	stop("Please convert data to delay vectors")
    }
	
    # Get start time
    start.time <- Sys.time()

    max.num     <- 2^(maxlag)-1
    gtarray     <- double(max.num)
    
    for(i in 1 : max.num)
    {
	m <- create.binary(i, maxlag)
	gtarray[i] <- gammatest(data, mask=m, plot=FALSE, summary=FALSE, eps=annerror,...)$Gamma       }
    
    y <- sort(abs(gtarray), index.return=TRUE)
    r <- y$ix
    bins <- array(dim=c(max.num, maxlag))
    for(i in 1:max.num){bins[i,] <- create.binary(r[i], maxlag)}
    
    if(plot==TRUE)
    	gammahist(list(Gammas=y, mask.array=bins), col="red")
    
    # Get finish time
    finish.time <- Sys.time()
    time.per.test <- (finish.time - start.time)/max.num
    
    if(attr(finish.time - start.time,"units") == "secs")
    	per.sec <- max.num /as.double(finish.time - start.time)
    
    if(attr(finish.time - start.time,"units") == "mins")
    	per.sec <- (max.num /as.double(finish.time - start.time))/60

	
    cat("\n")
    cat("===================================================\n")
    cat("            Full Embedding Statistics\n")
    cat("===================================================\n")
    cat("Full embedding search time:      ", finish.time - start.time, attr(finish.time - start.time,"units"),
	"\n")
	cat("Number of Gamma tests:     ", length(bins[,1]))
	cat("\n")
    
    return(list(Gammas=y, mask.array=bins, input.names=n, time=finish.time - start.time))
}


#===========================================================
gammahist <- function(fe.results, ...)
#===========================================================
{
	hist(fe.results$Gammas$x, main="", 
	xlab=paste(expression(Gamma), "Bins"), 
	cex.lab=1.15, cex.axis=1.25, breaks=50, ...)
}

#===========================================================
durrantsmethod <- function(mask.array, percentage=10)
#===========================================================
{
    	num.of.inputs       <- length(mask.array[1,])
    	num.of.masks        <- length(mask.array[,1])
    	sample              <- as.integer((num.of.masks/100)*percentage)
    	LGRcounter          = double(num.of.inputs)
    	HGRcounter          = double(num.of.inputs)
    	HGRptr              <- num.of.masks
        
    	for(i in 1 : sample)
    	{
        	for(j in 1 : num.of.inputs)
        	{
                LGRcounter[j]   <- LGRcounter[j] + mask.array[i,j]
                if(mask.array[HGRptr, j] == 0)
                    HGRcounter[j]   <- HGRcounter[j] +1
        	}
        	HGRptr <- HGRptr-1
    	}
    
    	result           <- LGRcounter/sample
    	results          <- HGRcounter/sample
      	LGRresults       <- result
    	HGRresults   	 <- results
  
    	
    	return(data.frame(lgr=LGRresults, hgr=HGRresults))
}    


#===========================================================  
mask2input <- function(mask, timeseries, multiple=FALSE)
#===========================================================
{
	 createIO <- function(data,lag)
    {
		newData <- array(dim=c(length(data)-lag, lag+1))	
		inc <- 1
	
		for(i in 1: (length(data) - lag))
		{
		    jj <- lag
	    	for(j in 1: (lag+1))
	    	{
	    		if(j == (lag+1))
				    newData[i,(lag+1)] <- data[inc]
			else
		    	newData[i,jj] <- data[inc]
			
			inc <- inc+1
			jj <- jj-1
	    	}
	    	inc <- inc-lag
		}
	
		return(data.frame(newData))
    }

    timeseries <- as.matrix(timeseries)
    if(length(mask) != (dim(timeseries)[2] - 1))
    {
        timeseries <- createIO(as.matrix(timeseries), length(mask))
    }
        
    arraySize <- 0
    # Count the size of the returned input/output set
    for(i in 1 : length(mask))
    {
        if(mask[i] == 1)
        {
            arraySize <- arraySize + 1
        }
    }
    
    # set up the array dimensions
    if(multiple == FALSE)
    {
    	newIO <- array(dim=c(dim(timeseries)[1], arraySize+1))
    	newIO[,arraySize+1] <- timeseries[,dim(timeseries)[2]]
    }
    
    if(multiple == TRUE)
    	newIO <- array(dim=c(dim(timeseries)[1], arraySize)) 
    
    index <- 1  
    for(i in 1 : length(mask))
    {
        if(mask[i] == 1)
        {
            newIO[,index] <- timeseries[,i]
            index <- index +1
        }
    }
    
    n <- NULL
    d <- length(mask)
    ii <- 1
    for(i in 1 : d)
    {
	if(mask[i] == 1)
	{
	    n[ii] <- paste("lag.", i, sep="")
	    ii <- ii + 1
	}
    }
    n[length(newIO[1,])] <- "output"
	
    io <- data.frame(newIO)
    names(io) <- as.character(n)
    
    return(data.frame(io))
}

#===========================================================
dvec <- function(time.series, lag)
#===========================================================
{
	# Check that the time series is a double, if not coerce
	if(is.double(time.series) == FALSE)
		time.series <- as.double(time.series)
		
    createIO <- function(data,lag)
    {
		newData <- array(dim=c(length(data)-lag, lag+1))	
		inc <- 1
	
		for(i in 1: (length(data) - lag))
		{
		    jj <- lag
	    	for(j in 1: (lag+1))
	    	{
	    		if(j == (lag+1))
				    newData[i,(lag+1)] <- data[inc]
			else
		    	newData[i,jj] <- data[inc]
			
			inc <- inc+1
			jj <- jj-1
	    	}
	    	inc <- inc-lag
		}
	
		return(data.frame(newData))
    }
    
    time.series <- as.matrix(time.series)
    time.series <- createIO(time.series, lag)
	
    n <- NULL
    d <- lag
    for(i in 1 : d)
    {
	n[i] <- paste("lag.", i, sep="")
    }
	
    n[lag + 1] <- "output"
	    
    io <- data.frame(time.series)
    names(io) <- as.character(n)
	    
    return(data.frame(io))
}

#===========================================================
iesearch <- function(data, ...)
#===========================================================
{
    d <- length(data[1,]) - 1
    mask <- double(d)
    gtarray <- double(d)
    for(i in 1 : d)
    {
    	mask[i] <- 1
        io <- mask2input(mask, as.matrix(data))
        res <- gammatest(io, plot=FALSE, summary=FALSE, ...)$Gamma
        gtarray[i] <- abs(res)
    }
    
    plot(gtarray, type="l", xlab="Lag", ylab=expression(Gamma), 
    lwd=2, cex.lab=1.35, cex.axis=1.15)
    
    return(data.frame(Gammas=gtarray))
}
