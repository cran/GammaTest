loom <- function(data, ...)
{
	if(length(data[1,]) <= 1)
		stop("Please create an input/output dataset")
	
	maxlag		<- length(data[1,]) - 1
    mask		<- seq(from=1, to=1, length=maxlag)    
  	results 	<- NULL
  	
  	ptr <- maxlag
    for(i in 1 : maxlag)
    {
    	if(i == maxlag)
    		mask[maxlag] <- 1
    		
    	mask[ptr]		<- 0
        io          <- mask2input(mask, data)
        jj			<- gammatest(io, summary=FALSE, plot=FALSE, ...)$Gamma
        results[ptr]<- abs(jj)
        ptr 		<- ptr-1
    }
        
   	plot(results, type="l", xlab="Lag", ylab=expression(Gamma), lwd=2, cex.lab=1.35, cex.axis=1.15,
   	xlim=c(maxlag, 1))
    return(data.frame(results))
}
