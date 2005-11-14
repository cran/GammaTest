inspectgammahist <- function(fe.results, gamma.range, ones=TRUE)
{
    num.of.inputs       <- length(fe.results$mask.array[1,])
    num.of.masks        <- length(fe.results$mask.array[,1])
    
    counter		        <- double(num.of.inputs)
    sample	            <- 0
    bit					<- NULL
    
    if(ones == TRUE)
    	bit <- 1
    else
    	bit <- 0
        
    for(i in 1 : num.of.masks)
    {
    	if((fe.results$Gammas$x[i] >= gamma.range[1]) && (fe.results$Gammas$x[i] <= gamma.range[2]))
        {   
        	for(j in 1 : num.of.inputs)
        	{       	
        		if(fe.results$mask.array[i, j] == bit)
            		counter[j]   <- counter[j] +1  	
            }
            sample <- sample + 1
        }        
    }
    
    result           <- counter / sample
    
    plot(result, xlab="Lag", ylab="Frequency", ylim=c(0,1), type="h", lwd=2, cex.axis=1.2, cex.lab=1.15)
       
    return(result)
    
}
