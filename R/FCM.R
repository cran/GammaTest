# FCM.R
# File written by: Samuel E. Kemp 2006
# This file contains an R implementation of the Frequency Combination method
#Êas outlined in Kemp (2006).

FCM <- function(mask.array, percentage=12.5)
{
   	if(!is.matrix(mask.array))
   		mask.array <- as.matrix(mask.array)
   		
   	num.of.inputs   <- length(mask.array[1,])
   	num.of.masks    <- length(mask.array[,1])
   	sample          <- as.integer((num.of.masks/100)*percentage)
	counter         = double(num.of.inputs)
	HGRptr          <- num.of.masks
    
    for(i in 1 : sample)
    {
       	for(j in 1 : num.of.inputs)
       	{
           	if(mask.array[i, j] == 1 && mask.array[HGRptr, j] == 0)
               	counter[j] <- counter[j]+1
       	}
       	HGRptr <- HGRptr-1
    }
    
   	results          <- counter/sample
  	
  	# Get the upper limit, which in our case is the upper quartile as calculated
  	# in the box-whisker plots.
  	upper.limit <- boxplot.stats(results)$stats[4] 

    plot(results, ylim=c(0,1), xlab="Lag", ylab="Frequency",type="h")
       abline(h=upper.limit,lty="dashed", col="blue")
           
    return(data.frame(frequency=results))
}
