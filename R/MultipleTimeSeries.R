mdvec <- function(data, lag) # lags is a list of vectors
{
	num.of.ts 	<- (length(data[1,]) - 1)
	d <- NULL
	max.lag <- max(lag)
	
	
	for(i in 1 : num.of.ts)
	{
		Mask	<- seq(from=1, to=1, length=lag)
		ts.io	<- mask2input(Mask, as.double(data[,i]), multiple=TRUE)
		data.length <- length(ts.io[,1] - max.lag)
		temp	<- data.frame(ts.io[1:data.length,])
		
		if(i == 1)
			d	<- temp
		else
			d		<- cbind(d, temp)
	}
	
	# Generate some useful names i.e. inputname.lag1, inputname.lag2
	ts.names <- names(data)
	new.names <- NULL
	inc <- 1
	for(i in 1 : num.of.ts)
	{
		for(j in 1 : lag)
		{
			new.names[inc] <- paste(ts.names[i], paste(".lag", j, sep=""), sep="")
			inc <- inc +1
		}
	}
	new.names[inc] <- "output"
	
	output	<- data.frame(output=data[(lag +1):length(data[,1]),num.of.ts + 1])
	d <- cbind(d, output)
	names(d) <- new.names
	
	return(d)
}
