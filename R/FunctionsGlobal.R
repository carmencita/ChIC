####GLOBAL features ###

f_shortenFrame=function(smoothed.density)
{
	##shorten frame
	newSmoothedDensity=NULL
	chrl=names(smoothed.density)
	for (i in chrl)
	{
		#print(i)
		##appending tags and quality elements of all inputs to newControl
		newSmoothedDensity[[i]]$x=smoothed.density[[i]]$x[seq.int(1, length(smoothed.density[[i]]$x), 6)]
		newSmoothedDensity[[i]]$y=smoothed.density[[i]]$y[seq.int(1, length(smoothed.density[[i]]$y), 6)]
	}
	return(newSmoothedDensity)
}


f_sortAndBinning=function(shortframe)
{
	shortframeSorted=sort(shortframe)
	BINS_NUMBER<-1e4

	cumSum<-cumsum(shortframeSorted)
	cumSumBins<-quantile(cumSum, probs=seq(0,1,(1/BINS_NUMBER)))
	pj<-(cumSumBins/cumSum[length(cumSum)])
	normalizedCumSum=data.frame(x=seq(0,1,(1/BINS_NUMBER)),pj=pj)
	rownames(normalizedCumSum)<-NULL
	return(normalizedCumSum)
}


f_chancePlots=function(cumChip,cumInput,plotname="chancePlot.pdf")
{
	
	pdf(plotname)
	plot(cumChip,type="l",col="blue",lwd=2,xlab="Percentage of bins",ylab="Percentage of tags",main="Fingerprint: global read distribution")
	lines(cumInput,col="red",lwd=2)
	arrowx=cumChip[which.max(abs(cumInput$pj-cumChip$pj)),]$x
	abline(v =arrowx, col = "green",lty=2,lwd=2)
	#abline(h=schneidePunktY,col='cyan',lty=2,lwd=2)
	#abline(v=schneidePunktX,col='cyan',lty=2,lwd=2)
	legend("topleft", legend = c("Input","ChIP"), fill = c("red","blue"))
	dev.off()
}
