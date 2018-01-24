####################
##PRIVATE functions
####################

f_shortenFrame=function(smoothedDensity)
{
	##shorten frame with cumulative distribution
	newSmoothedDensity=NULL
	chrl=names(smoothedDensity)
	for (i in chrl)
	{
		#print(i)
		##appending tags and quality elements of all inputs to newControl
		newSmoothedDensity[[i]]$x=smoothedDensity[[i]]$x[seq.int(1, length(smoothedDensity[[i]]$x), 6)]
		newSmoothedDensity[[i]]$y=smoothedDensity[[i]]$y[seq.int(1, length(smoothedDensity[[i]]$y), 6)]
	}
	return(newSmoothedDensity)
}

##sorts and bins the dataframe
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

####################
##GLOBAL functions
####################

#'@title Finger Print Plot
#'
#' @description
#' The plot function creates the "Fingerprint plot" i.e. the cumulative distribution of the read counts across genomic bins for Input and ChIP.
#'
#' f_fingerPrintPlot
#'
#' @param cumChip The cumulative distribution of the read counts for the ChIP
#' @param cumInput The cumulative distribution of the read counts for the Input
#' @param plotname Name of the FingerPring plot (default is "chancePlot.pdf" and saved in the working directory )
#'
#' @return returnList
#' @examples
#'\dontrun{
#' print("create pdf with plot")
#' f_fingerPrintPlot(cumSumChip,cumSumInput,plotname=paste(file.path(getwd(),"FingerPrintPlot.pdf",sep="_")
#'}

f_fingerPrintPlot=function(cumChip,cumInput,plotname=file.path(getwd(),"FingerPrintPlot.pdf"))
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
