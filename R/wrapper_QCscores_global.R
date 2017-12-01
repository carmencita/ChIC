source(paste(path,"GlobalParameters.R",sep=""))
load(file=file.path(path,"TagDensityChip.RData"))
load(file=file.path(path,"TagDensityInput.RData"))



f_shortenFrame=function(smoothed.density)
{
	print("shorten frame")
	##shorten frame
	newSmoothedDensity=NULL
	chrl=names(smoothed.density)
	for (i in chrl)
	{
		print(i)
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


f_chancePlots=function(cumChip,cumInput,plotname="chance")
	plotname=file.path(getwd(),paste(plotname,".pdf",sep=""))
	pdf(plotname)
	plot(cumChip,type="l",col="blue",lwd=2,xlab="Percentage of bins",ylab="Percentage of tags",main=TFname)
	lines(cumInput,col="red",lwd=2)
	arrowx=cumChip[which.max(abs(cumInput$pj-cumChip$pj)),]$x
	abline(v =arrowx, col = "green",lty=2,lwd=2)
	#abline(h=schneidePunktY,col='cyan',lty=2,lwd=2)
	#abline(v=schneidePunktX,col='cyan',lty=2,lwd=2)
	legend("topleft", legend = c("Input","ChIP"), fill = c("red","blue"))
	dev.off()
}


QCscores_global=function(densityChip,densityInput)
{

	chip.smoothed.density=f_shortenFrame(densityChip)
	input.smoothed.density=f_shortenFrame(densityInput)
	

	#tag.shift <- round(strandShift/2)

	chip<-unlist(lapply(chip.smoothed.density, FUN=function(list_element) {return(list_element$y)}))
	input<-unlist(lapply(input.smoothed.density, FUN=function(list_element) {return(list_element$y)}))


	cumSumChip=f_sortAndBinning(chip)
	chipFractionOfReadsIntop1percentBins<-(1-cumSumChip[(which(cumSumChip$x>=0.99)[1]),"pj"])
	chipFractionOfBinsWithoutReads<-cumSumChip$x[(which(cumSumChip$pj>0)[1])]

	cumSumInput=f_sortAndBinning(input)
	inputFractionOfReadsIntop1percentBins<-(1-cumSumInput[(which(cumSumInput$x>=0.99)[1]),"pj"])
	inputFractionOfBinsWithoutReads<-cumSumInput$x[(which(cumSumInput$pj>0)[1])]


	sign_sign=sign(max(cumSumInput$pj-cumSumChip$pj)) ##input-chip
	#cross_pointX=cumSumInput[which((cumSumInput$pj==cumSumChip$pj)&(cumSumInput$pj!=0)),]$x
	#cross_pointY_input=cumSumInput[which((cumSumInput$pj==cumSumChip$pj)&(cumSumInput$pj!=0)),]$pj
	#cross_pointY_chip=cumSumChip[which((cumSumInput$pj==cumSumChip$pj)&(cumSumInput$pj!=0)),]$pj
	arrowx=cumSumChip[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$x

	#final=merge(cumSumInput,cumSumChip,by="x")
	#colnames(final)=c("x","y_input","y_chip")

	finalList=list("X-axis"=arrowx,
		"Y-Input"=cumSumInput[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$pj,
		"Y-Chip"=cumSumChip[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$pj,
		"Ch_sign_chipVSinput"=sign_sign,
		"FractionReadsTopbins_chip"=chipFractionOfReadsIntop1percentBins,
		"FractionReadsTopbins_input"=inputFractionOfReadsIntop1percentBins,
		"Fractions_without_reads_chip"=chipFractionOfBinsWithoutReads,
		"Fractions_without_reads_input"=inputFractionOfBinsWithoutReads)
		#"CrossPoint_X"=cross_pointX,
		#"CrossPoint_Y_chip"=cross_pointY_chip,
		#"CrossPoint_Y_input"=cross_pointY_input)
	
	write.table(finalList,file=paste(getwd(),"Chance.results",sep=""))

	return(finalList)
}