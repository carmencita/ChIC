###FOR DEVEL#########
debug=TRUE
path=getwd()
dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
chipName="ENCFF000BBB"
inputName="ENCFF000BAF"
source("Functions.R")
source(file.path(path,"GlobalParameters.R"))
load(file.path(getwd(), paste(chipName, inputName, "TagDensityChip.RData", sep="_")))
load(file.path(getwd(), paste(chipName, inputName, "TagDensityInput.RData", sep="_")))

#densityChip=smoothed.densityChip
#densityInput=smoothed.densityInput
####
##test=QCscores_global(smoothed.densityChip,smoothed.densityInput)

QCscores_global=function(densityChip,densityInput)
{
	print("shorten frame")
	chip.smoothed.density=f_shortenFrame(densityChip)
	input.smoothed.density=f_shortenFrame(densityInput)
	

	#tag.shift <- round(strandShift/2)
	chip<-unlist(lapply(chip.smoothed.density, FUN=function(list_element) {return(list_element$y)}))
	input<-unlist(lapply(input.smoothed.density, FUN=function(list_element) {return(list_element$y)}))

	print("Calculate cumsum")
	cumSumChip=f_sortAndBinning(chip)
	cumSumInput=f_sortAndBinning(input)

	chipFractionOfReadsIntop1percentBins<-(1-cumSumChip[(which(cumSumChip$x>=0.99)[1]),"pj"])
	chipFractionOfBinsWithoutReads<-cumSumChip$x[(which(cumSumChip$pj>0)[1])]

	
	inputFractionOfReadsIntop1percentBins<-(1-cumSumInput[(which(cumSumInput$x>=0.99)[1]),"pj"])
	inputFractionOfBinsWithoutReads<-cumSumInput$x[(which(cumSumInput$pj>0)[1])]


	sign_sign=sign(max(cumSumInput$pj-cumSumChip$pj)) ##input-chip
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

	f_chancePlots(cumSumChip,cumSumInput,plotname=file.path(getwd(),"chance.pdf"))
	if (debug)
	{
		write.table(finalList,file=file.path(getwd(),"Chance.results"))
	}

	return(finalList)
}