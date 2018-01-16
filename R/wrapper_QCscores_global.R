####FOR DEVEL#########
# debug=TRUE
# path=getwd()
# dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
# chipName="ENCFF000BBB"
# inputName="ENCFF000BAF"
# source("FunctionsGlobal.R")
# source(file.path(path,"GlobalParameters.R"))
# load(file.path(getwd(), paste(chipName, inputName, "TagDensityChip.RData", sep="_")))
# load(file.path(getwd(), paste(chipName, inputName, "TagDensityInput.RData", sep="_")))
####END DEVEL#########


#'@title Metrics taken from global read distribution
#'
#' @description
#' This set of values is based on the global read distribution along the genome for immunoprecipitation and input data (Diaz et al., 2012). 
#' The genome is binned and the read coverage counted for each bin. Then the function computes the cumulative distribution of reads density per genomic bin and plots the 
#' fraction of the coverage on the y-axis and the fraction of bins on the x-axis. Then different values can be sampled from the cumulative distribution: 
#' like the fraction of bins without reads for in immunoprecipitation and input,the point of the maximum distance between the ChIP and the input (x-axis, y-axis for 
#' immunoprecipitation and input, distance (as absolute difference), the sign of the differences), the fraction of reads in the top 1%bin for immunoprecipitation and input. 
#' Finally, the funciton returns 9 QC-measures
#'
#' f_QCscores_global
#'
#' @param densityChip A smoothed tag density object for the ChIP (returned by f_CrossCorrelation)
#' @param densityInput smoothed tag density object for the Inpt (returned by f_CrossCorrelation)
#' @param plotname Name of the FingerPring plot (default is "chancePlot.pdf" and saved in the working directory )
#' @param debug Only for debugging
#'
#' @return List with the following values: 
#' X-axis, Y-Input,Y-Chip (the point of the maximum distance between the ChIP and the input (x-axis, y-axis for immunoprecipitation and input) ,
#' Ch_sign_chipVSinput (the sign of the maximum distance) 
#' FractionReadsTopbins_chip, FractionReadsTopbins_input, Fractions_without_reads_chip, Fractions_without_reads_input
#'
#'
#' @examples
#'\dontrun{
#'CC_Result=f_CrossCorrelation(chipName, inputName, read_length=36, reads.aligner.type="bam", path=path, dataPath=dataPath, debug=debug,cluster=cluster,chrominfo_file=chrominfo_file,plotname=file.path(getwd(),paste(chipName,"CrossCorelationPlot.pdf",sep="_")))
#'smoothedDensityInput=CC_Result$TagDensityInput
#'smoothedDensityChip=CC_Result$TagDensityChip
#'plotname=file.path(getwd(),paste(chipName,"FingerPrintPlot.pdf",sep="_"))
#'Ch_Results=f_QCscores_global(densityChip=smoothedDensityChip,densityInput=smoothedDensityInput,plotname=plotname,debug=FALSE)
#' }



f_QCscores_global=function(densityChip,densityInput,plotname=file.path(getwd(),"FingerPrintPlot.pdf",sep="_"),debug=FALSE)
{
	##sourcing functions
	source("FunctionsGlobal.R")

	#shorten frame and reduce resolution
	print("shorten frame")
	chip.smoothed.density=f_shortenFrame(densityChip)
	input.smoothed.density=f_shortenFrame(densityInput)
	
	chip<-unlist(lapply(chip.smoothed.density, FUN=function(list_element) {return(list_element$y)}))
	input<-unlist(lapply(input.smoothed.density, FUN=function(list_element) {return(list_element$y)}))

	##create cumulative function for chip and input
	print("Calculate cumsum")
	cumSumChip=f_sortAndBinning(chip)
	cumSumInput=f_sortAndBinning(input)

	##caluclate QCvalues for chip
	chipFractionOfReadsIntop1percentBins<-(1-cumSumChip[(which(cumSumChip$x>=0.99)[1]),"pj"])
	chipFractionOfBinsWithoutReads<-cumSumChip$x[(which(cumSumChip$pj>0)[1])]

	##caluclate QCvalues for input
	inputFractionOfReadsIntop1percentBins<-(1-cumSumInput[(which(cumSumInput$x>=0.99)[1]),"pj"])
	inputFractionOfBinsWithoutReads<-cumSumInput$x[(which(cumSumInput$pj>0)[1])]

	##get the maximum distance between the two functions and the sign of the distance
	sign_sign=sign(max(cumSumInput$pj-cumSumChip$pj)) ##input-chip
	arrowx=cumSumChip[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$x

	##prepare list to be returned
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

	#create Fingerprint plot
	if (!is.null(plotname))
	{
		print("create pdf with plot")
		f_fingerPrintPlot(cumSumChip,cumSumInput,plotname=plotname)
	}
	#f_chancePlots(cumSumChip,cumSumInput,plotname=file.path(getwd(),paste(chipName,"chance.pdf",sep="_")))

	if (debug)
	{
		write.table(finalList,file=file.path(getwd(),"Chance.results"))
	}
	#return QC values
	return(finalList)
}
