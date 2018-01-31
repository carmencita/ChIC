#'@title Wrapper function for calculating cross-correlation analysis and metrics designed for TFs
#'
#' @description
#' We use cross-correlation analysis to obtain QC-metrics proposed for narrow-binding patterns. 
#' After calculating the strand cross-correlation coefficient (Kharchenko et al., 2008), we take the following values from the profile: the coordinates of the 
#' ChIP-peak (fragment length, height A), the coordinates at the phantom-peak (read length, height B) and the baseline (C), the strand-shift, the number of uniquely mapped 
#' reads (unique_tags), uniquely mapped reads corrected by the library size, the number of reads and the read lengths. We  calculate different values using the relative and 
#' absolute height of the cross-correlation peaks: the relative and normalized strand coefficient RSC and NSC  (Landt et al., 2012), and the quality control tag (Fig. 1B) 
#'(Marinov et al., 2013). Other values regarding the library complexity (Landt et al., 2012) like the fraction of non-redundant mapped reads 
#'(NRF; ratio between the number of uniquely mapped reads divided by the total number of reads), the NRF adjusted by library size and ignoring the strand direction 
#' (NRF_nostrand), and the PCR bottleneck coefficient PBC (number of genomic locations to which exactly one unique mapping read maps, divided by the number of unique 
#' mapping reads). Other measures we include in our analysis are the fraction of usable reads in the peak regions (FRiP) (Landt et al., 2012), for which the function calls sharp- and 
#' broad-binding peaks to obtain two types: the FRiP_sharpsPeak and the FRiP_broadPeak. The function takes the number of called of peaks using an FDR of 0.01 and an 
#' evalue of 10 (Kharchenko et al., 2008). And count the number of peaks called when using the sharp- and broad-binding option. 
# 'Finally 22 features are given back.
#'
#' crossCorrelation
#'
#' @param chipName String, filename (without extension) of the input file for ChIP
#' @param inputName String, filename (without extension) of the input file for ChIP
#' @param read_length Integer, length of the reads
#' @param dataPath Path were to find the input files chipName and inputName, default is working directory
#' @param debug Boolean value to enter in debugging mode (default= FALSE)
#' @param chrominfo_file Path to the chromatin information file
#'
#' @return returnList, containing (!!!DESCRIBE BETTER)
#' QCscores_ChIP List with Crosscorrelation values of the ChIP
#' QCscores_Input List with Crosscorrelation values of the Input
#' QCscores_binding List with QCscores from called peaks
#' TagDensityChip (!!!DESCRIBE BETTER)
#' TagDensityInput (!!!DESCRIBE BETTER)
#'
#' @examples
#'\{dontrun
#' CC_Result=crossCorrelation(chipName, inputName, read_length=36,  dataPath=dataPath, debug=debug,chrominfo_file=chrominfo_file)
#'}


crossCorrelation=function(chipName, inputName, read_length, dataPath=getwd(), annotationID="hg19", savePlotPath=NULL, debug=FALSE, mc=1)
{
	##load rngl
	filename=paste(annotationID,".chromInfo.RData",sep="")
	print(paste("load ",filename))
	load(filename)

	#chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
	#rownames(chrominfo)<-chrominfo$chrom
	#rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})

	###read files
	print(paste("reading bam files",sep=" "))
	chip.data=f_readFile(chipName,f_path=dataPath)
	input.data=f_readFile(inputName,f_path=dataPath)

	# ORDERED_data<-chip.data
	# for (chr in names(chip.data$tags)) {
 #  		orderingVector<-order(abs(chip.data$tags[[chr]]))
 #  		ORDERED_data$tags[[chr]]<-chip.data$tags[[chr]][orderingVector]
 #  		ORDERED_data$quality[[chr]]<-chip.data$quality[[chr]][orderingVector]
	# }
	# chip.data=ORDERED_data

	# ORDERED_data<-input.data
	# for (chr in names(input.data$tags)) {
 #  		orderingVector<-order(abs(input.data$tags[[chr]]))
 #  		ORDERED_data$tags[[chr]]<-input.data$tags[[chr]][orderingVector]
 #  		ORDERED_data$quality[[chr]]<-input.data$quality[[chr]][orderingVector]
	# }
	# input.data=ORDERED_data



	#plot and calculate cross correlation and phantom characteristics for the ChIP
	print("calculate binding characteristics ChIP")
	## cross_correlation parameters
	estimating_fragment_length_range<-c(0,500)
	estimating_fragment_length_bin<-5

	#chip_binding.characteristics<-get.binding.characteristics(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, accept.all.tags=T)

	chip_binding.characteristics<-get.binding.characteristicsMy(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin,accept.all.tags=T)
	print("calculate cross correlation QC-metrics for the Chip")
	crossvalues_Chip<-calculateCrossCorrelation(chip.data,chip_binding.characteristics,read_length=read_length,savePlotPath=savePlotPath,plotname="ChIP")
	##save the tag.shift
	final.tag.shift<-crossvalues_Chip$tag.shift


	#plot and calculate cross correlation and phantom characteristics for the input
	print("calculate binding characteristics Input")
	
	input_binding.characteristics<-get.binding.characteristicsMy(input.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, accept.all.tags=T)
	print("calculate cross correlation QC-metrics for the Input")
	crossvalues_Input=calculateCrossCorrelation(input.data,input_binding.characteristics,read_length=read_length,savePlotPath=savePlotPath,plotname="Input")


	
	##get chromosome information and order chip and input by it
	chrl_final=intersect(names(chip.data$tags),names(input.data$tags))
	chip.data$tags=chip.data$tags[chrl_final]
	chip.data$quality=chip.data$quality[chrl_final]
	input.data$tags=input.data$tags[chrl_final]
	input.data$quality=input.data$quality[chrl_final]


	##select informative tags
	selectInformativeTags=f_selectInformativeTag(chip.data,input.data,chip_binding.characteristics,input_binding.characteristics)

	input.dataSelected=selectInformativeTags$input.dataSelected
	chip.dataSelected=selectInformativeTags$chip.dataSelected

	#get QC-values from peak calling
	bindingScores=getBindingRegionsScores(chip.data,input.data, chip.dataSelected,input.dataSelected,final.tag.shift)
	
	
	##objects of smoothed tag density for ChIP and Input
	smoothed.densityChip=f_tagDensity(chip.dataSelected,final.tag.shift,rngl=rngl,mc=mc)
	smoothed.densityInput=f_tagDensity(input.dataSelected,final.tag.shift,rngl=rngl,mc=mc)
	

	returnList=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores,
		"TagDensityChip"=smoothed.densityChip,
		"TagDensityInput"=smoothed.densityInput)


	if (debug)
	{
		writeout=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores)
		write.table(writeout,file=file.path(getwd(),"CC.results"))
	}

	return(returnList)

}

