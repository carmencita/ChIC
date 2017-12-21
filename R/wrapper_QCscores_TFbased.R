

#MAIN

#library("snow")
#library("girafe")
#library("caTools")


#########################
##### FOR DEVEL ONLY

require("spp")
##neds caTools
library("girafe")
#private variables

source("Functions.R")
path=getwd()
dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
source(paste(path,"GlobalParameters.R",sep="/"))
chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
rownames(chrominfo)<-chrominfo$chrom
rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})
chipName="ENCFF000BBB"
inputName="ENCFF000BAF"
debug=TRUE
cluster=NULL
dataPath="/lustre//data/FF/Carmen/BitBucket/chic/data"
#test=f_CrossCorrelation(chipName, inputName, 36, "bam", path, dataPath, debug=TRUE)
##### FOR DEVEL ONLY END
#########################




f_CrossCorrelation=function(chipName, inputName, read_length=36, reads.aligner.type="bam", path=getwd(), dataPath=getwd(), debug=FALSE,cluster=NULL)
{

	###read files
	chip.data=f_readFile(chipName,f_path=dataPath)
	input.data=f_readFile(inputName,f_path=dataPath)

	#plot and calculate cross correlation and phantom
	#chip_tagdistribution<-sapply(chip.data$tags, length)
	#input_tagdistribution<-sapply(input.data$tags, length)
	
	
	chip_binding.characteristics<-get.binding.characteristics(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)
	crossvalues_Chip=f_calculateCrossCorrelation(chip.data,chip_binding.characteristics)
	final.tag.shift= crossvalues_Chip$tag.shift

	input_binding.characteristics<-get.binding.characteristics(input.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)
	crossvalues_Input=f_calculateCrossCorrelation(input.data,input_binding.characteristics)

	chrl_final=intersect(names(chip.data$tags),names(input.data$tags))
	chip.data$tags=chip.data$tags[chrl_final]
	chip.data$quality=chip.data$quality[chrl_final]
	input.data$tags=input.data$tags[chrl_final]
	input.data$quality=input.data$quality[chrl_final]



	selectInformativeTags=f_selectInformativeTag(chip.data,input.data,chip_binding.characteristics,input_binding.characteristics)

		##Save data object to be read afterwards for the densityt distribution

	input.dataSelected=selectInformativeTags$input.dataSelected
	chip.dataSelected=selectInformativeTags$chip.dataSelected

	if (debug) {
		save(input.dataSelected,final.tag.shift,chip.dataSelected,file=file.path(getwd(), paste(chipName, inputName, "dataSelected.RData", sep="_")))
	}


	bindingAnalysis=f_getBindingRegionsScores(chip.data,input.data,	chip.dataSelected,input.dataSelected,final.tag.shift,custom_chrorder)

	bindingScores=bindingAnalysis
	#input.dataSelected=bindingAnalysis$input.dataSelected
	#chip.dataSelected=bindingAnalysis$chip.dataSelected


	smoothed.densityChip=f_tagDensity(chip.dataSelected)
	smoothed.densityInput=f_tagDensity(input.dataSelected)
	

	returnList=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores,
		"TagDensityChip"=smoothed.densityChip,
		"TagDensityInput"=smoothed.densityInput)

	return(returnList)

################RESULTsS save for different purposes

	if (debug) {
		save(smoothed.densityChip,file=file.path(getwd(), paste(chipName, inputName, "TagDensityChip.RData", sep="_")))
		save(smoothed.densityInput,file=file.path(getwd(), paste( chipName, inputName, "TagDensityInput.RData", sep="_")))

		write.table(crossvalues_Chip,file=file.path(getwd(), paste(chipName, inputName, "chip.results", sep="_")))
		write.table(crossvalues_Input,file=file.path(getwd(), paste(chipName, inputName, "Input.results", sep="_")))
		write.table(bindingScores,file=file.path(getwd(), paste(chipName, inputName, "BindingScores.results", sep="_")))
	}




}



