

#MAIN

#library("snow")
#library("girafe")
#library("caTools")

require("spp")
##neds caTools
library("girafe")
#private variables

source("Functions.R")
path=getwd()
source(paste(path,"GlobalParameters.R",sep="/"))
chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
rownames(chrominfo)<-chrominfo$chrom
rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})


CrossCorrelationInput=function(chipName="ENCFF000BBB", inputName="ENCFF000BAF",read_length=36,reads.aligner.type<-"bam",path=getwd())
{

	###read files
	chip.data=f_readFile(chipName,f_path="/lustre/data/FF/Carmen/BitBucket/chic/data/")
	input.data=f_readFile(inputName,f_path="/lustre/data/FF/Carmen/BitBucket/chic/data/")

	#plot and calculate cross correlation and phantom
	#chip_tagdistribution<-sapply(chip.data$tags, length)
	#input_tagdistribution<-sapply(input.data$tags, length)
	cluster=NULL
	
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
	save(input.dataSelected,final.tag.shift,chip.dataSelected,file=paste(getwd(),"dataSelected.RData",sep=""))


	bindingAnalysis=f_getBindingRegionsScores(chip.data,input.data,	chip.dataSelected,input.dataSelected,final.tag.shift,custom_chrorder)

	bindingScores=bindingAnalysis$QCscoreList
	input.dataSelected=bindingAnalysis$input.dataSelected
	chip.dataSelected=bindingAnalysis$chip.dataSelected


	smoothed.densityChip=f_tagDensity(chip.dataSelected)
	smoothed.densityInput=f_tagDensity(input.dataSelected)
	

	returnList=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores,
		"TagDensityChip"=smoothed.densityChip,
		"TagDensityInput"=smoothed.densityInput)

	return(returnList)

################RESULTS save for different purposes

save(smoothed.densityChip,file=file.path(path,"TagDensityChip.RData"))
save(smoothed.densityInput,file=file.path(path,"TagDensityInput.RData"))

write.table(crossvalues_Chip,file=paste(getwd(),"CC_chip.results",sep=""))
write.table(crossvalues_Input,file=paste(getwd(),"CC_Input.results",sep=""))
write.table(bindingScores,file=paste(getwd(),"CC_BindingScores.results",sep=""))

}