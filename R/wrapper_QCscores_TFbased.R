

# #MAIN

# #library("snow")
# #library("girafe")
# #library("caTools")


# #########################
# ##### FOR DEVEL ONLY

# require("spp")
# ##neds caTools
# library("girafe")
# #private variables

# source("Functions.R")
# path=getwd()
# dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
# source(paste(path,"GlobalParameters.R",sep="/"))
# chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
# rownames(chrominfo)<-chrominfo$chrom
# rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})
# chipName="ENCFF000BBB"
# inputName="ENCFF000BAF"
# debug=TRUE
# cluster=NULL
# dataPath="/lustre//data/FF/Carmen/BitBucket/chic/data"
# #test=f_CrossCorrelation(chipName, inputName, 36, "bam", path, dataPath, debug=TRUE)
# ##### FOR DEVEL ONLY END
# #########################




f_CrossCorrelation=function(chipName, inputName, read_length=36, reads.aligner.type="bam", path=getwd(), dataPath=getwd(), debug=FALSE,cluster=NULL,chrominfo_file=pwd())
{
	source("Functions.R")

	chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
	rownames(chrominfo)<-chrominfo$chrom
	rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})

	###read files
	chip.data=f_readFile(chipName,f_path=dataPath)
	input.data=f_readFile(inputName,f_path=dataPath)

	#plot and calculate cross correlation and phantom
	#chip_tagdistribution<-sapply(chip.data$tags, length)
	#input_tagdistribution<-sapply(input.data$tags, length)
	
	
	chip_binding.characteristics<-get.binding.characteristics(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)
	crossvalues_Chip=f_calculateCrossCorrelation(chip.data,chip_binding.characteristics,plotname=paste("phantomCrossCorrelation",chipName,sep="_"))
	final.tag.shift= crossvalues_Chip$tag.shift

	input_binding.characteristics<-get.binding.characteristics(input.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)
	crossvalues_Input=f_calculateCrossCorrelation(input.data,input_binding.characteristics,plotname=paste("phantomCrossCorrelation",inputName,sep="_"))

	chrl_final=intersect(names(chip.data$tags),names(input.data$tags))
	chip.data$tags=chip.data$tags[chrl_final]
	chip.data$quality=chip.data$quality[chrl_final]
	input.data$tags=input.data$tags[chrl_final]
	input.data$quality=input.data$quality[chrl_final]



	selectInformativeTags=f_selectInformativeTag(chip.data,input.data,chip_binding.characteristics,input_binding.characteristics)

		##Save data object to be read afterwards for the densityt distribution

	input.dataSelected=selectInformativeTags$input.dataSelected
	chip.dataSelected=selectInformativeTags$chip.dataSelected

	#if (debug) {
	#	save(input.dataSelected,final.tag.shift,chip.dataSelected,file=file.path(getwd(), paste(chipName, inputName, "dataSelected.RData", sep="_")))
	#}


	bindingAnalysis=f_getBindingRegionsScores(chip.data,input.data,	chip.dataSelected,input.dataSelected,final.tag.shift)#,custom_chrorder)

	bindingScores=bindingAnalysis
	#input.dataSelected=bindingAnalysis$input.dataSelected
	#chip.dataSelected=bindingAnalysis$chip.dataSelected


	smoothed.densityChip=f_tagDensity(chip.dataSelected,final.tag.shift,rngl=rngl)
	smoothed.densityInput=f_tagDensity(input.dataSelected,final.tag.shift,rngl=rngl)
	

	returnList=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores,
		"TagDensityChip"=smoothed.densityChip,
		"TagDensityInput"=smoothed.densityInput)


	return(returnList)

}






################RESULTsS save for different purposes

	# if (debug) {
	# 	save(smoothed.densityChip,file=file.path(path, paste(chipName, inputName, "TagDensityChip.RData", sep="_")))
	# 	save(smoothed.densityInput,file=file.path(path, paste( chipName, inputName, "TagDensityInput.RData", sep="_")))

	# 	write.table(crossvalues_Chip,file=file.path(path, paste(chipName, inputName, "chip.results", sep="_")))
	# 	write.table(crossvalues_Input,file=file.path(path, paste(chipName, inputName, "Input.results", sep="_")))
	# 	write.table(bindingScores,file=file.path(path, paste(chipName, inputName, "BindingScores.results", sep="_")))


	# 	print("write results in file")
	# 	################RESULTS
	# 	outname=file.path(path,paste(chipName,"TF.results",sep="_"))
	# 	file.remove(outname)
	# 	#sampleIndex datafilename
	# 	#readCount read_length
		
	# 	write(paste("Input: ",chipName,sep=" "),file=outname,append=T)
	# 	write(paste("Input: ",chipName ,sep=" "),file=outname,append=T)

	# 	write(paste("ReadCount",readCount,sep=" "),file=outname,append=T)
	# 	write(paste("Read_length",read_length,sep=" "),file=outname,append=T)
	# 	#strandShift<-binding.characteristics$peak$x
	# 	write(paste("StrandShift",strandShift,sep=" "),file=outname,append=T)
	# 	write(paste("Substitution of StrandShift from ",oldShift," to ",strandShift,sep=" "),file=outname,append=T)

	# 	#binding.characteristics$peak$whs
	# 	write(paste("bindingCharacteristicsPeak (x,y,whs)",binding.characteristics$peak$x,binding.characteristics$peak$y,binding.characteristics$whs,sep=" "),file=outname,append=T)

	# 	#phantom.characteristics$peak$x 
	# 	#phantom.characteristics$peak$y
	# 	#phantom.characteristics$peak$whs

	# 	write(paste("phantomCharacteristicsPeak (x,y,whs)",phantom.characteristics$peak$x,phantom.characteristics$peak$y,phantom.characteristics$whs,sep=" "),file=outname,append=T)
	# 	#phantomScores

	# 	write(paste("ALL_TAGS",ALL_TAGS,sep=" "),file=outname,append=T)
	# 	write(paste("NSC",round(NSC, 2),sep=" "),file=outname,append=T)
	# 	write(paste("RSC",round(RSC, 2),sep=" "),file=outname,append=T)
	# 	write(paste("Quality flag: ", qflag,sep=" "),file=outname,append=T)
	# 	write(paste("shift: ", round(as.double(phantomScores["shift"]),2),sep=" "),file=outname,append=T)

	# 	write(paste("read length",round(as.double(phantomScores["read_length"]),2),sep=" "),file=outname,append=T)
	# 	write(paste("A: ", round(as.double(phantomScores["A"]),2),sep=" "),file=outname,append=T)
	# 	write(paste("B: ", round(as.double(phantomScores["B"]),2),sep=" "),file=outname,append=T)
	# 	write(paste("C: ", round(as.double(phantomScores["C"]),2),sep=" "),file=outname,append=T)

	# 	write(paste("FDR detected",sum(unlist(lapply(bp_FDR$npl,function(d) length(d$x)))),"peaks",sep=" "),file=outname,append=T)
	# 	write(paste("eval detected",sum(unlist(lapply(bp_eval$npl,function(d) length(d$x)))),"peaks",sep=" "),file=outname,append=T)

	# 	#STATS_NRF

	# 	write(paste("UNIQUE_TAGS_LibSizeadjusted",UNIQUE_TAGS_LibSizeadjusted,sep=" "),file=outname,append=T)
	# 	write(paste("NRF_LibSizeadjusted",NRF_LibSizeadjusted,sep=" "),file=outname,append=T)
	# 	write(paste("ALL_TAGS",ALL_TAGS,sep=" "),file=outname,append=T)
	# 	write(paste("UNIQUE_TAGS",UNIQUE_TAGS,sep=" "),file=outname,append=T)
	# 	write(paste("UNIQUE_TAGS_nostrand",UNIQUE_TAGS_nostrand,sep=" "),file=outname,append=T)
	# 	write(paste("NRF",NRF,sep=" "),file=outname,append=T)
	# 	write(paste("NRF_LibSizeadjusted",NRF_LibSizeadjusted,sep=" "),file=outname,append=T)
	# 	write(paste("NRF_nostrand",NRF_nostrand,sep=" "),file=outname,append=T)
	# 	write(paste("PBC",PBC,sep=" "),file=outname,append=T)
	# 	write(paste("N1",N1,sep=" "),file=outname,append=T)
	# 	write(paste("Nd",Nd,sep=" "),file=outname,append=T)
	# 	#Frip
	# 	write(paste("Total_reads",TOTAL_reads,"FRiP_broadPeak",round(FRiP_broadPeak, 2),"outcountsBroadPeak", outcountsBroadPeak, "FRiP_sharpPeak", FRiP_sharpPeak, "outcountsSharpPeak", outcountsSharpPeak, sep=" "),file=outname,append=T)


	# }
