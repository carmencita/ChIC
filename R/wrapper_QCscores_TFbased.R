

#MAIN

#library("snow")
#library("girafe")
#library("caTools")

require("spp")
##neds caTools
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

	bindingAnalysis=f_getBindingRegionsScores(chip.data,input.data,	final.tag.shift)

	bindingScores=bindingAnalysis$QCscoreList
	input.dataSelected=bindingAnalysis$input.dataSelected
	chip.dataSelected=bindingAnalysis$chip.dataSelected


	f_tagDensity=function(data,parallel.mc=NULL)
	##takes dataSelected as input, parallel is the number of CPUs used for parallelization
	{
		## density distribution for data
		print("Smooth tag density")
		ts <- sum(unlist(lapply(data,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
		##parallelisation
		chromosomes_list<-names(data)
		##creates a list of lists
		data<-lapply(chromosomes_list, FUN=function(x) {
			return(data[x])
		})
		if (parallel!=NULL)
		{
			smoothed.density<-mclapply(data, FUN=function(current_chr_list)
			{
			    current_chr<-names(current_chr_list)
			    str(current_chr_list)
			    if (length(current_chr) != 1) 
			    {
			        stop("unexpected input.dataSelected structure")
			    }
			    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
			}, mc.preschedule = FALSE,mc.cores=parallel.mc)
		}else{
			smoothed.density<-lapply(data, FUN=function(current_chr_list)
			{
			    current_chr<-names(current_chr_list)
			    str(current_chr_list)
			    if (length(current_chr) != 1) 
			    {
			        stop("unexpected input.dataSelected structure")
			    }
			    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
			})
		}
		smoothed.density=(unlist(smoothed.density,recursive=FALSE))
		#normalizing smoothed tag density by library size
		smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })
		return(smoothed.density)
	}

	smoothed.densityChip=f_tagDensity(chip.dataSelected)
	smoothed.densityInput=f_tagDensity(input.dataSelected)
	

	save(smoothed.densityChip,file=file.path(path,"TagDensityChip.RData"))
	save(smoothed.densityInput,file=file.path(path,"TagDensityInput.RData"))

	returnList=list("QCscores_ChIP"=crossvalues_Chip,
		"QCscores_Input"=crossvalues_Input,
		"QCscores_binding"=bindingScores,
		"TagDensityChip"=smoothed.densityChip,
		"TagDensityInput"=smoothed.densityInput)

	return(returnList)

################RESULTS

file.remove(outname)
#sampleIndex datafilename
#readCount read_length
write.table(crossvalues_Chip,file=paste(getwd(),"CC_chip.results",sep=""))
write.table(crossvalues_Input,file=paste(getwd(),"CC_Input.results",sep=""))
write.table(bindingScores,file=paste(getwd(),"CC_BindingScores.results",sep=""))

}