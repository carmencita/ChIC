
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
#require(parallel)
#library("caTools")

path=getwd()
dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
#source(paste(path,"GlobalParameters.R",sep="/"))
chrominfo_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/hg19.chromInfo.txt"



# chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
# rownames(chrominfo)<-chrominfo$chrom
# rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})
# #chipName="ENCFF000BBB"
#inputName="ENCFF000BAF"
chipName="ENCFF000BLL"
inputName="ENCFF000BKA"
#debug=TRUE
cluster=NULL

##### FOR DEVEL ONLY END
#########################

frameOfAllValues=NULL

source("wrapper_QCscores_TFbased.R")

CC_Result=f_CrossCorrelation(chipName, inputName, read_length=36, reads.aligner.type="bam", path=path, dataPath=dataPath, debug=debug,cluster=cluster,chrominfo_file=chrominfo_file)

##add values to completeListOfValues
completeListOfValues=rbind(cbind(unlist(CC_Result$QCscores_ChIP)),cbind(unlist(CC_Result$QCscores_binding)))

##CC_Result$TagDensityChip    #CC_Result$TagDensityInput
#CC_Result$QCscores_Input not needed but are calculated

##save the tagshift as it is needed later
tag.shift=CC_Result$QCscores_ChIP$tag.shift
##smoothed tagdensities used as input for the next steps
smoothedDensityInput=CC_Result$TagDensityInput
smoothedDensityChip=CC_Result$TagDensityChip

###GLOBAL features#########

source("wrapper_QCscores_global.R")
Ch_Results=f_QCscores_global(densityChip=smoothedDensityChip,densityInput=smoothedDensityInput,plotname=plotname=file.path(getwd(),paste(chipName,"chance.pdf",sep="_")),debug=FALSE)
completeListOfValues=append(completeListOfValues,Ch_Results)
###LOCAL


source("wrapper_QCscores_local.R")
require(parallel)
mc=1

geneAnnotations_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/RefSeqAllGenesFiltered.RData"
#CreateMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,path=NULL,debug=FALSE)
Meta_Result=f_CreateMetageneProfile(smoothedDensityChip,smoothedDensityInput,tag.shift,path=path,geneAnnotations_file=geneAnnotations_file,debug=FALSE)

source("wrapper_plot_TSS_TES_allGenes.R")
TSS_Plot=f_plotMetageneProfile_onePoint(Meta_Result$TSS$chip,Meta_Result$TSS$input,tag="TSS",path=getwd(),plotName=chipName,debug=FALSE)
completeListOfValues=append(completeListOfValues,TSS_Plot)

TES_Plot=f_plotMetageneProfile_onePoint(Meta_Result$TES$chip,Meta_Result$TES$input,tag="TES",path=getwd(),plotName=chipName,debug=FALSE)
completeListOfValues=append(completeListOfValues,TES_Plot)

source("wrapper_twopoint_allGenes_plot.R")
geneBody_Plot=f_plotMetageneProfile(Meta_Result$twopoint$chip,Meta_Result$twopoint$input,path=getwd(),plotName=chipName,debug=FALSE)
completeListOfValues=append(completeListOfValues,geneBody_Plot)



##additional plots
source("createMetaGeneProfilePlots_ForComparison.R")

f_metagenePlotsForComparison(chrommark="H3K4me1",Meta_Result$twopoint, Meta_Result$TSS, Meta_Result$TES,plotName=chipName,profilePath="/lustre/data/FF/Carmen/BitBucket/chic/data/Profiles")
