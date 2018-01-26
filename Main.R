
#MAIN

#library("snow")
#library("girafe")
#library("caTools")


#########################
##### FOR DEVEL ONLY

require("spp")
##neds caTools
require("girafe")

#library(parallel)

#private variables
#require(parallel)
#library("caTools")

path=getwd()
dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
#source(paste(path,"GlobalParameters.R",sep="/"))
chrominfo_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/hg19.chromInfo.txt"

#chipName="ENCFF000BBB"
#inputName="ENCFF000BAF"
chipName="ENCFF000BLL"
inputName="ENCFF000BKA"
debug=TRUE

##### FOR DEVEL ONLY END
#########################

print(chipName)
print(inputName)

cluster=TRUE
mc=1
if (cluster)
{mc=5}




##calculate values
source("wrapper_QCscores_TFbased.R")

CC_Result=f_CrossCorrelation(chipName=chipName,inputName=inputName, read_length=36, reads.aligner.type="bam", dataPath=dataPath, debug=debug,cluster=cluster,mc=mc,chrominfo_file=chrominfo_file,savePlotPath=getwd())

##add values to completeListOfValues
completeListOfValues=NULL
completeListOfValues=rbind(cbind(unlist(CC_Result$QCscores_ChIP)),cbind(unlist(CC_Result$QCscores_binding)))
colnames(completeListOfValues)=c("Value")

##CC_Result$TagDensityChip    #CC_Result$TagDensityInput
#CC_Result$QCscores_Input not needed but are calculated

##save the tagshift as it is needed later
tag.shift=CC_Result$QCscores_ChIP$tag.shift
##smoothed tagdensities used as input for the next steps
smoothedDensityInput=CC_Result$TagDensityInput
smoothedDensityChip=CC_Result$TagDensityChip

###GLOBAL features#########

source("wrapper_QCscores_global.R")
Ch_Results=f_QCscores_global(densityChip=smoothedDensityChip,densityInput=smoothedDensityInput,savePlotPath=getwd(),debug=debug)
completeListOfValues=rbind(completeListOfValues,cbind(unlist(Ch_Results)))

###LOCAL features########


source("wrapper_QCscores_local.R")
#require(parallel)


#geneAnnotations_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/RefSeqAllGenesFiltered.RData"
#CreateMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,path=NULL,debug=FALSE)
Meta_Result=f_CreateMetageneProfile(smoothedDensityChip,smoothedDensityInput,tag.shift,annotationID="hg19",debug=debug)

source("wrapper_plot_TSS_TES_allGenes.R")
TSS_Plot=f_plotMetageneProfile_onePoint(Meta_Result$TSS$chip,Meta_Result$TSS$input,tag="TSS",savePlotPath=getwd(),debug=debug)
helper=TSS_Plot
helper$Feature=NULL
rownames(helper)=TSS_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)

TES_Plot=f_plotMetageneProfile_onePoint(Meta_Result$TES$chip,Meta_Result$TES$input,tag="TES",savePlotPath=getwd(),debug=debug)
helper=TES_Plot
helper$Feature=NULL
rownames(helper)=TES_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)

source("wrapper_twopoint_allGenes_plot.R")
geneBody_Plot=f_plotMetageneProfile(Meta_Result$twopoint$chip,Meta_Result$twopoint$input,savePlotPath=getwd(),debug=debug)
helper=geneBody_Plot
helper$Feature=NULL
rownames(helper)=geneBody_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)



##additional plots
source("createMetaGeneProfilePlots_ForComparison.R")

f_metagenePlotsForComparison(chrommark="H3K36me3",Meta_Result$twopoint, Meta_Result$TSS, Meta_Result$TES,savePlotPath=getwd())

f_plotReferenceDistribution(chrommark="H3K36me3",metricToBePlotted="RSC",currentValue=CC_Result$QCscores_ChIP$CC_RSC,savePlotPath=getwd())

f_plotPredictionScore(chrommark="H3K36me3",featureVector=completeListOfValues,savePlotPath=getwd())
