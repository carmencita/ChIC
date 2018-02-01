#########################
##### FOR DEVEL ONLY

require("spp")

## load gene annotations
require(girafe)
require(snow)
require(parallel)
requre(caret)
path=getwd()

dataDirectory<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
#source(paste(path,"GlobalParameters.R",sep="/"))
#chrominfo_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/hg19.chromInfo.txt"

#chipName="ENCFF000BBB"
#inputName="ENCFF000BAF"
chipName="ENCFF000BLL"
inputName="ENCFF000BKA"
debug=TRUE

##### FOR DEVEL ONLY END
#########################

print(chipName)
print(inputName)

mc=5


#annotationPath="/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/"



##calculate values

source("PrivateFunctions.R")
source("GlobalFunctions.R")
source("wrapper_QCscores_TFbased.R")
source("wrapper_plot_TSS_TES_allGenes.R")
source("wrapper_twopoint_allGenes_plot.R")

source("createMetaGeneProfilePlots_ForComparison.R")

CC_Result=crossCorrelation(chipName=chipName,inputName=inputName, read_length=36, dataPath=dataDirectory, debug=debug,mc=mc,annotationID="hg19",savePlotPath=getwd())

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

Ch_Results=QCscores_global(densityChip=smoothedDensityChip,densityInput=smoothedDensityInput,savePlotPath=getwd(),debug=debug)
completeListOfValues=rbind(completeListOfValues,cbind(unlist(Ch_Results)))

###LOCAL features########
#geneAnnotations_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/RefSeqAllGenesFiltered.RData"
#CreateMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,path=NULL,debug=FALSE)
Meta_Result=createMetageneProfile(smoothedDensityChip,smoothedDensityInput,tag.shift,annotationID="hg19",debug=debug,mc=mc)

TSS_Plot=nonScaledMetageneProfile(Meta_Result$TSS$chip,Meta_Result$TSS$input,tag="TSS",savePlotPath=getwd(),debug=debug)
helper=TSS_Plot
helper$Feature=NULL
rownames(helper)=TSS_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)

TES_Plot=nonScaledMetageneProfile(Meta_Result$TES$chip,Meta_Result$TES$input,tag="TES",savePlotPath=getwd(),debug=debug)
helper=TES_Plot
helper$Feature=NULL
rownames(helper)=TES_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)

geneBody_Plot=scaledMetageneProfile(Meta_Result$twopoint$chip,Meta_Result$twopoint$input,savePlotPath=getwd(),debug=debug)
helper=geneBody_Plot
helper$Feature=NULL
rownames(helper)=geneBody_Plot$Feature
completeListOfValues=rbind(completeListOfValues,helper)

##additional plots

load("ENCFF000BLL_ENCFF000BKA_OnePointTES.RData")
load("ENCFF000BLL_ENCFF000BKA_OnePointTSS.RData")
load("ENCFF000BLL_ENCFF000BKA_Twopoint.RData")

metagenePlotsForComparison(chrommark="H3K36me3",Meta_Result$twopoint, Meta_Result$TSS, Meta_Result$TES,savePlotPath=getwd())
plotReferenceDistribution(chrommark="H3K36me3",metricToBePlotted="RSC",currentValue=CC_Result$QCscores_ChIP$CC_RSC,savePlotPath=getwd())

plotPredictionScore(chrommark="H3K36me3",featureVector=completeListOfValues,savePlotPath=getwd())
