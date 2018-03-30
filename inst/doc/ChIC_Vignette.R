### R code from vignette source '/lustre/data/FF/Carmen/BitBucket/Bioc_submission/FinalVignette/vignettes/ChIC_Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ChIC_Vignette.Rnw:58-66
###################################################
##load ChIC
library(ChIC)
## set path for working directory
filepath=tempdir()
setwd(filepath)

data("chipBam", package = "ChIC.data", envir = environment())
data("inputBam", package = "ChIC.data", envir = environment())


###################################################
### code chunk number 2: ChIC_Vignette.Rnw:77-93 (eval = FALSE)
###################################################
## ##caluclate first set of QC-metrics: EM 
## mc=4
## filepath=tempdir()
## setwd(filepath)
## 
## system("wget 
## https://www.encodeproject.org/files/ENCFF000BLL/@@download/ENCFF000BLL.bam")
## system("wget 
## https://www.encodeproject.org/files/ENCFF000BKA/@@download/ENCFF000BKA.bam")
## 
## chipName=file.path(filepath,"ENCFF000BLL")
## inputName=file.path(filepath,"ENCFF000BKA")
## 
## CC_Result=qualityScores_EM(chipName=chipName, inputName=inputName, 
## read_length=36, mc=mc)
## finalTagShift=CC_Result$QCscores_ChIP$tag.shift


###################################################
### code chunk number 3: ChIC_Vignette.Rnw:120-124 (eval = FALSE)
###################################################
## chipName=file.path(filepath,"ENCFF000BLL")
## inputName=file.path(filepath,"ENCFF000BKA")
## chipBam=readBamFile(chipName)
## inputBam=readBamFile(inputName)


###################################################
### code chunk number 4: ChIC_Vignette.Rnw:142-143
###################################################
mc=36


###################################################
### code chunk number 5: ChIC_Vignette.Rnw:147-159
###################################################
## calculate binding characteristics 

chip_binding.characteristics<-spp::get.binding.characteristics(
    chipBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)

input_binding.characteristics<-spp::get.binding.characteristics(
    inputBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)

## calculate cross correlation QC-metrics
crossvalues_Chip<-getCrossCorrelationScores(chipBam, 
chip_binding.characteristics, read_length=36,savePlotPath=filepath)
finalTagShift=crossvalues_Chip$tag.shift


###################################################
### code chunk number 6: ChIC_Vignette.Rnw:165-168 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics for input
## crossvalues_input<-getCrossCorrelationScores(inputBam, 
## chip_binding.characteristics, read_length=36,savePlotPath=filepath)


###################################################
### code chunk number 7: ChIC_Vignette.Rnw:193-208
###################################################
##get chromosome information and order chip and input by it
chrl_final=intersect(names(chipBam$tags),
    names(inputBam$tags))
chipBam$tags=chipBam$tags[chrl_final]
chipBam$quality=chipBam$quality[chrl_final]
inputBam$tags=inputBam$tags[chrl_final]
inputBam$quality=inputBam$quality[chrl_final]

##remove sigular positions with extremely high read counts with 
##respect to the neighbourhood
selectedTags=removeLocalTagAnomalies(chipBam, inputBam, 
chip_binding.characteristics, input_binding.characteristics)

inputBamSelected=selectedTags$input.dataSelected
chipBamSelected=selectedTags$chip.dataSelected


###################################################
### code chunk number 8: ChIC_Vignette.Rnw:216-221
###################################################
##Finally run function
bindingScores=getPeakCallingScores(chip=chipBam, 
input=inputBam, chip.dataSelected=chipBamSelected, 
input.dataSelected=inputBamSelected, 
tag.shift=finalTagShift)


###################################################
### code chunk number 9: ChIC_Vignette.Rnw:233-237
###################################################
smoothedChip=tagDensity(chipBamSelected, 
    tag.shift=finalTagShift)
smoothedInput=tagDensity(inputBamSelected, 
    tag.shift=finalTagShift)


###################################################
### code chunk number 10: ChIC_Vignette.Rnw:253-255
###################################################
Ch_Results=qualityScores_GM(densityChip=smoothedChip,
densityInput=smoothedInput,savePlotPath=filepath)


###################################################
### code chunk number 11: ChIC_Vignette.Rnw:288-291
###################################################
Meta_Result=createMetageneProfile(
    smoothed.densityChip=smoothedChip, 
    smoothedInput,tag.shift=finalTagShift, mc=mc)


###################################################
### code chunk number 12: ChIC_Vignette.Rnw:297-301
###################################################
TSS_Scores=qualityScores_LM(data=Meta_Result$TSS, tag="TSS",
savePlotPath=filepath)
TES_Scores=qualityScores_LM(data=Meta_Result$TES, tag="TES",
savePlotPath=filepath)


###################################################
### code chunk number 13: ChIC_Vignette.Rnw:307-310
###################################################
#create scaled metagene profile
geneBody_Scores=qualityScores_LMgenebody(Meta_Result$geneBody,
savePlotPath=filepath)


###################################################
### code chunk number 14: ChIC_Vignette.Rnw:338-342
###################################################
metagenePlotsForComparison(data=Meta_Result$geneBody,
    chrommark="H3K36me3", tag="geneBody", savePlotPath=filepath)
metagenePlotsForComparison(data=Meta_Result$TSS,
    chrommark="H3K36me3", tag="TSS", savePlotPath=filepath)


###################################################
### code chunk number 15: ChIC_Vignette.Rnw:353-355
###################################################
plotReferenceDistribution(chrommark="H3K4me1", 
  metricToBePlotted="RSC", currentValue=0.49, savePlotPath=filepath)


###################################################
### code chunk number 16: ChIC_Vignette.Rnw:377-379
###################################################
data("EM_scores", package = "ChIC.data", envir = environment())
CC_Result=EM_scores


###################################################
### code chunk number 17: ChIC_Vignette.Rnw:382-389
###################################################
te=predictionScore(chrommark="H3K36me3", 
    features_cc=CC_Result,
    features_global=Ch_Results,
    features_TSS=TSS_Scores, 
    features_TES=TES_Scores, 
    features_scaled=geneBody_Scores)
print(te)


