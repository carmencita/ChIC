### R code from vignette source '/lustre/data/FF/Carmen/BitBucket/Bioc_submission/ChIC/vignettes/ChIC_Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ChIC_Vignette.Rnw:59-73
###################################################
##load ChIC
library(ChIC)

## set path for working directory
#filepath=tempdir()
#setwd(filepath)
filepath=getwd()
#load tag-list with reads aligned to a subset of chromosomes


data("chipSubset", package = "ChIC.data", envir = environment())
chipBam=chipSubset
data("inputSubset", package = "ChIC.data", envir = environment())
inputBam=inputSubset


###################################################
### code chunk number 2: ChIC_Vignette.Rnw:84-101 (eval = FALSE)
###################################################
## ##caluclate first set of QC-metrics: EM 
## mc=4
## 
## filepath=tempdir()
## setwd(filepath)
## 
## system("wget 
## https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")
## system("wget 
## https://www.encodeproject.org/files/ENCFF000BDQ/@@download/ENCFF000BDQ.bam")
## 
## chipName=file.path(filepath,"ENCFF000BFX")
## inputName=file.path(filepath,"ENCFF000BDQ")
## 
## CC_Result=qualityScores_EM(chipName=chipName, inputName=inputName, 
## read_length=36, mc=mc)
## finalTagShift=CC_Result$QCscores_ChIP$tag.shift


###################################################
### code chunk number 3: ChIC_Vignette.Rnw:127-131 (eval = FALSE)
###################################################
## chipName=file.path(filepath,"ENCFF000BFX")
## inputName=file.path(filepath,"ENCFF000BDQ")
## chipBam=readBamFile(chipName)
## inputBam=readBamFile(inputName)


###################################################
### code chunk number 4: ChIC_Vignette.Rnw:149-165
###################################################

data("chipSubset", package = "ChIC.data", envir = environment())
new=list(tags=list(chr17=chipSubset$tags$chr17),
    quality= list(chr17=chipSubset$quality$chr17))
chipBam=new

data("inputSubset", package = "ChIC.data", envir = environment())
new=list(tags=list(chr17=inputSubset$tags$chr17),
    quality= list(chr17=inputSubset$quality$chr17))

inputBam=new


mc=1
data("crossvalues_Chip", package = "ChIC.data", envir = environment())
#tagshift=98


###################################################
### code chunk number 5: ChIC_Vignette.Rnw:169-183
###################################################
cluster <- parallel::makeCluster( mc )

## calculate binding characteristics 

chip_binding.characteristics<-spp::get.binding.characteristics(
    chipBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE, 
    cluster = cluster)

input_binding.characteristics<-spp::get.binding.characteristics(
    inputBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE,
    cluster = cluster)

parallel::stopCluster( cluster )



###################################################
### code chunk number 6: ChIC_Vignette.Rnw:186-190 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics
## crossvalues_Chip<-getCrossCorrelationScores( chipBam , 
##     chip_binding.characteristics, read_length = 36, 
##     savePlotPath = filepath, mc = mc)


###################################################
### code chunk number 7: ChIC_Vignette.Rnw:196-199
###################################################
str(crossvalues_Chip)

finalTagShift <- crossvalues_Chip$tag.shift


###################################################
### code chunk number 8: ChIC_Vignette.Rnw:205-209 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics for input
## crossvalues_input <- getCrossCorrelationScores(inputBam, 
##     chip_binding.characteristics, read_length = 36, 
##     savePlotPath = filepath, mc = mc)


###################################################
### code chunk number 9: ChIC_Vignette.Rnw:234-240
###################################################
##get chromosome information and order chip and input by it
chrl_final <- intersect(names(chipBam$tags), names(inputBam$tags))
chipBam$tags <- chipBam$tags[chrl_final]
chipBam$quality <- chipBam$quality[chrl_final]
inputBam$tags <- inputBam$tags[chrl_final]
inputBam$quality <- inputBam$quality[chrl_final]


###################################################
### code chunk number 10: ChIC_Vignette.Rnw:243-250
###################################################
##remove sigular positions with extremely high read counts with 
##respect to the neighbourhood
selectedTags <- removeLocalTagAnomalies(chipBam, inputBam, 
    chip_binding.characteristics, input_binding.characteristics)

inputBamSelected <- selectedTags$input.dataSelected
chipBamSelected <- selectedTags$chip.dataSelected


###################################################
### code chunk number 11: ChIC_Vignette.Rnw:259-264
###################################################
##Finally run function
bindingScores <- getPeakCallingScores(chip = chipBam, 
    input = inputBam, chip.dataSelected = chipBamSelected, 
    input.dataSelected = inputBamSelected, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 12: ChIC_Vignette.Rnw:276-280
###################################################
smoothedChip <- tagDensity(chipBamSelected, 
    tag.shift = finalTagShift, mc = mc)
smoothedInput <- tagDensity(inputBamSelected, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 13: ChIC_Vignette.Rnw:296-298 (eval = FALSE)
###################################################
## Ch_Results <- qualityScores_GM(densityChip = smoothedChip,
##     densityInput = smoothedInput, savePlotPath = filepath)


###################################################
### code chunk number 14: Fingerprint (eval = FALSE)
###################################################
## Ch_Results=qualityScores_GM(densityChip=smoothedChip,
## densityInput=smoothedInput)


###################################################
### code chunk number 15: Fingerprint
###################################################
Ch_Results=qualityScores_GM(densityChip=smoothedChip,
densityInput=smoothedInput)


###################################################
### code chunk number 16: ChIC_Vignette.Rnw:339-343
###################################################
Meta_Result <- createMetageneProfile(
    smoothed.densityChip = smoothedChip, 
    smoothed.densityInput = smoothedInput, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 17: ChIC_Vignette.Rnw:350-354 (eval = FALSE)
###################################################
## TSS_Scores <- qualityScores_LM(data = Meta_Result$TSS, tag = "TSS",
##     savePlotPath = filepath)
## TES_Scores <- qualityScores_LM(data = Meta_Result$TES, tag = "TES",
##     savePlotPath = filepath)


###################################################
### code chunk number 18: ChIC_Vignette.Rnw:357-359
###################################################
TSS_Scores=qualityScores_LM(data=Meta_Result$TSS, tag="TSS")
TES_Scores=qualityScores_LM(data=Meta_Result$TES, tag="TES")


###################################################
### code chunk number 19: ChIC_Vignette.Rnw:365-368 (eval = FALSE)
###################################################
## #create scaled metagene profile
## geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody,
##     savePlotPath = filepath)


###################################################
### code chunk number 20: geneBody (eval = FALSE)
###################################################
## #create scaled metagene profile
## geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody)


###################################################
### code chunk number 21: geneBody
###################################################
#create scaled metagene profile
geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody)


###################################################
### code chunk number 22: ChIC_Vignette.Rnw:422-431 (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     chrommark = "H3K4me3", 
##     tag = "geneBody", 
##     savePlotPath = filepath)
## 
## metagenePlotsForComparison(data = Meta_Result$TSS,
##     chrommark = "H3K4me3", 
##     tag = "TSS", 
##     savePlotPath = filepath)


###################################################
### code chunk number 23: ComparisonInput (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     chrommark = "H3K4me3", tag = "geneBody")


###################################################
### code chunk number 24: ComparisonInput
###################################################
metagenePlotsForComparison(data = Meta_Result$geneBody,
    chrommark = "H3K4me3", tag = "geneBody")


###################################################
### code chunk number 25: ChIC_Vignette.Rnw:478-482 (eval = FALSE)
###################################################
## plotReferenceDistribution(chrommark = "H3K4me3", 
##     metricToBePlotted = "RSC", 
##     currentValue = crossvalues_Chip$CC_RSC, 
##     savePlotPath = filepath)


###################################################
### code chunk number 26: ReferenceDistr (eval = FALSE)
###################################################
## plotReferenceDistribution(chrommark = "H3K4me3", 
##     metricToBePlotted = "RSC", currentValue = crossvalues_Chip$CC_RSC )


###################################################
### code chunk number 27: ReferenceDistr
###################################################
plotReferenceDistribution(chrommark = "H3K4me3", 
    metricToBePlotted = "RSC", currentValue = crossvalues_Chip$CC_RSC )


###################################################
### code chunk number 28: ChIC_Vignette.Rnw:507-516
###################################################
EM_scoresNew=NULL

EM_scoresNew$QCscores_ChIP=crossvalues_Chip
EM_scoresNew$QCscores_binding=bindingScores
EM_scoresNew$TagDensityInput=list()
EM_scoresNew$TagDensityChip=list()

CC_Result=EM_scoresNew



###################################################
### code chunk number 29: ChIC_Vignette.Rnw:519-526
###################################################
te <- predictionScore(chrommark = "H3K4me3", 
    features_cc = CC_Result,
    features_global = Ch_Results,
    features_TSS = TSS_Scores, 
    features_TES = TES_Scores, 
    features_scaled = geneBody_Scores)
print(te)


