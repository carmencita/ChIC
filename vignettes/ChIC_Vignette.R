### R code from vignette source '/lustre/data/FF/Carmen/BitBucket/Bioc_submission/ChIC/vignettes/ChIC_Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ChIC_Vignette.Rnw:66-89 (eval = FALSE)
###################################################
## ##calculating EM
## 
## mc=4 #for parallelization
## filepath <- tempdir()
## setwd(filepath)
## 
## system("wget 
## https://www.encodeproject.org/files/
##     ENCFF000BFX/@@download/ENCFF000BFX.bam")
## system("wget 
## https://www.encodeproject.org/files/
##     ENCFF000BDQ/@@download/ENCFF000BDQ.bam")
## 
## chipName <- file.path(filepath,"ENCFF000BFX")
## inputName <- file.path(filepath,"ENCFF000BDQ")
## 
## CC_Result <- qualityScores_EM(chipName=chipName, 
##     inputName=inputName, 
##     annotationID="hg19",
##     read_length=36, 
##     mc=mc)
## 
## finalTagShift <- CC_Result$QCscores_ChIP$tag.shift


###################################################
### code chunk number 2: ChIC_Vignette.Rnw:128-132 (eval = FALSE)
###################################################
## chipName <- file.path(filepath,"ENCFF000BFX")
## inputName <- file.path(filepath,"ENCFF000BDQ")
## chipBam <- readBamFile(chipName)
## inputBam <- readBamFile(inputName)


###################################################
### code chunk number 3: ChIC_Vignette.Rnw:139-148
###################################################
library(ChIC)
## set path for working directory
filepath=getwd()

#load tag-list with reads aligned to a subset of chromosomes
data("chipSubset", package = "ChIC.data", envir = environment())
chipBam <- chipSubset
data("inputSubset", package = "ChIC.data", envir = environment())
inputBam <- inputSubset


###################################################
### code chunk number 4: ChIC_Vignette.Rnw:162-164
###################################################
mc=5
data("crossvalues_Chip", package = "ChIC.data", envir = environment())


###################################################
### code chunk number 5: ChIC_Vignette.Rnw:167-180
###################################################
cluster <- parallel::makeCluster( mc )

## calculate binding characteristics 

chip_binding.characteristics <- spp::get.binding.characteristics(
    chipBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE, 
    cluster = cluster)

input_binding.characteristics <- spp::get.binding.characteristics(
    inputBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE,
    cluster = cluster)

parallel::stopCluster( cluster )


###################################################
### code chunk number 6: ChIC_Vignette.Rnw:183-190 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics
## crossvalues_Chip <- getCrossCorrelationScores( chipBam , 
##     chip_binding.characteristics, 
##     read_length = 36, 
##     annotationID="hg19",
##     savePlotPath = filepath, 
##     mc = mc)


###################################################
### code chunk number 7: ChIC_Vignette.Rnw:202-204
###################################################
str(crossvalues_Chip)
finalTagShift <- crossvalues_Chip$tag.shift


###################################################
### code chunk number 8: ChIC_Vignette.Rnw:222-228
###################################################
##get chromosome information and order chip and input by it
chrl_final <- intersect(names(chipBam$tags), names(inputBam$tags))
chipBam$tags <- chipBam$tags[chrl_final]
chipBam$quality <- chipBam$quality[chrl_final]
inputBam$tags <- inputBam$tags[chrl_final]
inputBam$quality <- inputBam$quality[chrl_final]


###################################################
### code chunk number 9: ChIC_Vignette.Rnw:231-240
###################################################
##remove sigular positions with extremely high read counts with 
##respect to the neighbourhood
selectedTags <- removeLocalTagAnomalies(chipBam, 
    inputBam, 
    chip_binding.characteristics, 
    input_binding.characteristics)

inputBamSelected <- selectedTags$input.dataSelected
chipBamSelected <- selectedTags$chip.dataSelected


###################################################
### code chunk number 10: ChIC_Vignette.Rnw:249-256
###################################################
bindingScores <- getPeakCallingScores(chip = chipBam, 
    input = inputBam, 
    chip.dataSelected = chipBamSelected, 
    input.dataSelected = inputBamSelected, 
    annotationID="hg19",
    tag.shift = finalTagShift, 
    mc = mc)


###################################################
### code chunk number 11: ChIC_Vignette.Rnw:270-276
###################################################
smoothedChip <- tagDensity(chipBamSelected,
    annotationID = "hg19", 
    tag.shift = finalTagShift, mc = mc)
smoothedInput <- tagDensity(inputBamSelected, 
    annotationID = "hg19",
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 12: ChIC_Vignette.Rnw:289-290
###################################################
listAvailableElements(target="H3K36me3")


###################################################
### code chunk number 13: ChIC_Vignette.Rnw:307-310 (eval = FALSE)
###################################################
## Ch_Results <- qualityScores_GM(densityChip = smoothedChip,
##     densityInput = smoothedInput, 
##     savePlotPath = filepath)


###################################################
### code chunk number 14: Fingerprint (eval = FALSE)
###################################################
## Ch_Results <- qualityScores_GM(densityChip=smoothedChip,
## densityInput=smoothedInput)


###################################################
### code chunk number 15: Fingerprint
###################################################
Ch_Results <- qualityScores_GM(densityChip=smoothedChip,
densityInput=smoothedInput)


###################################################
### code chunk number 16: ChIC_Vignette.Rnw:350-356
###################################################
Meta_Result <- createMetageneProfile(
    smoothed.densityChip = smoothedChip, 
    smoothed.densityInput = smoothedInput, 
    annotationID="hg19",
    tag.shift = finalTagShift, 
    mc = mc)


###################################################
### code chunk number 17: ChIC_Vignette.Rnw:363-369 (eval = FALSE)
###################################################
## TSS_Scores <- qualityScores_LM(data = Meta_Result$TSS, 
##     tag = "TSS",
##     savePlotPath = filepath)
## TES_Scores <- qualityScores_LM(data = Meta_Result$TES, 
##     tag = "TES",
##     savePlotPath = filepath)


###################################################
### code chunk number 18: ChIC_Vignette.Rnw:372-376
###################################################
TSS_Scores=qualityScores_LM(data=Meta_Result$TSS, 
    tag="TSS")
TES_Scores=qualityScores_LM(data=Meta_Result$TES, 
    tag="TES")


###################################################
### code chunk number 19: ChIC_Vignette.Rnw:382-385 (eval = FALSE)
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
### code chunk number 22: ChIC_Vignette.Rnw:439-448 (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     target = "H3K4me3", 
##     tag = "geneBody", 
##     savePlotPath = filepath)
## 
## metagenePlotsForComparison(data = Meta_Result$TSS,
##     target = "H3K4me3", 
##     tag = "TSS", 
##     savePlotPath = filepath)


###################################################
### code chunk number 23: ComparisonInput (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     target = "H3K4me3", tag = "geneBody")


###################################################
### code chunk number 24: ComparisonInput
###################################################
metagenePlotsForComparison(data = Meta_Result$geneBody,
    target = "H3K4me3", tag = "geneBody")


###################################################
### code chunk number 25: ChIC_Vignette.Rnw:489-493 (eval = FALSE)
###################################################
## plotReferenceDistribution(target = "H3K4me3", 
##     metricToBePlotted = "RSC", 
##     currentValue = crossvalues_Chip$CC_RSC, 
##     savePlotPath = filepath)


###################################################
### code chunk number 26: ReferenceDistr (eval = FALSE)
###################################################
## plotReferenceDistribution(target = "H3K4me3", 
##     metricToBePlotted = "RSC", currentValue = crossvalues_Chip$CC_RSC )


###################################################
### code chunk number 27: ReferenceDistr
###################################################
plotReferenceDistribution(target = "H3K4me3", 
    metricToBePlotted = "RSC", currentValue = crossvalues_Chip$CC_RSC )


###################################################
### code chunk number 28: ChIC_Vignette.Rnw:519-527
###################################################
EM_scoresNew=NULL

EM_scoresNew$QCscores_ChIP <- crossvalues_Chip
EM_scoresNew$QCscores_binding <- bindingScores
EM_scoresNew$TagDensityInput <- list()
EM_scoresNew$TagDensityChip <-list()

CC_Result <- EM_scoresNew


###################################################
### code chunk number 29: ChIC_Vignette.Rnw:530-537
###################################################
te <- predictionScore(target = "H3K4me3", 
    features_cc = CC_Result,
    features_global = Ch_Results,
    features_TSS = TSS_Scores, 
    features_TES = TES_Scores, 
    features_scaled = geneBody_Scores)
print(te)


