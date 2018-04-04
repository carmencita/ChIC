### R code from vignette source '/lustre/data/FF/Carmen/BitBucket/Bioc_submission/FinalVignette/vignettes/ChIC_Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ChIC_Vignette.Rnw:59-71
###################################################
##load ChIC
library(ChIC)

## set path for working directory
#filepath=tempdir()
#setwd(filepath)
filepath=getwd()
#load tag-list with reads aligned to chromosome 1 
data("chr1Bam", package = "ChIC.data", envir = environment())
chipBam=chr1Bam
data("chr1BamInput", package = "ChIC.data", envir = environment())
inputBam=chr1BamInput


###################################################
### code chunk number 2: ChIC_Vignette.Rnw:82-99 (eval = FALSE)
###################################################
## ##caluclate first set of QC-metrics: EM 
## mc=4
## 
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
### code chunk number 3: ChIC_Vignette.Rnw:125-129 (eval = FALSE)
###################################################
## chipName=file.path(filepath,"ENCFF000BLL")
## inputName=file.path(filepath,"ENCFF000BKA")
## chipBam=readBamFile(chipName)
## inputBam=readBamFile(inputName)


###################################################
### code chunk number 4: ChIC_Vignette.Rnw:147-150
###################################################
start <- proc.time()
mc=20
finalTagShift=85


###################################################
### code chunk number 5: ChIC_Vignette.Rnw:154-168
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
### code chunk number 6: ChIC_Vignette.Rnw:171-176 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics
## crossvalues_Chip<-getCrossCorrelationScores( chipBam , 
##     chip_binding.characteristics, read_length = 36, 
##     savePlotPath = filepath, mc = mc)
## finalTagShift <- crossvalues_Chip$tag.shift


###################################################
### code chunk number 7: ChIC_Vignette.Rnw:182-186 (eval = FALSE)
###################################################
## ## calculate cross correlation QC-metrics for input
## crossvalues_input <- getCrossCorrelationScores(inputBam, 
##     chip_binding.characteristics, read_length = 36, 
##     savePlotPath = filepath, mc = mc)


###################################################
### code chunk number 8: ChIC_Vignette.Rnw:211-225
###################################################
##get chromosome information and order chip and input by it
chrl_final <- intersect(names(chipBam$tags), names(inputBam$tags))
chipBam$tags <- chipBam$tags[chrl_final]
chipBam$quality <- chipBam$quality[chrl_final]
inputBam$tags <- inputBam$tags[chrl_final]
inputBam$quality <- inputBam$quality[chrl_final]

##remove sigular positions with extremely high read counts with 
##respect to the neighbourhood
selectedTags <- removeLocalTagAnomalies(chipBam, inputBam, 
    chip_binding.characteristics, input_binding.characteristics)

inputBamSelected <- selectedTags$input.dataSelected
chipBamSelected <- selectedTags$chip.dataSelected


###################################################
### code chunk number 9: ChIC_Vignette.Rnw:233-238 (eval = FALSE)
###################################################
## ##Finally run function
## bindingScores <- getPeakCallingScores(chip = chipBam, 
##     input = inputBam, chip.dataSelected = chipBamSelected, 
##     input.dataSelected = inputBamSelected, 
##     tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 10: ChIC_Vignette.Rnw:242-247
###################################################
##Finally run function
bindingScores <- getPeakCallingScores(chip = chipBam, 
    input = inputBam, chip.dataSelected = chipBamSelected, 
    input.dataSelected = inputBamSelected, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 11: ChIC_Vignette.Rnw:259-263
###################################################
smoothedChip <- tagDensity(chipBamSelected, 
    tag.shift = finalTagShift, mc = mc)
smoothedInput <- tagDensity(inputBamSelected, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 12: ChIC_Vignette.Rnw:279-281 (eval = FALSE)
###################################################
## Ch_Results <- qualityScores_GM(densityChip = smoothedChip,
##     densityInput = smoothedInput, savePlotPath = filepath)


###################################################
### code chunk number 13: Fingerprint (eval = FALSE)
###################################################
## Ch_Results=qualityScores_GM(densityChip=smoothedChip,
## densityInput=smoothedInput)


###################################################
### code chunk number 14: Fingerprint
###################################################
Ch_Results=qualityScores_GM(densityChip=smoothedChip,
densityInput=smoothedInput)


###################################################
### code chunk number 15: ChIC_Vignette.Rnw:322-326
###################################################
Meta_Result <- createMetageneProfile(
    smoothed.densityChip = smoothedChip, 
    smoothed.densityInput = smoothedInput, 
    tag.shift = finalTagShift, mc = mc)


###################################################
### code chunk number 16: ChIC_Vignette.Rnw:332-336 (eval = FALSE)
###################################################
## TSS_Scores <- qualityScores_LM(data = Meta_Result$TSS, tag = "TSS",
##     savePlotPath = filepath)
## TES_Scores <- qualityScores_LM(data = Meta_Result$TES, tag = "TES",
##     savePlotPath = filepath)


###################################################
### code chunk number 17: ChIC_Vignette.Rnw:339-341
###################################################
TSS_Scores=qualityScores_LM(data=Meta_Result$TSS, tag="TSS")
TES_Scores=qualityScores_LM(data=Meta_Result$TES, tag="TES")


###################################################
### code chunk number 18: ChIC_Vignette.Rnw:347-350 (eval = FALSE)
###################################################
## #create scaled metagene profile
## geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody,
##     savePlotPath = filepath)


###################################################
### code chunk number 19: geneBody (eval = FALSE)
###################################################
## #create scaled metagene profile
## geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody)


###################################################
### code chunk number 20: geneBody
###################################################
#create scaled metagene profile
geneBody_Scores <- qualityScores_LMgenebody(Meta_Result$geneBody)


###################################################
### code chunk number 21: ChIC_Vignette.Rnw:393-397 (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     chrommark = "H3K36me3", tag = "geneBody", savePlotPath = filepath)
## metagenePlotsForComparison(data = Meta_Result$TSS,
##     chrommark = "H3K36me3", tag = "TSS", savePlotPath = filepath)


###################################################
### code chunk number 22: Comparison (eval = FALSE)
###################################################
## metagenePlotsForComparison(data = Meta_Result$geneBody,
##     chrommark = "H3K36me3", tag = "geneBody")


###################################################
### code chunk number 23: Comparison
###################################################
metagenePlotsForComparison(data = Meta_Result$geneBody,
    chrommark = "H3K36me3", tag = "geneBody")


###################################################
### code chunk number 24: ChIC_Vignette.Rnw:431-433 (eval = FALSE)
###################################################
## plotReferenceDistribution(chrommark = "H3K36me3", 
##     metricToBePlotted = "RSC", currentValue = 0.49, savePlotPath = filepath)


###################################################
### code chunk number 25: ReferenceDistr (eval = FALSE)
###################################################
## plotReferenceDistribution(chrommark = "H3K4me1", 
##     metricToBePlotted = "RSC", currentValue = 0.49)


###################################################
### code chunk number 26: ReferenceDistr
###################################################
plotReferenceDistribution(chrommark = "H3K4me1", 
    metricToBePlotted = "RSC", currentValue = 0.49)


###################################################
### code chunk number 27: ChIC_Vignette.Rnw:456-469
###################################################
data("EM_scores", package = "ChIC.data", envir = environment())

EM_scoresNew=NULL

EM_scoresNew$QCscores_ChIP=EM_scores$QCscores_ChIP
EM_scoresNew$QCscores_binding=bindingScores
EM_scoresNew$TagDensityInput=list()
EM_scoresNew$TagDensityChip=list()

CC_Result=EM_scoresNew
end <- proc.time()
print(start)
print(end)


###################################################
### code chunk number 28: ChIC_Vignette.Rnw:472-479
###################################################
te <- predictionScore(chrommark = "H3K36me3", 
    features_cc = CC_Result,
    features_global = Ch_Results,
    features_TSS = TSS_Scores, 
    features_TES = TES_Scores, 
    features_scaled = geneBody_Scores)
print(te)


