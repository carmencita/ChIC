#'@title Function to create metage plots for comparison
#'
#' @description
#' QC-metrics of newly analysed ChIP-seq samples can be compared with 
#' the reference values of the compendium and enrichment profiles 
#' can be plotted against pre-computed profiles of published 
#' datasets. The metagene profiles show the problematic samples signal 
#' (red line) for the ChIP, for the input and their relative 
#' enrichment when compared to the compendium’s mean signal 
#' (black line) and its 2 x standard error (blue shadow). 
#' Additionally the function plots the desired QC-metric
#' as a red dashed line for the sample plotted against the 
#' reference distribution (density plots) of the 
#' compendium values stratified by chromatin marks.
#'
#' metagenePlotsForComparison
#' @param chrommark String, chromatin mark to be analysed. 
#' Has to be one of the following: H3K36me3, H3K4me3, 
#' H3K79me2, H4K20me1,H2AFZ,H3K27me3, H3K9me3,H3K27ac,
#' POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1, H3K4me1,
#' H3K79me1, H3K4ac, H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac,
#' H4K91ac, H4K8ac, H3K18ac, H2BK12ac, H3K56ac, 
#' H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, 
#' H2A.Z, H3K23me2, H2AK9ac, H3T11ph. 
#' For RNAPOL2 different variants are available: POLR2A 
#' (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param data metagene-object of metagene profile by 
#' createMetageneProfile() containing 
#' input and chip profile
#' @param tag indicating the kind of profile to plot. 
#' Can be either: 'geneBody', 'TES' or 'TSS'.
#' @param savePlotPath if set the plot will be saved under 
#' 'savePlotPath'. Default=NULL and plot will be forwarded to stdout. 
#'
#' @return nothing, creates a figure under 'savePlotPath'
#'
#' @export
#'
#' @examples
#' print('Example Code')
#' ## This command is time intensive to run
#' ##To run the example code the user must provide 2 bam files: 
#' ##one for ChIP and one for the input'. Here we used ChIP-seq 
#' ##data from ENCODE. Two example files can be downloaded as follows
#' ## get bam file
#' setwd(tempdir())
#'
#' \dontrun{
#' ##chip data
#' system('wget 
#' https://www.encodeproject.org/files/ENCFF000BLL/@@download/ENCFF000BLL.bam')
#' chipName='ENCFF000BLL'
#' chip.data=readBamFile(chipName)
#'
#' ##input data
#' system('wget 
#' https://www.encodeproject.org/files/ENCFF000BLL/@@download/ENCFF000BKA.bam')
#' inputName='ENCFF000BKA'
#' input.data=readBamFile(inputName)
#'
#' ## calculate binding characteristics 
#' chip_binding.characteristics<-spp::get.binding.characteristics(chip.data, 
#' srange=c(0,500), bin=5,accept.all.tags=TRUE)
#'
#' ## calculate binding characteristics 
#' input_binding.characteristics<-spp::get.binding.characteristics(input.data, 
#' srange=c(0,500), bin=5,accept.all.tags=TRUE)
#'
#' ##get chromosome information and order chip and input by it
#' chrl_final=intersect(names(chip.data$tags),names(input.data$tags))
#' chip.data$tags=chip.datatags[chrl_final]
#' chip.data$quality=chip.data$quality[chrl_final]
#' input.data$tags=input.data$tags[chrl_final]
#' input.data$quality=input.data$quality[chrl_final]
#'
#' ##remove sigular positions with extremely high tag counts with 
#' ##respect to the neighbourhood
#' selectedTags=removeLocalTagAnomalies(chip.data, input.data, 
#' chip_binding.characteristics, input_binding.characteristics)
#' input.dataSelected=selectedTags$input.dataSelected
#' chip.dataSelected=selectedTags$chip.dataSelected
#'
#' ##get smoothed tagdensity 
#' smoothedChip=tagDensity(chip.dataSelected, 
#' tag.shift=82)
#' smoothedInput=tagDensity(input.dataSelected, 
#' tag.shift=82)
#'
#' Meta_Result=createMetageneProfile(smoothed.densityChip=smoothedChip, 
#' smoothedInput,tag.shift=82)
#'
#' metagenePlotsForComparison(data=Meta_Result$geneBody,chrommark='H3K36me3',
#' tag='geneBody',savePlotPath=getwd())
#' metagenePlotsForComparison(data=Meta_Result$TSS,chrommark='H3K36me3',
#' tag='TSS',savePlotPath=getwd())
#'}

metagenePlotsForComparison <- function(data, chrommark, tag, savePlotPath=NULL)
{
    
    Hlist <- f_metaGeneDefinition("Hlist")
    # pseudocount, required to avoid log2 of 0
    psc <- 1
    
    ########## check if input format is ok 
    ##stopifnot(tag %in% c('geneBody','TES','TSS'))
    ########## stopifnot(length(data) == 2L)
    if (!(is.list(data) & (length(data) == 2L))) 
        stop("Invalid format for data")
    # stopifnot(chrommark %in% Hlist)
    if (!(chrommark %in% Hlist)) 
        stop("Chromatin mark not valid. Check manual for valid options.")
    if (!(tag %in% c("geneBody", "TES", "TSS"))) 
        stop("tag not valid! Please use: geneBody, TES or TSS")
    ########## 
    
    iframe <- log2(do.call(rbind, data$input) + psc)
    cframe <- log2(do.call(rbind, data$chip) + psc)
    
    # load average dataframe normalized
    n_mean <- f_loadDataCompendium(endung = "norm", 
        chrommark = chrommark, tag = tag)
    normMin <- min(n_mean$mean - n_mean$sderr)
    normMax <- max(n_mean$mean + n_mean$sderr)
    ## load average dataframe chip
    c_mean <- f_loadDataCompendium(endung = "chip", 
        chrommark = chrommark, tag = tag)
    absoluteMin <- min(c_mean$mean - c_mean$sderr)
    absoluteMax <- max(c_mean$mean + c_mean$sderr)
    ## load average dataframe input
    i_mean <- f_loadDataCompendium(endung = "input", 
        chrommark = chrommark, tag = tag)
    ## get the range for x-y axis fro the final plot
    if ((min(i_mean$mean - i_mean$sderr)) < absoluteMin) {
        absoluteMin <- min(i_mean$mean - i_mean$sderr)
    }
    if ((max(i_mean$mean + i_mean$sderr)) > absoluteMax) {
        absoluteMax <- max(i_mean$mean + i_mean$sderr)
    }
    
    nframe <- colMeans(t(t(cframe) - t(iframe)), na.rm = TRUE)
    iframe <- colMeans(iframe, na.rm = TRUE)
    cframe <- colMeans(cframe, na.rm = TRUE)
    
    iframe <- f_prepareData(c_mean, iframe)
    cframe <- f_prepareData(c_mean, cframe)
    nframe <- f_prepareData(c_mean, nframe)
    
    ## get max and min for same y-axis values for chip and input
    newMin <- min(cframe$mean, absoluteMin, iframe$mean)
    newMax <- max(cframe$mean, absoluteMax, iframe$mean)
    
    ## get max and min for y-axis values for input
    normMin <- min(nframe$mean, normMin)
    normMax <- max(nframe$mean, normMax)
    
    # create comparison plots
    f_plotProfiles(i_mean, iframe, tag, c(newMin - 0.001, newMax + 0.001), 
        maintitel = paste(chrommark, tag, "Input", sep = "_"), 
        savePlotPath = savePlotPath)
    
    f_plotProfiles(c_mean, cframe, tag, c(newMin - 0.001, newMax + 0.001), 
        maintitel = paste(chrommark, tag, "Chip", sep = "_"), 
        savePlotPath = savePlotPath)
    
    f_plotProfiles(n_mean, nframe, tag, c(normMin - 0.001, normMax + 0.001), 
        maintitel = paste(chrommark, tag, "norm", sep = "_"), 
        ylab = "mean log2 enrichment (signal/input)", 
        savePlotPath = savePlotPath)
}


#'@title Function to create reference distribution plot for comparison
#'
#' @description
#' Creates a density plot (in pdf) for the sample against the reference 
#' distribution (density plots) of the compendium values stratified by 
#' chromatin marks.
#'
#' plotReferenceDistribution
#'
#' @param chrommark String, chromatin mark to be analysed. Has to be one 
#' of the following: H3K36me3, H3K4me3, H3K79me2, H4K20me1,H2AFZ,H3K27me3, 
#' H3K9me3,H3K27ac, 
#' POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1, H3K4me1, H3K79me1, H3K4ac, 
#' H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac, H4K91ac, H4K8ac, H3K18ac, 
#' H2BK12ac, H3K56ac, H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, H2A.Z, 
#' H3K23me2,
#' H2AK9ac, H3T11ph. For RNAPOL2 different variants are available: POLR2A 
#' (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param metricToBePlotted The metric to be plotted (Default='RSC')
#' @param currentValue The value of the current sample
#' @param savePlotPath if set the plot will be saved under 
#' 'savePlotPath'. Default=NULL and plot will be forwarded to stdout. 
#'
#' @export
#'
#' @return nothing, creates a figure under 'savePlotPath'
#'@examples
#' print ('Plot distribution of RSC')
#' plotReferenceDistribution(chrommark='H3K4me1',metricToBePlotted='RSC',
#' currentValue=0.8)
#'

plotReferenceDistribution <- function(chrommark, metricToBePlotted = "RSC", 
    currentValue, savePlotPath = NULL) {
    Hlist <- f_metaGeneDefinition("Hlist")
    
    ########## check if input format is ok stopifnot(chrommark %in% Hlist)
    if (!(chrommark %in% Hlist)) 
        stop("Chromatin mark not valid! (Check manual for valid options)")
    if (!is.numeric(currentValue)) 
        stop("currentValue is not numeric!")
    
    ########## 
    
    allChrom <- f_metaGeneDefinition("Classes")
    ## reading compendium compendium_db=NULL
    data("compendium_db", package = "ChIC.data", envir = environment())
    # compendium_db=ChIC.data::compendium_db select the class for the 
    # respective chromatin mark load('data/compendium_db.rda')
    profileInfo <- f_getBindingClass(chrommark)
    
    message(profileInfo$tag)
    # get values for chrommark alias=paste('CC',metricToBePlotted,sep='_')
    alias <- NULL
    if (paste("CC", metricToBePlotted, sep = "_") %in% colnames(compendium_db))
    {
        alias <- paste("CC", metricToBePlotted, sep = "_")
    } else {
        helpi <- paste("Ch", metricToBePlotted, sep = "_")
        if (helpi %in% colnames(compendium_db)) {
            alias <- paste("Ch", metricToBePlotted, sep = "_")
        }
    }
    ## get the values of respective set
    subset <- compendium_db[
        which(compendium_db$CC_TF %in% profileInfo$profileSet), alias]
    ## plot distribution
    f_plotValueDistribution(subset, 
        title = paste(metricToBePlotted, "\n", 
            chrommark, profileInfo$tag, set = " "),
        currentValue, savePlotPath)
}


#'@title Predict score
#'
#' @description
#'  
#' predictionScore
#'
#' @param chrommark String, chromatin mark to be analysed. Has to be 
#' one of the following: H3K36me3, H3K4me3, H3K79me2, H4K20me1,H2AFZ,
#' H3K27me3, H3K9me3,H3K27ac, POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1,
#' H3K4me1, H3K79me1, H3K4ac, 
#' H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac, H4K91ac, H4K8ac, H3K18ac, 
#' H2BK12ac, H3K56ac, H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, 
#' H2A.Z, H3K23me2,
#' H2AK9ac, H3T11ph. For RNAPOL2 different variants are available: 
#' POLR2A (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param features_cc list, with QC-metrics returned from 
#' qualityScores_EM()
#' @param features_global list, list with QC-metrics returned from 
#' qualityScores_GM()
#' @param features_TSS list, list with QC-metrics returned from 
#' qualityScores_LM() with option TSS
#' @param features_TES list,  list with QC-metrics returned from 
#' qualityScores_LM() with option TES
#' @param features_scaled list, list with QC-metrics returned from 
#' qualityScores_LMgenebody()
#'
#' @export
#'
#' @return something something
#'
#' @examples
#' ## This command is time intensive to run
#' ##To run the example code the user must provide 2 bam files: 
#' ##one for ChIP and one for the input'. Here we used ChIP-seq 
#' ##data from ENCODE. Two example files can be downloaded as follows
#' ## get bam file
#' setwd(tempdir())
#'
#' chipName='ENCFF000BLL'
#' \dontrun{
#' ##chip data
#' system('wget 
#' https://www.encodeproject.org/files/ENCFF000BLL/@@download/ENCFF000BLL.bam')
#' chip.data=readBamFile(chipName)
#'
#' ##input data
#' system('wget 
#' https://www.encodeproject.org/files/ENCFF000BLL/@@download/ENCFF000BKA.bam')
#' inputName='ENCFF000BKA'
#' input.data=readBamFile(inputName)
#'
#' ##caluclate first set of QC-metrics
#' CC_Result=qualityScores_EM(chipName=chipName, inputName=inputName, 
#' read_length=36, mc=5, annotationID='hg19')
#'
#' ##save the tagshift as it is needed later
#' tag.shift=CC_Result$QCscores_ChIP$tag.shift
#'
#' ##smoothed tagdensities used as input for the next steps
#' smoothedDensityInput=CC_Result$TagDensityInput
#' smoothedChip=CC_Result$TagDensityChip
#'
#' ###GLOBAL features#########
#' ##caluclate second set of QC-metrics
#' Ch_Results=qualityScores_GM(densityChip=smoothedChip,
#' densityInput=smoothedInput)
#'
#' ###LOCAL features########
#' ##caluclate third set of QC-metrics
#' Meta_Result=createMetageneProfile(smoothed.densityChip=smoothedChip, 
#' smoothed.densityInput=smoothedInput,tag.shift=tag.shift,annotationID='hg19')
#'
#' ##create plots and get values
#' TSS_scores=qualityScores_LM(Meta_Result$TSS,tag='TSS')
#' TES_scores=qualityScores_LM(Meta_Result$TES,tag='TES')
#' geneBody_scores=qualityScores_LMgenebody(Meta_Result$geneBody)
#'
#' te=predictionScore(chrommark='H3K36me3', features_cc=CC_Result,
#' features_global=Ch_Results,features_TSS=TSS_scores, features_TES=TES_scores,
#' features_scaled=geneBody_scores)
#'}


predictionScore <- function(chrommark, features_cc, features_global, 
    features_TSS, features_TES, features_scaled) 
{
    Hlist <- f_metaGeneDefinition("Hlist")
    # stopifnot((chrommark %in% Hlist)&)
    
    ########## check if input format is ok
    if (!(chrommark %in% Hlist)) 
        stop("Chromatin mark not valid. Check manual for valid options.")
    if (!(is.list(features_cc) & (length(features_cc) == 5L))) 
        stop("Invalid format for features_cc")
    if (!(is.list(features_global) & (length(features_global) == 9L))) 
        stop("Invalid format for features_global")
    if (!(is.list(features_TSS) & (length(features_TSS) == 3L))) 
        stop("Invalid format for features_TSS")
    if (!(is.list(features_TES) & (length(features_TES) == 3L))) 
        stop("Invalid format for features_TES")
    if (!(is.list(features_scaled) & (length(features_scaled) == 3L))) 
        stop("Invalid format for features_scaled")
    ########## 
    
    message("load prediction models...")
    
    ## loaad prediction model
    pmodel <- f_getPredictionModel(chrommark = chrommark, histList = Hlist)
    ## get features
    selectedFeat <- c(colnames(pmodel$trainingData))
    selectedFeat <- selectedFeat[-(which(selectedFeat == ".outcome"))]
    
    TSS <- f_convertframe(features_TSS)
    TES <- f_convertframe(features_TES)
    geneBody <- f_convertframe(features_scaled)
    
    p1 <- c(features_cc$QCscores_ChIP, 
        features_cc$QCscores_binding, 
        features_global)
    fframe <- data.frame(matrix(p1), nrow = 35, byrow = TRUE)
    rownames(fframe) <- names(p1)
    
    fframe$nrow <- NULL
    fframe$byrow <- NULL
    colnames(fframe) <- "values"
    
    helper <- rbind(TSS, TES, geneBody, fframe)
    
    a <- lapply(as.list(rownames(helper)), function(element) {
        word <- strsplit(element, "-")[[1]]
        if (length(word) > 1) {
            new <- paste(word[1], word[2], sep = ".")
            word <- new
        }
        
        if (length(grep("%", word)) > 0) {
            word <- strsplit(word, "%")[[1]]
            new <- paste(word, ".", sep = "")
            word <- new
        }
        
        ## Remove this part once featuresname 'twopoint' is 
        ## substitued by 'geneBody' in
        ## the prediction modeds
        if (length(grep("geneBody", element)) > 0) {
            word <- strsplit(element, "geneBody")[[1]]
            new <- paste(word[1], "twopoint", word[2], sep = "")
            word <- new
            if (length(grep("norm_localMax", word)) > 0) {
                new <- strsplit(word, "twopoint")[[1]]
                word <- paste(new[1], "twopoints", new[2], sep = "")
            }
        }
        return(word)
    })
    rownames(helper) <- a
    
    selectedFeat[!(selectedFeat %in% rownames(helper))]
    fVector <- data.frame(values = unlist(
        helper[rownames(helper) %in% selectedFeat, ]))
    rownames(fVector) <- rownames(helper)[rownames(helper) %in% selectedFeat]
    fVector <- data.frame(t(fVector))
    fVector$Class <- as.factor("P")
    
    prediction <- predict(pmodel, newdata = fVector, type = "prob")
    return(prediction$P)
}
