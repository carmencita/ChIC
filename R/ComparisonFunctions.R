#' @title Function to create metage plots for comparison
#'
#' @description
#' QC-metrics of newly analysed ChIP-seq samples can be compared with the 
#' reference values of the compendium and enrichment profiles can be 
#' plotted against pre-computed profiles of published datasets. The metagene
#' profiles show the problematic samples signal (red line) for the ChIP, for
#' the input and their relative enrichment when compared to the compendium’s 
#' mean signal (black line) and its 2 x standard error (blue shadow). 
#' Additionally the function plots the desired QC-metric as a red dashed line 
#' for the sample plotted against the reference distribution (density plots) 
#' of the compendium values stratified by chromatin marks.
#'
#' metagenePlotsForComparison
#' @param target String, chromatin mark or transcription factor to be analysed. 
#' Use listAvailableElements() function to check availability.
#' @param data metagene-object of metagene profile by createMetageneProfile() 
#' containing input and chip profile
#' @param tag indicating the kind of profile to plot. Can be either: 
#' geneBody, TES or TSS.
#' @param savePlotPath if set the plot will be saved under 'savePlotPath'. 
#' Default=NULL and plot will be forwarded to stdout. 
#'
#' @return Creates a pdf figure under 'savePlotPath'
#'
#' @export
#'
#' @examples
#'
#'
#' ## This command is time intensive to run
#' ##
#' ## To run the example code the user must provide two bam files for the ChIP
#' ## and the input and read them with the readBamFile() function. To make it
#' ## easier for the user to run the example code we provide tow bam examples 
#' ## (chip and input) in our ChIC.data package that have already been loaded 
#' ## with the readBamFile() function.
#'
#' mc=4
#' finalTagShift=98
#'
#' \dontrun{
#'
#' filepath=tempdir()
#' setwd(filepath)
#'
#' data("chipSubset", package = "ChIC.data", envir = environment())
#' chipBam=chipSubset
#' data("inputSubset", package = "ChIC.data", envir = environment())
#' inputBam=inputSubset
#'
#' ## calculate binding characteristics 
#'
#' chip_binding.characteristics<-spp::get.binding.characteristics(
#'     chipBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
#' input_binding.characteristics<-spp::get.binding.characteristics(
#'     inputBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
#'
#' ##get chromosome information and order chip and input by it
#' chrl_final=intersect(names(chipBam$tags),names(inputBam$tags))
#' chipBam$tags=chipBam$tags[chrl_final]
#' chipBam$quality=chipBam$quality[chrl_final]
#' inputBam$tags=inputBam$tags[chrl_final]
#' inputBam$quality=inputBam$quality[chrl_final]
#' 
#' ##remove sigular positions with extremely high read counts with 
#' ##respect to the neighbourhood
#' selectedTags=removeLocalTagAnomalies(chipBam, inputBam, 
#' chip_binding.characteristics, input_binding.characteristics)
#' 
#' inputBamSelected=selectedTags$input.dataSelected
#' chipBamSelected=selectedTags$chip.dataSelected
#'
#' ##smooth input and chip tags
#' smoothedChip <- tagDensity(chipBamSelected, 
#'     tag.shift = finalTagShift, mc = mc)
#' smoothedInput <- tagDensity(inputBamSelected, 
#'     tag.shift = finalTagShift, mc = mc)
#'
#' ##calculate metagene profiles
#' Meta_Result <- createMetageneProfile(
#'     smoothed.densityChip = smoothedChip, 
#'     smoothed.densityInput = smoothedInput, 
#'     tag.shift = finalTagShift, mc = mc)
#'
#' ##compare metagene features of the geneBody with the compendium
#' metagenePlotsForComparison(data = Meta_Result$geneBody,
#'     target = "H3K4me3", tag = "geneBody")
#'}

metagenePlotsForComparison <- function(data, target, tag, 
    savePlotPath = NULL)
{
    # pseudocount, required to avoid log2 of 0
    psc <- 1
    
    ########## check if input format is ok 
    if (!(is.list(data) & (length(data) == 2L))) 
        stop("Invalid format for data")
    # stopifnot(target %in% Hlist)
    if ((!(target %in% f_metaGeneDefinition("Hlist"))) &  
        (!(target %in% f_metaGeneDefinition("TFlist"))))
        stop("Chromatin mark or TF not valid. Check manual for valid options.")
    if (!(tag %in% c("geneBody", "TES", "TSS"))) 
        stop("tag not valid! Please use: geneBody, TES or TSS")
    ########## 
    
    message("Calculating metagene profiles...")
    iframe <- log2(do.call(rbind, data$input) + psc)
    cframe <- log2(do.call(rbind, data$chip) + psc)
    
    # load average dataframe normalized
    n_mean <- f_loadDataCompendium(endung = "norm", 
        target = target, tag = tag)
    normMin <- min(n_mean$mean - n_mean$sderr)
    normMax <- max(n_mean$mean + n_mean$sderr)
    ## load average dataframe chip
    c_mean <- f_loadDataCompendium(endung = "chip", 
        target = target, tag = tag)
    absoluteMin <- min(c_mean$mean - c_mean$sderr)
    absoluteMax <- max(c_mean$mean + c_mean$sderr)
    ## load average dataframe input
    i_mean <- f_loadDataCompendium(endung = "input", 
        target = target, tag = tag)
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
    message("Creating comparison plots...")
    f_plotProfiles(i_mean, iframe, tag, c(newMin - 0.001, newMax + 0.001), 
        maintitel = paste(target, tag, "Input", sep = "_"), 
        savePlotPath = savePlotPath)
    
    f_plotProfiles(c_mean, cframe, tag, c(newMin - 0.001, newMax + 0.001), 
        maintitel = paste(target, tag, "Chip", sep = "_"), 
        savePlotPath = savePlotPath)
    
    f_plotProfiles(n_mean, nframe, tag, c(normMin - 0.001, normMax + 0.001), 
        maintitel = paste(target, tag, "norm", sep = "_"), 
        ylab = "mean log2 enrichment (signal/input)", 
        savePlotPath = savePlotPath)
    
    if (!is.null(savePlotPath))
    {
        message("Plots have been saved under ",savePlotPath)
    }
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
#' @param target String, chromatin mark or transcription factor to be analysed.
#' Use listAvailableElements() function to check availability.
#' @param metricToBePlotted The metric to be plotted (Default='RSC')
#' @param currentValue The value of the current sample
#' @param savePlotPath if set the plot will be saved under 'savePlotPath'. 
#' Default=NULL and plot will be forwarded to stdout. 
#'
#' @export
#'
#' @return nothing, creates a figure under 'savePlotPath'
#' 
#' @examples
#' print ('Plot distribution of RSC')
#' \dontrun{
#' filepath=tempdir()
#' setwd(filepath)
#'
#' plotReferenceDistribution(target="H3K4me1", 
#'     metricToBePlotted="RSC", currentValue=0.49, savePlotPath=filepath)
#'}

plotReferenceDistribution <- function(target, metricToBePlotted = "RSC", 
    currentValue, savePlotPath = NULL) 
{
    ########## check if input format is ok stopifnot(target %in% Hlist)
    if ((!(target %in% f_metaGeneDefinition("Hlist"))) &  
        (!(target %in% f_metaGeneDefinition("TFlist"))))
        stop("Chromatin mark or transcription factor not valid. 
            Check manual for valid options.")
    if (!is.numeric(currentValue)) 
        stop("currentValue is not numeric!")

    profileInfo=NULL
    ########## 
    if (target %in% f_metaGeneDefinition("TFlist"))
    {
        data("compendium_db_tf", package = "ChIC.data", envir = environment())
        finalCompendium=compendium_db_tf
    }else{
        #allChrom <- f_metaGeneDefinition("Classes")
        ## reading compendium compendium_db=NULL
        data("compendium_db", package = "ChIC.data", envir = environment())
        finalCompendium=compendium_db
        profileInfo <- f_getBindingClass(target) 
        message(profileInfo$tag)
        tag=profileInfo$tag
    }

    # get values for target alias=paste('CC',metricToBePlotted,sep='_')
    alias <- NULL
    if (paste("CC", metricToBePlotted, sep = "_") %in% 
        colnames(finalCompendium))
    {
        alias <- paste("CC", metricToBePlotted, sep = "_")
    } else {
        helpi <- paste("Ch", metricToBePlotted, sep = "_")
        if (helpi %in% colnames(finalCompendium)) {
            alias <- paste("Ch", metricToBePlotted, sep = "_")
        }
    }

    if(target %in% f_metaGeneDefinition("TFlist")){
        subset=finalCompendium[,alias]
    }else{
        ## get the values of respective set
        subset <- finalCompendium[
            which(finalCompendium$CC_TF %in% profileInfo$profileSet), alias]
    }
    ## plot distribution
    f_plotValueDistribution(subset, 
        title = paste(metricToBePlotted, "\n", 
            target, profileInfo$tag, set = " "),
        currentValue, savePlotPath)
    
    if (!is.null(savePlotPath))
    {
        message("Plots have been saved under ",savePlotPath)
    }
}


#'@title Predict score
#'
#' @description
#'  
#' predictionScore
#'
#' @param target String, chromatin mark or transcription factor 
#' to be analysed. Use listAvailableElements() function to check availability.
#' If the specific transcription factor is not available the 
#' keyword "TF" can be used to call the TF-model.
#' @param features_cc list, with QC-metrics returned from qualityScores_EM()
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
#' @return prediction
#'
#' @examples
#'
#' ## To execute this command the user has to run the entire pipeline
#' ## (time intensive to run)
#'
#' ## To run this example code the user MUST provide 2 bam files: one for ChIP 
#' ## and one for the input". Here we used ChIP-seq data from ENCODE. Two 
#' ## example files can be downloaded using the following link:
#' ## https://www.encodeproject.org/files/ENCFF000BFX/
#' ## https://www.encodeproject.org/files/ENCFF000BDQ/
#' ## and save them in the working directory (here given in the temporary 
#' ## directory "filepath"
#'
#' mc=4
#'
#' \dontrun{
#' filepath=tempdir()
#' setwd(filepath)
#'
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BFX/@@@download/ENCFF000BFX.bam")
#' 
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BDQ/@@@download/ENCFF000BDQ.bam")
#'
#' chipName=file.path(filepath,"ENCFF000BFX")
#' inputName=file.path(filepath,"ENCFF000BDQ")
#'
#' CC_Result=qualityScores_EM(chipName=chipName, inputName=inputName, 
#' read_length=36, mc=mc,savePlotPath=filepath)
#' 
#' ##save tag.shift value
#' finalTagShift=CC_Result$QCscores_ChIP$tag.shift
#'
#' ##save the smooted profile
#' smoothedDensityInput=CC_Result$TagDensityInput
#' smoothedDensityChip=CC_Result$TagDensityChip
#'
#' ##caluclate GM QC-metrics
#' Ch_Results=qualityScores_GM(densityChip=smoothedDensityChip,
#' densityInput=smoothedDensityInput,savePlotPath=filepath)
#' 
#' ##caluclate metagene profiles
#' Meta_Result=createMetageneProfile(smoothedDensityChip,smoothedDensityInput,
#' finalTagShift,annotationID="hg19",mc=mc)
#'
#' ##get LM QC-values
#' TSSProfile=qualityScores_LM(Meta_Result$TSS,tag="TSS",savePlotPath=filepath)
#' TESProfile=qualityScores_LM(Meta_Result$TES,tag="TES",savePlotPath=filepath)
#' geneBody_Plot=qualityScores_LMgenebody(Meta_Result$geneBody, 
#' savePlotPath=filepath)
#'
#' ##Finally use all calculated QC-metrics to predict the final score
#' ##example for chromatin mark H3K4me3
#' predictionScore(target="H3K4me3", features_cc=CC_Result,
#' features_global=Ch_Results,features_TSS=TSSProfile, features_TES=TESProfile,
#' features_scaled=geneBody_Plot)
#'
#' ##example for TF
#' predictionScore(target="TF", features_cc=CC_Result,
#' features_global=Ch_Results,features_TSS=TSSProfile, features_TES=TESProfile,
#' features_scaled=geneBody_Plot)
#'}

predictionScore <- function(target, features_cc, features_global, 
    features_TSS, features_TES, features_scaled) 
{

    ########## check if input format is ok
    if ((!(target %in% f_metaGeneDefinition("Hlist"))) &
        (!(target %in% f_metaGeneDefinition("TFlist"))) &
        (!(target == "TF"))) 
            stop("Chromatin mark or transcription factor not valid.
                Check \"listAvailableElements()\" function for valid options.")

    if (!(is.list(features_cc) & (length(features_cc) == 4L)))
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

    message("loading prediction model for ", target)

    ## loaad prediction model
    pmodel <- f_getPredictionModel(id = target)
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
        #if (length(grep("geneBody", element)) > 0) {
        #    word <- strsplit(element, "geneBody")[[1]]
        #    new <- paste(word[1], "twopoint", word[2], sep = "")
        #    word <- new
        #    if (length(grep("norm_localMax", word)) > 0) {
        #        new <- strsplit(word, "twopoint")[[1]]
        #        word <- paste(new[1], "twopoints", new[2], sep = "")
        #    }
        #    if (length(grep("norm_auc", word)) > 0) {
        #        new <- strsplit(word, "twopoint")[[1]]
        #        word <- paste(new[1], "twopoints", new[2], sep = "")
        #    }
        #}

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
    message("Predicted QC score is ", prediction$P )
    return(prediction$P)
}
