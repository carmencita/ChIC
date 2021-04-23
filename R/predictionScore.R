
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
#' @return predictions for positive and negative class
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
#' ##example for TF not available in compendium
#' predictionScore(target="TF", features_cc=CC_Result,
#' features_global=Ch_Results,features_TSS=TSSProfile, features_TES=TESProfile,
#' features_scaled=geneBody_Plot)
#'
#' ##example for CTCF
#' predictionScore(target="CTCF", features_cc=CC_Result,
#' features_global=Ch_Results,features_TSS=TSSProfile, features_TES=TESProfile,
#' features_scaled=geneBody_Plot)
#'}

predictionScore <- function(target, features_cc, features_global, 
    features_TSS, features_TES, features_scaled) 
{

    ########## check if input format is ok
    # adding support for generic "broad", "sharp", "RNAPol2" models
    if ((!(target %in% c( f_metaGeneDefinition("Hlist"), f_metaGeneDefinition("TFlist"),  "TF", "broad", "sharp", "RNAPol2"  ))))
           stop("Target model selection not valid not valid.
                Check \"listAvailableElements()\" function for valid options.")

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

    message("loading prediction model for ", target)

    ## loaad prediction model
    pmodel <- f_getPredictionModel(id = target)
    ## get features
    selectedFeat <- names(pmodel$forest$xlevels)

    TSS <- f_convertframe(features_TSS)
    TES <- f_convertframe(features_TES)
    geneBody <- f_convertframe(features_scaled)

    p1 <- c(features_cc$QCscores_ChIP, 
        features_cc$QCscores_binding, 
        features_global)

    fframe <- data.frame(matrix(unlist(p1)), nrow = 35, byrow = TRUE)
    rownames(fframe) <- names(p1)
    message("collecting features...")
    fframe$nrow <- NULL
    fframe$byrow <- NULL
    colnames(fframe) <- "values"

    helper <- rbind(TSS, TES, geneBody, fframe)

    a1 <- lapply(rownames(helper), function(element) {
        word=element
        ##delete % , otherwise it is not matching 
        ##the features of the model
        if (length(grep("%", element)) > 0) {
            word <- strsplit(element, "%")[[1]]
        }
        return(word)
    })
    a <- lapply(a1, function(element) {
        word=element
        ##substitue - with ., otherwise it is not matching 
        ##the features of the model
        if (length(grep("-", element)) > 0) { 
            word <- paste(strsplit(element, "-")[[1]][1],
                strsplit(element, "-")[[1]][2],sep=".")
        }
        return(word)
    })
    rownames(helper) = unlist(a)
    #selectedFeat[!(which(selectedFeat %in% rownames(helper)))]
    fVector <- data.frame(values = unlist(
        helper[which(rownames(helper) %in% selectedFeat), ]))
    rownames(fVector) <- rownames(helper)[rownames(helper) %in% selectedFeat]
    fVector <- data.frame(t(fVector))
    fVector$Class <- as.factor("P")
    message("computing prediction score...")
    #getFeat=rownames(helper)[which(rownames(helper) %in% selectedFeat)]
    prediction <- predict(pmodel, newdata = fVector, type = "prob")
    print("Predicted RF score for positive class is ", round(prediction[[1]],4) )
    print("Predicted RF score for negative class is ", round(prediction[[2]],4) )
    return(prediction)
}

