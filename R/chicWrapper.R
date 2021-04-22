


############################################################################
####                                                      ##################
#### FUNCTION to produce  all plots in one pdf            ##################
#### Author: Ilario Tagliaferri                           ##################
############################################################################


#'@title ChIC analysis in one command  
#'
#' @description
#' This function creates a single document (in pdf format) containing all 
#' the analysis plots produced by ChIC: the CrossCorrelation profile of the 
#' immunoprecipitation, the fingerprint plot and the different metagene 
#' profiles and returns the ChIC RF score.
#' 
#' chicWrapper
#'
#' @param chipName String, filename and path to the ChIP bam file 
#' (without ".bam" extension)
#' @param inputName String, filename and path to the Input bam file
#' (without ".bam" extension)
#' @param read_length Integer, length of the sequencing reads
#' @param target String, chromatin mark or transcription factor to be 
#' analysed. Using the function "listAvailableElements" with the keywords 
#' "mark" and "TF" shows a list with the available elements.
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param savePlotPath path, needs to be set to save the summary plot 
#' under "savePlotPath" (default=getwd()). In the current version of the code we MUST save the output summary plots in a PDF file.
#' If the end users really needs to suppress the PDF output, they can use the "/dev/null" as output savePlotPath 
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
#' @return ChIC RF score, returns the rnadom forest prediction score of the 
#' predictionScore() function, produces the summary report and saves it 
#' under "savePlotPath".
#'
#' @export
#'
#' @examples
#'
#' ## This command is time intensive to run
#'
#' ## To run this example code the user MUST provide 2 bam files: one for 
#' ## ChIP and one for the input. Here we used ChIP-seq data from ENCODE. 
#' ## Two example files can be downloaded using the following link:
#' ## https://www.encodeproject.org/files/ENCFF000BFX/
#' ## https://www.encodeproject.org/files/ENCFF000BDQ/
#' ## and save them in the working directory (here given in the temporary 
#' ## directory "filepath"
#'
#' mc=5
#' \dontrun{
#' 
#' filepath=tempdir()
#' setwd(filepath)
#' 
#' target= "H3K4me3"
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BDQ/@@download/ENCFF000BDQ.bam")
#' 
#' chipName=file.path(filepath,"ENCFF000BFX")
#' inputName=file.path(filepath,"ENCFF000BDQ")
#' 
#' prediction=chicWrapper(chipName=chipName, inputName=inputName, 
#' target= "H3K4me3", read_length=36, mc=mc , savePlotPath=filepath)
#' print(prediction)
#'}


chicWrapper<-function(chipName, inputName, read_length, 
    savePlotPath=getwd(), target, annotationID="hg19", 
    mc=1, debug=FALSE) {

    ### check input data
    checkTargetForPredictor <-suppressMessages(target %in% c(listAvailableElements("mark"), listAvailableElements("TF") , "TF", "broad", "sharp", "RNAPol2" ))
    if (!checkTargetForPredictor) {
message("
###################
###################
####
####    The ChIC RF score will not be computed as 
####    the target is not in the reference compendium.
####    You can use one of the more generic machine learning models
####    by specifying a different target parameter):
####    possible options are \"TF\", \"broad\", \"sharp\", \"RNAPol2\". 
####
###################
###################
")
    }


    checkTargetForComparisonplots<-(target %in% c( f_metaGeneDefinition("Hlist"), f_metaGeneDefinition("TFlist")))
    if (!checkTargetForComparisonplots) {
message("
###################
###################
####")
message(paste("####    The comparison plots will not be produced
####    for the selected target:", target)
)
message("####    as it is not among the available ones in the reference compendium.
####    See the package vignette for more details.
####
###################
###################
")
    }


    # in the current version of the code we MUST save the output summary plots in a PDF file
    # if the end users really needs to suppress the PDF output, they can use the "/dev/null" as output savePlotPath
    summaryPlotsFilename<-paste(chipName, inputName, "ChIC_report.pdf", sep="_")
    if ((!is.null(savePlotPath)) && is.finite(savePlotPath)) { #is.finite is to avoid NAs errors
        if (dir.exists(savePlotPath)) {
            # filepath is a directory, just add the filename
            savePlotPath <- file.path(savePlotPath,summaryPlotsFilename)
        } else {
            # note that if filepath is a full file path already existing, it will be overwritten
            savePlotPath <- file.path(savePlotPath)
        }
    } else {
        savePlotPath <- file.path(getwd(),summaryPlotsFilename)
    }
     message(paste("Creating summary pdf report in", savePlotPath))
     pdf(savePlotPath, onefile=TRUE)


    ## Read BAM files
    message("Reading BAM files")
        chip.data <- readBamFile(chipName)
        input.data <- readBamFile(inputName)

                
    ## Compute ENCODE Metrics (EM)
    EM_Results=qualityScores_EM(
        chipName=chipName,
        inputName=inputName,
        chip.data=chip.data,
        input.data=input.data,
        read_length=read_length, 
        debug=debug,
        mc=mc,
        annotationID=annotationID,
        savePlotPath=NULL
    )
    
    
    ## Compute Global Enrichment profile Metrics (GM)
    GM_Results=qualityScores_GM(
        selectedTagsChip=EM_Results$SelectedTagsChip,
        selectedTagsInput=EM_Results$SelectedTagsInput,
        tag.shift=EM_Results$QCscores_ChIP$tag.shift,
        savePlotPath=NULL,
        mc=mc,
        returnDensities=TRUE
    )
    
    smoothed.densityChip<-GM_Results$densities$densityChip
    smoothed.densityInput<-GM_Results$densities$densityInput
    GM_Results<-GM_Results[(-1*which(names(GM_Results)=="densities"))]  ## drop densities

    ## Compute Local Enrichment profile Metrics (LM)
    Meta_Results=createMetageneProfile(
        selectedTagsChip=EM_Results$SelectedTagsChip,
        selectedTagsInput=EM_Results$SelectedTagsInput,
        smoothed.densityChip=smoothed.densityChip,
        smoothed.densityInput=smoothed.densityInput,
        tag.shift=EM_Results$QCscores_ChIP$tag.shift,
        annotationID=annotationID,
        debug=debug,
        mc=mc
    )
    
    ##create plots and get values
    TSSProfile=qualityScores_LM(
        data=Meta_Results,
        tag="TSS",
        savePlotPath=NULL
    )
    
    TESProfile=qualityScores_LM(
        data=Meta_Results,
        tag="TES",
        savePlotPath=NULL
    )
    
    geneBody_Plot=qualityScores_LM(
        data=Meta_Results,
        tag="geneBody",
        savePlotPath=NULL
    )


    if (checkTargetForComparisonplots) {
        ##additional plots        
        plotReferenceDistribution(
            target=target,
            metricToBePlotted="RSC",
            currentValue=EM_Results$QCscores_ChIP$CC_RSC,
            savePlotPath=NULL
        )

        metagenePlotsForComparison(
            target=target,
            data=Meta_Results,
            tag="TSS",
            savePlotPath=NULL
        )
        
        metagenePlotsForComparison(
            target=target,
            data=Meta_Results,
            tag="TES",
            savePlotPath=NULL
        )
        metagenePlotsForComparison(
            target=target,
            data=Meta_Results,
            tag="geneBody",
            savePlotPath=NULL
        )
           
    } 

    ## close the summary PDF
    message(paste("Closing output summary in", savePlotPath))
    dev.off()


    # adding support for generic "broad", "sharp", "RNAPol2" models
 
    if (checkTargetForPredictor) {
        message("Calculating the prediction score...")
        
        predictedScore=predictionScore(
            target=target,
            features_cc=EM_Results,
            features_global=GM_Results,
            features_TSS=TSSProfile,
            features_TES=TESProfile,
            features_scaled=geneBody_Plot
        )

        print("prediction")
        print(predictedScore)
        return(predictedScore)
    } else {
        return(NA)
    }
}

