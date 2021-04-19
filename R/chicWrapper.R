


############################################################################
####                                                      ##################
#### FUNCTION to produce  all plots in one pdf            ##################
#### Author: Ilario Tagliaferri                           ##################
############################################################################


#'@title Wrapper function to create a single document (pdf) that summarizes  
#' all plots produced by ChIC  
#'
#' @description
#' This function creates a single document (in pdf format) containing all 
#' the analysis plots produced by ChIC: the CrossCorrelation profile of the 
#' immunoprecipitation, the fingerprint plot and the different metagene 
#' profiles.
#' 
#' chicWrapper
#'
#' @param chipName String, filename and path to the ChIP bam file 
#' (without extension)
#' @param inputName String, filename and path to the Input bam file
#' (without extension)
#' @param read_length Integer, length of the reads
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
#'
#' @return predictedScore, returns the prediction score of the 
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
    mc=1, debug=FALSE)
{

    # in the current version of the code we MUST save the output summary plots in a PDF file
    # if the end users really needs to suppress the PDF output, they can use the "/dev/null" as output savePlotPath
    summaryPlotsFilename<-paste(chipName, inputName, "ChIC_report.pdf", sep="_")
    
    if (!is.null(savePlotPath) && is.finite(savePlotPath)) { #is.finite is to avoid NAs errors
        if (dir.exists(savePlotPath) {
            # filepath is a directory, just add the filename
            savePlotPath <- file.path(savePlotPath,summaryPlotsFilename)
        } else {
            # note that if filepath is a full file path already existing, it will be overwritten
            savePlotPath <- file.path(savePlotPath)
        }
    } else {
        savePlotPath <- file.path(getwd(),summaryPlotsFilename)
    }
     message(paste("Creating summary pdf report in", savePlotPath)
     pdf(savePlotPath, onefile=TRUE)

                
    ## get Encode Metrics
    CC_Result=qualityScores_EM(
        chipName=chipName,
        inputName=inputName,
        read_length=read_length, 
        debug=debug,
        mc=mc,
        annotationID=annotationID,
        savePlotPath=NULL
    )
    
    ##save the tagshift as it is needed later
    
    tag.shift=CC_Result$QCscores_ChIP$tag.shift
    
    ##smoothed tagdensities used as input for the next steps

    ##Old

    ##smoothedDensityInput=CC_Result$TagDensityInput
    ##smoothedDensityChip=CC_Result$TagDensityChip

    ##New

    smoothedDensityInput=CC_Result$SelectedTagsInput
    smoothedDensityChip=CC_Result$SelectedTagsChip

    ##GLOBAL features#########
    ##caluclate second set of QC-metrics
    Ch_Results=qualityScores_GM(
        selectedTagsChip=smoothedDensityChip,
        selectedTagsInput=smoothedDensityInput,
        tag.shift=tag.shift,
        savePlotPath=NULL,
        mc=mc
        
    )
    
    
    ##LOCAL features########
    ##caluclate third set of QC-metrics
    Meta_Result=createMetageneProfile(
        selectedTagsChip=smoothedDensityChip,
        selectedTagsInput=smoothedDensityInput,
        tag.shift=tag.shift,
        annotationID=annotationID,
        debug=debug,
        mc=mc
    )
    
    ##create plots and get values
    TSSProfile=qualityScores_LM(
        Meta_Result$TSS,
        tag="TSS",
        savePlotPath=NULL
    )
    
    TESProfile=qualityScores_LM(
        Meta_Result$TES,
        tag="TES",
        savePlotPath=NULL
    )
    
    geneBody_Plot=qualityScores_LMgenebody(
        Meta_Result$geneBody,
        savePlotPath=NULL
    )


    if(target %in% c( f_metaGeneDefinition("Hlist"), f_metaGeneDefinition("TFlist"))) {
        ##additional plots
        
        metagenePlotsForComparison(
            target=target,
            data=Meta_Result$geneBody,
            tag="geneBody",
            savePlotPath=NULL
        )
        
        metagenePlotsForComparison(
            target=target,
            Meta_Result$TSS,
            tag="TSS",
            savePlotPath=NULL
        )
        
        metagenePlotsForComparison(
            target=target,
            Meta_Result$TES,
            tag="TES",
            savePlotPath=NULL
        )
        
        plotReferenceDistribution(
            target=target,
            metricToBePlotted="RSC",
            currentValue=CC_Result$QCscores_ChIP$CC_RSC,
            savePlotPath=NULL
        )
    }else{
        message(paste("The comparison plots will not be produced for the selected target:", target))
    }

    # adding support for generic "broad", "sharp", "RNAPol2" models
    if (target %in% c(listAvailableElements("mark"), listAvailableElements("TF") , "TF", "broad", "sharp", "RNAPol2" )) {
        message("Calculating the prediction score...")
        
        predictedScore=predictionScore(
            target=target,
            features_cc=CC_Result,
            features_global=Ch_Results,
            features_TSS=TSSProfile,
            features_TES=TESProfile,
            features_scaled=geneBody_Plot
        )

        print("prediction")
        print(predictedScore)

    } else {
        stop( "Histone mark or TF not found. 
            Could not calculate the prediction score 
            using chicWrapper(). You might try the 
            predictionScore() function wihtout the wrapper.
            Alternatively, you can use one of the more generic models (taget parameter):
            possible options are \"TF\", \"broad\", \"sharp\", \"RNAPol2\". ")
    }
    
    dev.off()
    return(predictedScore)
}

