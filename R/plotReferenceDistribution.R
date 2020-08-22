
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
        finalCompendium <- compendium_db_tf
    }else{
        #allChrom <- f_metaGeneDefinition("Classes")
        ## reading compendium compendium_db=NULL
        data("compendium_db", package = "ChIC.data", envir = environment())
        finalCompendium <- compendium_db
        profileInfo <- f_getBindingClass( target ) 
        message( profileInfo$tag )
        tag <- profileInfo$tag
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
        subset <- finalCompendium[ ,alias]
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

