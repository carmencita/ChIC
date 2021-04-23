
#' @title Shows available chromatin marks and factors
#'
#' @description
#' Lists chromatin marks and transcription factors that are available 
#' to be used for the comparison analysis by the fuctions 
#' metagenePlotsForComparison(), plotReferenceDistribution() 
#' and predictionScore().
#' 
#' listAvailableElements
#'
#' @param target String, chromatin mark or transcription factor to be analysed.
#' With the keywords "mark" and "TF" the respective lists with the available 
#' elements are listed.
#'
#' @return List of elements or single string
#'
#' @export
#'
#' @examples
#'
#' listAvailableElements(target="CTCF")
#' listAvailableElements(target="H3K36me3")
#' listAvailableElements(target="TF")
#' listAvailableElements(target="mark")
#'
listAvailableElements <- function( target )
{
    if ( target == "TF" ) {
        message( "The following TF are available:" )
        f_metaGeneDefinition( "TFlist" )
    } else if ( target == "mark" ) { 
        message( "The following chromatin marks are available:" )
        f_metaGeneDefinition( "Hlist" )
    } else if ( target == "broad" ) { 
        message( "The following chromatin marks are available for the category broad:" )
        f_metaGeneDefinition("Classes")$allBroad
    } else if ( target == "sharp" ) { 
        message( "The following chromatin marks are available for the category sharp:" )
        f_metaGeneDefinition("Classes")$allSharp
    } else if ( target == "RNAPol2" ) { 
        message( "The following chromatin marks are available for the category RNAPol2:" )
        f_metaGeneDefinition("Classes")$RNAPol2
    } else if ( target %in% f_metaGeneDefinition("Hlist") ) { 
        message( "Chromatin mark available for comparison analysis" )
    } else if ( target %in% f_metaGeneDefinition( "TFlist" )) { 
        message( "transcription factor available for comparison analysis" )
    } else {
        stop( "No match found." )
    }
}


#' @title Lists the IDs of samples included in the compendium
#'
#' @description
#' Shows the IDs of all analysed ChIP-seq samples 
#' included in the compendium from ENCODE and Roadmap. 
#' 
#' listDatasets
#'
#' @param dataset String, to specify the dataset for which the IDs 
#' have to be returned. Valid keywords are "ENCODE" and "Roadmap".
#'
#' @return "ENCODE" returns a vecor of transcription factor, chromatin mark  
#' and RNAPol2 sample IDs from ENCODE, "Roadmap" returns a vector of 
#' chromatin mark IDs from Roadmap that have been included in the compendium.
#'
#' @export
#'
#' @examples
#'
#' listDatasets(dataset="ENCODE")
#' listDatasets(dataset="Roadmap")
#'

listDatasets <- function( dataset )
{
    ALLIDs_histones<-as.character(ChIC.data::compendium_db$ID)
    EIDs <- ALLIDs_histones[grep("ENC", ALLIDs_histones)]
    ALLIDs_TFs<-as.character(ChIC.data::compendium_db_tf$ID)

    if ( dataset == "ENCODE" ) {
        message( "ENCODE IDs: " )
        return(c( EIDs,  ALLIDs_TFs))

    } else if ( dataset == "Roadmap" ) { 

        message( "Roadmap IDs: " )
        return(ALLIDs_histones[(!(ALLIDs_histones %in% EIDs))])

    } else {
        stop( "No match found. Please use the keywords: 
            \"ENCODE\" or \"Roadmap\"" )
    }
}


#' @title Lists the metrics available in the compendium 
#'
#' @description
#' Lists the metrics available in the compendium that can be used with plotReferenceDistribution()
#' 
#' listMetrics
#'
#' @param category String, to specify the category for which the list of metrics 
#' should to be returned. Valid keywords are "all", "EM", "GM", "LM".
#'
#' @return returns a character vecor of metrics available in the
#' compendium that can be used with plotReferenceDistribution() as "metricToBePlotted" parameter
#'
#' @export
#'
#' @examples
#'
#' listMetrics(category="all")
#' listMetrics(category="EM")
#'

listMetrics <- function( category = "all" )
{

    if (!(category %in% c("all","EM","GM","LM"))) {
        stop("Invalid value for category. Possible values are \"all\",\"EM\",\"GM\",\"LM\"")
    }

    data("compendium_db_tf", package = "ChIC.data", envir = environment())
    compendium_db_tf<-get("compendium_db_tf") #just to solve the warning on no visible binding for variable loaded from data pacakge

    data("compendium_db", package = "ChIC.data", envir = environment())
    compendium_db<-get("compendium_db") #just to solve the warning on no visible binding for variable loaded from data pacakge
    
    if (!all(names(compendium_db_tf) == names(compendium_db))) {
        stop("compendium_db_tf and compendium_dd have different featrue names")
    }

    finalCompendium <- compendium_db_tf
    allAvailableMetricsForPlot<-colnames(finalCompendium)
    allAvailableMetricsForPlot<-allAvailableMetricsForPlot[(!(allAvailableMetricsForPlot %in% c("ID", "target")))]
    
    if (category == "all") {
        allAvailableMetricsForPlot<-allAvailableMetricsForPlot
    } else if (category == "EM") {
        allAvailableMetricsForPlot<-allAvailableMetricsForPlot[c(which(allAvailableMetricsForPlot %in% c("tag.shift", "N1", "Nd")), grep(pattern="^CC_", perl=TRUE,  x=allAvailableMetricsForPlot))]
        allAvailableMetricsForPlot<-gsub(pattern="^CC_", perl=TRUE,  replacement="", x=allAvailableMetricsForPlot)
    } else if (category == "GM") {
        allAvailableMetricsForPlot<-grep(pattern="^Ch_", perl=TRUE, value=TRUE,  x=allAvailableMetricsForPlot)
        allAvailableMetricsForPlot<-gsub(pattern="^Ch_", perl=TRUE,  replacement="", x=allAvailableMetricsForPlot)
    } else if (category == "LM") {
        allAvailableMetricsForPlot<-allAvailableMetricsForPlot[which(!(allAvailableMetricsForPlot %in% c("tag.shift", "N1", "Nd")))]
        allAvailableMetricsForPlot<-grep(pattern="^CC_", perl=TRUE, value=TRUE, invert=TRUE, x=allAvailableMetricsForPlot)
        allAvailableMetricsForPlot<-grep(pattern="^Ch_", perl=TRUE, value=TRUE, invert=TRUE, x=allAvailableMetricsForPlot)
    } 
    return(allAvailableMetricsForPlot)
}        


