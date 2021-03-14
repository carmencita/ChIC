
#'@title Shows available chromatin marks and factors
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


#'@title Lists the IDs of samples included in the compendium
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

    EIDs <- rownames( ChIC.data::compendium_db )[
        grep( "ENC", rownames( ChIC.data::compendium_db ))]

    if ( dataset == "ENCODE" ) {
        message( "ENCODE IDs: " )
        c( EIDs, rownames(ChIC.data::compendium_db_tf) )

    } else if ( dataset == "Roadmap" ) { 

        message( "Roadmap IDs: " )
        rownames( ChIC.data::compendium_db ) [
            ! rownames (ChIC.data::compendium_db) %in% EIDs ]

    } else {
        stop( "No match found. Please use the keywords: 
            \"ENCODE\" or \"Roadmap\"" )
    }
}


