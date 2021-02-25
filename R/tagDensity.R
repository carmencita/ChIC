
#'@title Smoothed tag density
#'
#' @description
#' Calculates the smoothed tag density using spp::get.smoothed.tag.density
#' 
#' tagDensity
#'
#' @param data, data-structure with tag information (see readBamFile())
#' @param tag.shift, Integer containing the value of the tag shif, calculated 
#' by getCrossCorrelationScores()
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#'
#' @return smoothed.density A list of lists, each list corresponds to a 
#' chromosome and contains a vector of coordinates of the 5' ends of the 
#' aligned tags
#'
#' @export
#'
#' @examples
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
#' chipBamSelected=selectedTags$chip.dataSelected
#'
#' smoothedChip=tagDensity(chipBamSelected, tag.shift=finalTagShift)
#' }

tagDensity <- function(data, tag.shift, annotationID = "hg19", mc = 1) {
    ########## check if input format is ok
    if ( !is.list(data) ) 
        stop( "Invalid format for data" )
    
    if ( !is.numeric(tag.shift) ) 
        stop( "tag.shift must be numeric" )
    if ( tag.shift < 1 ) 
        stop( "tag.shift must be > 0" )
    
    annotationID <- f_annotationCheck( annotationID )
    
    if ( !is.numeric(mc) ) {
        warning( "mc must be numeric" )
        mc <- 1
    }
    
    if ( mc < 1 ) {
        warning( "mc set to 1" )
        mc <- 1
    }
    ########## 
    message( "load chrom_info" )
    chromInfo <- f_chromInfoLoad( annotationID )
    data <- f_clearChromStructure( data,annotationID )
    smoothed.density <- f_tagDensity( data = data, 
        tag.shift = tag.shift, 
        chromDef = chromInfo, 
        mc = mc)
    return( smoothed.density )
}
