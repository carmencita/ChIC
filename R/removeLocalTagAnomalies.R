

#'@title Removes loval anomalies
#'
#' @description
#' The removeLocalTagAnomalies function removes tags from regions with 
#' extremely high tag counts compared to the neighbourhood.
#'
#' removeLocalTagAnomalies
#'
#' @param chip, data-structure with tag information for the ChIP (see 
#' readBamFile())
#' @param input, data-structure with tag information for the Input (see 
#' readBamFile())
#' @param chip_b.characteristics binding.characteristics of the ChIP. Is a 
#' data-structure containing binding information for binding peak separation 
#' distance and cross-correlation profile (see get.binding.characteristics)
#' @param input_b.characteristics, binding.characteristics of the Input. Is a 
#' data-structure containing binding information for binding preak separation 
#' distance and cross-correlation profile (see get.binding.characteristics)
#'
#' @return result A list containing filtered data structure for ChIP and Input
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
#' ##with the readBamFile() function.
#'
#' mc=4
#' \dontrun{
#'
#' filepath=tempdir()
#' setwd(filepath)
#'
#' ##load the data
#' data("chipSubset", package = "ChIC.data", envir = environment())
#' chipBam=chipSubset
#' data("inputSubset", package = "ChIC.data", envir = environment())
#' inputBam=inputSubset
#'
#' ## calculate binding characteristics 
#'
#' chip_binding.characteristics<-spp::get.binding.characteristics(
#' chipBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
#' 
#' input_binding.characteristics<-spp::get.binding.characteristics(
#' inputBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
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
#' }

removeLocalTagAnomalies <- function(chip, input, chip_b.characteristics, 
    input_b.characteristics) 
{
    ########## check if input format is ok
    if (!(is.list(chip) & (length(chip) == 2L))) 
        stop("Invalid format for chip")
    if (!(is.list(chip_b.characteristics) & 
        (length(chip_b.characteristics) == 3L))) 
        stop("Invalid format for chip_b.characteristics")
    
    if (!(is.list(input) & (length(input) == 2L))) 
        stop("Invalid format for input")
    if (!(is.list(input_b.characteristics) & 
        (length(input_b.characteristics) == 3L))) 
        stop("Invalid format for input_b.characteristics")
    ######### 
    
    result <- f_removeLocalTagAnomalies(chip, input, chip_b.characteristics, 
        input_b.characteristics, 
        remove.local.tag.anomalies = TRUE, select.informative.tags = FALSE)
    return(result)
}

