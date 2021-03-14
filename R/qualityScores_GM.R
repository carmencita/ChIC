

###############################################################################
#####                                                      ####################
##### FUNCTIONS QC-metrics for global read distribution    ####################
#####                                                      ####################
###############################################################################

#'@title Wrapper function to calculate GM metrics from global read distribution
#'
#' @description
#' This set of values is based on the global read distribution along the genome 
#' for immunoprecipitation and input data (Diaz et al., 2012). The genome is 
#' binned and the read coverage counted for each bin. Then the function 
#' computes the cumulative distribution of reads density per genomic bin and 
#' plots the fraction of the coverage on the y-axis and the fraction of bins 
#' on the x-axis. Then different values can be sampled from the cumulative 
#' distribution: like the fraction of bins without reads for in 
#' immunoprecipitation and input,the point of the maximum distance between the 
#' ChIP and the input (x-axis, y-axis for immunoprecipitation and input, 
#' distance (as absolute difference), the sign of the differences), the 
#' fraction of reads in the top 1 percent bin for immunoprecipitation and 
#' input. Finally, the funciton returns 9 QC-measures.
#'

#' qualityScores_GM
#'
#' @param selectedTagsChip Data-structure with selected tag information for ChIP 
#' (returned by qualityScores_EM). 
#' @param selectedTagsInput Data-structure with selected tag information for Input 
#' (returned by qualityScores_EM)
#' @param tag.shift, Integer containing the value of the tag shif, calculated 
#' by getCrossCorrelationScores()
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param savePlotPath if set the plot will be saved under "savePlotPath". 
#' Default=NULL and plot will be forwarded to stdout.
#' @param mc
#' @return finalList List with 9 QC-values
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
#' Ch_Results <- qualityScores_GM(densityChip = smoothedChip,
#'     densityInput = smoothedInput, savePlotPath = filepath)
#'}

qualityScores_GM <- function(selectedTagsChip, selectedTagsInput, tag.shift,
    annotationID="hg19", savePlotPath = NULL, mc=1) 
{
    message("***Calculating GM...***")

    pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 5, 
        clear = FALSE, width = 60)
    pb$tick()

    ########## check if input format is ok
    if (!is.list(selectedTagsChip)) 
        stop("Invalid format for selectedChip")
    if (!is.list(selectedTagsInput)) 
        stop("Invalid format for selectedInput")

    if ( !is.numeric(tag.shift) ) 
        stop( "tag.shift must be numeric" )
    if ( tag.shift < 1 ) 
        stop( "tag.shift must be > 0" )

    ########## 
    
    ## objects of smoothed tag density for ChIP and Input
    densityChip <- tagDensity(selectedTagsChip, tag.shift, 
        annotationID = annotationID, mc = mc)

    densityInput <- tagDensity(selectedTagsInput, tag.shift, 
        annotationID = annotationID, mc = mc)
    
    ## shorten frame and reduce resolution
    message("shorten frame...")
    chip.smoothed.density <- f_shortenFrame(densityChip)
    input.smoothed.density <- f_shortenFrame(densityInput)
    
    chip <- unlist(lapply(chip.smoothed.density, 
        FUN = function(list_element) {
        return(list_element$y)
    }))

    input <- unlist(lapply(input.smoothed.density, 
        FUN = function(list_element) {
        return(list_element$y)
    }))
    pb$tick()
    
    ## create cumulative function for chip and input
    cumSumChip <- f_sortAndBinning(chip)
    cumSumInput <- f_sortAndBinning(input)
    
    ## caluclate QCvalues for chip
    message("\nCalculate GM for ChIP ...")
    chipFracTopPercent <- (1-cumSumChip[(which(cumSumChip$x >= 0.99)[1]),"pj"])
    chipFracOfBinsWithoutReads <- cumSumChip$x[(which(cumSumChip$pj > 0)[1])]
    pb$tick()

    ## caluclate QCvalues for input
    message("\nCalculate GM for Input ...")
    inputFracTopPercent <- (1-
        cumSumInput[(which(cumSumInput$x >= 0.99)[1]),"pj"])
    inputFracWithoutReads <- cumSumInput$x[(which(cumSumInput$pj > 0)[1])]
    
    ## get the point of maximum distance between the two functions
    ##  and the sign of the distance
    sign_sign <- sign(max(cumSumInput$pj - cumSumChip$pj))  ###input-chip
    arrowx <- cumSumChip[which.max(abs(cumSumInput$pj - cumSumChip$pj)), ]$x
    pb$tick()

    ## get the distance between the curves
    yIn <- round(
        cumSumInput[which.max(abs(cumSumInput$pj - cumSumChip$pj)), ]$pj, 3)
    yCh <- round(
        cumSumChip[which.max(abs(cumSumInput$pj - cumSumChip$pj)), ]$pj, 3)
    dist <- abs(yIn - yCh)
    ## prepare list to be returned
    finalList <- list(Ch_X.axis = round(arrowx, 3), 
        Ch_Y.Input = yIn, Ch_Y.Chip = yCh, 
        Ch_sign_chipVSinput = sign_sign, 
        Ch_FractionReadsTopbins_chip = round(chipFracTopPercent, 3), 
        Ch_FractionReadsTopbins_input = round(inputFracTopPercent, 3), 
        Ch_Fractions_without_reads_chip = round(chipFracOfBinsWithoutReads, 3),
        Ch_Fractions_without_reads_input = round(inputFracWithoutReads, 3), 
        Ch_DistanceInputChip = dist)

    ## create Fingerprint plot
    f_fingerPrintPlot(cumSumChip, cumSumInput, savePlotPath = savePlotPath)
    if (!is.null(savePlotPath))    
        message("pdf saved under ", savePlotPath)

    
    ## return QC values
    pb$tick()

    message("Calculation of GM done!")
    return(finalList)
}


