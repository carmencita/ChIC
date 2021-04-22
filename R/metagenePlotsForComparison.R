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
    if (!(is.list(data) & (length(data) == 3L))) 
        stop("Invalid format for data")
    # stopifnot(target %in% Hlist)
    if ((!(target %in% f_metaGeneDefinition("Hlist"))) &  
        (!(target %in% f_metaGeneDefinition("TFlist"))))
        stop("Chromatin mark or TF not valid. Check manual for valid options.")
    if (!(tag %in% c("geneBody", "TES", "TSS"))) 
        stop("tag not valid! Please use: geneBody, TES or TSS")
    ########## 
    
    message("Calculating metagene profiles...")
    iframe <- log2(do.call(rbind, data[[tag]]$input) + psc)
    cframe <- log2(do.call(rbind, data[[tag]]$chip) + psc)
    norm <- do.call(rbind, data[[tag]]$norm)

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
    
    #nframe <- colMeans(t(t(cframe) - t(iframe)), na.rm = TRUE)
    ##average columns
    nframeB <- colMeans(norm, na.rm = TRUE)
    iframeB <- colMeans(iframe, na.rm = TRUE)
    cframeB <- colMeans(cframe, na.rm = TRUE)
    
    #format data frame
    iframeC= data.frame("x"=as.numeric(names(iframeB)),"mean"=iframeB)
    cframeC= data.frame("x"=as.numeric(names(cframeB)),"mean"=cframeB)
    nframeC= data.frame("x"=as.numeric(names(nframeB)),"mean"=nframeB)

    ## get max and min for same y-axis values for chip and input
    newMin <- min(cframeC$mean, absoluteMin, iframeC$mean)
    newMax <- max(cframeC$mean, absoluteMax, iframeC$mean)
    
    ## get max and min for y-axis values for norm
    normMin <- min(nframeC$mean, normMin)
    normMax <- max(nframeC$mean, normMax)
    
    # create comparison plots
    message ("Creating comparison plots...")

    f_plotProfiles(c_mean, cframeC, tag, c(newMin - (newMin/10), newMax + (newMax/10)), 
        maintitel = paste(target, tag, "Chip", sep = "_"), 
        savePlotPath = savePlotPath)


    f_plotProfiles(i_mean, iframeC, tag, c(newMin - (newMin/10), newMax + (newMax/10)), 
        maintitel = paste(target, tag, "Input", sep = "_"), 
        savePlotPath = savePlotPath)


    f_plotProfiles(n_mean, nframeC, tag, c(normMin - 0.001, normMax + 0.001), 
        maintitel = paste(target, tag, "norm", sep = "_"), 
        ylab = "mean log2 enrichment (signal/input)", 
        savePlotPath = savePlotPath)
    
    if ( !is.null(savePlotPath) )
    {
        message("Plots have been saved under ", savePlotPath)
    }
}

