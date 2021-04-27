
############################################################################
####                                                      ##################
#### FUNCTIONS QC-metrics for local ENRICHMENT            ##################
####                                                      ##################
############################################################################


#'@title Wrapper function to create scaled and non-scaled metageneprofiles 
#'
#' @description
#' Metagene plots show the signal enrichment around a region of interest like 
#' the TSS or over a predefined set of genes. The tag density of the 
#' immunoprecipitation is taken over all RefSeg annotated human genes, averaged
#' and log2 transformed. The same is done for the input. The normalized profile
#' is calculated as the signal enrichment (immunoprecipitation over the input).
#' Two objects are created: a non-scaled profile for the TSS  and TES, and a 
#' scaled profile for the entire gene, including the gene body. The non-scaled 
#' profile is constructed around the TSS/TES, with 2KB up- and downstream 
#' regions respectively. 
#' 
#' CreateMetageneProfile
#'
#' @param selectedTagsChip Data-structure with selected tag information for ChIP 
#' (returned by qualityScores_EM). 
#' @param selectedTagsInput Data-structure with selected tag information for Input 
#' (returned by qualityScores_EM)
#' @param smoothed.densityChip Optional, output of tagDensity on the ChIP taglist object.
#' It is just used as part of the ChIC_wrapper workflow to skip
#' the call to tagDensity function as it is already performed in the qualityScores_GM function.
#' @param smoothed.densityInput Optional, output of tagDensity on the Input taglist object.
#' It is just used as part of the ChIC_wrapper workflow to skip
#' the call to tagDensity function as it is already performed in the qualityScores_GM function.
#' @param tag.shift Integer containing the value of the tag shif, calculated by
#' getCrossCorrelationScores()
#' @param annotationID String indicating the genome assembly
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#'
#' @return list with 3 objects: scaled profile ("geneBody"), non-scaled profile
#' for TSS (TSS) and TES (TES). Each object is made of a list containing the 
#' chip and the input profile
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
#' ##calculate metagene profiles
#' Meta_Result <- createMetageneProfile(
#'     smoothed.densityChip = smoothedChip, 
#'     smoothed.densityInput = smoothedInput, 
#'     tag.shift = finalTagShift, mc = mc)
#'}

createMetageneProfile <- function( selectedTagsChip, selectedTagsInput,
    smoothed.densityChip=NULL, smoothed.densityInput=NULL,
    tag.shift, annotationID , debug = FALSE, mc = 1) 
{
    ########## check if input format is ok
    if (!is.list(selectedTagsChip)) 
        stop("Invalid format for smoothed.densityChip")
    if (!is.list(selectedTagsInput)) 
        stop("Invalid format for smoothed.densityInput")
    
    if (!is.numeric(tag.shift)) 
        stop("tag.shift must be numeric")
    if (tag.shift < 1) 
        stop("tag.shift must be > 0")
    
    annotationID=f_annotationCheck(annotationID)
    
    if (!is.numeric(mc)) {
        warning("mc must be numeric")
        mc <- 1
    }
    
    if (mc < 1) {
        warning("mc set to 1")
        mc <- 1
    }
    ########## 
    
    message("***Computing metagene profiles...***")
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 6, 
        clear = FALSE, width = 60)

    annotObject <- f_annotationLoad(annotationID)
        
    annotObjectNew <- data.frame(annotObject@.Data, annotObject@annotation, 
        stringsAsFactors = FALSE)
    annotObjectNew$interval_starts <-as.integer(annotObjectNew$interval_starts)
    annotObjectNew$interval_ends <- as.integer(annotObjectNew$interval_ends)
    annotObjectNew$seq_name <- as.character(annotObjectNew$seq_name)
    annotatedGenesPerChr <- split(annotObjectNew, f = annotObjectNew$seq_name)
    

    ## two.point.scaling create scaled metageneprofile input
    message("***Calculating scaled metageneprofile ...***")
    ## objects of smoothed tag density for ChIP and Input

    if (is.null(smoothed.densityChip)) {
    message("Computing reads density profile for ChIP ...")
    smoothed.densityChip <- tagDensity(selectedTagsChip, tag.shift, 
        annotationID = annotationID, mc = mc)
    }

    if (is.null(smoothed.densityInput)) {
    message("Computing reads density profile for Input control ...")
    smoothed.densityInput <- tagDensity(selectedTagsInput, tag.shift, 
        annotationID = annotationID, mc = mc)
    }

    smoothingBandwidth <- 50
    smoothingStep <- 20 
    normalizedProfile= spp::get.smoothed.enrichment.mle(
        selectedTagsChip, selectedTagsInput,background.density.scaling = T, bandwidth=smoothingBandwidth, 
        step=smoothingStep,tag.shift=tag.shift)

    #smoothed.densityInput <- list(td = smoothed.densityInput)
    
    ## Chip
    #smoothed.densityChip <- list(td = smoothed.densityChip)
    message("\nprocess ChIP")
    pb$tick()

    binned_Chip <- masked_t.get.gene.av.density(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc)
    pb$tick()

    message("process input")
    pb$tick()
    
    binned_Input <- masked_t.get.gene.av.density(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc)
    pb$tick()
        

    message("\nprocess ChIP over Input")
   

    binned_Norm <- masked_t.get.gene.av.density(normalizedProfile, 
        gdl = annotatedGenesPerChr, 
        mc = mc)

    geneBody <- list(chip = binned_Chip, input = binned_Input, norm=binned_Norm)
    
    ## one.point.scaling create non-scaled metageneprofile for TSS
    message("\nCreating non-scaled metageneprofiles...")
    message("...TSS")
    
    binnedInput_TSS <- masked_getGeneAvDensity_TES_TSS(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TSS")
    binnedChip_TSS <- masked_getGeneAvDensity_TES_TSS(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TSS")
    binnedNorm_TSS <- masked_getGeneAvDensity_TES_TSS(normalizedProfile, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TSS")

    onepointTSS <- list(chip = binnedChip_TSS, input = binnedInput_TSS,
        norm=binnedNorm_TSS)
    pb$tick()

    ## one.point.scaling create non-scaled metageneprofile for TES
    message("...TES")
    binnedInput_TES <- masked_getGeneAvDensity_TES_TSS(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TES")
    binnedChip_TES <- masked_getGeneAvDensity_TES_TSS(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TES")
    binnedNorm_TES <- masked_getGeneAvDensity_TES_TSS(normalizedProfile, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TES")
    onepointTES <- list(chip = binnedChip_TES, input = binnedInput_TES, 
        norm=binnedNorm_TES)
    pb$tick()

    if ( debug ) {
        message("Debuggin mode ON...")
        message("writing metageneprofiles Rdata objects")
        save(binned_Chip, binned_Input, binned_Norm, 
            file = file.path(getwd(), "geneBody.RData"))
        save(binnedChip_TSS, binnedInput_TSS, binnedNorm_TSS, 
            file = file.path(getwd(), "OnePointTSS.RData"))
        save(binnedChip_TES, binnedInput_TES, binnedNorm_TSS,
            file = file.path(getwd(), "OnePointTES.RData"))
    }
    message("Metageneprofile objects created!")
    return(list(geneBody = geneBody, TSS = onepointTSS, TES = onepointTES))
}
