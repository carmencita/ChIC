#' @import ChIC.data
#' @import caret
# @rawNamespace import(girafe, except = c(plot,reduce))f
#' @import spp 
#' @import GenomicRanges
#' @importFrom graphics abline axis legend lines matplot par
#' plot polygon text
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom BiocGenerics strand
#' @importFrom S4Vectors Rle queryHits subjectHits
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom stats density na.omit predict quantile sd var
#' @importFrom utils str write.table data
# @importFrom BiocParallel bplapply
#' @importFrom methods new
#' @importFrom parallel makeCluster stopCluster mclapply
#' @importFrom progress progress_bar


#######################################################################
###############                                         ###############  
############### FUNCTIONS QC-metrics for narrow binding PROFILES ######
###############                                         ###############
#######################################################################



#' @title Wrapper function to calculate EM metrics
#'
#' @description
#' Wrapper that reads bam files and provides EM QC-metrics from 
#' cross-correlation analysis, peak calling and general metrics like 
#' for example the read-length or NRF. In total 22 features are calculated.
#'
#' qualityScores_EM
#'
#' @param chipName String, filename and path to the ChIP bam file 
#' (without extension)
#' @param inputName String, filename and path to the Input bam file
#' (without extension)
#' @param read_length Integer, length of the reads
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param savePlotPath, set if Cross-correlation plot should be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
#'
#' @return returnList, contains
#' QCscores_ChIP List of QC-metrics with crosscorrelation values for the ChIP
#' QCscores_binding List of QCscores from peak calls
#' TagDensityChip Tag-density profile, smoothed by the Gaussian kernel 
#' (for further details see "spp" package)
#' TagDensityInput Tag density-profile, smoothed by the Gaussian kernel 
#' (for further details see "spp" package)
#'
#' @export
#'
#' @examples
#'
#' ## This command is time intensive to run
#'
#' ## To run this example code the user MUST provide 2 bam files: one for ChIP 
#' ## and one for the input". Here we used ChIP-seq data from ENCODE. Two 
#' ## example files can be downloaded using the following link:
#' ## https://www.encodeproject.org/files/ENCFF000BFX/
#' ## https://www.encodeproject.org/files/ENCFF000BDQ/
#' ## and save them in the working directory (here given in the temporary 
#' ## directory "filepath"
#'
#' mc=4
#' \dontrun{
#' 
#' filepath=tempdir()
#' setwd(filepath)
#' 
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BDQ/@@download/ENCFF000BDQ.bam")
#' 
#' chipName=file.path(filepath,"ENCFF000BFX")
#' inputName=file.path(filepath,"ENCFF000BDQ")
#' 
#' CC_Result=qualityScores_EM(chipName=chipName, inputName=inputName, 
#' read_length=36, mc=mc, annotationID = "hg19")
#'}



qualityScores_EM <- function(chipName, inputName, read_length, 
    annotationID = "hg19", mc = 1, savePlotPath = NULL, debug = FALSE) 
{
    start_time <- Sys.time()
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 6, 
        clear = FALSE, width = 60)
    pb$tick()

    ########## check if input format is ok
    if (!(is.character(chipName) & is.character(inputName))) 
        stop("Invalid chipName or inputName (String required)")
    
    if (!is.numeric(read_length)) 
        stop("read_length must be numeric")
    if (read_length < 1) 
        stop("read_length must be > 0")
    
    annotationID=f_annotationCheck(annotationID)
    
    if (!is.numeric(mc)) {
        warning("mc must be numeric")
        mc <- 1
    }
    
    if (mc < 1) {
        warning("mc set to 1")
        mc <- 1
    }
    cluster=NULL
    ########## 
    pb$tick()

    message("reading bam files")
    message("...for ChIP")
    chip.data <- readBamFile(chipName)
    
    pb$tick()

    message("\n...for Input")
    input.data <- readBamFile(inputName)


    if ( debug ) {
        message("Debugging mode ON")
        save(chip.data, input.data, 
            file = file.path(getwd(), "bamFiles.RData"))
    }

    pb$tick()

    ## plot and calculate cross correlation and phantom 
    ## characteristics for the ChIP
    message("\ncalculating binding characteristics for ChIP... ")

    estimating_fragment_length_range <- c(0, 500)
    estimating_fragment_length_bin <- 5
    
    #switch cluster on
    if (mc > 1) {
        cluster <- parallel::makeCluster( mc )
    }
    pb$tick()

    chip_binding.characteristics<-spp::get.binding.characteristics(chip.data, 
        srange = estimating_fragment_length_range, 
        bin = estimating_fragment_length_bin, 
        accept.all.tags = TRUE, cluster = cluster)
    
    pb$tick()

    #switch cluster off   
    if (mc > 1) {
        parallel::stopCluster( cluster )
    }

    message("\n***Calculating cross correlation QC-metrics for Chip...***")

    crossvalues_Chip <- getCrossCorrelationScores(chip.data, 
        chip_binding.characteristics, 
        read_length = read_length, 
        savePlotPath = savePlotPath, 
        mc = mc,
        annotationID = annotationID)
    ## save the tag.shift
    final.tag.shift <- crossvalues_Chip$tag.shift
    
    ## plot and calculate cross correlation and phantom 
    ## characteristics for the input
    # message("calculating binding characteristics for Input...")
    # pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 3, 
    #     clear = FALSE, width = 60)
    # pb$tick()

    # #switch cluster on
    # if (mc > 1) {
    #     cluster <- parallel::makeCluster( mc )
    # }

    # input_binding.characteristics<-spp::get.binding.characteristics(input.data,
    #     srange = estimating_fragment_length_range, 
    #     bin = estimating_fragment_length_bin, 
    #     accept.all.tags = TRUE, cluster = cluster)
    # pb$tick()

    # #switch cluster off   
    # if (mc > 1) {
    #     parallel::stopCluster( cluster )
    # }


    # if ( debug ) {
    #     save(chip_binding.characteristics, input_binding.characteristics, 
    #         file = file.path(getwd(),"bindingCharacteristics.RData"))
    # }

    if ( debug ) {
        save(chip_binding.characteristics, 
            file = file.path(getwd(),"bindingCharacteristics.RData"))
    }


    ## get chromosome information and order chip and input by it
    chrl_final <- intersect(names(chip.data$tags), names(input.data$tags))
    chip.data$tags <- chip.data$tags[chrl_final]
    chip.data$quality <- chip.data$quality[chrl_final]
    input.data$tags <- input.data$tags[chrl_final]
    input.data$quality <- input.data$quality[chrl_final]
    
    ## remove sigular positions with extremely high tag counts with respect 
    ## to the neighbourhood
    message("\nremoving loval tag anomalies...")
    selectedTags <- removeLocalTagAnomalies(chip.data, input.data, 
        chip_binding.characteristics)
        #input_binding.characteristics)
    pb$tick()

    #cleaning up memory space
    remove(chip_binding.characteristics)
    #remove(input_binding.characteristics)

    input.dataSelected <- selectedTags$input.dataSelected
    chip.dataSelected <- selectedTags$chip.dataSelected
    
    if ( debug) {
        save(chip.dataSelected, input.dataSelected, 
            file = file.path(getwd(),"dataSelected.RData"))
    }
    
    ## get QC-values from peak calling
    message("\n***Calculating QC-metrics from peak-calling...***")
    bindingScores <- getPeakCallingScores(chip.data, 
        input.data, 
        chip.dataSelected, 
        input.dataSelected, 
        final.tag.shift, 
        mc=mc,
        annotationID=annotationID,
        debug=debug)
    
    message("\n***Tag smooting...***")
    
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 3, 
        clear = FALSE, width = 60)
    pb$tick()

    ## objects of smoothed tag density for ChIP and Input
    smoothed.densityChip <- tagDensity(chip.dataSelected, final.tag.shift, 
        annotationID = annotationID, mc = mc)
    pb$tick()

    smoothed.densityInput <- tagDensity(input.dataSelected, final.tag.shift, 
        annotationID = annotationID, mc = mc)
    pb$tick()
    
    if ( debug )
    {
        message("saving tracks as wig...")
        f_writewig(smoothed.densityChip, 
            file.path(getwd(), "chip.wig"),"track chip")
        f_writewig(smoothed.densityInput, 
            file.path(getwd(),"input.wig"),"track input")
    }

    returnList <- list(QCscores_ChIP = crossvalues_Chip, 
        #QCscores_Input = crossvalues_Input, 
        QCscores_binding = bindingScores, 
        TagDensityChip = smoothed.densityChip, 
        TagDensityInput = smoothed.densityInput)
    
    if ( debug ) {
        writeout <- list(QCscores_ChIP = crossvalues_Chip, 
            #QCscores_Input = crossvalues_Input, 
            QCscores_binding = bindingScores)
        filename <- file.path(getwd(), "CC.results")
        write.table(writeout, file = filename,
            row.names = TRUE, col.names = TRUE, 
            append = FALSE, quote = FALSE)
        save(smoothed.densityChip, smoothed.densityInput, 
            file = file.path(getwd(), "smoothed.RData"))
    }
    message("Calculation of EM done!")
    end_time <- Sys.time()
    message("Time used: ")
    message(end_time - start_time)
    return(returnList)
    
}
