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
############### FUNCTION QC-metrics for narrow binding PROFILES ######
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
#' @param chipName Character, filename (and optional path) for the ChIP bam file (without the .bam extension)
#' @param inputName Character, filename (and optional path) for the Input control bam file (without the .bam extension)
#' @param read_length Integer, length of the reads
#' @param chip.data Optional, taglist object for ChIP reads as returned by spp or readBamFile() function. If not set (NULL) the data will be read from 
#' @param input.data Optional, taglist object for Input control reads as returned by spp or readBamFile() function
#' @param readAlignerType string, bam (default) tagAlign file format are supported
#' @param annotationID Character, indicating the genome assembly (Default="hg19")
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param crossCorrelation_Input Boolean, calculates cross-correlation and 
#' and EM metrics for the input. The default=FALSE as the running time 
#' increases and the metrics are not used in quality prediction. 
#' @param savePlotPath, set if Cross-correlation plot should be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout
#' @param downSamplingChIP Boolean, to be used to downsample reads within enrichment peaks.
#' This option was used for generating simulated low quality (low enrichment) profiles for testing the prediction models.
#' The default is FALSE and should generally not be used by end users.
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
#' @param writeWig Boolean, saves smoothed tag density in  
#' wig format in working directory for Input and ChIP
#'
#' @return returnList, contains
#' QCscores_ChIP List of QC-metrics with crosscorrelation values for the ChIP
#' QCscores_Input List of QC-metrics with crosscorrelation values for the Input 
#' if "crossCorrelation_Input" parameter was set to TRUE, NULL otherwise
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
    chip.data=NULL, input.data=NULL, readAlignerType = "bam",
    annotationID = "hg19",  mc = 1, crossCorrelation_Input=FALSE,
    downSamplingChIP=FALSE, writeWig=FALSE,
    savePlotPath = NULL, debug = FALSE) 
{
    start_time <- Sys.time()
   message("***Calculating EM metrics...***")
   pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 8, 
        clear = FALSE, width = 60)
    pb$tick()

    ########## check if input format is ok
    if (!(is.character(chipName) & is.character(inputName))) 
        stop("Invalid chipName or inputName (String required)")
    
    if (!is.numeric(read_length)) 
        stop("read_length must be numeric")
    if (read_length < 1) 
        stop("read_length must be > 0")
    if(crossCorrelation_Input)
    {
        message("Calculating the cross-correlation of the Input 
            will increase running time significantly! ")
    }
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
    pb$tick()

    # read input data from bamfile (unless they are passed as taglist object among input parameters)
    if (is.null(chip.data)) {
        message("reading bam file for ChIP")
        chip.data <- readBamFile(chipName,readAlignerType = readAlignerType)
        pb$tick()
    }
    if (is.null(input.data)) {
        message("reading bam file for Input control")
        input.data <- readBamFile(inputName,readAlignerType = readAlignerType)
        pb$tick()
    }

        if ( debug ) {
            message("Debugging mode ON, saving taglist objects to bamFiles.RData")
            save(chip.data, input.data, 
                file = file.path(getwd(), "bamFiles.RData"))
        }
        pb$tick()


    if (downSamplingChIP){
        message("downsampling ChIP data. This can take a while!")
        chip.dataNew=downsample_ChIPpeaks(chip.data=chip.data, input.data=input.data,
            read_length=read_length,
            annotationID=annotationID,mc=mc,debug=debug)
        chip.data=chip.dataNew

        if ( debug ) 
        {
            save(chip.dataNew, 
                file = file.path(getwd(),"ChIP_downsampledBam.RData"))
        }

    }

    ## plot and calculate cross correlation and phantom 
    ## characteristics for the ChIP
    message("\nCalculating binding characteristics for ChIP... ")

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
    
    
    if ( debug ) {
        save(chip_binding.characteristics, 
            file = file.path(getwd(),"bindingCharacteristics.RData"))
    }

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

    pb$tick()

    crossvalues_Input=NULL
    if (crossCorrelation_Input){
        message("\n***Calculating cross correlation QC-metrics for Input...***")

        # plot and calculate cross correlation and phantom 
        # characteristics for the input        
        #switch cluster on
        if (mc > 1) {
            cluster <- parallel::makeCluster( mc )
        }
        message("calculating binding characteristics...")
        input_binding.characteristics<-spp::get.binding.characteristics(input.data,
            srange = estimating_fragment_length_range, 
            bin = estimating_fragment_length_bin, 
            accept.all.tags = TRUE, cluster = cluster)

        #switch cluster off   
        if (mc > 1) {
            parallel::stopCluster( cluster )
        }
        message("calculating cross-correlation scores...")
        crossvalues_Input <- getCrossCorrelationScores(input.data, 
            input_binding.characteristics, 
            read_length = read_length, 
            savePlotPath = savePlotPath, 
            mc = mc,
            tag="Input",
            annotationID = annotationID)

        if ( debug ) {
            save(input_binding.characteristics, 
                file = file.path(getwd(),"bindingCharacteristicsInput.RData"))
        }
        message("calculating cross-correlation scores for Input... done!")

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

    input.tagsSelected <- selectedTags$input.dataSelected
    chip.tagsSelected <- selectedTags$chip.dataSelected
    
    if ( debug) {
        save(chip.tagsSelected, input.tagsSelected, 
            file = file.path(getwd(),"dataSelected.RData"))
    }
    
    ## get QC-values from peak calling
    message("\n***Calculating QC-metrics from peak-calling...***")
    bindingScores <- getPeakCallingScores(chip.data, 
        input.data, 
        chip.tagsSelected, 
        input.tagsSelected, 
        final.tag.shift, 
        mc=mc,
        annotationID=annotationID,
        debug=debug)
    
    returnList <- list(QCscores_ChIP = crossvalues_Chip, 
        QCscores_Input = crossvalues_Input, 
        QCscores_binding = bindingScores, 
        SelectedTagsChip = chip.tagsSelected, 
        SelectedTagsInput = input.tagsSelected
        )

    if ( writeWig )
    {
        message("compute track for ChIP")
        smoothed.densityChip<- tagDensity(data=chip.data, tag.shift=final.tag.shift, annotationID = annotationID, mc = mc)
        message("saving ChIP track as wig...")
        f_writewig(smoothed.densityChip, file.path(getwd(), "chip.wig"),"track chip")
        
        message("compute track for input")
        smoothed.densityInput<- tagDensity(data=input.data, tag.shift=final.tag.shift, annotationID = annotationID, mc = mc)
        message("saving input track as wig...")
        f_writewig(smoothed.densityInput, file.path(getwd(),"input.wig"),"track input")

    }
    

    message("Calculation of EM metrics done!")
    end_time <- Sys.time()
    message("Time used: ")
    message(end_time - start_time)
    return(returnList)
    
}
