#' @import ChIC.data
#' @import caret
# @rawNamespace import(girafe, except = c(plot,reduce))
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

#######################################################################
###############                                         ###############  
############### FUNCTIONS QC-metrics for narrow binding PROFILES ######
###############                                         ###############
#######################################################################



#' @title Wrapper function to calculate QC-metrics from cross-correlation 
#' analysis, QC-metrics designed for TFs and QC-metrics from peak-calls
#'
#' @description
#' Wrapper function that reads input bam files and provides QC-metrics from 
#' cross-correlation analysis, from peak calling and general metrics like 
#' read-length or NRF. In total 22 features are calculated.
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
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout.
#' @param debug Boolean, to enter debugging mode (default= FALSE)
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
#' mc=4
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
#' read_length=36, mc=mc)
#'}


qualityScores_EM <- function(chipName, inputName, read_length, 
    annotationID = "hg19", mc = 1, savePlotPath = NULL, debug = FALSE) 
{
    ########## check if input format is ok
    if (!(is.character(chipName) & is.character(inputName))) 
        stop("Invalid chipName or inputName (String required)")
    
    if (!is.numeric(read_length)) 
        stop("read_length must be numeric")
    if (read_length < 1) 
        stop("tag.shift must be > 0")
    
    if (!(is.character(annotationID) & (annotationID == "hg19"))) {
        warning("Invalid annotationID. Setting back to default value. 
        (Currently only hg19 is supported)")
        annotationID <- "hg19"
    }
    
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

    message("reading bam files")
    chip.data <- readBamFile(chipName)
    input.data <- readBamFile(inputName)
    
    if (debug) {
        save(chip.data, input.data, file = "bamFiles.RData")
    }
    
    ## plot and calculate cross correlation and phantom 
    ## characteristics for the ChIP
    message("calculate binding characteristics ChIP")
    ## cross_correlation parameters
    estimating_fragment_length_range <- c(0, 500)
    estimating_fragment_length_bin <- 5
    
    ## chip_binding.characteristics<-get.binding.characteristics(chip.data,
    ## srange=estimating_fragment_length_range, 
    ## bin=estimating_fragment_length_bin,
    ## accept.all.tags=T)

    #switch cluster on
    if (mc > 1) {
        cluster <- makeCluster( mc )
    }

    chip_binding.characteristics<-spp::get.binding.characteristics(chip.data, 
        srange = estimating_fragment_length_range, 
        bin = estimating_fragment_length_bin, 
        accept.all.tags = TRUE, cluster = cluster)
    
    #switch cluster off   
    if (mc > 1) {
        stopCluster( cluster )
    }
    message("calculate cross correlation QC-metrics for the Chip")
    crossvalues_Chip <- getCrossCorrelationScores(chip.data, 
        chip_binding.characteristics, 
        read_length = read_length, 
        savePlotPath = savePlotPath, 
        mc = mc,
        plotname = "ChIP")
    ## save the tag.shift
    final.tag.shift <- crossvalues_Chip$tag.shift
    
    ## plot and calculate cross correlation and phantom 
    ## characteristics for the input
    message("calculate binding characteristics Input")
        
    #switch cluster on
    if (mc > 1) {
        cluster <- makeCluster( mc )
    }

    input_binding.characteristics<-spp::get.binding.characteristics(input.data,
        srange = estimating_fragment_length_range, 
        bin = estimating_fragment_length_bin, 
        accept.all.tags = TRUE, cluster = cluster)
    
    #switch cluster off   
    if (mc > 1) {
        stopCluster( cluster )
    }

    #message("calculate cross correlation QC-metrics for the Input")
    #crossvalues_Input <- getCrossCorrelationScores(input.data, 
    #    input_binding.characteristics, 
    #    read_length = read_length, 
    #    savePlotPath = savePlotPath, plotname = "Input")
    
    if (debug) {
        save(chip_binding.characteristics, input_binding.characteristics, 
            file = "bindingCharacteristics.RData")
    }
    
    ## get chromosome information and order chip and input by it
    chrl_final <- intersect(names(chip.data$tags), names(input.data$tags))
    chip.data$tags <- chip.data$tags[chrl_final]
    chip.data$quality <- chip.data$quality[chrl_final]
    input.data$tags <- input.data$tags[chrl_final]
    input.data$quality <- input.data$quality[chrl_final]
    
    ## remove sigular positions with extremely high tag counts with respect 
    ## to the neighbourhood
    selectedTags <- removeLocalTagAnomalies(chip.data, input.data, 
        chip_binding.characteristics, 
        input_binding.characteristics)
    
    input.dataSelected <- selectedTags$input.dataSelected
    chip.dataSelected <- selectedTags$chip.dataSelected
    
    if (debug) {
        save(chip.dataSelected, input.dataSelected, 
            file = "dataSelected.RData")
    }
    
    ## get QC-values from peak calling
    bindingScores <- getPeakCallingScores(chip.data, 
        input.data, 
        chip.dataSelected, 
        input.dataSelected, 
        final.tag.shift, mc=mc)
    
    ## objects of smoothed tag density for ChIP and Input
    smoothed.densityChip <- tagDensity(chip.dataSelected, final.tag.shift, 
        annotationID = annotationID, mc = mc)
    smoothed.densityInput <- tagDensity(input.dataSelected, final.tag.shift, 
        annotationID = annotationID, mc = mc)
    
    returnList <- list(QCscores_ChIP = crossvalues_Chip, 
        #QCscores_Input = crossvalues_Input, 
        QCscores_binding = bindingScores, 
        TagDensityChip = smoothed.densityChip, 
        TagDensityInput = smoothed.densityInput)
    
    if (debug) {
        writeout <- list(QCscores_ChIP = crossvalues_Chip, 
            #QCscores_Input = crossvalues_Input, 
            QCscores_binding = bindingScores)
        write.table(writeout, file = file.path(getwd(), "CC.results"))
        save(smoothed.densityChip, smoothed.densityInput, 
            file = "smoothed.RData")
    }
    return(returnList)
    
}



#' @title QC-metrics from cross-correlation profile, phantom peak and 
#' general QC-metrics
#'
#' @description 
#' We use cross-correlation analysis to obtain QC-metrics proposed for 
#' narrow-binding patterns. After calculating the strand cross-correlation
#' coefficient (Kharchenko et al., 2008), we take the following values from
#' the profile: 
#' coordinates of the ChIP-peak (fragment length, height A), coordinates at 
#' the phantom-peak (read length, height B) and the baseline (C), the 
#' strand-shift, the number of uniquely mapped reads (unique_tags), uniquely 
#' mapped reads corrected by the library size, the number of reads and the 
#' read lengths. We calculate different values using the relative and absolute 
#' height of the cross-correlation peaks: the relative and normalized strand 
#' coefficient RSC and NSC  (Landt et al., 2012), and the quality control tag 
#' (Marinov et al., 2013). Other values regarding the library complexity 
#' (Landt et al., 2012) like the fraction of non-redundant mapped reads 
#' (NRF; ratio between the number of uniquely mapped reads divided by the 
#' total number of reads), the NRF adjusted by library size and ignoring the 
#' strand direction (NRF_nostrand), and the PCR bottleneck coefficient PBC 
#' (number of genomic locations to which exactly one unique mapping read maps, 
#' divided by the number of unique mapping reads).
#'
#' getCrossCorrelationScores
#'
#' @param data data-structure with tag information read from bam file 
#' (see readBamFile())
#' @param bchar binding.characteristics is a data-structure containing binding 
#' information for binding preak separation distance and cross-correlation 
#' profile (see spp::get.binding.characteristics).
#' @param read_length Integer, read length of "data" (Defaul="36") 
#' @param savePlotPath if set the plot will be saved under "savePlotPath". 
#' Default=NULL and plot will be forwarded to stdout. 
#' @param plotname Name of the crossCorrelation plot (pdf). Available only 
#' if savePlotPath is set
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#'
#' @return finalList List with QC-metrics 
#'
#' @export
#'
#' @examples
#'
#' ## This command is time intensive to run
#'
#' ## To run the example code the user must provide a bam file and read it 
#' ## with the readBamFile() function. To make it easier for the user to run 
#' ## the example code we provide a bam file in our ChIC.data package that has 
#' ## already been loaded with the readBamFile() function.
#' 
#' mc=4
#' print("Cross-correlation for ChIP")
#' \dontrun{
#' filepath=tempdir()
#' setwd(filepath)
#' data("chipSubset", package = "ChIC.data", envir = environment())
#' chipBam=chipSubset
#'
#' ## calculate binding characteristics 
#'
#' chip_binding.characteristics<-spp::get.binding.characteristics( chipBam, 
#'     srange=c(0,500), bin = 5, accept.all.tags = TRUE)
#'
#' crossvalues_Chip<-getCrossCorrelationScores( chipBam , 
#'     chip_binding.characteristics, read_length = 36, 
#'     savePlotPath = filepath, mc = mc)
#'}

getCrossCorrelationScores <- function(data, bchar, read_length, 
    savePlotPath = NULL, plotname = "phantom", mc=1) 
{
    ########## check if input format is ok
    if (!(is.list(data) & (length(data) == 2L))) 
        stop("Invalid format for data")
    
    if (!(is.list(bchar) & (length(bchar) == 3L))) 
        stop("Invalid format for bchar")
    
    if (!is.numeric(read_length)) 
        stop("read_length must be numeric")
    if (read_length < 1) 
        stop("tag.shift must be > 0")
    
    ########## 
    cluster=NULL
    
    ## cross_correlation_customShift_withSmoothing parameters
    ## ccRangeSubset=cross_correlation_range_subset
    ccRangeSubset <- 100:500
    cross_correlation_smoothing_k <- 10
    ## Phantom peaks
    PhantomPeak_range <- c(-500, 1500)
    PhPeakbin <- 5
    
    ## step 1.2: Phantom peak and cross-correlation
    if (mc > 1) {
        cluster <- makeCluster( mc )
    }

    message("Phantom peak and cross-correlation")
    phChar <- spp::get.binding.characteristics(data, 
        srange = PhantomPeak_range, 
        bin = PhPeakbin, accept.all.tags = TRUE,
        cluster = cluster)

    if (mc>1) {
        stopCluster( cluster )
    }

    ph_peakidx <- which(
        (phChar$cross.correlation$x >= (read_length - round(2 * PhPeakbin))) & 
        (phChar$cross.correlation$x <= (read_length + round(1.5 * PhPeakbin))))
    
    ph_peakidx <- ph_peakidx[which.max(phChar$cross.correlation$y[ph_peakidx])]
    phantom_cc <- phChar$cross.correlation[ph_peakidx, ]
    
    ## Minimum value of cross correlation in srange
    min_cc <- phChar$cross.correlation[which.min(phChar$cross.correlation$y), ]
    ## Normalized Strand cross-correlation coefficient (NSC)
    NSC <- phChar$peak$y/min_cc$y
    ## Relative Strand Cross correlation Coefficient (RSC)
    RSC <- (phChar$peak$y - min_cc$y)/(phantom_cc$y - min_cc$y)
    
    ## Quality flag based on RSC value
    qflag <- f_qflag(RSC)
    ## phScores=phantom_peak.scores
    phScores <- list(phantom_cc = phantom_cc, 
        NSC = NSC, RSC = RSC, quality_flag = qflag, 
        min_cc = min_cc, peak = phChar$peak, read_length = read_length)
    
    message("smooting...")
    ## 2.0 smoothed cross correlation ss is subset selection
    ss <- which(bchar$cross.correlation$x %in% ccRangeSubset)
    bchar$cross.correlation <- bchar$cross.correlation[ss, ]
    ## add smoothing
    bindCharCC_Ysmoothed <- caTools::runmean(bchar$cross.correlation$y, 
        k = cross_correlation_smoothing_k)
    ## assign the new maximum coordinates
    bchar$peak$y <- max(bindCharCC_Ysmoothed)
    iindex <- (which(bindCharCC_Ysmoothed == bchar$peak$y))
    bchar$peak$x <- bchar$cross.correlation$x[iindex]
    
    ## plot cross correlation curve with smoothing
    
    strandShift <- bchar$peak$x
    message("Check strandshift...")
    
    a <- tryCatch({
        deriv <- diff(bindCharCC_Ysmoothed)/diff(bchar$cross.correlation$x)
        deriv <- append(0, deriv)
        iindex <- which.max( abs(( diff(sign(deriv) ))))
        newShift <- bchar$cross.correlation$x[iindex]
    }, warning = function(w) {
        newShift <- strandShift
        message("strandshift problems...")
    })
    message(a)
    
    # # newShift=f_getCustomStrandShift(x=bchar$cross.correlation$x, #
    # y=bindCharCC_Ysmoothed) message('newShift is ',newShift) oldShift=NULL if
    # (newShift!='ERROR') {
    if (newShift != strandShift) {
        oldShift <- strandShift
        strandShift <- newShift
        message("Strandshift is substituted")
    } else {
        message("strandshift remains the same...")
    }
    
    ## 2.2 phantom peak with smoothing
    message("Phantom peak with smoothing")
    ## phantom.characteristics<-phantom.characteristics select a subset of 
    ## cross correlation profile where we expect the peak
    ss_forPeakcheck <- which(phChar$cross.correlation$x %in% ccRangeSubset)
    ## phantomSmoothing=phantom.characteristics_cross.correlation_y_smoothed
    phantomSmoothing <- caTools::runmean(phChar$cross.correlation$y, 
        k = cross_correlation_smoothing_k)
    ## assign the new maximum coordinates
    max_y_peakcheck <- max(phantomSmoothing[ss_forPeakcheck])
    ## indexMaxPeak=index_for_maxpeak
    indexMaxPeak <- (which(phantomSmoothing == max_y_peakcheck))
    ## use this additional control just in case there 
    ## is another cross correlation
    ## bin with exactly the ame height outside of the desired 
    ## search region (e.g.
    ## important if the phantom peak is as high or higher than the main cross
    ## correlation peak)
    indexMaxPeak <- indexMaxPeak[(indexMaxPeak %in% ss_forPeakcheck)]
    ## force index 1 just in case there are two points with 
    ## exactly the same cc value
    max_x_peakcheck <- phChar$cross.correlation$x[(indexMaxPeak[1])]
    
    if (max_x_peakcheck > phChar$peak$x) {
        whatever <- which(phChar$cross.correlation$x == max_x_peakcheck)
        phChar$peak$y <- phChar$cross.correlation$y[whatever]
        phChar$peak$x <- max_x_peakcheck
        phScores$peak <- phChar$peak
        ## Normalized Strand cross-correlation coefficient (NSC)
        NSC <- phChar$peak$y/phScores$min_cc$y
        phScores$NSC <- NSC
        ## Relative Strand Cross correlation Coefficient (RSC)
        divisor <- (phScores$phantom_cc$y - phScores$min_cc$y)
        RSC <- (phChar$peak$y - phScores$min_cc$y)/divisor
        phScores$RSC <- RSC
        ## Quality flag based on RSC value
        qflag <- f_qflag(RSC)
        phScores$quality_flag <- qflag
    } else {
        phChar$peak$y
        ## <-max(phantom.characteristics_cross.correlation_y_smoothed)
        phChar$peak$x
    }
    
    if (!is.null(savePlotPath)) {
        filename <- file.path(savePlotPath, 
            paste(plotname, "CrossCorrelation.pdf", sep = "_"))
        pdf(file = filename)
    }
    
    ## plot cross correlation curve with smoothing
    message("plot cross correlation curve with smoothing")
    par(mar = c(3.5, 3.5, 1, 0.5), mgp = c(2, 0.65, 0), cex = 0.8)
    plot(phChar$cross.correlation, type = "l", xlab = "strand shift", 
        ylab = "cross-correlation", 
        main = "CrossCorrelation Profile")
    lines(x = phChar$cross.correlation$x, y = phantomSmoothing, 
        lwd = 2, col = "blue")
    lines(x = rep(phScores$peak$x, times = 2), 
        y = c(0, phScores$peak$y), lty = 2, 
        lwd = 2, col = "red")
    lines(x = rep(phScores$phantom_cc$x, times = 2), 
        y = c(0, phScores$phantom_cc$y), 
        lty = 2, lwd = 2, col = "orange")
    abline(h = phScores$min_cc$y, lty = 2, lwd = 2, col = "grey")
    text(x = phScores$peak$x, y = phScores$peak$y, 
        labels = paste("A =", signif(phScores$peak$y, 3)), 
        col = "red", pos = 3)
    text(x = phScores$phantom_cc$x, y = phScores$phantom_cc$y, 
        labels = paste("B =", signif(phScores$phantom_cc$y, 3)), 
        col = "orange", pos = 2)
    text(x = min(phChar$cross.correlation$x), 
        y = phScores$min_cc$y, 
        labels = paste("C =", signif(phScores$min_cc$y, 3)), 
        col = "grey", adj = c(0, -1))
    legend(x = "topright", 
        legend = c(paste("NSC = A/C =", signif(phScores$NSC, 3)), 
        paste("RSC = (A-C)/(B-C) =", signif(phScores$RSC, 3)), 
        paste("Quality flag =", phScores$quality_flag), "", 
        paste("Shift =", (phScores$peak$x)), 
        paste("Read length =", (read_length))))
    
    if (!is.null(savePlotPath)) {
        dev.off()
        message("pdf saved under", filename)
    }
    
    phantomScores <- list(CC_NSC = round(phScores$NSC, 3), 
        CC_RSC = round(phScores$RSC, 3), 
        CC_QualityFlag = phScores$quality_flag, 
        CC_shift. = phScores$peak$x, 
        CC_A. = round(phScores$peak$y, 3), 
        CC_B. = round(phScores$phantom_cc$y, 3), 
        CC_C. = round(phScores$min_cc$y, 3))
    
    ## 4 NRF calculation
    message("NRF calculation")
    
    ALL_TAGS <- sum(lengths(data$tags))
    UNIQUE_TAGS <- sum(lengths(lapply(data$tags, unique)))
    UNIQUE_TAGS_nostrand <- sum(lengths(lapply(data$tags, FUN = function(x) {
        unique(abs(x))
    })))
    
    NRF <- UNIQUE_TAGS/ALL_TAGS
    NRF_nostrand <- UNIQUE_TAGS_nostrand/ALL_TAGS
    
    ## to compensate for lib size differences
    message("compensate for lib size differences")
    ## nomi<-rep(names(data$tags), sapply(data$tags, length))
    nomi <- rep(names(data$tags), lapply(data$tags, length))
    
    dataNRF <- unlist(data$tags)
    names(dataNRF) <- NULL
    dataNRF <- paste(nomi, dataNRF, sep = "")
    
    if (ALL_TAGS > 1e+07) {
        
        UNIQUE_TAGS_LibSizeadjusted <- round(mean(sapply(1:100, 
            FUN = function(x) {
            return(length(unique(sample(dataNRF, size = 1e+07))))
        })))
        
    } else {
        UNIQUE_TAGS_LibSizeadjusted <- round(mean(sapply(1:100, #6053517
            FUN = function(x) {
            return(length(unique(sample(dataNRF, 
                size = 1e+07, replace = TRUE))))
        })))
    }

    NRF_LibSizeadjusted <- UNIQUE_TAGS_LibSizeadjusted/1e+07
    STATS_NRF <- list(CC_ALL_TAGS = ALL_TAGS, 
        CC_UNIQUE_TAGS = UNIQUE_TAGS, 
        CC_UNIQUE_TAGS_nostrand = UNIQUE_TAGS_nostrand, 
        CC_NRF = NRF, CC_NRF_nostrand = NRF_nostrand, 
        CC_NRF_LibSizeadjusted = NRF_LibSizeadjusted)
    
    ## N1= number of genomic locations to which EXACTLY one 
    ## unique mapping read maps
    ## Nd = the number of genomic locations to which AT LEAST 
    ## one unique mapping read
    ## maps, i.e. the number of non-redundant, unique mapping reads
    ## N1<-sum(sapply(data$tags, FUN=function(x) {checkDuplicate<-duplicated(x)
    ## duplicated_positions<-unique(x[checkDuplicate]) OUT<-x[!(x %in%
    ## duplicated_positions)] return(length(OUT)) }))
    
    N1 <- sum(unlist(lapply(data$tags, FUN = function(x) {
        checkDuplicate <- duplicated(x)
        duplicated_positions <- unique(x[checkDuplicate])
        OUT <- x[!(x %in% duplicated_positions)]
        return(length(OUT))
    })))
    
    ## Nd<-sum(sapply(lapply(data$tags, unique), length))
    Nd <- sum(unlist(lapply(data$tags, FUN = function(x) {
        length(unique(x))
    })))
    
    PBC <- N1/Nd
    tag.shift <- round(strandShift/2)
    
    finalList <- append(append(list(CC_StrandShift = strandShift, 
        tag.shift = tag.shift, 
        N1 = round(N1, 3), 
        Nd = round(Nd, 3), 
        CC_PBC = round(PBC, 3), 
        CC_readLength = read_length, 
        CC_UNIQUE_TAGS_LibSizeadjusted = UNIQUE_TAGS_LibSizeadjusted), 
        phantomScores), 
        STATS_NRF)
    return(finalList)
}


#' @title Calculating QC-values from peak calling procedure
#'
#' @description
#' QC-metrics based on the peak calling are the fraction of usable reads in 
#' the peak regions (FRiP) (Landt et al., 2012), for which the function calls 
#' sharp- and broad-binding peaks to obtain two types: the FRiP_sharpsPeak and 
#' the FRiP_broadPeak. The function takes the number of called of peaks using 
#' an FDR of 0.01 and an evalue of 10 (Kharchenko et al., 2008). And count 
#' the number of peaks called when using the sharp- and broad-binding option. 
#'
#' getPeakCallingScores
#'
#' @param chip data-structure with tag information for the ChIP 
#' (see readBamFile())
#' @param input data-structure with tag information for the Input
#' (see readBamFile())
#' @param chip.dataSelected selected ChIP tags after running 
#' removeLocalTagAnomalies() which removes local tag anomalies
#' @param input.dataSelected selected Input tags after running 
#' removeLocalTagAnomalies() which removes local tag anomalies
#' @param tag.shift Integer containing the value of the tag shift, 
#' calculated by getCrossCorrelationScores(). Default=75
#' @param chrorder chromosome order (default=NULL) 
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#'
#' @return QCscoreList List with 6 QC-values
#'
#' @export
#' @examples
#'
#' mc=4
#' finalTagShift=98
#' print("Cross-correlation for ChIP")
#'
#' \dontrun{
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
#'     chipBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE)
#'
#' input_binding.characteristics<-spp::get.binding.characteristics(
#'     inputBam, srange=c(0,500), bin = 5, accept.all.tags = TRUE)
#'
#' ##get chromosome information and order chip and input by it
#' chrl_final <- intersect(names(chipBam$tags), names(inputBam$tags))
#' chipBam$tags <- chipBam$tags[chrl_final]
#' chipBam$quality <- chipBam$quality[chrl_final]
#' inputBam$tags <- inputBam$tags[chrl_final]
#' inputBam$quality <- inputBam$quality[chrl_final]
#'
#' ##remove sigular positions with extremely high read counts with 
#' ##respect to the neighbourhood
#' selectedTags <- removeLocalTagAnomalies(chipBam, inputBam, 
#'     chip_binding.characteristics, input_binding.characteristics)
#'
#' inputBamSelected <- selectedTags$input.dataSelected
#' chipBamSelected <- selectedTags$chip.dataSelected
#'
#' ##Finally run function
#' bindingScores <- getPeakCallingScores(chip = chipBam, 
#'     input = inputBam, chip.dataSelected = chipBamSelected, 
#'     input.dataSelected = inputBamSelected, 
#'     tag.shift = finalTagShift, mc = mc)
#'}


getPeakCallingScores <- function(chip, input, chip.dataSelected, 
    input.dataSelected, tag.shift = 75, mc=1, chrorder = NULL) 
{
    ########## check if input format is ok
    if (!(is.list(chip) & (length(chip) == 2L))) 
        stop("Invalid format for data")
    if (!(is.list(input) & (length(input) == 2L))) 
        stop("Invalid format for data")
    
    if (!is.list(chip.dataSelected)) 
        stop("Invalid format for chip.dataSelected")
    if (!is.list(input.dataSelected)) 
        stop("Invalid format for input.dataSelected")
    
    if (!is.numeric(tag.shift)) 
        stop("tag.shift must be numeric")
    if (tag.shift < 1) 
        stop("tag.shift must be > 0")
    ########## 
    cluster=NULL
    ## for simplicity we use currently a shorter chromosome 
    ## frame to avoid problems
    #chrorder <- paste("chr", c(1:19, "X", "Y"), sep = "")
    chrorder <- paste("chr", c(seq_len(19), "X", "Y"), sep = "")
    
    ## 5 broadRegions 6 enrichment broad regions zthresh_list<-c(3,4)
    current_window_size <- 1000
    message("Broad regions of enrichment")
    ## for (current_window_size in window_sizes_list) { for (current_zthresh in
    ## zthresh_list) {
    current_zthresh <- 4
    broad.clusters <- spp::get.broad.enrichment.clusters(chip.dataSelected, 
        input.dataSelected, 
        window.size = current_window_size, 
        z.thr = current_zthresh, 
        tag.shift = tag.shift)
    ### start end logE znrichment write.broadpeak.info(broad.clusters,
    ### paste('broadRegions', chip.data.samplename,
    ## input,input.data.samplename, 'window', current_window_size, 'zthresh',
    ## current_zthresh,'broadPeak', sep='.'))
    md <- f_convertFormatBroadPeak(broad.clusters)
    broadPeakRangesObject <- MakeGRangesObject(Chrom = md$chr, 
        Start = md$start, 
        End = md$end)
    
    ## 12 get binding sites with FDR and eval
    chip.data12 <- chip.dataSelected[(names(chip.dataSelected) %in% chrorder)]
    input.data12<-input.dataSelected[(names(input.dataSelected) %in% chrorder)]
    
    if (mc > 1) {
        cluster <- makeCluster( mc )
    }

    message("Binding sites detection fdr")
    fdr <- 0.01
    detection.window.halfsize <- tag.shift
    message("Window Tag Density method - WTD")
    bp_FDR <- spp::find.binding.positions(signal.data = chip.data12, 
        control.data = input.data12, 
        fdr = fdr, whs = detection.window.halfsize * 2, 
        tag.count.whs = detection.window.halfsize, 
        cluster = cluster)
    FDR_detect <- sum(unlist(lapply(bp_FDR$npl, function(d) length(d$x))))
    
    message("Binding sites detection evalue")
    eval <- 10
    bp_eval <- spp::find.binding.positions(signal.data = chip.data12, 
        control.data = input.data12, 
        e.value = eval, whs = detection.window.halfsize * 2, cluster = cluster)
    eval_detect <- sum(unlist(lapply(bp_eval$npl, function(d) length(d$x))))
    
    if (mc > 1) {
        stopCluster( cluster )
    }
    
    ## output detected binding positions output.binding.results(results=bp,
    ## filename=paste('WTC.binding.positions.evalue',
    ## chip.data.samplename,'input',input.data.samplename, 'txt', sep='.'))
    if (length(bp_eval$npl) > 1) {
        ## 14 get precise binding position using escore LARGE PEAKS
        bp_broadpeak <- spp::add.broad.peak.regions(chip.data12, 
            input.data12, 
            bp_eval, 
            window.size = 1000, z.thr = 3)

        md <- f_converNarrowPeakFormat(bp_broadpeak)
        sharpPeakRangesObject <- MakeGRangesObject(Chrom = md[, 1], 
            Start = md[, 2], End = md[, 3])
        
        ## FRiP
        chip.test <- lapply(chip$tags, FUN = function(x) {
            x + tag.shift
        })
        TOTAL_reads <- sum(lengths(chip.test))
        
        ## Frip broad binding sites (histones)
        broadPeakRangesObject<-ReduceOverlappingRegions(broadPeakRangesObject)

        regions_data_list <- split(as.data.frame(broadPeakRangesObject), 
            f = seqnames(broadPeakRangesObject))
        
        chrl <- names(regions_data_list)
        names(chrl) <- chrl
        
        outcountsBroadPeak <- sum(unlist(lapply(chrl, FUN = function(chr) {
            sum(spp::points_withinFunction(x = abs(chip.test[[chr]]), 
                fs = regions_data_list[[chr]]$start, 
                fe = regions_data_list[[chr]]$end, 
                return.point.counts = TRUE))
        })))
        FRiP_broadPeak <- outcountsBroadPeak/TOTAL_reads
        
        ## Frip sharp peaks 14
        sharpPeakRangesObject<-ReduceOverlappingRegions(sharpPeakRangesObject)
        
        regions_data_list <- split(as.data.frame(sharpPeakRangesObject), 
            f = seqnames(sharpPeakRangesObject))
        
        chrl <- names(regions_data_list)
        names(chrl) <- chrl
        
        outcountsSharpPeak <- sum(unlist(lapply(chrl, FUN = function(chr) {
            sum(spp::points_withinFunction(x = abs(chip.test[[chr]]), 
                fs = regions_data_list[[chr]]$start, 
                fe = regions_data_list[[chr]]$end, 
                return.point.counts = TRUE))
        })))
        
        FRiP_sharpPeak <- outcountsSharpPeak/TOTAL_reads
    } else {
        TOTAL_reads <- 0
        FRiP_broadPeak <- 0
        outcountsBroadPeak <- 0
        FRiP_sharpPeak <- 0
        outcountsSharpPeak <- 0
    }
    
    QCscoreList <- list(CC_FDRpeaks = round(FDR_detect, 3), 
        CC_evalpeaks = round(eval_detect, 3), 
        CC_FRiP_broadPeak = round(FRiP_broadPeak, 3), 
        CC_FRiP_sharpPeak = round(FRiP_sharpPeak, 3), 
        CC_outcountsBroadPeak = outcountsBroadPeak, 
        CC_outcountsSharpPeak = outcountsSharpPeak)
    
    return(QCscoreList)
}



###############################################################################
#####                                                      ####################
##### FUNCTIONS QC-metrics for global read distribution    ####################
#####                                                      ####################
###############################################################################

#'@title Metrics taken from global read distribution
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
#' @param densityChip Smoothed tag-density object for ChIP (returned by
#' qualityScores_EM). 
#' @param densityInput Smoothed tag density object for Input (returned by
#' qualityScores_EM)
#' @param savePlotPath if set the plot will be saved under "savePlotPath". 
#' Default=NULL and plot will be forwarded to stdout.
#' @param debug Boolean to enter debugging mode (default= FALSE)
#'
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

qualityScores_GM <- function(densityChip, densityInput, savePlotPath = NULL, 
    debug = FALSE) 
{
    ########## check if input format is ok
    if (!is.list(densityChip)) 
        stop("Invalid format for densityChip")
    if (!is.list(densityInput)) 
        stop("Invalid format for densityInput")
    ########## 
    
    ## shorten frame and reduce resolution
    message("shorten frame")
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
    
    ## create cumulative function for chip and input
    message("Calculate cumsum")
    cumSumChip <- f_sortAndBinning(chip)
    cumSumInput <- f_sortAndBinning(input)
    
    ## caluclate QCvalues for chip
    chipFracTopPercent <- (1-cumSumChip[(which(cumSumChip$x >= 0.99)[1]),"pj"])
    chipFracOfBinsWithoutReads <- cumSumChip$x[(which(cumSumChip$pj > 0)[1])]
    
    ## caluclate QCvalues for input
    inputFracTopPercent <- (1-
        cumSumInput[(which(cumSumInput$x >= 0.99)[1]),"pj"])
    inputFracWithoutReads <- cumSumInput$x[(which(cumSumInput$pj > 0)[1])]
    
    ## get the point of maximum distance between the two functions
    ##  and the sign of the distance
    sign_sign <- sign(max(cumSumInput$pj - cumSumChip$pj))  ###input-chip
    arrowx <- cumSumChip[which.max(abs(cumSumInput$pj - cumSumChip$pj)), ]$x
    
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
    
    if (debug) {
        write.table(finalList, file = file.path(getwd(), "Chance.results"))
    }
    ## return QC values
    return(finalList)
}



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
#' @param smoothed.densityChip Smoothed tag-density object for ChIP (returned 
#' by qualityScores_EM)
#' @param smoothed.densityInput Smoothed tag-density object for Input (returned
#' by qualityScores_EM)
#' @param tag.shift Integer containing the value of the tag shif, calculated by
#' getCrossCorrelationScores()
#' @param annotationID String indicating the genome assembly (Default="hg19")
#' @param debug Boolean to enter debugging mode (default= FALSE)
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

createMetageneProfile <- function(smoothed.densityChip, smoothed.densityInput, 
    tag.shift, annotationID = "hg19", debug = FALSE, mc = 1) 
{
    ########## check if input format is ok
    if (!is.list(smoothed.densityChip)) 
        stop("Invalid format for smoothed.densityChip")
    if (!is.list(smoothed.densityInput)) 
        stop("Invalid format for smoothed.densityInput")
    
    if (!is.numeric(tag.shift)) 
        stop("tag.shift must be numeric")
    if (tag.shift < 1) 
        stop("tag.shift must be > 0")
    
    if (!(is.character(annotationID) & (annotationID == "hg19"))) {
        warning("Invalid annotationID. Setting back to default value. 
        (Currently only hg19 is supported)")
        annotationID <- "hg19"
    }
    
    if (!is.numeric(mc)) {
        warning("mc must be numeric")
        mc <- 1
    }
    
    if (mc < 1) {
        warning("mc set to 1")
        mc <- 1
    }
    ########## 
    
    message("Load gene annotation")
    ## require(ChIC.data)
    if (annotationID == "hg19") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("hg19_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- hg19_refseq_genes_filtered_granges
    }
    ## geneAnnotations_file=
    ## paste(annotationID,'_refseq_genes_filtered_granges',sep='')
    ## annotObject=current_annotations_object 
    ##annotObject=get(geneAnnotations_file)
    ## annotObject=get(hg19_refseq_genes_filtered_granges)
    ## current_annotations_object=refseq_genes_filtered_granges 
    ## format annotations as chromosome lists
    
    # annotObject<-data.frame(annotObject@.Data,
    # genomeIntervals::annotation(annotObject), stringsAsFactors=FALSE)
    # annotObject$interval_starts<-as.integer(annotObject$interval_starts)
    # annotObject$interval_ends<-as.integer(annotObject$interval_ends)
    # annotObject$seq_name<-as.character(annotObject$seq_name) 
    # annotatedGenesPerChr
    # <-split(annotObject, f=annotObject$seq_name)
    
    annotObjectNew <- data.frame(annotObject@.Data, annotObject@annotation, 
        stringsAsFactors = FALSE)
    annotObjectNew$interval_starts <-as.integer(annotObjectNew$interval_starts)
    annotObjectNew$interval_ends <- as.integer(annotObjectNew$interval_ends)
    annotObjectNew$seq_name <- as.character(annotObjectNew$seq_name)
    
    annotatedGenesPerChr <- split(annotObjectNew, f = annotObjectNew$seq_name)
    
    ## two.point.scaling create scaled metageneprofile input
    message("Calculating scaled metageneprofile ...")
    smoothed.densityInput <- list(td = smoothed.densityInput)
    message("process input")
    
    binned_Input <- masked_t.get.gene.av.density(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc)
    
    ## Chip
    smoothed.densityChip <- list(td = smoothed.densityChip)
    message("process ChIP")
    binned_Chip <- masked_t.get.gene.av.density(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc)
    
    geneBody <- list(chip = binned_Chip, input = binned_Input)
    
    ## one.point.scaling create non-scaled metageneprofile for TSS
    message("Creating non-scaled metageneprofiles...")
    message("...TSS")
    
    binnedInput_TSS <- masked_getGeneAvDensity_TES_TSS(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TSS")
    binnedChip_TSS <- masked_getGeneAvDensity_TES_TSS(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TSS")
    onepointTSS <- list(chip = binnedChip_TSS, input = binnedInput_TSS)
    
    ## one.point.scaling create non-scaled metageneprofile for TES
    message("...TES")
    binnedInput_TES <- masked_getGeneAvDensity_TES_TSS(smoothed.densityInput, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TES")
    binnedChip_TES <- masked_getGeneAvDensity_TES_TSS(smoothed.densityChip, 
        gdl = annotatedGenesPerChr, 
        mc = mc, tag = "TES")
    onepointTES <- list(chip = binnedChip_TES, input = binnedInput_TES)
    
    if (debug) {
        save(binned_Chip, binned_Input, file = file.path(getwd(), 
            paste("geneBody.RData", sep = "_")))
        save(binnedChip_TSS, binnedInput_TSS, file = file.path(getwd(), 
            paste("OnePointTSS.RData", sep = "_")))
        save(binnedChip_TES, binnedInput_TES, file = file.path(getwd(), 
            paste("OnePointTES.RData", sep = "_")))
    }
    
    return(list(geneBody = geneBody, TSS = onepointTSS, TES = onepointTES))
}

#'@title Read bam file
#'
#' @description
#' Reading bam file format
#' 
#' readBamFile
#'
#' @param filename, name/path of the bam file to be read (without extension)
#'
#' @return result list of lists, every list corresponds to a chromosome and 
#' contains a vector of coordinates of the 5' ends of the aligned tags.
#'
#' @export
#'
#' @examples
#'
#' ## To run this example code the user MUST provide a bam file: The user can 
#' ## download a ChIP-seq bam file for example from ENCODE:
#' ## https://www.encodeproject.org/files/ENCFF000BFX/
#' ## and save it in the working directory 
#'
#' bamID="ENCFF000BFX"
#' \dontrun{
#' filepath=tempdir()
#' setwd(filepath)
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")
#'
#' bamName=file.path(filepath,bamID)
#' chipBam=readBamFile(bamName)
#' }

readBamFile <- function(filename) {
    ########## check if input format is ok
    if (!is.character(filename)) 
        stop("Invalid filename (String required)")
    ######### 
    result <- f_readFile(filename = filename, reads.aligner.type = "bam")
    return(result)
}



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
    if (!is.list(data)) 
        stop("Invalid format for data")
    
    if (!is.numeric(tag.shift)) 
        stop("tag.shift must be numeric")
    if (tag.shift < 1) 
        stop("tag.shift must be > 0")
    
    if (!(is.character(annotationID) & (annotationID == "hg19"))) {
        warning("Invalid annotationID. Setting back to default value. 
        (Currently only hg19 is supported)")
        annotationID <- "hg19"
    }
    
    if (!is.numeric(mc)) {
        warning("mc must be numeric")
        mc <- 1
    }
    
    if (mc < 1) {
        warning("mc set to 1")
        mc <- 1
    }
    ########## 
    
    ## load chrom_info
    if (annotationID == "hg19") {
        # hg19_chrom_info=NULL
        data("hg19_chrom_info", package = "ChIC.data", envir = environment())
        chromInfo <- hg19_chrom_info
    }
    smoothed.density <- f_tagDensity(data = data, 
        tag.shift = tag.shift, 
        chromDef = chromInfo, 
        mc = mc)
    return(smoothed.density)
}
