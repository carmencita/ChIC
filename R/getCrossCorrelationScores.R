
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
#' Default=NULL and plot will be omitted. 
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param tag String,can be used to personalize the prefix of the filename 
#' for the cross- correltion plot (default="ChIP" and "Input" in case of 
#' cross-correlation plot for the input')
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
#'     annotationID="hg19",
#'     savePlotPath = filepath, mc = mc)
#'}

getCrossCorrelationScores <- function(data, bchar, annotationID="hg19", 
    read_length, savePlotPath = NULL, mc=1,tag="ChIP") 
{
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 12, 
        clear = FALSE, width = 60)

    ########## check if input format is ok
    if (!(is.list(data) & (length(data) == 2L))) 
        stop("Invalid format for data")
    
    if (!(is.list(bchar) & (length(bchar) == 3L))) 
        stop("Invalid format for bchar")
    
    if (!is.numeric(read_length)) 
        stop("read_length must be numeric")
    if (read_length < 1) 
        stop("read_length must be > 0")

    annotationID=f_annotationCheck(annotationID)
    ########## 
    cluster=NULL
        
    #chromInfo=f_chromInfoLoad(annotationID)
    
    ##filter for standard chromosomes
    #standardChroms = names(data$tags)[names(data$tags) %in%  names(chromInfo)]
    #data$tags=data$tags[standardChroms]
    #data$quality=data$quality[standardChroms]
    data=f_clearChromStructure(data,annotationID)
    pb$tick()

    ## cross_correlation_customShift_withSmoothing parameters
    ## ccRangeSubset=cross_correlation_range_subset
    ccRangeSubset <- 100:500
    cross_correlation_smoothing_k <- 10
    ## Phantom peaks
    PhantomPeak_range <- c(-500, 1500)
    PhPeakbin <- 5
    
    ## step 1.2: Phantom peak and cross-correlation
    if (mc > 1) {
        cluster <- parallel::makeCluster( mc )
    }

    pb$tick()
    message("\nPhantom peak and cross-correlation...")

    phChar <- spp::get.binding.characteristics(data, 
        srange = PhantomPeak_range, 
        bin = PhPeakbin, accept.all.tags = TRUE,
        cluster = cluster)

    if (mc>1) {
        parallel::stopCluster( cluster )
    }

    pb$tick()
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
    
    pb$tick()
    message("\nsmoothing...")
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
        message("Strandshift is adjusted")
    } else {
        message("strandshift not adjusted...")
    }
    pb$tick()
    ## 2.2 phantom peak with smoothing
    message("\nPhantom peak with smoothing...")
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
        message("Crosscorrelation plot saved under ",savePlotPath)
        filename <- file.path(savePlotPath, paste0(tag,"CrossCorrelation.pdf"))
        pdf(file = filename)
    
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
        
        
        dev.off()
        message("pdf saved under", filename)
    }
    
    phantomScores <- list(CC_NSC = round(phScores$NSC, 3), 
        CC_RSC = round(phScores$RSC, 3), 
        CC_QualityFlag = phScores$quality_flag, 
        CC_shift = phScores$peak$x, 
        CC_A = round(phScores$peak$y, 3), 
        CC_B = round(phScores$phantom_cc$y, 3), 
        CC_C = round(phScores$min_cc$y, 3))
    
    pb$tick()
    ## 4 NRF calculation
    message("\nNRF calculation")
    
    ALL_TAGS <- sum(lengths(data$tags))
    UNIQUE_TAGS <- sum(lengths(lapply(data$tags, unique)))
    UNIQUE_TAGS_nostrand <- sum(lengths(lapply(data$tags, FUN = function(x) {
        unique(abs(x))
    })))
    
    NRF <- UNIQUE_TAGS/ALL_TAGS
    NRF_nostrand <- UNIQUE_TAGS_nostrand/ALL_TAGS
    
    ## to compensate for lib size differences
    pb$tick()
    message("\ncalculate different QC values...")
    ## nomi<-rep(names(data$tags), sapply(data$tags, length))
    nomi <- rep(names(data$tags), lapply(data$tags, length))
    
    dataNRF <- unlist(data$tags)
    names(dataNRF) <- NULL
    dataNRF <- paste(nomi, dataNRF, sep = "")
    pb$tick()
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
    pb$tick()
    NRF_LibSizeadjusted <- UNIQUE_TAGS_LibSizeadjusted/1e+07

    STATS_NRF <- list(CC_ALL_TAGS = ALL_TAGS, 
        CC_UNIQUE_TAGS = UNIQUE_TAGS, 
        CC_UNIQUE_TAGS_nostrand = UNIQUE_TAGS_nostrand, 
        CC_NRF = round(NRF,3), 
        CC_NRF_nostrand = round(NRF_nostrand,3), 
        CC_NRF_LibSizeadjusted = round(NRF_LibSizeadjusted,3))
    
    pb$tick()
    
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

    pb$tick()
    
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
    pb$tick()
    return(finalList)
}

