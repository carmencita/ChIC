
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
#' @param annotationID String indicating the genome assembly (Default="hg19")
#' @param chrorder chromosome order (default=NULL) 
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
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
#'     annotationID="hg19",
#'     tag.shift = finalTagShift, mc = mc)
#'}


getPeakCallingScores <- function(chip, input, chip.dataSelected, 
    input.dataSelected, annotationID="hg19",
    tag.shift = 75, mc=1, chrorder = NULL,
    debug = FALSE) 
{
    ## pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 6, 
    ##     clear = FALSE, width = 60)
    ## pb$tick()

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
    if (annotationID=="dm3")
    {   
        chrorder=c("chr2L","chr2R","chr3L","chr3R","chr4")
        ##chrorder <- names(ChIC.data::dm3_chrom_info)
        ##without M!! remove M from

    }else{
        #chrorder <- paste("chr", c(seq_len(19), "X", "Y"), sep = "")
        chrorder <- paste("chr", c(seq_len(22)), sep = "")
    }
    
    ## 5 broadRegions 6 enrichment broad regions zthresh_list<-c(3,4)
    current_window_size <- 1000
    ## pb$tick()

    message("\nBroad regions of enrichment")
    ## for (current_window_size in window_sizes_list) { for (current_zthresh in
    ## zthresh_list) {
    current_zthresh <- 4

    message("checking chromosomes")
    chip.dataSelected  <- f_clearChromStructure(chip.dataSelected,annotationID)
    input.dataSelected <-f_clearChromStructure(input.dataSelected,annotationID)


    broad.clusters <- spp::get.broad.enrichment.clusters(chip.dataSelected, 
        input.dataSelected, 
        window.size = current_window_size, 
        z.thr = current_zthresh, 
        tag.shift = tag.shift)

    if (debug)
    {
        message("Debugging mode ON")
        filename=file.path(getwd(),"broadEnrichmentCluster.broadPeak")
        spp::write.broadpeak.info(broad.clusters,filename)
    }
    ## pb$tick()

    ### start end logE znrichment write.broadpeak.info(broad.clusters,
    ### paste('broadRegions', chip.data.samplename,
    ## input,input.data.samplename, 'window', current_window_size, 'zthresh',
    ## current_zthresh,'broadPeak', sep='.'))
    md <- f_convertFormatBroadPeak(broad.clusters)
    broadPeakRangesObject <- f_makeGRangesObject(Chrom = md$chr, 
        Start = md$start, 
        End = md$end)
    
    ## 12 get binding sites with FDR and eval
    chip.data12 <- chip.dataSelected[(names(chip.dataSelected) %in% chrorder)]
    input.data12<-input.dataSelected[(names(input.dataSelected) %in% chrorder)]

    if (mc > 1) {
        cluster <- makeCluster( mc )
    }
    ## pb$tick()

    message("\nBinding sites detection fdr")
    fdr <- 0.01
    detection.window.halfsize <- tag.shift
    message("Window Tag Density method - WTD")
    bp_FDR <- spp::find.binding.positions(signal.data = chip.data12, 
        control.data = input.data12, 
        fdr = fdr, whs = detection.window.halfsize * 2, 
        tag.count.whs = detection.window.halfsize, 
        cluster = NULL)
    FDR_detect <- sum(unlist(lapply(bp_FDR$npl, function(d) length(d$x))))
    
    ## pb$tick()
    message("\nBinding sites detection evalue")
    eval <- 10
    bp_eval <- spp::find.binding.positions(signal.data = chip.data12, 
        control.data = input.data12, 
        e.value = eval, whs = detection.window.halfsize * 2, cluster = NULL)
    eval_detect <- sum(unlist(lapply(bp_eval$npl, function(d) length(d$x))))
    
    if (mc > 1) {
        stopCluster( cluster )
    }

    ## pb$tick()
    if ( debug )
    {
        ## output detected binding positions
        message("writing detected binding positions to workingdirectory")
        spp::output.binding.results(bp_FDR, 
            file.path(getwd(),"FDR_bindingPositions.txt"))
        spp::output.binding.results(bp_eval, 
            file.path(getwd(),"eval_bindingPositions.txt"))
        save(bp_eval, file = file.path(getwd(), 
            paste("bp_eval.RData", sep = "_")))
    }
    
    ## output detected binding positions output.binding.results(results=bp,
    ## filename=paste('WTC.binding.positions.evalue',
    ## chip.data.samplename,'input',input.data.samplename, 'txt', sep='.'))
    if (length(bp_eval$npl) > 1) {
        pb <- progress_bar$new(format = "(:spin) [:bar] :percent",total = 5, 
            clear = FALSE, width = 60)
        ## pb$tick()
        ## 14 get precise binding position using escore LARGE PEAKS
        bp_broadpeak <- spp::add.broad.peak.regions(chip.data12, 
            input.data12, 
            bp_eval, 
            window.size = 1000, z.thr = 3)
        
        if (debug)
        {
            message("writing output in narrowpeak format...")
            ## output using narrowpeak format
            spp::write.narrowpeak.binding(bp_broadpeak, 
                file.path(getwd(),".peaks.narrowPeak"))
        }

        ## pb$tick()
        md <- f_converNarrowPeakFormat(bp_broadpeak)
        sharpPeakRangesObject <- f_makeGRangesObject(Chrom = md[, 1], 
            Start = md[, 2], End = md[, 3])
        
        ## FRiP
        chip.test <- lapply(chip$tags, FUN = function(x) {
            x + tag.shift
        })
        TOTAL_reads <- sum(lengths(chip.test))
        
        ## Frip broad binding sites (histones)
        broadPeakRangesObject<-f_reduceOverlappingRegions(broadPeakRangesObject)

        regions_data_list <- split(as.data.frame(broadPeakRangesObject), 
            f = seqnames(broadPeakRangesObject))
        
        if ( debug )
        {        
            if ( file.exists( file.path (getwd(),"broadPeakRanges.bed" )))
            {
                file.remove( file.path (getwd(),"broadPeakRanges.bed" ))
            }
            list=lapply(regions_data_list,function(chr){
                #print(chr)
                text <- cbind(as.character(chr$seqnames), 
                    as.numeric(chr$start),
                    as.numeric(chr$end))
                write.table(text, 
                    file = file.path(getwd(),"broadPeakRanges.bed"),
                append =TRUE, 
                quote = FALSE, 
                row.names = FALSE,
                col.names = FALSE)
            })
        }

        chrl <- names(regions_data_list)
        names(chrl) <- chrl

        ## pb$tick()



        outcountsBroadPeak <- sum(unlist(lapply(chrl, FUN = function(chr) {
            sum(spp::points_within(x = abs(chip.test[[chr]]), 
                fs = regions_data_list[[chr]]$start, 
                fe = regions_data_list[[chr]]$end, 
                return.point.counts = TRUE))
        })))
        
        FRiP_broadPeak <- outcountsBroadPeak/TOTAL_reads



        ## Frip sharp peaks 14

        sharpPeakRangesObject<-f_reduceOverlappingRegions(sharpPeakRangesObject)
    
        regions_data_list <- split(as.data.frame(sharpPeakRangesObject), 
            f = seqnames(sharpPeakRangesObject))

        if ( debug )
        {        
            if ( file.exists( file.path (getwd(),"sharpPeakRanges.bed" )))
            {
                file.remove( file.path (getwd(),"sharpPeakRanges.bed" ))
            }
            list <- lapply(regions_data_list, function(chr){
                #print(chr)
                text <- cbind(as.character(chr$seqnames), 
                    as.numeric(chr$start),
                    as.numeric(chr$end))

                write.table(text, 
                    file = file.path(getwd(),"sharpPeakRanges.bed"), 
                    append =TRUE, 
                    quote = FALSE, 
                    row.names = FALSE,
                    col.names = FALSE)
            })
        }

        ## pb$tick()

        chrl <- names(regions_data_list)
        names(chrl) <- chrl
        

        outcountsSharpPeak <- sum(unlist(lapply(chrl, FUN = function(chr) {
            sum(spp::points_within(x = abs(chip.test[[chr]]), 
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

    ## pb$tick()

    QCscoreList <- list(CC_FDRpeaks = round(FDR_detect, 3), 
        CC_evalpeaks = round(eval_detect, 3), 
        CC_FRiP_broadPeak = round(FRiP_broadPeak, 3), 
        CC_FRiP_sharpPeak = round(FRiP_sharpPeak, 3), 
        CC_outcountsBroadPeak = outcountsBroadPeak, 
        CC_outcountsSharpPeak = outcountsSharpPeak)
    
    return(QCscoreList)
}



