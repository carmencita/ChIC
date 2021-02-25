#' @title Function for downsampling peaks in ChIP-data 
#' to simulate failed antibody
#'
#' @description
#' Function to downsample bam files: reads ChIP and Input 
#' bam files, calles peaks and removes randomly 60percent 
#' of the reads from the peaks. Returns a downsampled ChIP dataframe.
#'
#' downsample_ChIPpeaks
#'
#' @param chip.data data-structure with tag information reads from bam file 
#' (see readBamFile())
#' @param input.data data-structure with tag information reads from bam file 
#' (see readBamFile())
#' @param read_length Integer, length of the reads
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param mc Integer, the number of CPUs for parallelization (default=1)
#' @param debug Boolean, to enter debugging mode. Intermediate files are 
#' saved in working directory
#'
#' @return chip.dataDownSampeld, data-structure with downsampled tags 
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
#' data("chipSubset", package = "ChIC.data", envir = environment())
#' chipBam=chipSubset
#' data("inputSubset", package = "ChIC.data", envir = environment())
#' inputBam=inputSubset
#' chip.dataNew=downsample_ChIPpeaks(chip.data=chipBam, 
#'    input.data=inputBa, read_length=read_length,
#'    annotationID="hg19", mc=mc, debug=FALSE)
#'
#' message(" downsampling from", sum(unlist(lapply(chipBam$tags,length)))," to ",
#' sum(unlist(lapply(chip.dataNew$tags,length))))'
#'}



downsample_ChIPpeaks <- function(chip.data, input.data,read_length, 
    annotationID = "hg19", mc = 1, debug = FALSE) 
{
    start_time <- Sys.time()
    
    ########## check if input format is ok
    if (!(is.list(chip.data) & (length(chip.data) == 2L))) 
        stop("Invalid format for data")
    if (!(is.list(input.data) & (length(input.data) == 2L))) 
        stop("Invalid format for data")
    
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
    
    ## plot and calculate cross correlation and phantom 
    ## characteristics for the ChIP
    message("\ncalculating binding characteristics for ChIP... ")

    estimating_fragment_length_range <- c(0, 500)
    estimating_fragment_length_bin <- 5
    
    #switch cluster on
    if (mc > 1) {
        cluster <- parallel::makeCluster( mc )
    }
    chip_binding.characteristics<-spp::get.binding.characteristics(chip.data, 
        srange = estimating_fragment_length_range, 
        bin = estimating_fragment_length_bin, 
        accept.all.tags = TRUE, cluster = cluster)
    
    #switch cluster off   
    if (mc > 1) {
        parallel::stopCluster( cluster )
    }

    message("\n***Calculating cross correlation QC-metrics for Chip...***")

    crossvalues_Chip <- getCrossCorrelationScores(chip.data, 
        chip_binding.characteristics, 
        read_length = read_length, 
        savePlotPath = NULL, 
        mc = mc,
        annotationID = annotationID)
    ## save the tag.shift
    final.tag.shift <- crossvalues_Chip$tag.shift
    
    ## get chromosome information and order chip and input by it
    chrl_final <- intersect(names(chip.data$tags), names(input.data$tags))
    chip.data$tags <- chip.data$tags[chrl_final]
    chip.data$quality <- chip.data$quality[chrl_final]
    input.data$tags <- input.data$tags[chrl_final]
    input.data$quality <- input.data$quality[chrl_final]
    
    ## remove singular positions with extremely high tag counts with respect 
    ## to the neighbourhood
    message("\nremoving loval tag anomalies...")
    selectedTags <- removeLocalTagAnomalies(chip.data, input.data, 
        chip_binding.characteristics)
    #cleaning up memory space
    remove(chip_binding.characteristics)
    input.dataSelected <- selectedTags$input.dataSelected
    chip.dataSelected <- selectedTags$chip.dataSelected
    
    if ( debug) {
        save(chip.dataSelected, input.dataSelected, 
            file = file.path(getwd(),"dataSelected.RData"))
    }
    
	#chrorder <- paste("chr", c(seq_len(22), "X", "Y"), sep = "")

 	## call broad regions 
    current_window_size <- 1000
    ## pb$tick()

    message("\nBroad regions of enrichment")
    current_zthresh <- 4

    message("checking chromosomes")
    chip.dataSelected  <- f_clearChromStructure(chip.dataSelected,annotationID)
    input.dataSelected <-f_clearChromStructure(input.dataSelected,annotationID)

    message("calling broad peak enrichments... ")
    broad.clusters <- spp::get.broad.enrichment.clusters(chip.dataSelected, 
        input.dataSelected, 
        window.size = current_window_size, 
        z.thr = current_zthresh, 
        tag.shift = final.tag.shift)

    if (debug)
    {
        message("Debugging mode ON")
        filename=file.path(getwd(),"broadEnrichmentCluster.broadPeak")
        spp::write.broadpeak.info(broad.clusters,filename)
    }
   
    ### start end logE znrichment write.broadpeak.info(broad.clusters,
    ### paste('broadRegions', chip.data.samplename,
    ## input,input.data.samplename, 'window', current_window_size, 'zthresh',
    ## current_zthresh,'broadPeak', sep='.'))
    md <- f_convertFormatBroadPeak(broad.clusters)
    broadPeakRangesObject <- f_makeGRangesObject(Chrom = md$chr, 
        Start = md$start, 
        End = md$end)


    message("Starting with downsampling.... ")
	totDownsampled=0
	chip.dataDownSampeld=chip.data

	#iToRemove=lapply(unique(broadPeakRangesObject@seqnames),function(chrID){
	for (chrID in unique(broadPeakRangesObject@seqnames)){
	    removeIndexes=NULL
	    message(chrID)
  	    select=data.frame(broadPeakRangesObject[which(chrID==broadPeakRangesObject@seqnames),])

		for (i in seq(1,nrow(select)))
	    { 
	    	peakToDS=select[i,]
  	    	#print(peakToDS)
  	    	start=peakToDS$start
  	    	end=peakToDS$end

			##Total region to be downsampeled (INDEX)
			frameToDel_indexed=which((abs(chip.data$tags[[chrID]])<=end) & (abs(chip.data$tags[[chrID]])>=start))

	        ##take a random number from 1 to 9, take it as percentage and downsample it
	        #length(frameToDel_indexed) 
	        percentageToDelete=60
	        numbersToDownSample=as.integer(length(frameToDel_indexed)*percentageToDelete/100)
	        totDownsampled=totDownsampled+numbersToDownSample
	        ##deletes final elements
	        removeIndexes=c(removeIndexes,sample(frameToDel_indexed, numbersToDownSample))
    	}
    	#return(removeIndexes)
     	chip.dataDownSampeld$tags[[chrID]]=chip.dataDownSampeld$tags[[chrID]][-c(removeIndexes)]
		 chip.dataDownSampeld$quality[[chrID]]=chip.dataDownSampeld$quality[[chrID]][-c(removeIndexes)]
	 
	}
	message(" downsampling from", sum(unlist(lapply(chip.data$tags,length)))," to ",
		sum(unlist(lapply(chip.dataDownSampeld$tags,length))))

	
	return(chip.dataDownSampeld)

}
