#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for narrow binding PROFILES 	###########################
#######														###########################
#######################################################################################


#'@title Creating cross-correlation profile, phantom peak and calcualting derived QC-metrics
#'
#' @description 
#' We use cross-correlation analysis to obtain QC-metrics proposed for narrow-binding patterns. 
#' After calculating the strand cross-correlation coefficient (Kharchenko et al., 2008), we take the following values from the profile: the coordinates of the 
#' ChIP-peak (fragment length, height A), the coordinates at the phantom-peak (read length, height B) and the baseline (C), the strand-shift, the number of uniquely mapped 
#' reads (unique_tags), uniquely mapped reads corrected by the library size, the number of reads and the read lengths. We  calculate different values using the relative and 
#' absolute height of the cross-correlation peaks: the relative and normalized strand coefficient RSC and NSC  (Landt et al., 2012), and the quality control tag (Fig. 1B) 
#' (Marinov et al., 2013). Other values regarding the library complexity (Landt et al., 2012) like the fraction of non-redundant mapped reads 
#' (NRF; ratio between the number of uniquely mapped reads divided by the total number of reads), the NRF adjusted by library size and ignoring the strand direction 
#' (NRF_nostrand), and the PCR bottleneck coefficient PBC (number of genomic locations to which exactly one unique mapping read maps, divided by the number of unique 
#' mapping reads). 
#'
#' calculateCrossCorrelation
#'
#' @param data spp data structure, structure with tag information from bam file
#' @param binding.characteristics spp data structure, containing binding information for binding preak separation and cross-correlation profile.
#' @param read_length Integer, read length of "data" (Defaul="36")
#' @param savePlotPath, set if Cross-correlation plot should be saved under "savePlotPath". Default=NULL and plot will be shown on screen.
#' @param plotname Name and path were the CrossCorrelation plot (pdf) should be stored (by DEFAULT stored as "phantomCrossCorrelation.pdf" under the working directory)
#'
#' @return finalList, list with different QC-metrics
#' @export
#'
#' @examples
#'\dontrun{
#'	crossvalues_Chip<-calculateCrossCorrelation(chip.data,chip_binding.characteristics,read_length=read_length,savePlotPath=savePlotPath,plotname="ChIP")
#'}

calculateCrossCorrelation=function(data,binding.characteristics,read_length=70,savePlotPath=NULL,plotname="name")
{

	## cross_correlation_customShift_withSmoothing parameters
	cross_correlation_range_subset<-100:500
	cross_correlation_smoothing_k<-10	
	### Phantom peaks
	PhantomPeak_range<-c(-500, 1500)
	PhantomPeak_bin<-5

	###step 1.2: Phantom peak and cross-correlation
	print("Phantom peak and cross-correlation")
	phantom.characteristics<-get.binding.characteristics(data, srange=PhantomPeak_range, bin=PhantomPeak_bin, accept.all.tags=T)
	
	
	ph_peakidx <- which( ( phantom.characteristics$cross.correlation$x >= ( read_length - round(2*PhantomPeak_bin) ) ) & ( phantom.characteristics$cross.correlation$x <= ( read_length + round(1.5*PhantomPeak_bin) ) ) )
	ph_peakidx <- ph_peakidx[ which.max(phantom.characteristics$cross.correlation$y[ph_peakidx]) ]
	phantom_cc <- phantom.characteristics$cross.correlation[ph_peakidx,]

	# Minimum value of cross correlation in srange
	min_cc <- phantom.characteristics$cross.correlation[ which.min(phantom.characteristics$cross.correlation$y) , ]
	# Normalized Strand cross-correlation coefficient (NSC)
	NSC <- phantom.characteristics$peak$y / min_cc$y
	# Relative Strand Cross correlation Coefficient (RSC)
	RSC <- (phantom.characteristics$peak$y - min_cc$y) / (phantom_cc$y - min_cc$y)

	# Quality flag based on RSC value
	qflag <- f_qflag(RSC)
	phantom_peak.scores <- list(phantom_cc=phantom_cc, NSC=NSC, RSC=RSC, quality_flag=qflag, min_cc=min_cc, peak=phantom.characteristics$peak, read_length=read_length)
	
	print("smooting...")
	###2.0 smoothed cross correlation
	subset_selection<- which(binding.characteristics$cross.correlation$x %in% cross_correlation_range_subset)
	binding.characteristics$cross.correlation<-binding.characteristics$cross.correlation[subset_selection,]
	# add smoothing
	binding.characteristics_cross.correlation_y_smoothed<-caTools::runmean(binding.characteristics$cross.correlation$y, k=cross_correlation_smoothing_k)
	# assign the new maximum coordinates
	binding.characteristics$peak$y<-max(binding.characteristics_cross.correlation_y_smoothed)
	binding.characteristics$peak$x<-binding.characteristics$cross.correlation$x[(which(binding.characteristics_cross.correlation_y_smoothed==binding.characteristics$peak$y))]

	# plot cross correlation curve with smoothing

	strandShift<-binding.characteristics$peak$x
	print("Check strandshift...")
	newShift=f_getCustomStrandShift(x=binding.characteristics$cross.correlation$x, y=binding.characteristics_cross.correlation_y_smoothed)
	print(paste("newShift  is ",newShift,sep=""))
	oldShift=NULL
	if (newShift!="ERROR")
	{
		if (newShift!=strandShift)
		{
			oldShift=strandShift
			strandShift=newShift
			print("Strandshift is substituted")
		}
	}else{
		print(paste("strandshift remains the same...",sep=""))
	}
		###2.2 phantom peak with smoothing
	print("Phantom peak with smoothing")
	# phantom.characteristics<-phantom.characteristics
	# select a subset of cross correlation profile where we expect the peak
	subset_selection_forPeakcheck<- which(phantom.characteristics$cross.correlation$x %in% cross_correlation_range_subset)
	phantom.characteristics_cross.correlation_y_smoothed<-caTools::runmean(phantom.characteristics$cross.correlation$y, k=cross_correlation_smoothing_k)
	# assign the new maximum coordinates
	max_y_peakcheck<-max(phantom.characteristics_cross.correlation_y_smoothed[subset_selection_forPeakcheck])
	index_for_maxpeak<-(which(phantom.characteristics_cross.correlation_y_smoothed==max_y_peakcheck))
	## %%% use this additional control just in case there is another cross correlation bin with exactly the ame height outside of the desired search region
	## %%% (e.g. important if the phantom peak is as high or higher than the main cross correlation peak)
	index_for_maxpeak<-index_for_maxpeak[(index_for_maxpeak %in% subset_selection_forPeakcheck)]
	## %% force index 1 just in case there are two points with exactly the same cc value
	max_x_peakcheck<-phantom.characteristics$cross.correlation$x[(index_for_maxpeak[1])]

	if (max_x_peakcheck>phantom.characteristics$peak$x) 
	{
		phantom.characteristics$peak$y<-phantom.characteristics$cross.correlation$y[which(phantom.characteristics$cross.correlation$x==max_x_peakcheck)]
		phantom.characteristics$peak$x<-max_x_peakcheck
		phantom_peak.scores$peak<-phantom.characteristics$peak
		# Normalized Strand cross-correlation coefficient (NSC)
		NSC <- phantom.characteristics$peak$y / phantom_peak.scores$min_cc$y
		phantom_peak.scores$NSC<-NSC
		# Relative Strand Cross correlation Coefficient (RSC)
		RSC <- (phantom.characteristics$peak$y - phantom_peak.scores$min_cc$y) / (phantom_peak.scores$phantom_cc$y - phantom_peak.scores$min_cc$y)
		phantom_peak.scores$RSC<-RSC
		# Quality flag based on RSC value
		qflag <- f_qflag(RSC)
		phantom_peak.scores$quality_flag<-qflag	
	}else{
		phantom.characteristics$peak$y
		#<-max(phantom.characteristics_cross.correlation_y_smoothed)
		phantom.characteristics$peak$x
		#<-phantom.characteristics$cross.correlation$x[(which(phantom.characteristics_cross.correlation_y_smoothed==phantom.characteristics$peak$y))]
	}
	
	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,paste(plotname,"CrossCorrelation.pdf",sep="_"))
		pdf(file=filename)
	}

	# plot cross correlation curve with smoothing
	print("plot cross correlation curve with smoothing")
	par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
	plot(phantom.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation",main="CrossCorrelation Profile")
	lines(x=phantom.characteristics$cross.correlation$x, y=phantom.characteristics_cross.correlation_y_smoothed, lwd=2, col="blue")
	lines(x=rep(phantom_peak.scores$peak$x, times=2), y=c(0,phantom_peak.scores$peak$y), lty=2,lwd=2, col="red")
	lines(x=rep(phantom_peak.scores$phantom_cc$x, times=2), y=c(0,phantom_peak.scores$phantom_cc$y), lty=2,lwd=2, col="orange")
	abline(h=phantom_peak.scores$min_cc$y, lty=2,lwd=2, col="grey")
	text(x=phantom_peak.scores$peak$x, y=phantom_peak.scores$peak$y, labels=paste("A =",signif(phantom_peak.scores$peak$y,3)), col="red", pos=3)
	text(x=phantom_peak.scores$phantom_cc$x, y=phantom_peak.scores$phantom_cc$y, labels=paste("B =",signif(phantom_peak.scores$phantom_cc$y,3)), col="orange", pos=2)
	text(x=min(phantom.characteristics$cross.correlation$x), y=phantom_peak.scores$min_cc$y, labels=paste("C =",signif(phantom_peak.scores$min_cc$y,3)), col="grey", adj=c(0,-1))
	legend(x="topright", legend=c(
		paste("NSC = A/C =", signif(phantom_peak.scores$NSC,3)),
			paste("RSC = (A-C)/(B-C) =", signif(phantom_peak.scores$RSC,3)),
			paste("Quality flag =", phantom_peak.scores$quality_flag),
			"",
			paste("Shift =", (phantom_peak.scores$peak$x)),
			paste("Read length =", (read_length))))
	
	if (!is.null(savePlotPath))
	{
		dev.off()
    	print(paste("pdf saved under",filename,sep=" "))
	}

	phantomScores= list(
		"CC_NSC"=round(phantom_peak.scores$NSC,3),
		"CC_RSC"=round(phantom_peak.scores$RSC,3),
		"CC_QualityFlag"=phantom_peak.scores$quality_flag,
		"CC_shift."=phantom_peak.scores$peak$x,
		#"read_length"=phantom_peak.scores$read_length,
		"CC_A."=round(phantom_peak.scores$peak$y,3),
		"CC_B."=round(phantom_peak.scores$phantom_cc$y,3),
		"CC_C."=round(phantom_peak.scores$min_cc$y,3))
	
	#4 NRF calculation

	print("NRF calculation")
	ALL_TAGS<-sum(sapply(data$tags, length))
	UNIQUE_TAGS<-sum(sapply(lapply(data$tags, unique), length))
	UNIQUE_TAGS_nostrand<-sum(sapply(lapply(data$tags, FUN=function(x) {unique(abs(x))}), length))

	NRF<-UNIQUE_TAGS/ALL_TAGS
	NRF_nostrand<-UNIQUE_TAGS_nostrand/ALL_TAGS

	## to compensate for lib size differences
	print("compensate for lib size differences")
	nomi<-rep(names(data$tags), sapply(data$tags, length))
	dataNRF<-unlist(data$tags)
	names(dataNRF)<-NULL
	dataNRF<-paste(nomi, dataNRF, sep="")

	if (ALL_TAGS > 10e6) {
		UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:100, FUN=function(x) {
			return(length(unique(sample(dataNRF, size=10e6))))
		})))
	} else {
		UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:100, FUN=function(x) {
			return(length(unique(sample(dataNRF, size=10e6, replace=TRUE))))
		})))
	}
	
	NRF_LibSizeadjusted<-UNIQUE_TAGS_LibSizeadjusted/10e6
	STATS_NRF=list("CC_ALL_TAGS"=ALL_TAGS, 
		"CC_UNIQUE_TAGS"=UNIQUE_TAGS, 
		"CC_UNIQUE_TAGS_nostrand"=UNIQUE_TAGS_nostrand, 
		"CC_NRF"=NRF, 
		"CC_NRF_nostrand"=NRF_nostrand, 
		"CC_NRF_LibSizeadjusted"=NRF_LibSizeadjusted)

	#N1= number of genomic locations to which EXACTLY one unique mapping read maps
	#Nd = the number of genomic locations to which AT LEAST one unique mapping read maps, i.e. the number of non-redundant, unique mapping reads
	N1<-sum(sapply(data$tags, FUN=function(x) {
		checkDuplicate<-duplicated(x)
		duplicated_positions<-unique(x[checkDuplicate])
		OUT<-x[!(x %in% duplicated_positions)]
		return(length(OUT))
	}))
	Nd<-sum(sapply(lapply(data$tags, unique), length))
	PBC = N1/Nd
	tag.shift <- round(strandShift/2)
		
	finalList <- append(append(
		list("CC_StrandShift"=strandShift,
			"tag.shift"=tag.shift,
			"N1"=round(N1,3),
			"Nd"=round(Nd,3),
			"CC_PBC"=round(PBC,3),
			"CC_readLength"=read_length,
			"CC_UNIQUE_TAGS_LibSizeadjusted"=UNIQUE_TAGS_LibSizeadjusted),
		phantomScores),STATS_NRF)
	return(finalList)
}





#'@title Calculating QC-values from peak calling procedure
#'
#' @description
#' Other measures we include in our analysis are the fraction of usable reads in the peak regions (FRiP) (Landt et al., 2012), for which the function calls sharp- and 
#' broad-binding peaks to obtain two types: the FRiP_sharpsPeak and the FRiP_broadPeak. The function takes the number of called of peaks using an FDR of 0.01 and an 
#' evalue of 10 (Kharchenko et al., 2008). And count the number of peaks called when using the sharp- and broad-binding option. 
#' getBindingRegionsScores
#'
#' @param chip spp data structure, structure with tag information from ChIP file
#' @param input spp data structure, structure with tag information from Input file
#' @param chip.dataSelected,
#' @param input.dataSelected,
#' @param tag.shift integer, value from calculateCrossCorrelation list
#' @param chrorder, chromosome order (default=NULL) 
#' @return QCscoreList, containing 6 QC-values
#'
#' @export
#'
#' @examples
#'\{dontrun
#'	bindingScores=getBindingRegionsScores(chip.data,input.data, chip.dataSelected,input.dataSelected,final.tag.shift)
#'}



getBindingRegionsScores=function(chip,input,chip.dataSelected,input.dataSelected,tag.shift=75,chrorder=NULL)
{

	#load("Settings.RData")
	##for simplicity we use currently a shorter chromosome frame to avoid problems with M
	chrorder<-paste("chr", c(1:19, "X","Y"), sep="")
	#custom_chrorder<-paste("chr", c(1:19, "X","Y"), sep="")
	#custom_chrorder<-paste("chr", c(1:22, "X","Y"), sep="")


	###5 broadRegions
	###6 enrichment broad regions
	#zthresh_list<-c(3,4)
	current_window_size<-1000
	print("Broad regions of enrichment")
	#for (current_window_size in window_sizes_list) 
	#{
	#for (current_zthresh in zthresh_list) 
	#{
	current_zthresh=4
	broad.clusters <- get.broad.enrichment.clusters(chip.dataSelected, input.dataSelected, window.size=current_window_size, z.thr=current_zthresh, tag.shift=tag.shift) ##start end logEnrichment
	#write.broadpeak.info(broad.clusters, paste("broadRegions", chip.data.samplename,"input",input.data.samplename, "window", current_window_size, "zthresh", current_zthresh,"broadPeak", sep="."))
	md=f_convertFormatBroadPeak(broad.clusters)
	broadPeak_selected_genes_genomeIntervals_object<-new("Genome_intervals",
		.Data=(cbind(
			as.integer(as.character(md[,2])), 
			as.integer(as.character(md[,3]))) ), 
			closed=TRUE,
		annotation = data.frame(
			seq_name = factor(md[,1]), 
				inter_base=FALSE, 
				start = as.integer(as.character(md[,2])),
				end =as.integer(as.character(md[,3]))))
		#}
	#}

	###12 get binding sites with FDR and eval
	chip.data12<-chip.dataSelected[(names(chip.dataSelected) %in% chrorder)]
	input.data12<-input.dataSelected[(names(input.dataSelected) %in% chrorder)]

	print("Binding sites detection fdr")
	fdr <- 1e-2
	detection.window.halfsize <- tag.shift
	print("Window Tag Density method - WTD")
	bp_FDR <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,fdr=fdr,whs=detection.window.halfsize*2, tag.count.whs=detection.window.halfsize, cluster=NULL)
	FDR_detect=	sum(unlist(lapply(bp_FDR$npl,function(d) length(d$x))))
		
	print("Binding sites detection evalue")
	eval<-10
	bp_eval <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,e.value=eval,whs=detection.window.halfsize*2,cluster=NULL)
	eval_detect=sum(unlist(lapply(bp_eval$npl,function(d) length(d$x))))


	# output detected binding positions
	#output.binding.results(results=bp,filename=paste("WTC.binding.positions.evalue", chip.data.samplename,"input",input.data.samplename, "txt", sep="."))
	if (length(bp_eval$npl)>1)
	{
		###14 get precise binding position using escore LARGE PEAKS
		bp_broadpeak <- add.broad.peak.regions(chip.data12,input.data12,bp_eval, window.size=1000, z.thr=3)
		md=f_converNarrowPeakFormat(bp_broadpeak)
		sharpPeak_selected_genes_genomeIntervals_object<-new("Genome_intervals",
			.Data=(cbind(
				as.integer(as.character(md[,2])), 
				as.integer(as.character(md[,3]))) ), 
				closed=TRUE,
			annotation = data.frame(
				seq_name = factor(md[,1]), 
					inter_base=FALSE, 
					start = as.integer(as.character(md[,2])),
					end =as.integer(as.character(md[,3]))))

		###13 is only looping and writing into file
		###15 FRiP
		chip.test<-lapply(chip$tags, FUN=function(x) {x+tag.shift}) ###SORTS the tags for each chrom
		TOTAL_reads<-sum(sapply(chip.test, length))

		##Frip broad binding sites (histones)
		broadPeak_selected_genes_genomeIntervals_object<-close_intervals(interval_union(broadPeak_selected_genes_genomeIntervals_object))
		regions_data_list<-split(as.data.frame(broadPeak_selected_genes_genomeIntervals_object), f=seqnames(broadPeak_selected_genes_genomeIntervals_object))
		chrl<-names(regions_data_list)
		names(chrl)<-chrl
		outcountsBroadPeak<-sum(sapply(chrl, FUN=function(chr) {
			sum(spp:::points.within(x=abs(chip.test[[chr]]), fs=((regions_data_list[[chr]])[,1]), fe=((regions_data_list[[chr]])[,2]), return.point.counts = TRUE))
	    }))
		FRiP_broadPeak<-outcountsBroadPeak/TOTAL_reads

		###Frip sharp peaks 14
		sharpPeak_selected_genes_genomeIntervals_object<-close_intervals(interval_union(sharpPeak_selected_genes_genomeIntervals_object))
		regions_data_list<-split(as.data.frame(sharpPeak_selected_genes_genomeIntervals_object), f=seqnames(sharpPeak_selected_genes_genomeIntervals_object))
		chrl<-names(regions_data_list)
		names(chrl)<-chrl
		outcountsSharpPeak<-sum(sapply(chrl, FUN=function(chr) {
			sum(spp:::points.within(x=abs(chip.test[[chr]]), fs=((regions_data_list[[chr]])[,1]), fe=((regions_data_list[[chr]])[,2]), return.point.counts = TRUE))
		}))
		FRiP_sharpPeak<-outcountsSharpPeak/TOTAL_reads
	}else{
		TOTAL_reads=0
		FRiP_broadPeak=0
		outcountsBroadPeak=0
		FRiP_sharpPeak=0
		outcountsSharpPeak=0
	}

	QCscoreList=list("CC_FDRpeaks"=round(FDR_detect,3),
		"CC_evalpeaks"=round(eval_detect,3),
		"CC_FRiP_broadPeak"=round(FRiP_broadPeak,3),  
		"CC_FRiP_sharpPeak"=round(FRiP_sharpPeak,3), 
		"CC_outcountsBroadPeak"=outcountsBroadPeak,
		"CC_outcountsSharpPeak"= outcountsSharpPeak)

	return(QCscoreList)
}



#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for global read distribution 	###########################
#######														###########################
#######################################################################################

#'@title Metrics taken from global read distribution
#'
#' @description
#' This set of values is based on the global read distribution along the genome for immunoprecipitation and input data (Diaz et al., 2012). 
#' The genome is binned and the read coverage counted for each bin. Then the function computes the cumulative distribution of reads density per genomic bin and plots the 
#' fraction of the coverage on the y-axis and the fraction of bins on the x-axis. Then different values can be sampled from the cumulative distribution: 
#' like the fraction of bins without reads for in immunoprecipitation and input,the point of the maximum distance between the ChIP and the input (x-axis, y-axis for 
#' immunoprecipitation and input, distance (as absolute difference), the sign of the differences), the fraction of reads in the top 1%bin for immunoprecipitation and input. 
#' Finally, the funciton returns 9 QC-measures
#'
#' QCscores_global
#'
#' @param densityChip Smoothed tag density object for the ChIP (returned by f_CrossCorrelation). 
#' @param densityInput Smoothed tag density object for the Inpt (returned by f_CrossCorrelation)
#' @param savePlotPath Path, path set forces the fingerprint plot to be saved under the respective directory. Default=NULL, plot not saved but shown on screen
#' @param debug Boolean, to enter debugging mode (default= FALSE)
#'
#' @return finalList, list containing 9 QC-values
#'
#' @export
#'
#' @examples
#'\dontrun{
#' CC_Result=crossCorrelation(chipName=chipName,
#' inputName=inputName, read_length=36, 
#' dataPath=dataDirectory, annotationID="hg19",savePlotPath=getwd())
#' tag.shift=CC_Result$QCscores_ChIP$tag.shift
#' smoothedDensityInput=CC_Result$TagDensityInput
#' smoothedDensityChip=CC_Result$TagDensityChip
#' Ch_Results=QCscores_global(densityChip=smoothedDensityChip,
#' densityInput=smoothedDensityInput)
#' }

QCscores_global=function(densityChip,densityInput,savePlotPath=NULL,debug=FALSE)
{

	#shorten frame and reduce resolution
	print("shorten frame")
	chip.smoothed.density=f_shortenFrame(densityChip)
	input.smoothed.density=f_shortenFrame(densityInput)
	
	chip<-unlist(lapply(chip.smoothed.density, FUN=function(list_element) {return(list_element$y)}))
	input<-unlist(lapply(input.smoothed.density, FUN=function(list_element) {return(list_element$y)}))

	##create cumulative function for chip and input
	print("Calculate cumsum")
	cumSumChip=f_sortAndBinning(chip)
	cumSumInput=f_sortAndBinning(input)

	##caluclate QCvalues for chip
	chipFractionOfReadsIntop1percentBins<-(1-cumSumChip[(which(cumSumChip$x>=0.99)[1]),"pj"])
	chipFractionOfBinsWithoutReads<-cumSumChip$x[(which(cumSumChip$pj>0)[1])]

	##caluclate QCvalues for input
	inputFractionOfReadsIntop1percentBins<-(1-cumSumInput[(which(cumSumInput$x>=0.99)[1]),"pj"])
	inputFractionOfBinsWithoutReads<-cumSumInput$x[(which(cumSumInput$pj>0)[1])]

	##get the point of maximum distance between the two functions and the sign of the distance
	sign_sign=sign(max(cumSumInput$pj-cumSumChip$pj)) ##input-chip
	arrowx=cumSumChip[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$x

	#get the distance between the curves
	y_input=round(cumSumInput[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$pj,3)
	y_chip=round(cumSumChip[which.max(abs(cumSumInput$pj-cumSumChip$pj)),]$pj,3)
	dist=abs(y_input-y_chip)
	##prepare list to be returned
	finalList=list("Ch_X.axis"=round(arrowx,3),
		"Ch_Y.Input"=y_input,
		"Ch_Y.Chip"=y_chip,
		"Ch_sign_chipVSinput"=sign_sign,
		"Ch_FractionReadsTopbins_chip"=round(chipFractionOfReadsIntop1percentBins,3),
		"Ch_FractionReadsTopbins_input"=round(inputFractionOfReadsIntop1percentBins,3),
		"Ch_Fractions_without_reads_chip"=round(chipFractionOfBinsWithoutReads,3),
		"Ch_Fractions_without_reads_input"=round(inputFractionOfBinsWithoutReads,3),
		"Ch_DistanceInputChip"=dist)
		#"CrossPoint_X"=cross_pointX,
		#"CrossPoint_Y_chip"=cross_pointY_chip,
		#"CrossPoint_Y_input"=cross_pointY_input)


	#create Fingerprint plot
	f_fingerPrintPlot(cumSumChip,cumSumInput,savePlotPath=savePlotPath)
	
	if (debug)
	{
		write.table(finalList,file=file.path(getwd(),"Chance.results"))
	}
	#return QC values
	return(finalList)
}



#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for local ENRICHMENT 			###########################
#######														###########################
#######################################################################################


#'@title Wrapper function to create scaled and non-scaled metageneprofiles 
#'
#' @description
#' Metagene plots show the signal enrichment around a region of interest like the TSS or over a predefined set of genes.
#' The tag density of the immunoprecipitation is taken over all RefSeg annotated human genes, averaged and log2 transformed. The same is done for the 
#' input. The normalized profile is calculated as the signal enrichment (immunoprecipitation over the input). 
#' We created two types of metagene profiles: a non-scaled profile for the TSS  and TES, and a scaled profile 
#' for the entire gene, including the gene body. 
#' The non-scaled profile is constructed around the TSS/TES, with 2KB up- and downstream regions respectively. 
#' 
#' CreateMetageneProfile
#'
#' @param smoothed.densityChip Smoothed tag density object for the ChIP (returned by f_CrossCorrelation)
#' @param smoothed.densityInput, Smoothed tag density object for the input(returned by f_CrossCorrelation)
#' @param tag.shift Integer, tag shift returned by f_CrossCorrelation()
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @return list with 3 objects: scaled profile ("twopoint"), non-scaled profile for TSS (TSS) and TES (TES). Each object is made of a list containing 
#' the chip and the input profile
#'
#' @export
#'
#' @examples
#'\{dontrun
#' Meta_Result=createMetageneProfile(smoothedDensityChip,smoothedDensityInput,tag.shift,annotationID="hg19")
#'}

createMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,annotationID="hg19",debug=FALSE,mc=1)
{
  annotation=paste(annotationID,"RefSeqAllGenesFiltered.RData",sep=".")
  geneAnnotations_file<-annotation
  print("Load geneannotation")
  load(geneAnnotations_file) #RefSeqGenesAll_object

  #data(paste(annotationID,"RefSeqAllGenesFiltered",sep="."))

  current_annotations_type<-gsub(pattern=".RData", replacement="", fixed=TRUE, x=basename(geneAnnotations_file))
  current_annotations_object=RefSeqGenes_annotated_filteredByOverlap_geneLength
  # format annotations as chromosome lists
  current_annotations_object<-data.frame(current_annotations_object@.Data, annotation(current_annotations_object), stringsAsFactors=FALSE)
  current_annotations_object$interval_starts<-as.integer(current_annotations_object$interval_starts)
  current_annotations_object$interval_ends<-as.integer(current_annotations_object$interval_ends)
  current_annotations_object$seq_name<-as.character(current_annotations_object$seq_name)
  annotatedGenesPerChr <-split(current_annotations_object, f=current_annotations_object$seq_name)

  ##two.point.scaling
  #create scaled metageneprofile
  #input
  print("Calculating scaled metageneprofile ...")
  smoothed.densityInput=list(td=smoothed.densityInput)
  print("process input")
	

  binned_Input = f_t.get.gene.av.density(smoothed.densityInput,gdl=annotatedGenesPerChr,mc=mc)
  #Chip
  smoothed.densityChip=list(td=smoothed.densityChip)
  print("process ChIP")
  binned_Chip= f_t.get.gene.av.density(smoothed.densityChip,gdl=annotatedGenesPerChr,mc=mc)

  twopoint=list(chip=binned_Chip,input= binned_Input)

  ##one.point.scaling
  #create non-scaled metageneprofile for TSS
  print("Creating non-scaled metageneprofiles...")

  print("...TSS")
  binnedInput_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityInput,gdl=annotatedGenesPerChr,mc=mc)
  binnedChip_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityChip,gdl=annotatedGenesPerChr,mc=mc)

  onepointTSS=list(chip=binnedChip_TSS,input= binnedInput_TSS)
  
  ##one.point.scaling
  #create non-scaled metageneprofile for TES
  print("...TES")
  binnedInput_TES <- f_t.get.gene.av.density_TES(smoothed.densityInput,gdl=annotatedGenesPerChr,mc=mc)
  binnedChip_TES <- f_t.get.gene.av.density_TES(smoothed.densityChip,gdl=annotatedGenesPerChr,mc=mc)

  onepointTES=list(chip=binnedChip_TES,input= binnedInput_TES)
  
  if (debug)
  {
    save(binned_Chip, binned_Input,file=file.path(getwd(), paste(chipName,inputName,"Twopoint.RData",sep="_")))
    save(binnedChip_TSS, binnedInput_TSS,file=file.path(getwd(), paste(chipName,inputName,"OnePointTSS.RData",sep="_")))
    save(binnedChip_TES, binnedInput_TES,file=file.path(getwd(), paste(chipName,inputName,"OnePointTES.RData",sep="_")))
  }

  return(list("twopoint"=twopoint,"TSS"=onepointTSS,"TES"=onepointTES))
}

