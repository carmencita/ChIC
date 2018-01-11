##private variables


###
### cross_correlation parameters
###
estimating_fragment_length_range<-c(0,500)
estimating_fragment_length_bin<-5


###
### Phantom peaks
###

PhantomPeak_range<-c(-500, 1500)
PhantomPeak_bin<-5
#PhantomPeak_exclude.min<-10


###
### cross_correlation_customShift_withSmoothing parameters
###
cross_correlation_range_subset<-100:500
cross_correlation_smoothing_k<-10

##Smoothing parameters
smoothingBandwidth<-50 ##is slidingwindow shift 
smoothingStep<-20	##is size of sw ##before it was 10


###
### broadRegionsOfEnrichment_loop parameters
###
window_sizes_list<-1000 #c(1000, 500, 200)[1]
remove.local.tag.anomalies_filter<-TRUE
select.informative.tags_filter<-FALSE
####

f_convertFormatBroadPeak <- function(given.clusters)
{
	chrl <- names(given.clusters)
	names(chrl) <- chrl
	chrl <- chrl[unlist(lapply(given.clusters, function(d) length(d$s))) > 0]
	md <- do.call(rbind, lapply(chrl, function(chr)		
	{ 
		df <- given.clusters[[chr]]
		cbind(chr, df$s, df$e, ".", "0", ".", df$rv, -1, -1)
	}
	))
	md <- md[order(as.numeric(md[, 7]), decreasing = T), ]
	md=data.frame(md)
	return(md)
}


f_converNarrowPeakFormat =function(bd, margin = bd$whs) 
{
	if (is.null(margin)) 
	{
		margin <- 50
	}
	chrl <- names(bd$npl)
	names(chrl) <- chrl
	md <- do.call(rbind, lapply(chrl, function(chr) 
	{
		df <- bd$npl[[chr]]
		x <- df$x
		rs <- df$rs
		if (is.null(rs)) {rs <- rep(NA, length(x))}
			re <- df$re
		if (is.null(re)) {re <- rep(NA, length(x))}
			ivi <- which(is.na(rs))
		if (any(ivi)) {rs[ivi] <- x[ivi] - margin}
			ivi <- which(is.na(re))
		if (any(ivi)) {re[ivi] <- x[ivi] + margin}
			cbind(chr, rs, re, ".", "0", ".", df$y, -1, format(df$fdr, scientific = T, digits = 3), x - rs)
	}))
	md <- md[order(as.numeric(md[, 7]), decreasing = T), ]
	md=data.frame(md)
	#write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
	return(md)
}


f_getCustomStrandShift= function(x,y){
	#x=binding.characteristics$cross.correlation$x
	#y=binding.characteristics_cross.correlation_y_smoothed
	deriv=diff(y)/diff(x)
	#globalminY=min(abs(deriv))
	deriv=append(0,deriv) ##to catch up the right index, because in deriv I loose one array-field
	#globalminX=x[which(abs(deriv)==globalminY)] ###shift all um 5 nach rechts
	#check from the right side and pick the points with the x (closest to zero) and y largest 
	##means: search for regions with Vorzeichen change
	startVorzeichen=sign(deriv[length(deriv)])
	field=NULL
	for (index in rev(seq(2,length(deriv))))
	{
		derivpoint=deriv[index]
		xpoint=x[index]
		ypoint=y[index]
		#print(paste(xpoint,ypoint,sep=" "))
		if (startVorzeichen!=sign(derivpoint))
		{ 
			##Vorzeichenwechsel
			field=rbind(field,c(xpoint,ypoint,derivpoint))
			startVorzeichen=sign(derivpoint)
		}
	}
	if (is.null(field) )
	{
		newShift="ERROR"
	}else{
		field=data.frame(field)
		colnames(field)=c("x","y","deriv")
		newShift=field[which(max(field$y)==field$y),]$x
	}
	return(newShift)
}


f_readFile=function(f_filename,f_reads.aligner.type="bam",f_path=getwd())
{
	read.tags.current_function<-get(paste("read", f_reads.aligner.type , "tags", sep="."))
	if (f_reads.aligner.type=="bam")
	{
	    data<-read.tags.current_function(file.path(f_path,paste(f_filename,".bam",sep="")))
	}
	if (f_reads.aligner.type=="tagalign")
	{
	    data<-read.tags.current_function(file.path(f_path,paste(f_filename,".tagAlign",sep="")))
	}
	return(data)
}

f_qflag=function(RSC)
{
	qflag=NA
	if ( (RSC >= 0) & (RSC < 0.25) ){
		qflag <- -2
	} else if ( (RSC >= 0.25) & (RSC < 0.5) ) {
		qflag <- -1
	} else if ( (RSC >= 0.5) & (RSC < 1) ) {
		qflag <- 0
	} else if ( (RSC >= 1) & (RSC < 1.5) ) {
		qflag <- 1
	} else if ( (RSC >= 1.5) ) {
		qflag <- 2
	}
	return(qflag)
}

f_calculateCrossCorrelation=function(data,binding.characteristics,read_length=36,plotname="phantomCrossCorrelation",TFname="NA")
{
	#readCount=sum(sapply(data$tags, length))
	###step 1.2: Phantom peak and cross-correlation
	print("Phantom peak and cross-correlation")
	phantom.characteristics<-get.binding.characteristics(data, srange=PhantomPeak_range, bin=PhantomPeak_bin, cluster=cluster)
	
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
	#f_plot(binding.characteristics$cross.correlation,datafilename,maintitel=TFname,plotname="customCrossCorrelation",xlab="strand shift",ylab="cross-correlation",line=binding.characteristics$peak$x, lineplotX=binding.characteristics$cross.correlation$x,lineplotY=binding.characteristics_cross.correlation_y_smoothed)
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
	print("Phantom peak with smooting")
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

	# plot cross correlation curve with smoothing
	print("plot cross correlation curve with smoothing")
	pdf((filename=file.path(path, paste(plotname, "pdf", sep="."))))
		par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
		plot(phantom.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation",main=TFname)
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
		title(TFname)
	dev.off()

	phantomScores= list(
		NSC=round(phantom_peak.scores$NSC,3),
		RSC=round(phantom_peak.scores$RSC,3),
		quality_flag=phantom_peak.scores$quality_flag,
		shift=phantom_peak.scores$peak$x,
		read_length=phantom_peak.scores$read_length,
		A=round(phantom_peak.scores$peak$y,3),
		B=round(phantom_peak.scores$phantom_cc$y,3),
		C=round(phantom_peak.scores$min_cc$y,3))
	
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
	STATS_NRF=list(ALL_TAGS=ALL_TAGS, 
		UNIQUE_TAGS=UNIQUE_TAGS, 
		UNIQUE_TAGS_nostrand=UNIQUE_TAGS_nostrand, 
		NRF=NRF, NRF_nostrand=NRF_nostrand, NRF_LibSizeadjusted=NRF_LibSizeadjusted)
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
		list("strandShift"=strandShift,
			"tag.shift"=tag.shift,
			"N1"=round(N1,3),"Nd"=round(Nd,3),"PBC"=round(PBC,3),
			"read_length"=read_length,
			"UNIQUE_TAGS_LibSizeadjusted"=UNIQUE_TAGS_LibSizeadjusted),
		phantomScores),STATS_NRF)
	return(finalList)
}



f_selectInformativeTag=function(chip,input,chip_b.characteristics,input_b.characteristics)
{
	print("Filter tags")

	if (select.informative.tags_filter) {
	      print("select.informative.tags filter")
	     #load(paste("sppdata", "binding", chip.data.samplename, "RData", sep="."))
	      chip.dataSelected <- select.informative.tags(chip, chip_b.characteristics)
	      #load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
	      input.dataSelected <- select.informative.tags(input.data, input_b.characteristics)
	} else {
	      print("SKIP select.informative.tags filter")
	      chip.dataSelected<-chip$tags
	      input.dataSelected<-input$tags
	}
	if (remove.local.tag.anomalies_filter) {
		print("remove.local.tag.anomalies filter")
	      # restrict or remove singular positions with very high tag counts
		chip.dataSelected <- remove.local.tag.anomalies(chip.dataSelected)
		input.dataSelected <- remove.local.tag.anomalies(input.dataSelected)
	} else {
		print("SKIP remove.local.tag.anomalies filter")
		input.dataSelected<-input$tags
		chip.dataSelected=chip$tags
	}

	finalList=list("input.dataSelected"=input.dataSelected,
		"chip.dataSelected"=chip.dataSelected)
	return(finalList)
}

#> bindingAnalysis=f_getBindingRegionsScores(chip.data,input.data,chip.dataSelected,input.dataSelected,final.tag.shift,custom_chrorder)
f_getBindingRegionsScores=function(chip,input,chip.dataSelected,input.dataSelected,tag.shift=75)#,chrorder=NULL)
{
	chrorder<-paste("chr", c(1:19, "X","Y"), sep="")
	#custom_chrorder<-paste("chr", c(1:19, "X","Y"), sep="")
	#custom_chrorder<-paste("chr", c(1:22, "X","Y"), sep="")

	###5 broadRegions
	###6 enrichment broad regions
	zthresh_list<-c(3,4)
	extension_vector_list<-list()
	print("Broad regions of enrichment")
	for (current_window_size in window_sizes_list) 
	{
		print(paste("Window size", current_window_size))
		for (current_zthresh in zthresh_list) 
		{
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
		}
	}

	###12 get binding sites with FDR and eval
	chip.data12<-chip.dataSelected[(names(chip.dataSelected) %in% chrorder)]
	input.data12<-input.dataSelected[(names(input.dataSelected) %in% chrorder)]

	print("Binding sites detection fdr")
	fdr <- 1e-2
	detection.window.halfsize <- tag.shift
	print("Window Tag Density method - WTD")
	bp_FDR <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,fdr=fdr,whs=detection.window.halfsize*2, tag.count.whs=detection.window.halfsize, cluster=cluster)
	FDR_detect=	sum(unlist(lapply(bp_FDR$npl,function(d) length(d$x))))
		
	print("Binding sites detection evalue")
	eval<-10
	bp_eval <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,e.value=eval,whs=detection.window.halfsize*2,cluster=cluster)
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
		regions_data_list<-split(as.data.frame(broadPeak_selected_genes_genomeIntervals_object), f=seq_name(broadPeak_selected_genes_genomeIntervals_object))
		chrl<-names(regions_data_list)
		names(chrl)<-chrl
		outcountsBroadPeak<-sum(sapply(chrl, FUN=function(chr) {
			sum(points.within(x=abs(chip.test[[chr]]), fs=((regions_data_list[[chr]])[,1]), fe=((regions_data_list[[chr]])[,2]), return.point.counts = TRUE))
	    }))
		FRiP_broadPeak<-outcountsBroadPeak/TOTAL_reads

		###Frip sharp peaks 14
		sharpPeak_selected_genes_genomeIntervals_object<-close_intervals(interval_union(sharpPeak_selected_genes_genomeIntervals_object))
		regions_data_list<-split(as.data.frame(sharpPeak_selected_genes_genomeIntervals_object), f=seq_name(sharpPeak_selected_genes_genomeIntervals_object))
		chrl<-names(regions_data_list)
		names(chrl)<-chrl
		outcountsSharpPeak<-sum(sapply(chrl, FUN=function(chr) {
			sum(points.within(x=abs(chip.test[[chr]]), fs=((regions_data_list[[chr]])[,1]), fe=((regions_data_list[[chr]])[,2]), return.point.counts = TRUE))
		}))
		FRiP_sharpPeak<-outcountsSharpPeak/TOTAL_reads
	}else{
		TOTAL_reads=0
		FRiP_broadPeak=0
		outcountsBroadPeak=0
		FRiP_sharpPeak=0
		outcountsSharpPeak=0
	}

	QCscoreList=list("FDR_detected"=round(FDR_detect,3),
		"eval_detected"=round(eval_detect,3),
		"FRiP_broadPeak"=round(FRiP_broadPeak,3),  
		"FRiP_sharpPeak"=round(FRiP_sharpPeak,3), 
		"outcountsBroadPeak"=outcountsBroadPeak,
		"outcountsSharpPeak"= outcountsSharpPeak)

	return(QCscoreList)
}



f_tagDensity=function(data,tag.shift,rngl,parallel.mc=1)
##takes dataSelected as input, parallel is the number of CPUs used for parallelization
{
	## density distribution for data
	print("Smooth tag density")
	ts <- sum(unlist(lapply(data,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
	##parallelisation
	chromosomes_list<-names(data)
	##creates a list of lists
	data<-lapply(chromosomes_list, FUN=function(x) {
		return(data[x])
	})
	if (parallel.mc>1)
	{
		smoothed.density<-mclapply(data, FUN=function(current_chr_list)
		{
		    current_chr<-names(current_chr_list)
		    str(current_chr_list)
		    if (length(current_chr) != 1) 
		    {
		        stop("unexpected input.dataSelected structure")
		    }
		    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
		}, mc.preschedule = FALSE,mc.cores=parallel.mc)
	}else{
		smoothed.density<-lapply(data, FUN=function(current_chr_list)
		{
		    current_chr<-names(current_chr_list)
		    str(current_chr_list)
		    if (length(current_chr) != 1) 
		    {
		        stop("unexpected input.dataSelected structure")
		    }
		    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
		})
	}
	smoothed.density=(unlist(smoothed.density,recursive=FALSE))
	#normalizing smoothed tag density by library size
	smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })
	return(smoothed.density)
}
