library("spp")
library("snow")
library("girafe")
library("Rmpi")
library("caTools")



f_plot <- function(datax_y, chipname, maintitel="title", plotname="plot",xlabname="x-axis",ylabname="y-axis",line=NULL, lineplotX=NULL,lineplotY=NULL) 
{
	#options(bitmapType='cairo')
	filename=file.path(plotsdir, paste(plotname, chipname, "pdf", sep="."))
	print(filename)
	pdf(filename)	
	#print(file.path(outputdir, paste(plotname, chipname, "pdf", sep=".")))
	#bitmap(filename,"png16m")
	#pdf(filename=file.path(outputdir, paste(plotname, chipname, "pdf", sep=".")))
	par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
	plot(datax_y,type='l',xlab=xlabname,ylab=ylabname)
	abline(v=line,lty=2,lwd=2, col="red")
	title(maintitel)
	!is.null(lineplotX)
	{
		lines(x=lineplotX, y=lineplotY, lwd=2, col="blue")
	}
	
	dev.off()
	return(TRUE)
}


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
	if (is.null(margin)) {
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
	x=binding.characteristics$cross.correlation$x
	y=binding.characteristics_cross.correlation_y_smoothed
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

#path="/gpfs/work/IscrC_CONCEPt/QC_pipeline/"
#datafilename="ENCFF000OYO"
#read_length=50


t0<-proc.time()[3]


arg <- commandArgs()
print(arg)
#3sampleIndex=as.integer(arg[7])
path=arg[6]
datafilename=arg[7]
print(datafilename)
read_length=as.integer(arg[8])


#print (sampleIndex)
print(path)
print(read_length)

source(paste(path,"GlobalParameters.R",sep=""))
sampleinfo_file<-paste(path,"matchlist.txt",sep="")
sampleinfo<-read.table(sampleinfo_file,  header=TRUE, quote="", stringsAsFactors=FALSE)


workingdir<-paste(path,"CrossCorrelation/",sep="")
outputdir<-paste(workingdir,"out/",sep="")
timedir<-paste(workingdir,"time_stamps/",sep="")
plotsdir<-paste(workingdir,"plots/",sep="")

mc <- getOption("mc.cores",mpi.universe.size() )
if (cluster_ON_OFF==TRUE)
{
	cluster=makeCluster(mc, type="MPI")
}else{cluster=NULL}

chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
rownames(chrominfo)<-chrominfo$chrom
rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})


##step 00: binding characteristics
#get dataset
sampleIndex=which(sampleinfo$Filename==datafilename)
print(sampleIndex)

#datafilename=sampleinfo$Filename[sampleIndex]
#print(datafilename)
TFname=sampleinfo$IG[sampleIndex]


#get readcount using specific aligner 
read.tags.current_function<-get(paste("read", reads.aligner.type , "tags", sep="."))

if (reads.aligner.type=="bam")
{
    chip.data<-read.tags.current_function(file.path(bamdir,paste(sampleinfo$Filename[sampleIndex],"bam",sep=".")))
}
if (reads.aligner.type=="tagalign")
{
    chip.data<-read.tags.current_function(file.path(bamdir,paste(sampleinfo$Filename[sampleIndex],"tagAlign",sep=".")))
}

readCount=sum(sapply(chip.data$tags, length))
cat("READS count", "\n", datafilename, readCount, "\n")



###step 1: tag distribution 
#get binding characteristics 
chip.data_tagdistribution<-sapply(chip.data$tags, length)
binding.characteristics<-get.binding.characteristics(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)
#plot cross correlation curve
print("Cross correlation plot")
f_plot(binding.characteristics$cross.correlation,datafilename,maintitel=TFname,plotname="CrossCorrelation",xlab="strand shift",ylab="cross-correlation",line=binding.characteristics$peak$x)


if (reads.aligner.type=="bam")
{
    file.remove(paste(bamdir,datafilename,".bam",sep=""))
    
}
if (reads.aligner.type=="tagalign")
{
    file.remove(paste(bamdir,datafilename,".tagAlign",sep=""))
    file.remove(paste(bamdir,datafilename,".tagAlign.gz",sep=""))
}



print(binding.characteristics$peak$x)
#save(input.data,file=paste("/data/FF/Carmen/Pipeline/",datafilename,".RData",sep=""))
#load(file.path(datadir, paste(datafilename,".RData",sep="")))


###step 1.2: Phantom peak and cross-correlation
print("Estimating fragment lengths")
phantom.characteristics<-get.binding.characteristics(chip.data, srange=PhantomPeak_range, bin=PhantomPeak_bin, cluster=cluster)
print(datafilename)
print(phantom.characteristics$peak$x)

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
qflag <- NA
if ( (RSC >= 0) & (RSC < 0.25) ) {
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

phantom_peak.scores <- list(phantom_cc=phantom_cc, NSC=NSC, RSC=RSC, quality_flag=qflag, min_cc=min_cc, peak=phantom.characteristics$peak, read_length=read_length)
cat("NSC=", round(NSC, 2), ", RSC=", round(RSC, 2), ", Quality flag: ", qflag, "\n", sep="")


###2.0 smoothed cross correlation
subset_selection<- which(binding.characteristics$cross.correlation$x %in% cross_correlation_range_subset)
binding.characteristics$cross.correlation<-binding.characteristics$cross.correlation[subset_selection,]
# add smoothing
binding.characteristics_cross.correlation_y_smoothed<-caTools::runmean(binding.characteristics$cross.correlation$y, k=cross_correlation_smoothing_k)
# assign the new maximum coordinates
binding.characteristics$peak$y<-max(binding.characteristics_cross.correlation_y_smoothed)
binding.characteristics$peak$x<-binding.characteristics$cross.correlation$x[(which(binding.characteristics_cross.correlation_y_smoothed==binding.characteristics$peak$y))]


# plot cross correlation curve with smoothing
f_plot(binding.characteristics$cross.correlation,datafilename,maintitel=TFname,plotname="customCrossCorrelation",xlab="strand shift",ylab="cross-correlation",line=binding.characteristics$peak$x, lineplotX=binding.characteristics$cross.correlation$x,lineplotY=binding.characteristics_cross.correlation_y_smoothed)

strandShift<-binding.characteristics$peak$x



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


if (max_x_peakcheck>phantom.characteristics$peak$x) {
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
    qflag <- NA
    if ( (RSC >= 0) & (RSC < 0.25) ) {
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
    phantom_peak.scores$quality_flag<-qflag


} else {
  phantom.characteristics$peak$y#<-max(phantom.characteristics_cross.correlation_y_smoothed)
  phantom.characteristics$peak$x#<-phantom.characteristics$cross.correlation$x[(which(phantom.characteristics_cross.correlation_y_smoothed==phantom.characteristics$peak$y))]
}


  # plot cross correlation curve with smoothing
  pdf((filename=file.path(plotsdir, paste("phantomCrossCorrelation", datafilename, "pdf", sep="."))))
  par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
  plot(phantom.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation",main=TFname)
  lines(x=phantom.characteristics$cross.correlation$x, y=phantom.characteristics_cross.correlation_y_smoothed, lwd=2, col=cross_correlation_smoothing_color)
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
    paste("Read length =", (phantom_peak.scores$read_length))
    ))
  title(datafilename)
  dev.off()

phantomScores= c(
    Sample=datafilename,
    NSC=phantom_peak.scores$NSC,
    RSC=phantom_peak.scores$RSC,
    quality_flag=phantom_peak.scores$quality_flag,
    shift=phantom_peak.scores$peak$x,
    read_length=phantom_peak.scores$read_length,
    A=phantom_peak.scores$peak$y,
    B=phantom_peak.scores$phantom_cc$y,
    C=phantom_peak.scores$min_cc$y)


###4 NRF calculation
ALL_TAGS<-sum(sapply(chip.data$tags, length))
UNIQUE_TAGS<-sum(sapply(lapply(chip.data$tags, unique), length))
UNIQUE_TAGS_nostrand<-sum(sapply(lapply(chip.data$tags, FUN=function(x) {unique(abs(x))}), length))

NRF<-UNIQUE_TAGS/ALL_TAGS
NRF_nostrand<-UNIQUE_TAGS_nostrand/ALL_TAGS



## to compensate for lib size differences
nomi<-rep(names(chip.data$tags), sapply(chip.data$tags, length))
chip.dataNRF<-unlist(chip.data$tags)
names(chip.dataNRF)<-NULL
chip.dataNRF<-paste(nomi, chip.dataNRF, sep="")

if (ALL_TAGS > 10e6) {
   UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:100, FUN=function(x) {
     return(length(unique(sample(chip.dataNRF, size=10e6))))
   })))
} else {
   UNIQUE_TAGS_LibSizeadjusted<-round(mean(sapply(1:100, FUN=function(x) {
     return(length(unique(sample(chip.dataNRF, size=10e6, replace=TRUE))))
   })))
}

NRF_LibSizeadjusted<-UNIQUE_TAGS_LibSizeadjusted/10e6

STATS_NRF=c(Sample=datafilename, ALL_TAGS=ALL_TAGS, UNIQUE_TAGS=UNIQUE_TAGS, UNIQUE_TAGS_nostrand=UNIQUE_TAGS_nostrand, NRF=NRF, NRF_nostrand=NRF_nostrand, NRF_LibSizeadjusted=NRF_LibSizeadjusted)



#N1= number of genomic locations to which EXACTLY one unique mapping read maps
#Nd = the number of genomic locations to which AT LEAST one unique mapping read maps, i.e. the number of non-redundant, unique mapping reads

#PBC = N1/Nd


N1<-sum(sapply(chip.data$tags, FUN=function(x) {
	checkDuplicate<-duplicated(x)
	duplicated_positions<-unique(x[checkDuplicate])
	OUT<-x[!(x %in% duplicated_positions)]
	return(length(OUT))
}))

Nd<-sum(sapply(lapply(chip.data$tags, unique), length))
PBC = N1/Nd



#### join input and chip
###3.0 save wig smoothed and 
###3.1 CHIP-INPUT enrichment
tag.shift <- round(strandShift/2)

#load inputdata
inputIndex=sampleIndex
inputname=sampleinfo$ControlName[inputIndex]



print("Load control file...")
newControl=NULL
if (length(unlist(strsplit(inputname,"\\.")))>1)
{
	print("Concat multiple inputs...")
	for (element in unlist(strsplit(inputname,"\\.")))
	{	
		if (element!="")
		{
			if (is.null(newControl))
			{	
				##if newControl is empty, fill it with the first input-data
				print(element)
				load(file.path(path, paste(element,".RData",sep="")))
				newControl<-input.data
				print("original size")
			
				for (i in seq(length(newControl$tags)))
				{
					print(length(newControl$tags[[i]]))
				}
				#newControl.tags=(sapply(newControl$tags,c))
				#newControl.q=(sapply(newControl$quality,c))
				#print(length(newControl.tags[[25]]))
				#print(length(newControl.q[[25]]))
			}else{
				##otherwise add the input-data and sort
				print(element)
				load(file.path(path, paste(element,".RData",sep="")))
				temp<-input.data
				#temp.tags=(sapply(temp$tags,c))
				#temp.q=(sapply(temp$quality,c))
				chrl=names(input.data$tags)
				for (i in chrl)
				{
					print(i)
					##appending tags and quality elements of all inputs to newControl
					newControl$tags[[i]]=append(newControl$tags[[i]],temp$tags[[i]])
					newControl$quality[[i]]=append(newControl$quality[[i]],temp$quality[[i]])
				}
				}
		}
	}
	print("final size")
	for (i in seq(length(newControl$tags)))
	{
		print(length(newControl$tags[[i]]))
	}
	##sort tags and quality of the final NewControl to get an increasing list independently from the sign:  -57 -95 -112 151 159 166 169 -217...
	print("Sort tags and quality flags...")
	chrl=names(newControl$tags)
	for (i in chrl)
	{
		print(i)
		index=sort.list(abs(newControl$tags[[i]]))
		newControl$tags[[i]]=newControl$tags[[i]][index]	
		newControl$quality[[i]]=newControl$quality[[i]][index]	
	}
	input.data=newControl
}else{
		#inputIndex=as.integer(arg[8])
	print("Single input...")
	load(file.path(path, paste(inputname,".RData",sep="")))
	input.data<-input.data
}


chrl_final=intersect(names(chip.data$tags),names(input.data$tags))


chip.data$tags=chip.data$tags[chrl_final]
chip.data$quality=chip.data$quality[chrl_final]

input.data$tags=input.data$tags[chrl_final]
input.data$quality=input.data$quality[chrl_final]


print("Filter tags")
if (select.informative.tags_filter) {
      print("select.informative.tags filter")
     #load(paste("sppdata", "binding", chip.data.samplename, "RData", sep="."))
      chip.dataSelected <- select.informative.tags(chip.data, binding.characteristics)
      #load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
      input.dataSelected <- select.informative.tags(input.data, binding.characteristics)
} else {
      print("SKIP select.informative.tags filter")
      chip.dataSelected<-chip.data$tags
      input.dataSelected<-input.data$tags
}

if (remove.local.tag.anomalies_filter) {
	print("remove.local.tag.anomalies filter")
      # restrict or remove singular positions with very high tag counts
	chip.dataSelected <- remove.local.tag.anomalies(chip.dataSelected)
	input.dataSelected <- remove.local.tag.anomalies(input.dataSelected)
} else {
	print("SKIP remove.local.tag.anomalies filter")
	input.dataSelected<-input.data$tags
	chip.dataSelected=chip.data$tags
}




##Save data object to be read afterwards for the densityt distribution
save(input.dataSelected,tag.shift,chip.dataSelected,file=paste(storagedir,"dataSelected_",datafilename,"_",inputname,".RData",sep=""))


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

		extension_vector<-size(broadPeak_selected_genes_genomeIntervals_object) ##saves the size of the peaks (intervals)

		STATS= c("Sample"=datafilename, "window"=current_window_size,"zthresh"=current_zthresh, "Regions number"=length(extension_vector), "Total extension"=sum(extension_vector), "Mean extension"=mean(extension_vector), "Median extension"=median(extension_vector))
		extension_vector_list<-c(extension_vector_list, list(extension_vector))
	}
}


###12 get binding sites with FDR and eval


chip.data12<-chip.dataSelected[(names(chip.dataSelected) %in% custom_chrorder)]
input.data12<-input.dataSelected[(names(input.dataSelected) %in% custom_chrorder)]

print("Binding sites detection fdr")
fdr <- 1e-2
detection.window.halfsize <- tag.shift
print("Window Tag Density method - WTD")
bp_FDR <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,fdr=fdr,whs=detection.window.halfsize*2, tag.count.whs=detection.window.halfsize, cluster=cluster)
print(paste("detected",sum(unlist(lapply(bp_FDR$npl,function(d) length(d$x)))),"peaks"))
# output detected binding positions
#output.binding.results(results=bp, filename=paste("WTC.binding.positions", chip.data.samplename, "txt", sep="."))

print("Binding sites detection evalue")
eval<-10
bp_eval <- find.binding.positions(signal.data=chip.data12,control.data=input.data12,e.value=eval,whs=detection.window.halfsize*2,cluster=cluster)
print(paste("detected",sum(unlist(lapply(bp_eval$npl,function(d) length(d$x)))),"peaks"))
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
	#require(girafe)
	chip.test<-lapply(chip.data$tags, FUN=function(x) {x+tag.shift}) ###SORTS the tags for each chrom
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


spaceUsage=sum(sort(sapply(ls(),function(x){object.size(get(x))}))) ##in bytes

print("switch cluster off")
if (cluster_ON_OFF==TRUE)
{
	 stopCluster(cluster)
}


print("write results in file")
################RESULTS
outname=paste(outputdir,datafilename,".results",sep="")
file.remove(outname)
#sampleIndex datafilename
#readCount read_length
write(paste("TF: ",TFname,sep=" "),file=outname,append=T)
write(paste("Input: ",datafilename,sampleIndex,sep=" "),file=outname,append=T)
write(paste("Input: ",inputname,inputIndex,sep=" "),file=outname,append=T)

write(paste("ReadCount",readCount,sep=" "),file=outname,append=T)
write(paste("Read_length",read_length,sep=" "),file=outname,append=T)
#strandShift<-binding.characteristics$peak$x
write(paste("StrandShift",strandShift,sep=" "),file=outname,append=T)
write(paste("Substitution of StrandShift from ",oldShift," to ",strandShift,sep=" "),file=outname,append=T)

#binding.characteristics$peak$whs
write(paste("bindingCharacteristicsPeak (x,y,whs)",binding.characteristics$peak$x,binding.characteristics$peak$y,binding.characteristics$whs,sep=" "),file=outname,append=T)

#phantom.characteristics$peak$x 
#phantom.characteristics$peak$y
#phantom.characteristics$peak$whs

write(paste("phantomCharacteristicsPeak (x,y,whs)",phantom.characteristics$peak$x,phantom.characteristics$peak$y,phantom.characteristics$whs,sep=" "),file=outname,append=T)
#phantomScores

write(paste("ALL_TAGS",ALL_TAGS,sep=" "),file=outname,append=T)
write(paste("NSC",round(NSC, 2),sep=" "),file=outname,append=T)
write(paste("RSC",round(RSC, 2),sep=" "),file=outname,append=T)
write(paste("Quality flag: ", qflag,sep=" "),file=outname,append=T)
write(paste("shift: ", round(as.double(phantomScores["shift"]),2),sep=" "),file=outname,append=T)

write(paste("read length",round(as.double(phantomScores["read_length"]),2),sep=" "),file=outname,append=T)
write(paste("A: ", round(as.double(phantomScores["A"]),2),sep=" "),file=outname,append=T)
write(paste("B: ", round(as.double(phantomScores["B"]),2),sep=" "),file=outname,append=T)
write(paste("C: ", round(as.double(phantomScores["C"]),2),sep=" "),file=outname,append=T)

write(paste("FDR detected",sum(unlist(lapply(bp_FDR$npl,function(d) length(d$x)))),"peaks",sep=" "),file=outname,append=T)
write(paste("eval detected",sum(unlist(lapply(bp_eval$npl,function(d) length(d$x)))),"peaks",sep=" "),file=outname,append=T)

#STATS_NRF

write(paste("UNIQUE_TAGS_LibSizeadjusted",UNIQUE_TAGS_LibSizeadjusted,sep=" "),file=outname,append=T)
write(paste("NRF_LibSizeadjusted",NRF_LibSizeadjusted,sep=" "),file=outname,append=T)
write(paste("ALL_TAGS",ALL_TAGS,sep=" "),file=outname,append=T)
write(paste("UNIQUE_TAGS",UNIQUE_TAGS,sep=" "),file=outname,append=T)
write(paste("UNIQUE_TAGS_nostrand",UNIQUE_TAGS_nostrand,sep=" "),file=outname,append=T)
write(paste("NRF",NRF,sep=" "),file=outname,append=T)
write(paste("NRF_LibSizeadjusted",NRF_LibSizeadjusted,sep=" "),file=outname,append=T)
write(paste("NRF_nostrand",NRF_nostrand,sep=" "),file=outname,append=T)
write(paste("PBC",PBC,sep=" "),file=outname,append=T)
write(paste("N1",N1,sep=" "),file=outname,append=T)
write(paste("Nd",Nd,sep=" "),file=outname,append=T)
#Frip
write(paste("Total_reads",TOTAL_reads,"FRiP_broadPeak",round(FRiP_broadPeak, 2),"outcountsBroadPeak", outcountsBroadPeak, "FRiP_sharpPeak", FRiP_sharpPeak, "outcountsSharpPeak", outcountsSharpPeak, sep=" "),file=outname,append=T)

t1<-proc.time()[3]
deltat=t1-t0

write(paste("CrossCorrelation_Timing",datafilename,deltat,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""))
write(paste("Space_Usage",datafilename,spaceUsage,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""),append=T)
