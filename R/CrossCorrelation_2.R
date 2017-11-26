library("spp")
library("snow")
library("girafe")
library("caTools")
require("parallel")


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
	plot(datax_y,type='l',xlab=xlabname,ylab=ylabname,main=maintitel)
	abline(v=line,lty=2,lwd=2, col="red")
	title(chipname)
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


chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
rownames(chrominfo)<-chrominfo$chrom
rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(as.integer(c(1,x)))})


#sampleIndex=16
##step 00: binding characteristics
#get dataset
sampleIndex=which(sampleinfo$Filename==datafilename)
print(sampleIndex)

#datafilename=sampleinfo$Filename[sampleIndex]
#print(datafilename)
TFname=sampleinfo$IG[sampleIndex]

inputIndex=sampleIndex
inputname=sampleinfo$ControlName[inputIndex]
##load objects and tag.shift to calculate tagDensity
#load inputdata

load(file=paste(storagedir,"dataSelected_",datafilename,"_",inputname,".RData",sep=""))



####################################################################################
####################################################################################
####
#### calculate TagDensity
####



## density distribution for the input
print("Smooth tag density input")
ts <- sum(unlist(lapply(input.dataSelected,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
###original code
#smoothed.densityOrig <- get.smoothed.tag.density(input.dataSelected, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift,rngl=rngl)
##parallelisation
chromosomes_list<-names(input.dataSelected)

##creates a list of lists
input.dataSelected<-lapply(chromosomes_list, FUN=function(x) {
	return(input.dataSelected[x])
})

#names(input.dataSelected)=chromosomes_list

smoothed.density<-mclapply(input.dataSelected, 
    FUN=function(current_chr_list) {
    current_chr<-names(current_chr_list)
    str(current_chr_list)
    if (length(current_chr) != 1) {
        stop("unexpected input.dataSelected structure")
        }
    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
    }
, mc.preschedule = FALSE,mc.cores=mc)

rm(input.dataSelected)


smoothed.density=(unlist(smoothed.density,recursive=FALSE))
#normalizing smoothed tag density by library size
smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })


save(smoothed.density,file=paste(storagedir,"TagDensity_",datafilename,"_",inputname,".RData",sep=""))
str(smoothed.density)
rm(smoothed.density)


## density distribution for the chip
print("Smooth tag density Chip")
ts <- sum(unlist(lapply(chip.dataSelected,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
chromosomes_list<-names(chip.dataSelected)

chip.dataSelected<-lapply(chromosomes_list, FUN=function(x) {
	return(chip.dataSelected[x])
})

smoothed.density<-mclapply(chip.dataSelected, 
	FUN=function(current_chr_list) {
	current_chr<-names(current_chr_list)
	print(current_chr)
	if (length(current_chr) != 1) {
		stop("unexpected chip.dataSelected structure")
		}
	get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
	}
, mc.preschedule = FALSE,mc.cores=mc)

rm(chip.dataSelected)

smoothed.density=(unlist(smoothed.density,recursive=FALSE))
#names(smoothed.density)=chromosomes_list
#normalizing smoothed tag density by library size
smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })

save(smoothed.density,file=paste(storagedir,"TagDensity_",datafilename,".RData",sep=""))
str(smoothed.density)
rm(smoothed.density)

file.remove(paste(storagedir,"dataSelected_",datafilename,"_",inputname,".RData",sep=""))

