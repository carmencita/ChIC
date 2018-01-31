#'@title Wrapper to plot non-scaled profiles for TSS or TES and to collect feature values
#'@description The non-scaled profile is constructed around the TSS/TES, with 2KB up- and downstream regions respectively. 
#' Different values are taken at the TSS/TES and surroundings with +/-2KB, +/-1KB and +/-500 sizes. 
#' For all the genomic positions, we kept the values for the ChIP and the normalized profile, as the normalization already contains 
#' information from the input. Additionally, we calculated for all of the intervals between the predefined positions the area under the profile, 
#' the local maxima (x, y coordinates), the variance, the standard deviation and the quantiles at 0%, 25%, 50% and 75%. 
#' In total the function returns 43 QC-metrics
#'
#' nonScaledMetageneProfile
#'
#' @param binnedChip DESCRIBE!!
#' @param binnedInput DESCRiBE!!
#' @param tag String, can be "TSS" or "TES". Indicates if the TSS or the TES should be calcualted (Default="TSS")
#' @param savePlotPath Path in which plots (pdf format) should be saved. If NULL on screen (default=NULL) 
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @return returnList
#'
#' @examples
#'\{dontrun
#'source("wrapper_plot_TSS_TES_allGenes.R")
#'TSS_Plot=nonScaledMetageneProfile(Meta_Result$TSS$chip,Meta_Result$TSS$input,tag="TSS",path=getwd(),debug=TRUE)
#'completeListOfValues=append(completeListOfValues,TSS_Plot)
#'
#'TES_Plot=nonScaledMetageneProfile(Meta_Result$TES$chip,Meta_Result$TES$input,tag="TES",path=getwd(),debug=TRUE)
#'completeListOfValues=append(completeListOfValues,TES_Plot)
#'}

nonScaledMetageneProfile=function(binnedChip,binnedInput,tag="TSS",savePlotPath=NULL,debug=FALSE)
{	
	print("load metagene setting")
	load("Settings.RData")
	psc <- 1; # pseudocount # required to avoid log2 of 0
	chip <- log2(do.call(rbind,binnedChip)+psc)
	input<-log2(do.call(rbind,binnedInput)+psc)
	

	input.noNorm<-colMeans(input,na.rm=T)
	chip.noNorm <- colMeans(chip,na.rm=T)
	all.noNorm<-cbind(chip.noNorm, input.noNorm)
	colnames(all.noNorm)<-c("ChIP","Input")
	

	##values at specific predefined points
	hotSpotsValues=f_spotfunction(all.noNorm, break_points, estimated_bin_size_1P, tag=tag)
	##local maxima and area in all the predefined regions
	maxAucValues=f_maximaAucfunction(all.noNorm, break_points, estimated_bin_size_1P, tag=tag)

	# chip_dispersion_TES_-1000_0% 
	# chip_dispersion_TES_-1000_25% 
	# chip_dispersion_TES_-1000_50% 
	# chip_dispersion_TES_-1000_75% 
	# chip_dispersion_TES_-1000_sd 
	# chip_dispersion_TES_-1000_variance 

	
	variabilityValues=NULL
	midpoint=as.integer(length(break_points)/2)+1
	for (i in seq(as.integer(length(break_points)/2)))
	{
		start=which(as.integer(row.names(all.noNorm))==break_points[midpoint-i])
		end=which(as.integer(row.names(all.noNorm))==break_points[midpoint+i])
		interval=all.noNorm[start:end,]	

		#dispersion=paste(var(interval[,1]),sd(interval[,1]),quantile(interval[,1])[1],quantile(interval[,1])[2],quantile(interval[,1])[3],quantile(interval[,1])[4],quantile(interval[,1])[5],sep=" ")
		#write(paste("chip",break_points[midpoint-i],dispersion,0,"dispersion TSS",sep=" "),file=outname,append=TRUE)
		sdChip=cbind(paste("chip","dispersion",tag,break_points[midpoint-i],"variance",sep="_"), var(interval[,1]))
		varChip=cbind(paste("chip","dispersion",tag,break_points[midpoint-i],"sd",sep="_"), sd(interval[,1]))
		valueFrame=quantile(interval[,1])[1:4]		#colnames(t(valueFrame))	#t(valueFrame)
		quantilesChip=cbind(paste("chip","dispersion",tag,break_points[midpoint-i],colnames(t(valueFrame)),sep="_"), (valueFrame))
		rownames(quantilesChip)=NULL

		#dispersion=paste(var(interval[,2]),sd(interval[,2]),quantile(interval[,2])[1],quantile(interval[,2])[2],quantile(interval[,2])[3],quantile(interval[,2])[4],quantile(interval[,2])[5],sep=" ")
		#write(paste("input",break_points[midpoint-i],dispersion,0,"dispersion TSS",sep=" "),file=outname,append=TRUE)

		sdInput=cbind(paste("input","dispersion",tag,break_points[midpoint-i],"variance",sep="_"), var(interval[,1]))
		varInput=cbind(paste("input","dispersion",tag,break_points[midpoint-i],"sd",sep="_"), sd(interval[,1]))
		valueFrame=quantile(interval[,2])[1:4]
		quantilesInput=cbind(paste("input","dispersion",tag,break_points[midpoint-i],colnames(t(valueFrame)),sep="_"), (valueFrame))
		rownames(quantilesInput)=NULL

		variabilityValues=rbind(variabilityValues, sdChip,varChip,quantilesChip,sdInput,varInput, quantilesInput)
	}
	colnames(variabilityValues)=c("Feature","Value")

	##make plots

	colori<-c(rev(rainbow(ncol(all.noNorm)-1)), "black")
	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,paste("ChIP_Input_",tag,".pdf",sep=""))
		pdf(file=filename,width=10, height=7)
	}    
	par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	matplot(x=as.numeric(rownames(all.noNorm)),y=all.noNorm, type="l", lwd=2, lty=1,
	col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main=tag,xaxt='n')
	abline(v=0,lty=2,col="darkgrey", lwd=2) ##plot TSS	    
	plotPoints=c(-2000,-1000,-500,500,1000,2000) ##plot remaining be
	abline(v=plotPoints,lty=3,col="darkgrey", lwd=2)
	axis(side = 1, at = sort(c(plotPoints,0)), labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	if (!is.null(savePlotPath))
	{
		dev.off()
		print(paste("pdf saved under ",filename,sep=""))
	}


	##normalized plot and values
	common_genes<-rownames(input)[rownames(input) %in% rownames(chip)]
	#all.Norm<-colMeans(t(t(input[common_genes,])-t(chip[common_genes,])),na.rm=T)
	all.Norm<-colMeans(t(t(chip[common_genes,])-t(input[common_genes,])),na.rm=T)

	hotSpotsValuesNorm=f_spotfunctionNorm(all.Norm,  break_points, estimated_bin_size_1P, tag=tag)

	maxAucValuesNorm=f_maximaAucfunctionNorm(all.Norm,  break_points, estimated_bin_size_1P, tag=tag)

	variabilityValuesNorm=NULL
	midpoint=as.integer(length(break_points)/2)+1
	for (i in seq(as.integer(length(break_points)/2)))
	{
		start=which(as.integer(names(all.Norm))==break_points[midpoint-i])
		end=which(as.integer(names(all.Norm))==break_points[midpoint+i])
		interval=all.Norm[start:end]	
		#dispersion=paste(var(interval),sd(interval),quantile(interval)[1],quantile(interval)[2],quantile(interval)[3],quantile(interval)[4],quantile(interval)[5],sep=" ")
		#write(paste("norm",break_points[midpoint-i],dispersion,1,"dispersion TSS",sep=" "),file=outname,append=TRUE)
		sdNorm=cbind(paste("norm","dispersion",tag,break_points[midpoint-i],"variance",sep="_"), var(interval))
		varNorm=cbind(paste("norm","dispersion",tag,break_points[midpoint-i],"sd",sep="_"), sd(interval))
		valueFrame=quantile(interval)[1:4]		#colnames(t(valueFrame))	#t(valueFrame)
		quantilesNorm=cbind(paste("norm","dispersion",tag,break_points[midpoint-i],colnames(t(valueFrame)),sep="_"), (valueFrame))
		rownames(quantilesNorm)=NULL

		variabilityValuesNorm=rbind(variabilityValuesNorm, sdNorm,varNorm,quantilesNorm)

	}
	colnames(variabilityValuesNorm)=c("Feature","Value")

	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,paste("Normalized_",tag,".pdf",sep=""))
		pdf(file=filename,width=10, height=7)
	}
	par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	plot(x=as.numeric(names(all.Norm)),y=all.Norm, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="mean log2 enrichment (signal/input)",
		main=paste("normalized",tag,sep=" "), xaxt='n')#,cex.axis=1.3,cex.lab=1.3)
	abline(v=0,lty=2,col="darkgrey", lwd=2) ##plot TSS
	abline(v=plotPoints,lty=3,col="darkgrey", lwd=2)
	axis(side = 1, at = sort(c(plotPoints,0)), labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	if (!is.null(savePlotPath))
	{
		dev.off()
		print(paste("pdf saved under ",filename,sep=""))
	}

	result=NULL
	result=rbind(hotSpotsValues,maxAucValues,variabilityValues,hotSpotsValuesNorm,maxAucValuesNorm,variabilityValuesNorm)


	if (debug)
	{
		print("Debugging mode ON")
		outname=file.path(getwd(), paste(tag,"onepoints.result",sep="_"))
		file.remove(outname)
		write.table(result,file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	}


	return(result)

}