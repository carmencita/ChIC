#'@title Wrapper to plot scaled profile and to collect feature values
#'@description The scaled metagene profile that includes the gene body, the signal is 
#' captured on a real scale from the TSS and an upstream region of 2KB. From the TSS, 
#' the gene body is constructed with 0.5KB in real scale at the gene start 
#' (TSS + 0.5KB) and the gene end (TES - 0.5KB), whereas the remaining gene body is 
#' scaled to a virtual length of 2000. Considering the length of these regions, 
#' the minimum gene length required is 3KB and shorter genes are filtered out. 
#' From the profile, we take enrichment values at different coordinates: at 
#' -2KB, at the TSS, inner margin (0.5KB), gene body (2KB + 2 * inner margin), 
#' gene body+1KB. We collect in total 42 QC-metrics from the ChIP and 
#' normalized profile. 
#'
#' scaledMetageneProfile
#'
#' @param binnedChip DESCRIBE!!
#' @param binnedInput DESCRiBE!!
#' @param savePlotPath Path in which plots (pdf format) should be saved. 
#' If NULL on screen (default=NULL) 
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @return returnList
#'
#'@examples
#'\dontrun{
#' source("wrapper_twopoint_allGenes_plot.R")
#' geneBody_Plot=scaledMetageneProfile(Meta_Result$twopoint$chip,
#' Meta_Result$twopoint$input,path=getwd(),debug=TRUE)
#' completeListOfValues=append(completeListOfValues,geneBody_Plot)
#'}

scaledMetageneProfile=function(binnedChip,binnedInput,savePlotPath=NULL,debug=FALSE)
{
	print("load metagene setting")
	settings=f_metaGeneDefinition(selection="Settings")
	psc <- 1; # pseudocount # required to avoid log2 of 0
	break_points_2P=settings$break_points_2P
	estimated_bin_size_2P=settings$estimated_bin_size_2P

	chip <- log2(do.call(rbind,binnedChip)+psc)
	input<-log2(do.call(rbind,binnedInput)+psc)	

	input.noNorm<-colMeans(input,na.rm=T)
	chip.noNorm <- colMeans(chip,na.rm=T)
	all.noNorm<-cbind(chip.noNorm, input.noNorm)
	colnames(all.noNorm)<-c("ChIP","Input")
	

	##values at specific predefined points
	hotSpotsValues=f_spotfunction(all.noNorm, break_points_2P, estimated_bin_size_2P, tag="twopoints")
	##local maxima and area in all the predefined regions
	maxAucValues=f_maximaAucfunction(all.noNorm, break_points_2P, estimated_bin_size_2P, tag="twopoint")

	##make plots
	colori<-c(rev(rainbow(ncol(all.noNorm)-1)), "black")
	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,"ScaledMetaGene_ChIP_Input.pdf")
		pdf(file=filename,width=10, height=7)
	}    
	par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	matplot(x=as.numeric(rownames(all.noNorm)),y=all.noNorm, type="l", lwd=2, lty=1,
	col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main="metagene",xaxt='n')

	currBreak_points=break_points_2P[c(-2,-5)] #c(-2000,500,2500,4000)
	abline(v=c( break_points_2P[c(2,5)]),lty=2,col="darkgrey", lwd=3)
	abline(v=currBreak_points,lty=3,col="darkgrey", lwd=2)
	axis(side = 1, at = break_points_2P, labels = c("-2KB","TSS","TSS+500","TES-500","TES","+1KB"))


	#plotpoints=c(-2000,-1000,500,2500,4000)
	#abline(v=c(0,totalGeneLength),lty=2,col="darkgrey", lwd=3)
	#abline(v=plotpoints,lty=3,col="darkgrey", lwd=2)
	#axis(side = 1, at =sort(c(plotpoints,0,3000)), labels = c("-2KB","-1KB","TSS","500","500","TES","+1KB"))    
	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	if (!is.null(savePlotPath))
	{
		dev.off()
		print(paste("pdf saved under ",filename,sep=""))
	}

	###################
	##normalized plot and values
	###################

	common_genes<-rownames(input)[rownames(input) %in% rownames(chip)]
	#frameNormalized<-colMeans(t(t(input[common_genes,])-t(chip[common_genes,])),na.rm=T)
	frameNormalized<-colMeans(t(t(chip[common_genes,])-t(input[common_genes,])),na.rm=T)

	hotSpotsValuesNorm=f_spotfunctionNorm(frameNormalized, break_points_2P, estimated_bin_size_2P, tag="twopoints")

	maxAucValuesNorm=f_maximaAucfunctionNorm(frameNormalized, break_points_2P, estimated_bin_size_2P, tag="twopoints")

	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,"ScaledMetaGene_normalized.pdf")
		pdf(filename,width=10, height=7)
	}
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	plot(x=as.numeric(names(frameNormalized)),y=frameNormalized, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="mean log2 enrichment (signal/input)",
		main="normalized metagene", xaxt='n')#,cex.axis=1.3,cex.lab=1.3)
	
	currBreak_points=break_points_2P[c(-2,-5)] #c(-2000,500,2500,4000)
	abline(v=c( break_points_2P[c(2,5)]),lty=2,col="darkgrey", lwd=3)
	abline(v=currBreak_points,lty=3,col="darkgrey", lwd=2)
	axis(side = 1, at = break_points_2P, labels = c("-2KB","TSS","TSS+500","TES-500","TES","+1KB"))

#	abline(v=c(0,totalGeneLength),lty=2,col="darkgrey", lwd=3)
#	abline(v=plotpoints,lty=3,col="darkgrey", lwd=2)
 #   axis(side = 1, at =sort(c(plotpoints,0,3000)), labels = c("-2KB","-1KB","TSS","500","500","TES","+1KB"))    

	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	if (!is.null(savePlotPath))
	{
		dev.off()
		print(paste("pdf saved under ",filename,sep=""))
	}

	finalValues=rbind(hotSpotsValues,maxAucValues,hotSpotsValuesNorm,maxAucValuesNorm)
	if (debug)
	{
		print("Debugging mode ON")
		outname=file.path(path, "twopoints.result")
		print(outname)
		file.remove(outname)
		write.table(finalValues,file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	}

	return(finalValues)
	
}
