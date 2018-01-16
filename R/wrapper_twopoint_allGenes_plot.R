# path=getwd()
# source(file.path(path,"GlobalParameters.R"))
# source(file.path(path,"Functions_MetaGenePlots.R"))



# chipName="ENCFF000BBB"
# inputName="ENCFF000BAF"
# debug=TRUE
# cluster=NULL
# dataPath="/lustre//data/FF/Carmen/BitBucket/chic/data"


# load(file.path(path, paste(chipName,inputName,"Twopoint.RData",sep="_")))


# help=binned_Chip
# binned_Chip=binned_Input
# binned_Input=help

f_plotMetageneProfile=function(binnedChip,binnedInput,path=getwd(),debug=FALSE,plotName="NA")
{
	source("FunctionsLocal.R")
	psc <- 1; # pseudocount # required to avoid log2 of 0

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
	pdf(file=file.path(path,paste(plotName,"ChIP_Input_MetaGene.pdf",sep="_")),width=10, height=7)
	    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	    matplot(x=as.numeric(rownames(all.noNorm)),y=all.noNorm, type="l", lwd=2, lty=1,
	    col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main="metagene",xaxt='n')
		plotpoints=c(-2000,-1000,500,2500,4000)
		abline(v=c(0,totalGeneLength),lty=2,col="darkgrey", lwd=3)
		abline(v=plotpoints,lty=3,col="darkgrey", lwd=2)
	    axis(side = 1, at =sort(c(plotpoints,0,3000)), labels = c("-2KB","-1KB","TSS","500","500","TES","+1KB"))    
	    legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	dev.off()


	###################
	##normalized plot and values
	###################

	common_genes<-rownames(input)[rownames(input) %in% rownames(chip)]
	#frameNormalized<-colMeans(t(t(input[common_genes,])-t(chip[common_genes,])),na.rm=T)
	frameNormalized<-colMeans(t(t(chip[common_genes,])-t(input[common_genes,])),na.rm=T)

	hotSpotsValuesNorm=f_spotfunctionNorm(frameNormalized, break_points_2P, estimated_bin_size_2P, tag="norm")

	maxAucValuesNorm=f_maximaAucfunctionNorm(frameNormalized, break_points_2P, estimated_bin_size_2P, tag="norm")

	pdf(file.path(path, paste(plotName,"twopointsNormalized.pdf",sep="_")),width=10, height=7)
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	plot(x=as.numeric(names(frameNormalized)),y=frameNormalized, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="mean log2 enrichment (signal/input)",
		main="normalized metagene", xaxt='n')#,cex.axis=1.3,cex.lab=1.3)
	abline(v=c(0,totalGeneLength),lty=2,col="darkgrey", lwd=3)
	abline(v=plotpoints,lty=3,col="darkgrey", lwd=2)
	   axis(side = 1, at =sort(c(plotpoints,0,3000)), labels = c("-2KB","-1KB","TSS","500","500","TES","+1KB"))    
	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	dev.off()

	finalValues=rbind(hotSpotsValues,maxAucValues,hotSpotsValuesNorm,maxAucValuesNorm)
	if (debug)
	{
		outname=file.path(path, paste(plotName,"twopoints.result",sep="_"))

		file.remove(outname)

		write.table(finalValues,file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	}

	return(finalValues)
	
}



# chip_hotSpots_twopoints_-2000 
# chip_hotSpots_twopoints_0 
# chip_hotSpots_twopoints_2500 
# chip_hotSpots_twopoints_3000 
# chip_hotSpots_twopoints_4000 
# chip_hotSpots_twopoints_500
# chip_localMax_twopoint_1_x 
# chip_localMax_twopoint_1_y 
# chip_localMax_twopoint_2_x 
# chip_localMax_twopoint_2_y 
# chip_localMax_twopoint_3_x 
# chip_localMax_twopoint_3_y 
# chip_localMax_twopoint_4_x 
# chip_localMax_twopoint_4_y 
# chip_localMax_twopoint_5_x 
# chip_localMax_twopoint_5_y 
# chip_auc_twopoint_1 
# chip_auc_twopoint_2 
# chip_auc_twopoint_3 
# chip_auc_twopoint_4 
# chip_auc_twopoint_5  
# input_auc_twopoint_1 
# input_auc_twopoint_2 
# input_auc_twopoint_3 
# input_auc_twopoint_4 
# input_auc_twopoint_5 
# input_hotSpots_twopoints_-2000 
# input_hotSpots_twopoints_0 
# input_hotSpots_twopoints_2500 
# input_hotSpots_twopoints_3000 
# input_hotSpots_twopoints_4000 
# input_hotSpots_twopoints_500 
# input_localMax_twopoint_1_x 
# input_localMax_twopoint_1_y 
# input_localMax_twopoint_2_x 
# input_localMax_twopoint_2_y 
# input_localMax_twopoint_3_x 
# input_localMax_twopoint_3_y 
# input_localMax_twopoint_4_x 
# input_localMax_twopoint_4_y 
# input_localMax_twopoint_5_x 
# input_localMax_twopoint_5_y
# norm_auc_twopoints_1 
# norm_auc_twopoints_2 
# norm_auc_twopoints_3 
# norm_auc_twopoints_4 
# norm_auc_twopoints_5 
# norm_hotSpots_twopoints_-2000 
# norm_hotSpots_twopoints_0 
# norm_hotSpots_twopoints_2500 
# norm_hotSpots_twopoints_3000 
# norm_hotSpots_twopoints_4000 
# norm_hotSpots_twopoints_500
# norm_localMax_twopoints_1_x
# norm_localMax_twopoints_1_y
# norm_localMax_twopoints_2_x 
# norm_localMax_twopoints_2_y 
# norm_localMax_twopoints_3_x 
# norm_localMax_twopoints_3_y 
# norm_localMax_twopoints_4_x 
# norm_localMax_twopoints_4_y 
# norm_localMax_twopoints_5_x 
# norm_localMax_twopoints_5_y