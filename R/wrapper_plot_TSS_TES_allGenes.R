###FOR DEVEL ONLY#####
# path=getwd()
# source(file.path(path,"GlobalParameters.R"))

# chipName="ENCFF000BBB"
# inputName="ENCFF000BAF"
# debug=TRUE
# dataPath="/lustre//data/FF/Carmen/BitBucket/chic/data"


# load(file.path(path, paste(chipName,inputName,"TSS.RData",sep="_")))

# outname=file.path(path, paste(chipName,inputName,"TSS.result",sep="_"))

# binnedInput=binnedChip_TSS
# binnedChip=binnedInput_TSS
# tag="TSS"
###FOR DEVEL ONLY####


#' f_plotMetageneProfile_onePoint
#'
#' @param chipName
#' @param inputName
#' @param read_length
#' @param reads.aligner.type
#' @param path
#' @param dataPath
#' @param debug
#' @param cluster
#' @param chrominfo_file
#'
#' @return returnList
#' @export
#'
#' @examples

f_plotMetageneProfile_onePoint=function(binnedChip,binnedInput,tag="TSS",path=getwd(),debug=FALSE,plotName="NA")
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
	pdf(file=file.path(path,paste(plotName,"ChIP_Input_",tag,".pdf",sep="")),width=10, height=7)
	    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	    matplot(x=as.numeric(rownames(all.noNorm)),y=all.noNorm, type="l", lwd=2, lty=1,
	    col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main=tag,xaxt='n')
		abline(v=0,lty=2,col="darkgrey", lwd=2) ##plot TSS	    
		plotPoints=c(-2000,-1000,-500,500,1000,2000) ##plot remaining be
		abline(v=plotPoints,lty=3,col="darkgrey", lwd=2)
		axis(side = 1, at = sort(c(plotPoints,0)), labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
	    legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	dev.off()


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

	pdf(file=file.path(path,paste(plotName,"Normalized_",tag,".pdf",sep="")),width=10, height=7)
	par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	plot(x=as.numeric(names(all.Norm)),y=all.Norm, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="mean log2 enrichment (signal/input)",
		main=paste("normalized",tag,sep=" "), xaxt='n')#,cex.axis=1.3,cex.lab=1.3)
	abline(v=0,lty=2,col="darkgrey", lwd=2) ##plot TSS
	abline(v=plotPoints,lty=3,col="darkgrey", lwd=2)
	axis(side = 1, at = sort(c(plotPoints,0)), labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
	legend(x="topleft", fill=colori, legend=colnames(all.noNorm),bg="white",cex=0.8)
	dev.off()

	result=NULL
	result=rbind(hotSpotsValues,maxAucValues,variabilityValues,hotSpotsValuesNorm,maxAucValuesNorm,variabilityValuesNorm)
	return(result)

}


# chip_auc_TES_1 
# chip_auc_TES_2 
# chip_auc_TES_3 
# chip_auc_TES_4 
# chip_auc_TES_5 
# chip_auc_TES_6 
# chip_auc_TSS_1 
# chip_auc_TSS_2 
# chip_auc_TSS_3 
# chip_auc_TSS_4 
# chip_auc_TSS_5 
# chip_auc_TSS_6 
# chip_dispersion_TES_-1000_0% 
# chip_dispersion_TES_-1000_25% 
# chip_dispersion_TES_-1000_50% 
# chip_dispersion_TES_-1000_75% 
# chip_dispersion_TES_-1000_sd 
# chip_dispersion_TES_-1000_variance 
# chip_dispersion_TES_-2000_0% 
# chip_dispersion_TES_-2000_25% 
# chip_dispersion_TES_-2000_50% 
# chip_dispersion_TES_-2000_75% 
# chip_dispersion_TES_-2000_sd 
# chip_dispersion_TES_-2000_variance 
# chip_dispersion_TES_-500_0% 
# chip_dispersion_TES_-500_25% 
# chip_dispersion_TES_-500_50% 
# chip_dispersion_TES_-500_75% 
# chip_dispersion_TES_-500_sd 
# chip_dispersion_TES_-500_variance 
# chip_dispersion_TSS_-1000_0% 
# chip_dispersion_TSS_-1000_25% 
# chip_dispersion_TSS_-1000_50% 
# chip_dispersion_TSS_-1000_75% 
# chip_dispersion_TSS_-1000_sd 
# chip_dispersion_TSS_-1000_variance 
# chip_dispersion_TSS_-2000_0% chip_dispersion_TSS_-2000_25% chip_dispersion_TSS_-2000_50% 
# chip_dispersion_TSS_-2000_75% chip_dispersion_TSS_-2000_sd chip_dispersion_TSS_-2000_variance 
# chip_dispersion_TSS_-500_0% chip_dispersion_TSS_-500_25% chip_dispersion_TSS_-500_50% chip_dispersion_TSS_-500_75% 
# chip_dispersion_TSS_-500_sd chip_dispersion_TSS_-500_variance 
# chip_hotSpots_TES_-1000 
# chip_hotSpots_TES_-2000 
# chip_hotSpots_TES_-500 
# chip_hotSpots_TES_0 
# chip_hotSpots_TES_1000 
# chip_hotSpots_TES_2000 
# chip_hotSpots_TES_500 
# chip_hotSpots_TSS_-1000 
# chip_hotSpots_TSS_-2000 
# chip_hotSpots_TSS_-500 
# chip_hotSpots_TSS_0 
# chip_hotSpots_TSS_1000 
# chip_hotSpots_TSS_2000 
# chip_hotSpots_TSS_500 
# chip_localMax_TES_1_x 
# chip_localMax_TES_1_y 
# chip_localMax_TES_2_x 
# chip_localMax_TES_2_y chip_localMax_TES_3_x chip_localMax_TES_3_y chip_localMax_TES_4_x chip_localMax_TES_4_y chip_localMax_TES_5_x 
# chip_localMax_TES_5_y chip_localMax_TES_6_x chip_localMax_TES_6_y chip_localMax_TSS_1_x chip_localMax_TSS_1_y chip_localMax_TSS_2_x 
# chip_localMax_TSS_2_y chip_localMax_TSS_3_x chip_localMax_TSS_3_y chip_localMax_TSS_4_x chip_localMax_TSS_4_y chip_localMax_TSS_5_x 
# chip_localMax_TSS_5_y chip_localMax_TSS_6_x chip_localMax_TSS_6_y 
# input_auc_TES_1 input_auc_TES_2 input_auc_TES_3 input_auc_TES_4 input_auc_TES_5 input_auc_TES_6 input_auc_TSS_1 
# input_auc_TSS_2 input_auc_TSS_3 input_auc_TSS_4 input_auc_TSS_5 input_auc_TSS_6 
# input_dispersion_TES_-1000_0% input_dispersion_TES_-1000_25% input_dispersion_TES_-1000_50% input_dispersion_TES_-1000_75% 
# input_dispersion_TES_-1000_sd input_dispersion_TES_-1000_variance input_dispersion_TES_-2000_0% input_dispersion_TES_-2000_25% input_dispersion_TES_-2000_50% 
# input_dispersion_TES_-2000_75% input_dispersion_TES_-2000_sd input_dispersion_TES_-2000_variance input_dispersion_TES_-500_0% input_dispersion_TES_-500_25% 
# input_dispersion_TES_-500_50% input_dispersion_TES_-500_75% input_dispersion_TES_-500_sd input_dispersion_TES_-500_variance input_dispersion_TSS_-1000_0% 
# input_dispersion_TSS_-1000_25% input_dispersion_TSS_-1000_50% input_dispersion_TSS_-1000_75% input_dispersion_TSS_-1000_sd input_dispersion_TSS_-1000_variance 
# input_dispersion_TSS_-2000_0% input_dispersion_TSS_-2000_25% input_dispersion_TSS_-2000_50% input_dispersion_TSS_-2000_75% input_dispersion_TSS_-2000_sd 
# input_dispersion_TSS_-2000_variance input_dispersion_TSS_-500_0% input_dispersion_TSS_-500_25% input_dispersion_TSS_-500_50% input_dispersion_TSS_-500_75% 
# input_dispersion_TSS_-500_sd input_dispersion_TSS_-500_variance input_hotSpots_TES_-1000 input_hotSpots_TES_-2000 input_hotSpots_TES_-500 input_hotSpots_TES_0 
# input_hotSpots_TES_1000 input_hotSpots_TES_2000 input_hotSpots_TES_500 input_hotSpots_TSS_-1000 input_hotSpots_TSS_-2000 input_hotSpots_TSS_-500 
# input_hotSpots_TSS_0 input_hotSpots_TSS_1000 input_hotSpots_TSS_2000 input_hotSpots_TSS_500 
# input_localMax_TES_1_x input_localMax_TES_1_y input_localMax_TES_2_x input_localMax_TES_2_y input_localMax_TES_3_x input_localMax_TES_3_y 
# input_localMax_TES_4_x input_localMax_TES_4_y input_localMax_TES_5_x input_localMax_TES_5_y input_localMax_TES_6_x input_localMax_TES_6_y 
# input_localMax_TSS_1_x input_localMax_TSS_1_y input_localMax_TSS_2_x input_localMax_TSS_2_y input_localMax_TSS_3_x input_localMax_TSS_3_y 
# input_localMax_TSS_4_x input_localMax_TSS_4_y input_localMax_TSS_5_x input_localMax_TSS_5_y input_localMax_TSS_6_x input_localMax_TSS_6_y 
# norm_auc_TES_1 norm_auc_TES_2 norm_auc_TES_3 norm_auc_TES_4 norm_auc_TES_5 norm_auc_TES_6 norm_auc_TSS_1 norm_auc_TSS_2 norm_auc_TSS_3 
# norm_auc_TSS_4 norm_auc_TSS_5 norm_auc_TSS_6 
# norm_dispersion_TES_-1000_0% norm_dispersion_TES_-1000_25% norm_dispersion_TES_-1000_50% norm_dispersion_TES_-1000_75% norm_dispersion_TES_-1000_sd 
# norm_dispersion_TES_-1000_variance norm_dispersion_TES_-2000_0% norm_dispersion_TES_-2000_25% norm_dispersion_TES_-2000_50% norm_dispersion_TES_-2000_75% 
# norm_dispersion_TES_-2000_sd norm_dispersion_TES_-2000_variance norm_dispersion_TES_-500_0% norm_dispersion_TES_-500_25% norm_dispersion_TES_-500_50% 
# norm_dispersion_TES_-500_75% norm_dispersion_TES_-500_sd norm_dispersion_TES_-500_variance norm_dispersion_TSS_-1000_0% norm_dispersion_TSS_-1000_25% 
# norm_dispersion_TSS_-1000_50% norm_dispersion_TSS_-1000_75% norm_dispersion_TSS_-1000_sd norm_dispersion_TSS_-1000_variance norm_dispersion_TSS_-2000_0% 
# norm_dispersion_TSS_-2000_25% norm_dispersion_TSS_-2000_50% norm_dispersion_TSS_-2000_75% norm_dispersion_TSS_-2000_sd norm_dispersion_TSS_-2000_variance 
# norm_dispersion_TSS_-500_0% norm_dispersion_TSS_-500_25% norm_dispersion_TSS_-500_50% norm_dispersion_TSS_-500_75% norm_dispersion_TSS_-500_sd 
# norm_dispersion_TSS_-500_variance norm_hotSpots_TES_-1000 norm_hotSpots_TES_-2000 norm_hotSpots_TES_-500 norm_hotSpots_TES_0 norm_hotSpots_TES_1000 
# norm_hotSpots_TES_2000 norm_hotSpots_TES_500 norm_hotSpots_TSS_-1000 norm_hotSpots_TSS_-2000 norm_hotSpots_TSS_-500 norm_hotSpots_TSS_0 norm_hotSpots_TSS_1000 
# norm_hotSpots_TSS_2000 norm_hotSpots_TSS_500 
# norm_localMax_TES_1_x norm_localMax_TES_1_y norm_localMax_TES_2_x norm_localMax_TES_2_y norm_localMax_TES_3_x norm_localMax_TES_3_y 
# norm_localMax_TES_4_x norm_localMax_TES_4_y norm_localMax_TES_5_x norm_localMax_TES_5_y norm_localMax_TES_6_x norm_localMax_TES_6_y 
# norm_localMax_TSS_1_x norm_localMax_TSS_1_y norm_localMax_TSS_2_x norm_localMax_TSS_2_y norm_localMax_TSS_3_x norm_localMax_TSS_3_y 
# norm_localMax_TSS_4_x norm_localMax_TSS_4_y norm_localMax_TSS_5_x norm_localMax_TSS_5_y norm_localMax_TSS_6_x norm_localMax_TSS_6_y n