arg <- commandArgs()
print(arg)

path=arg[6]
datafilename=arg[7]

source(paste(path,"GlobalParameters.R",sep=""))
sampleinfo_file<-paste(path,"matchlist.txt",sep="")
sampleinfo<-read.table(sampleinfo_file,  header=TRUE, quote="", stringsAsFactors=FALSE)


sampleIndex=which(sampleinfo$Filename==datafilename)
print(sampleIndex)
TFname=sampleinfo$IG[sampleIndex]
#file_metaprofili<-metaprofiles.all.file_list[sampleIndex]
datafilename=sampleinfo$Filename[sampleIndex]
print(datafilename)
TFname=sampleinfo$IG[sampleIndex]



psc <- 1; # pseudocount # required to avoid log2 of 0

t0<-proc.time()[3]
workingdir<-paste(path,"/MetaGene",sep="")
outputdir<-paste(workingdir,"/out",sep="")
timedir<-paste(workingdir,"/time_stamps",sep="")
plotsdir<-paste(workingdir,"/plots",sep="")

outname=paste(outputdir,"/",datafilename,"_",TFname,".txt",sep="")
file.remove(outname)

inputIndex=sampleIndex
inputname=sampleinfo$ControlName[inputIndex]


print(inputname)



filterGenes=read.delim(paste(path,"HouseKeepingGenes_wikicell.org.txt",sep=""),header=TRUE) #X.Gene and geneSymbol
filterGenes=filterGenes$X.Gene

for (go in c("TSS","TES"))
{	
	if (go=="TSS")
	{
		file.metaprofili<-paste(outputdir,"/one.point.scalingTSS_",datafilename, ".RData", sep="")
		load(file.metaprofili)
		gd.gen=gd.gen_TSS
	}else{
		
		file.metaprofili<-paste(outputdir,"/one.point.scalingTES_",datafilename, ".RData", sep="")
		load(file.metaprofili)
		gd.gen=gd.gen_TES
	}

	print(go)
	chip_genelist_pre <- log2(do.call(rbind,gd.gen)+psc)
	chip_genelist= chip_genelist_pre[which((rownames(chip_genelist_pre) %in% filterGenes)==TRUE),]

	input_preMeans_pre<-log2(do.call(rbind,gd.gfp)+psc)
	input_preMeans= input_preMeans_pre[which((rownames(input_preMeans_pre) %in% filterGenes)==TRUE),]


	##
	input.noNorm<-colMeans(input_preMeans,na.rm=T)
	chip_genelist.noNorm <- colMeans(chip_genelist,na.rm=T)
	input.noNorm<-cbind(input.noNorm, chip_genelist.noNorm)
	colnames(input.noNorm)<-c(datafilename, inputname)

	hotSpots=NULL ##takes values at different predefined points
	for (i  in break_points){
		print(i)
		if (length(which(as.integer(row.names(input.noNorm))==i))<1){
			x_bin= which(as.integer(row.names(input.noNorm))==i+estimated_bin_size_1P/2)}else{
			x_bin= which(as.integer(row.names(input.noNorm))==i)
		}
		hotSpots=rbind(hotSpots,c(i, input.noNorm[x_bin,][1], input.noNorm[x_bin,][2]))

	}
	#colnames(hotSpots)=c("x","chip_y","input_y")
	if (go=="TSS")
	{
		write.table(cbind("chip",hotSpots[,1],hotSpots[,2],0,"hotSpots ",go),file=outname,row.names = FALSE,col.names=FALSE, quote = FALSE)
	}else{
		write.table(cbind("chip",hotSpots[,1],hotSpots[,2],0,"hotSpots ",go),file=outname,row.names = FALSE,col.names=FALSE, quote = FALSE,append=TRUE)	
	}
	write.table(cbind("input",hotSpots[,1],hotSpots[,3],0,"hotSpots ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	##local maxima and area in all the predefined regions
	localMaxima_auc=NULL
	for (i in seq(length(break_points)-1)){
		if (length(which(as.integer(row.names(input.noNorm))==break_points[i]))<1){
			xi= which(as.integer(row.names(input.noNorm))==break_points[i]+estimated_bin_size_1P/2)}else{
			xi= which(as.integer(row.names(input.noNorm))==break_points[i])
		}
		if (length(which(as.integer(row.names(input.noNorm))==break_points[i+1]))<1){
			xj= which(as.integer(row.names(input.noNorm))==break_points[i+1]+estimated_bin_size_1P/2)}else{
			xj= which(as.integer(row.names(input.noNorm))==break_points[i+1])
		}

		subsettt= input.noNorm[xi:xj,]
		l1=max(subsettt[,1])
		l2=max(subsettt[,2])
		print(break_points[i])
		localMaxima_auc=rbind(localMaxima_auc,c(rownames(subsettt)[which(subsettt[,1]==l1)], l1, sum((subsettt[,1])*estimated_bin_size_1P),rownames(subsettt)[which(subsettt[,2]==l2)],l2,sum((subsettt[,2])*estimated_bin_size_1P)))
	}
	#colnames(localMaxima_auc)=c("#noNorm_localMaxima_chip_x","y","auc","input_x","y","auc")
	write.table(cbind("chip",localMaxima_auc[,1],localMaxima_auc[,2],0,"localMax ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	write.table(cbind("input",localMaxima_auc[,4],localMaxima_auc[,5],0,"localMax ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)

	write.table(cbind("chip",localMaxima_auc[,3],0,"auc ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	write.table(cbind("input",localMaxima_auc[,6],0,"auc ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)


	write("#noNorm break_point var sd 0% 25% 50% 75%", file=outname,append=TRUE)
	dispersion=NULL
	midpoint=as.integer(length(break_points)/2)+1
	for (i in seq(as.integer(length(break_points)/2)))
	{
		start=which(as.integer(row.names(input.noNorm))==break_points[midpoint-i])
		end=which(as.integer(row.names(input.noNorm))==break_points[midpoint+i])
		interval=input.noNorm[start:end,]	

		dispersion=paste(var(interval[,1]),sd(interval[,1]),quantile(interval[,1])[1],quantile(interval[,1])[2],quantile(interval[,1])[3],quantile(interval[,1])[4],quantile(interval[,1])[5],sep=" ")
		write(paste("chip",break_points[midpoint-i],dispersion,0,"dispersion ",go,sep=" "),file=outname,append=TRUE)

		dispersion=paste(var(interval[,2]),sd(interval[,2]),quantile(interval[,2])[1],quantile(interval[,2])[2],quantile(interval[,2])[3],quantile(interval[,2])[4],quantile(interval[,2])[5],sep=" ")
		write(paste("input",break_points[midpoint-i],dispersion,0,"dispersion ",go,sep=" "),file=outname,append=TRUE)
	}

	##make plots

	colori<-c(rev(rainbow(ncol(input.noNorm)-1)), "black")
	pdf(file=paste(plotsdir,"/NOnormalizaztion_",go,"_",TFname,"_",datafilename , ".pdf", sep=""), width=1024)
	    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	    matplot(x=as.numeric(rownames(input.noNorm)),y=input.noNorm, type="l", lwd=2, lty=1, col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(TFname,"density ",go))
	#    matplot(as.numeric(names(gfp.noNorm)),gfp.noNorm,col="red",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(nam,"density"));
	#    plot(as.numeric(names(gfp.a)),gfp.a,col="orange",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 GRO normalized density",main=paste(nam,"density"));
	#     lines(as.numeric(names(gfp.noNorm)),gfp.noNorm,col="red",lwd=2)
	#    lines(as.numeric(names(gen.noNorm)),gen.noNorm,col="black",lwd=2)

	    abline(v=0,lty=2,col="black", lwd=2)
	    abline(v=break_points,lty=3,col="black", lwd=2)
	    legend(x="topleft", fill=colori, legend=colnames(input.noNorm),bg="white",cex=0.8)
	dev.off()


	##normalized plot and values
	common_genes<-rownames(input_preMeans)[rownames(input_preMeans) %in% rownames(chip_genelist)]
	input.Norm<-colMeans(t(t(input_preMeans[common_genes,])-t(chip_genelist[common_genes,])),na.rm=T)

	hotSpots=NULL ##takes values at different predefined points
	for (i  in break_points){
		if (length(which(as.integer(names(input.Norm))==i))<1){
			x_bin= which(as.integer(names(input.Norm))==i+estimated_bin_size_1P/2)}else{
			x_bin= which(as.integer(names(input.Norm))==i)
		}
		hotSpots=rbind(hotSpots,c(i, input.Norm[x_bin][1]))

	}
	write.table(cbind("norm",hotSpots[,1],hotSpots[,2],1,"hotSpots ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)


	##local maxima and area in all the predefined regions
	localMaxima_auc=NULL
	for (i in seq(length(break_points)-1)){
		if (length(which(as.integer(names(input.Norm))==break_points[i]))<1){
			xi= which(as.integer(names(input.Norm))==break_points[i]+estimated_bin_size_1P/2)}else{
			xi= which(as.integer(names(input.Norm))==break_points[i])
		}
		if (length(which(as.integer(names(input.Norm))==break_points[i+1]))<1){
			xj= which(as.integer(names(input.Norm))==break_points[i+1]+estimated_bin_size_1P/2)}else{
			xj= which(as.integer(names(input.Norm))==break_points[i+1])
		}

		subsettt= data.frame(input.Norm[xi:xj])
		l1=max(subsettt)
		localMaxima_auc=rbind(localMaxima_auc,c(rownames(subsettt)[which(subsettt[,1]==l1)], l1, sum((subsettt)*estimated_bin_size_1P)))
	}
	write.table(cbind("norm",localMaxima_auc[,1],localMaxima_auc[,2],1,"localMax ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
	write.table(cbind("norm",localMaxima_auc[,3],1,"auc ",go),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)


	dispersion=NULL
	midpoint=as.integer(length(break_points)/2)+1
	for (i in seq(as.integer(length(break_points)/2)))
	{
		start=which(as.integer(names(input.Norm))==break_points[midpoint-i])
		end=which(as.integer(names(input.Norm))==break_points[midpoint+i])
		interval=input.Norm[start:end]	
		dispersion=paste(var(interval),sd(interval),quantile(interval)[1],quantile(interval)[2],quantile(interval)[3],quantile(interval)[4],quantile(interval)[5],sep=" ")
		write(paste("norm",break_points[midpoint-i],dispersion,1,"dispersion ",go,sep=" "),file=outname,append=TRUE)
	}



	pdf(file=paste(plotsdir,"/Normalized_",go,"_",TFname,"_", datafilename,  ".pdf", sep=""), width=1024)
	    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
	    matplot(x=as.numeric(names(input.Norm)),y=input.Norm, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="tag density logration",main=paste(TFname,"normalized ",go))
	#    matplot(as.numeric(names(gfp.noNorm)),gfp.noNorm,col="red",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(nam,"density"));
	#    plot(as.numeric(names(gfp.a)),gfp.a,col="orange",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 GRO normalized density",main=paste(nam,"density"));
	#     lines(as.numeric(names(gfp.noNorm)),gfp.noNorm,col="red",lwd=2)
	#    lines(as.numeric(names(gen.noNorm)),gen.noNorm,col="black",lwd=2)
	    abline(v=0,lty=2,col="black", lwd=2)
	    abline(v=break_points,lty=3,col="black", lwd=2)
	    legend(x="topleft", fill=colori, legend=colnames(input.noNorm),bg="white",cex=0.8)

	dev.off()

}

write(paste("onePointScaling",datafilename,deltat,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""),append=TRUE)

t1<-proc.time()[3]
deltat=t1-t0


#file.remove((file.metaprofili))
