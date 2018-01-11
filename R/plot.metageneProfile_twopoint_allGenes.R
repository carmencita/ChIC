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
workingdir<-path #paste(path,"/MetaGene",sep="")
timedir<-paste(workingdir,"/time_stamps",sep="")
outputdir<-paste(workingdir,"/out",sep="")
plotsdir<-paste(workingdir,"/plots",sep="")



inputIndex=sampleIndex
inputname=sampleinfo$ControlName[inputIndex]



#filterGenes=read.delim(paste(path,"HouseKeepingGenes_wikicell.org.txt",sep=""),header=TRUE) #X.Gene and geneSymbol
#filterGenes=filterGenes$X.Gene


print(inputname)
file_metaprofili<-paste(storagedir,"/two.point.scaling_",datafilename, ".RData", sep="")
load(file_metaprofili)

outname=paste(outputdir,"/",datafilename,"_",TFname,".txt",sep="")


chip_genelist_pre <- log2(do.call(rbind,gd.gen)+psc)
#chip_genelist= chip_genelist_pre[which((rownames(chip_genelist_pre) %in% filterGenes)==TRUE),]
chip_genelist= chip_genelist_pre

#gfp_preMeans<-log2(do.call(rbind,gd.gfp)+psc)
input_preMeans_pre<-log2(do.call(rbind,gd.gfp)+psc)
#input_preMeans= input_preMeans_pre[which((rownames(input_preMeans_pre) %in% filterGenes)==TRUE),]
input_preMeans= input_preMeans_pre

input.noNorm<-colMeans(input_preMeans,na.rm=T)
chip_genelist.noNorm <- colMeans(chip_genelist,na.rm=T)
input.noNorm<-cbind(input.noNorm, chip_genelist.noNorm)
colnames(input.noNorm)<-c(datafilename, inputname)



hotSpots=NULL ##takes values at different predefined points
for (i  in break_points_2P){
	print(i)
	if (length(which(as.integer(row.names(input.noNorm))==i))<1){
		x_bin= which(as.integer(row.names(input.noNorm))==i+estimated_bin_size_2P/2)}else{
		print("ok")
		x_bin= which(as.integer(row.names(input.noNorm))==i)
	}
	hotSpots=rbind(hotSpots,c(i, input.noNorm[x_bin,][1], input.noNorm[x_bin,][2]))
}

write.table(cbind("chip",hotSpots[,1],hotSpots[,2],0,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
write.table(cbind("input",hotSpots[,1],hotSpots[,3],0,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)



##local maxima and area in all the predefined regions
##local maxima and area in all the predefined regions
localMaxima_auc=NULL
for (i in seq(length(break_points_2P)-1)){
	if (length(which(as.integer(row.names(input.noNorm))==break_points_2P[i]))<1){
		xi= which(as.integer(row.names(input.noNorm))==break_points_2P[i]+estimated_bin_size_2P/2)}else{
		xi= which(as.integer(row.names(input.noNorm))==break_points_2P[i])
	}
	if (length(which(as.integer(row.names(input.noNorm))==break_points_2P[i+1]))<1){
		xj= which(as.integer(row.names(input.noNorm))==break_points_2P[i+1]+estimated_bin_size_2P/2)}else{
		xj= which(as.integer(row.names(input.noNorm))==break_points_2P[i+1])
	}

	subsettt= input.noNorm[xi:xj,]
	l1=max(subsettt[,1])
	l2=max(subsettt[,2])
	print(break_points_2P[i])
	localMaxima_auc=rbind(localMaxima_auc,c(rownames(subsettt)[which(subsettt[,1]==l1)], l1, sum((subsettt[,1])*estimated_bin_size_2P),rownames(subsettt)[which(subsettt[,2]==l2)],l2,sum((subsettt[,2])*estimated_bin_size_2P)))
}

write.table(cbind("chip",localMaxima_auc[,1],localMaxima_auc[,2],0,"localMax twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
write.table(cbind("input",localMaxima_auc[,4],localMaxima_auc[,5],0,"localMax twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)

write.table(cbind("chip",localMaxima_auc[,3],0,"auc twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
write.table(cbind("input",localMaxima_auc[,6],0,"auc twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)


##make plots
colori<-c(rev(rainbow(ncol(input.noNorm)-1)), "black")
pdf(file=paste(plotsdir,"/NOnormalizaztion_",TFname,"_",datafilename ,".pdf", sep=""))
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    matplot(x=as.numeric(rownames(input.noNorm)),y=input.noNorm, type="l", lwd=2, lty=1, col=colori,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(TFname,"density"))
#    matplot(as.numeric(names(gfp.a_noNorm)),gfp.a_noNorm,col="red",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(nam,"density"));
#    plot(as.numeric(names(gfp.a)),gfp.a,col="orange",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 GRO normalized density",main=paste(nam,"density"));
#     lines(as.numeric(names(gfp.a_noNorm)),gfp.a_noNorm,col="red",lwd=2)
#    lines(as.numeric(names(gen.a_noNorm)),gen.a_noNorm,col="black",lwd=2)

    abline(v=c(0,totalGeneLength),lty=2,col="black", lwd=2)
    abline(v=c(inner_margin,totalGeneLength-inner_margin),lty=3,col="black", lwd=2)
    legend(x="topleft", fill=colori, legend=colnames(input.noNorm),bg="white",cex=0.8)

dev.off()


##normalized plot and values


common_genes<-rownames(input_preMeans)[rownames(input_preMeans) %in% rownames(chip_genelist)]
input.Norm<-colMeans(t(t(input_preMeans[common_genes,])-t(chip_genelist[common_genes,])),na.rm=T)

hotSpots=NULL ##takes values at different predefined points
for (i  in break_points_2P){
	if (length(which(as.integer(names(input.Norm))==i))<1){
		x_bin= which(as.integer(names(input.Norm))==i+estimated_bin_size_2P/2)}else{
		x_bin= which(as.integer(names(input.Norm))==i)
	}
	hotSpots=rbind(hotSpots,c(i, input.Norm[x_bin][1]))

}

write.table(cbind("norm",hotSpots[,1],hotSpots[,2],1,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)

##local maxima and area in all the predefined regions

localMaxima_auc=NULL
for (i in seq(length(break_points_2P)-1)){
	if (length(which(as.integer(names(input.Norm))==break_points_2P[i]))<1){
		xi= which(as.integer(names(input.Norm))==break_points_2P[i]+estimated_bin_size_2P/2)}else{
		xi= which(as.integer(names(input.Norm))==break_points_2P[i])
	}
	if (length(which(as.integer(names(input.Norm))==break_points_2P[i+1]))<1){
		xj= which(as.integer(names(input.Norm))==break_points_2P[i+1]+estimated_bin_size_2P/2)}else{
		xj= which(as.integer(names(input.Norm))==break_points_2P[i+1])
	}

	subsettt= data.frame(input.Norm[xi:xj])
	l1=max(subsettt)
	localMaxima_auc=rbind(localMaxima_auc,c(rownames(subsettt)[which(subsettt[,1]==l1)], l1, sum((subsettt)*estimated_bin_size_2P)))
}


write.table(cbind("norm",localMaxima_auc[,1],localMaxima_auc[,2],1,"localMax twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
write.table(cbind("norm",localMaxima_auc[,3],1,"auc twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)




pdf(file=paste(plotsdir,"/Normalized_",TFname,"_", datafilename, ".pdf", sep=""))
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    matplot(x=as.numeric(names(input.Norm)),y=input.Norm, type="l", lwd=2, lty=1, col="orange",xlab="metagene coordinates",ylab="tag density logration",main=paste(TFname,"normalized"))
#    matplot(as.numeric(names(gfp.a_noNorm)),gfp.a_noNorm,col="red",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 tag density",main=paste(nam,"density"));
#    plot(as.numeric(names(gfp.a)),gfp.a,col="orange",type='l',ylim=rng,lwd=2,xlab="metagene coordinates",ylab="mean of log2 GRO normalized density",main=paste(nam,"density"));
#     lines(as.numeric(names(gfp.a_noNorm)),gfp.a_noNorm,col="red",lwd=2)
#    lines(as.numeric(names(gen.a_noNorm)),gen.a_noNorm,col="black",lwd=2)

    abline(v=c(0,totalGeneLength),lty=2,col="black", lwd=2)
    abline(v=c(inner_margin,totalGeneLength-inner_margin),lty=3,col="black", lwd=2)
    legend(x="topleft", fill=colori, legend=colnames(input.noNorm),bg="white",cex=0.8)

dev.off()


t1<-proc.time()[3]
deltat=t1-t0
write(paste("twoPointScaling",datafilename,deltat,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""),append=TRUE)


