library("girafe")#, lib.loc="/home/clivi/programs/R-3.2.2/library/")
#library("Rmpi")#, lib.loc="/home/clivi/programs/R-3.2.2/library/")
library("caTools")#, lib.loc="/home/clivi/programs/R-3.2.2/library/")
library("spp")#, lib.loc="/home/clivi/programs/R-3.2.2/library/")
#library("snow")#, lib.loc="/home/clivi/programs/R-3.2.2/library/")

#cluster=makeCluster(8, type="MPI")
cluster=NULL


arg <- commandArgs()
print(arg)
#3sampleIndex=as.integer(arg[7])
path=arg[6]
datafilename=arg[7]
print(datafilename)
#strandShift=as.integer(arg[8])


#sampleIndex=as.integer(arg[7])
#path=arg[6]

print(path)



source(paste(path,"GlobalParameters.R",sep=""))
sampleinfo_file<-paste(path,"matchlist.txt",sep="")
sampleinfo<-read.table(sampleinfo_file,  header=TRUE, quote="", stringsAsFactors=FALSE)

t0<-proc.time()[3]

workingdir<-paste(path,"Chance/",sep="")
outputdir<-paste(workingdir,"out/",sep="")
timedir<-paste(workingdir,"time_stamps/",sep="")
plotsdir<-paste(workingdir,"plots/",sep="")

##step 00: binding characteristics
#get dataset
#datafilename=sampleinfo$Filename[sampleIndex]
print(datafilename)

#sampleIndex=16
##step 00: binding characteristics
#get dataset
sampleIndex=which(sampleinfo$Filename==datafilename)
print(sampleIndex)
TFname=sampleinfo$IG[sampleIndex]


inputIndex=sampleIndex
inputname=sampleinfo$ControlName[inputIndex]

#wig.bw <- 100;
#wig.step <- 50;
#wig.bw
#wig.step


#wig.bw<-MetaprofilesSmoothingBandwidth
#wig.step<-MetaprofilesSmoothingStep



#print("Load control file...")
#newControl=NULL
#if (length(unlist(strsplit(inputname,"\\.")))>1)
#{
#	print("Concat multiple inputs...")
#	for (element in unlist(strsplit(inputname,"\\.")))
#	{	
##		if (element!="")
#		{
#			if (is.null(newControl))
#			{	
				##if newControl is empty, fill it with the first input-data#
#				print(element)
#				load(file.path(path, paste(element,".RData",sep="")))
#				newControl<-input.data
#				print("original size")
#				for (i in seq(length(newControl$tags)))
#				{
#					print(length(newControl$tags[[i]]))
#				}
				#newControl.tags=(sapply(newControl$tags,c))
				#newControl.q=(sapply(newControl$quality,c))
				#print(length(newControl.tags[[25]]))
				#print(length(newControl.q[[25]]))
#			}else{
#				##otherwise add the input-data and sort
#				print(element)
#				load(file.path(path, paste(element,".RData",sep="")))
#				temp<-input.data
#				#temp.tags=(sapply(temp$tags,c))
				#temp.q=(sapply(temp$quality,c))
#
#				for (i in seq(length(newControl$tags)))
#				{
#					##appending tags and quality elements of all inputs to newControl
#					newControl$tags[[i]]=append(newControl$tags[[i]],temp$tags[[i]])
#					newControl$quality[[i]]=append(newControl$quality[[i]],temp$quality[[i]])
#				}
#
#			}
#		}
#	}
#	print("final size")
#	for (i in seq(length(newControl$tags)))
#	{
#		print(length(newControl$tags[[i]]))
#	}
#
	##sort tags and quality of the final NewControl to get an increasing list independently from the sign:  -57 -95 -112 151 159 166 169 -217...
#	print("Sort tags and quality flags...")
#	for (i in seq(length(newControl$tags)))
#	{
#		index=sort.list(abs(newControl$tags[[i]]))
#		newControl$tags[[i]]=newControl$tags[[i]][index]	
#		newControl$quality[[i]]=newControl$quality[[i]][index]	
#	}

#	input.data=newControl
#}else{
#
#	#inputIndex=as.integer(arg[8])
#	print("Single input...")
#	load(file.path(path, paste(inputname,".RData",sep="")))
#	input.data<-input.data
#}






#load inputdata

#print("Filter tags")
#if (select.informative.tags_filter) {
#      print("select.informative.tags filter")
#      #load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
#      input.dataSelected <- select.informative.tags(input.data, binding.characteristics)
#} else {
#      print("SKIP select.informative.tags filter")
#      input.dataSelected<-input.data$tags
#}



#print("Smooth tag density")
#ts <- sum(unlist(lapply(input.dataSelected,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
#smoothed.density <- get.smoothed.tag.density(input.dataSelected, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift)
#smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })

#save(smoothed.density,file=paste(path,"TagDensity_",inputname,".RData",sep=""))
load(paste(path,"TagDensity_",datafilename,"_",inputname,".RData",sep=""))

print("shorten frame")
##shorten frame
input.new=NULL
chrl=names(smoothed.density)
for (i in chrl)
{
	print(i)
	##appending tags and quality elements of all inputs to newControl
	input.new[[i]]$x=smoothed.density[[i]]$x[seq.int(1, length(smoothed.density[[i]]$x), 6)]
	input.new[[i]]$y=smoothed.density[[i]]$y[seq.int(1, length(smoothed.density[[i]]$y), 6)]
}

input.smoothed.density=input.new


load(paste(path,"TagDensity_",datafilename,".RData",sep=""))
new=NULL
chrl=names(smoothed.density)
for (i in chrl)
{
	print(i)
	##appending tags and quality elements of all inputs to newControl
	new[[i]]$x=smoothed.density[[i]]$x[seq.int(1, length(smoothed.density[[i]]$x), 6)]
	new[[i]]$y=smoothed.density[[i]]$y[seq.int(1, length(smoothed.density[[i]]$y), 6)]
}

smoothed.density=new





#tag.shift <- round(strandShift/2)

Y<-unlist(lapply(smoothed.density, FUN=function(list_element) {return(list_element$y)}))
X<-unlist(lapply(input.smoothed.density, FUN=function(list_element) {return(list_element$y)}))


Xsorted=sort(X)
Ysorted=sort(Y)

head(Xsorted)
tail(Xsorted)
head(Ysorted)
tail(Ysorted)
print("partial summing")

BINS_NUMBER<-1e4
print("input ...")
cumulative_sum_Y<-cumsum(Ysorted)
cumulative_sum_Y_bins<-quantile(cumulative_sum_Y, probs=seq(0,1,(1/BINS_NUMBER)))
pj<-(cumulative_sum_Y_bins/cumulative_sum_Y[length(cumulative_sum_Y)])
fcompY=data.frame(x=seq(0,1,(1/BINS_NUMBER)),pj=pj)
rownames(fcompY)<-NULL

Yfraction_of_reads_intop_1percent_bins<-(1-fcompY[(which(fcompY$x>=0.99)[1]),"pj"])
Yfraction_of_bins_without_reads<-fcompY$x[(which(fcompY$pj>0)[1])]


print("chip ...")
cumulative_sum_X<-cumsum(Xsorted)
cumulative_sum_X_bins<-quantile(cumulative_sum_X, probs=seq(0,1,(1/BINS_NUMBER)))
pj<-(cumulative_sum_X_bins/cumulative_sum_X[length(cumulative_sum_X)])
fcompX=data.frame(x=seq(0,1,(1/BINS_NUMBER)),pj=pj)
rownames(fcompX)<-NULL

Xfraction_of_reads_intop_1percent_bins<-(1-fcompX[(which(fcompX$x>=0.99)[1]),"pj"])
Xfraction_of_bins_without_reads<-fcompX$x[(which(fcompX$pj>0)[1])]


plotname=paste(plotsdir,paste(datafilename,inputname,sep="_"),".pdf",sep="")
pdf(plotname)
plot(fcompY,type="l",col="blue",lwd=2,xlab="Percentage of bins",ylab="Percentage of tags",main=TFname)
lines(fcompX,col="red",lwd=2)
arrowx=fcompY[which.max(abs(fcompX$pj-fcompY$pj)),]$x
abline(v =arrowx, col = "green",lty=2,lwd=2)
#abline(h=schneidePunktY,col='cyan',lty=2,lwd=2)
#abline(v=schneidePunktX,col='cyan',lty=2,lwd=2)
legend("topleft", legend = c("Input","ChIP"), fill = c("red","blue"))
dev.off()

sign_sign=sign(max(fcompX$pj-fcompY$pj)) ##chip-input
cross_pointX=fcompX[which((fcompX$pj==fcompY$pj)&(fcompX$pj!=0)),]$x
cross_pointY_input=fcompX[which((fcompX$pj==fcompY$pj)&(fcompX$pj!=0)),]$pj
cross_pointY_chip=fcompY[which((fcompX$pj==fcompY$pj)&(fcompX$pj!=0)),]$pj

final=merge(fcompX,fcompY,by="x")
colnames(final)=c("x","y_input","y_chip")


filename=paste(outputdir,datafilename,"_",TFname,".txt",sep="")
file.remove(filename)
write(TFname,file=filename,append=T)
write("#Max distance between curves",file=filename,append=T)
write(paste("X-axis",arrowx,sep=" "),file=filename,append=T)
write(paste("Y-Input",fcompX[which.max(abs(fcompX$pj-fcompY$pj)),]$pj,sep=" "),file=filename,append=T)
write(paste("Y-Chip",fcompY[which.max(abs(fcompX$pj-fcompY$pj)),]$pj,sep=" "),file=filename,append=T)
write(paste("sign: chip-input",sign_sign,sep=" "),file=filename,append=T)
write(paste("FractionReadsTopbins_chip",Yfraction_of_reads_intop_1percent_bins,sep=" "),file=filename,append=T)
write(paste("FractionReadsTopbins_input",Xfraction_of_reads_intop_1percent_bins,sep=" "),file=filename,append=T)
write(paste("Fractions_without_reads_chip",Yfraction_of_bins_without_reads,sep=" "),file=filename,append=T)
write(paste("Fractions_without_reads_input",Xfraction_of_bins_without_reads,sep=" "),file=filename,append=T)
write(paste("CrossPoint_X",cross_pointX,sep=" "),file=filename,append=T)
write(paste("CrossPoint_Y_chip",cross_pointY_chip,sep=" "),file=filename,append=T)
write(paste("CrossPoint_Y_input",cross_pointY_input,sep=" "),file=filename,append=T)
write.table(final,file=filename, append = T, sep = " ",row.names = F)


t1<-proc.time()[3]
deltat=t1-t0

write(paste(datafilename,deltat,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""))


