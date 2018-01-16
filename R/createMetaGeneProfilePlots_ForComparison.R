##########DEVEL##########


#source("/lustre/data/FF/Carmen/Pipeline_allEncodeBeta/CinecaScripts/QC_pipeline/GlobalParameters.R")


# chipName="ENCFF000BBB"
# ll="H3K4me3"
# dataPath="/lustre/data/FF/Carmen/BitBucket/chic/data/Profiles"


# load("two.point.scaling_ENCFF000BBB.RData")
# #input is a list: list("twopoint"=twopoint,"TSS"=TSS,"TES"=TES))
# twopointElements=NULL
# twopointElements$Input=gd.gfp
# twopointElements$Chip=gd.gen

# load("one.point.scalingTES_ENCFF000BBB.RData")
# TESElements=NULL
# TESElements$Input=gd.gfp
# TESElements$Chip=gd.gen_TES


# load("one.point.scalingTSS_ENCFF000BBB.RData")
# #input is a list: list("twopoint"=twopoint,"TSS"=TSS,"TES"=TES))
# TSSElements=NULL
# TSSElements$Input=gd.gfp
# TSSElements$Chip=gd.gen_TSS

##########DEVEL##########



##private function
f_plotProfiles <- function(meanFrame, currentFrame , endung="TWO", absoluteMinMax, maintitel="title",ylab="mean of log2 tag density",savePlotPath=getwd())#, color="green") 
{
	filename=paste(maintitel,".pdf", sep="")
	pdf(file=file.path(savePlotPath,filename),width=10, height=7)
	#The standard error of the mean (SEM) is the standard deviation of the sample-mean's estimate of a population mean. 
	#(It can also be viewed as the standard deviation of the error in the sample mean with respect to the true mean, 
	#since the sample mean is an unbiased estimator.) SEM is usually estimated by the sample estimate of the population 
	#standard deviation (sample standard deviation) divided by the square root of the sample size (assuming statistical independence of the values in the sample)
	plot(x=c(min(meanFrame$x),max(meanFrame$x)),y=c(absoluteMinMax[1],absoluteMinMax[2]),type="n",xlab="metagene coordinates",ylab=ylab,main=maintitel,xaxt='n')
	polygon(x=c(meanFrame$x,rev(meanFrame$x)),y=c(meanFrame$mean+2*meanFrame$sderr,rev(meanFrame$mean)),col="lightblue",border=NA)
	polygon(x=c(meanFrame$x,rev(meanFrame$x)),y=c(meanFrame$mean-2*meanFrame$sderr,rev(meanFrame$mean)),col="lightblue",border=NA)
	lines(x=meanFrame$x,y=meanFrame$mean,col="black",lwd=2)
	lines(x=currentFrame$x,y=currentFrame$mean,col="red",lwd=2)
	if (endung=="TWO")
	{
		break_points=c(-2000,-1000,500,2500,4000)
		abline(v=c(0,totalGeneLength),lty=2,col="darkgrey", lwd=3)
	   	#abline(v=c(inner_margin,totalGeneLength-inner_margin),lty=3,col="darkgrey", lwd=2)
	   	abline(v=break_points,lty=3,col="darkgrey", lwd=2)
	   	axis(side = 1, at = sort(c(break_points,0,3000)), labels = c("-2KB","-1KB","TSS","500","500","TES","+1KB"))
	}else{
	   	TSSbreak_points=c(-2000,-1000,-500,500,1000,2000)
	   	if (endung=="TSS")
		{
			axis(side = 1, at = sort(c(TSSbreak_points,0)), labels = c("-2KB","-1KB","-500","TSS","500","+1KB","+2KB"))#,cex.axis=1.3)
		}else{
			axis(side = 1, at = sort(c(TSSbreak_points,0)), labels = c("-2KB","-1KB","-500","TES","500","+1KB","+2KB"))#,cex.axis=1.3)
		}
		abline(v=0,lty=2,col="darkgrey", lwd=2)
	    abline(v=TSSbreak_points,lty=3,col="darkgrey", lwd=2)
	}
	legend("topleft",legend=c("mean","2*stdErr"),fill=c("black","lightblue"),bg="white")
	dev.off()
}

f_plotValueDistribution = function(compendium,title,coordinateLine,savePlotPath=NULL)
{   #png(paste(profileclass,label,".png",sep=""))
        #pdf(paste(name,label,".pdf",sep=""))
    if (!is.null(savePlotPath))
    {
    	filename=file.path(savePlotPath,"PlotValueDistribution.pdf")
		pdf(file=filename,width=10, height=7)
	}
    mycolors= c("lightsteelblue1","lightsteelblue3")
	maxi=max(compendium)
    mini=min(compendium)

	d=density(compendium)
	plot(d,main=title,xlim=c(mini,maxi))
       
	median=median(compendium)
	qq=quantile(compendium,probs = c(0.05, 0.25,0.50,0.75,0.95))
        
	i <-  d$x >= (qq[1]) &  d$x <= (qq[2])
	lines(d$x, d$y)
	polygon(c(qq[1],d$x[i],qq[2]), c(0,d$y[i],0), col=mycolors[1])

	i <-  d$x >= (qq[2]) &  d$x <= (qq[3])
    lines(d$x, d$y)
    polygon(c(qq[2],d$x[i],qq[3]), c(0,d$y[i],0), col=mycolors[2])

    i <-  d$x >= (qq[3]) &  d$x <= (qq[4])
    lines(d$x, d$y)
    polygon(c(qq[3],d$x[i],qq[4]), c(0,d$y[i],0), col=mycolors[2])
        
    i <-  d$x >= (qq[4]) &  d$x <= (qq[5])
    lines(d$x, d$y)
    polygon(c(qq[4],d$x[i],qq[5]), c(0,d$y[i],0), col=mycolors[1])

    abline(v =coordinateLine, ,col="red",lwd=3,lty=2)

    if (!is.null(savePlotPath))
    {
    	dev.off()
    	print(paste("pdf saved under",filename,sep=" "))
	}
    
}

#'@title Function to create metage plots for comparison
#'
#' @description
#' QC-metrics of newly analysed ChIP-seq samples can be compared with the reference values of the compendium and enrichment profiles can be 
#' plotted against pre-computed profiles of published datasets.  
#' The metagene profiles show the problematic samples signal (red line) for the ChIP, for the input and their relative enrichment when compared to the 
#' compendium’s mean signal (black line) and its 2 x standard error (blue shadow). Additionally the function plots the desired QC-metric
#' as a red dashed line for the sample plotted against the reference distribution (density plots) of the compendium values stratified by chromatin marks.
#'
#' f_metagenePlotsForComparison
#' @param chrommark String, chromatin mark to be analysed. Has to be one of the following: H3K36me3, H3K4me3, H3K79me2, H4K20me1,H2AFZ","H3K27me3","H3K9me3","H3K27ac","POLR2AphosphoS5","H3K9ac","H3K4me2",
#' H3K9me1, H3K4me1, H3K79me1,H3K4ac,H3K14ac,H2BK5ac,H2BK120ac,H2BK15ac,H4K91ac,H4K8ac,H3K18ac,H2BK12ac,H3K56ac
#' H3K23ac, H2AK5ac, H2BK20ac,H4K5ac,H4K12ac,H2A.Z,H3K23me2,H2AK9ac,H3T11ph. For RNAPOL2 different variants are available: POLR2A (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param twopointElements  DESCRIBE!!!
#' @param TSSElements DESCRIBE!!!
#' @param TESElements DESCRIBE!!!
#' @param profilePath Path and dirctory in which the precalculated stratified reference profiles are stored (default=)
#' @param savePlotPath Path in which the profiles (in pdf) should be saved (Default=working directory)
#'
#'
#' @examples
#'\{dontrun
#' f_metagenePlotsForComparison(chrommark="H3K4me1",Meta_Result$twopoint, Meta_Result$TSS, Meta_Result$TES,plotName=chipName,profilePath="/lustre/data/FF/Carmen/BitBucket/chic/data/Profiles")

#'}
##GLOBAL VARIABLES
dataDirectory="../data"
profilePath="../data/profile"
Hlist=c("H3K36me3","POLR2A","H3K4me3","POLR3G","H3K79me2","H4K20me1","H2AFZ","H3K27me3","H3K9me3","H3K27ac","POLR2AphosphoS5","H3K9ac","H3K4me2",
	"H3K9me1","H3K4me1","POLR2AphosphoS2","H3K79me1","H3K4ac","H3K14ac","H2BK5ac","H2BK120ac","H2BK15ac","H4K91ac","H4K8ac","H3K18ac","H2BK12ac","H3K56ac",
	"H3K23ac","H2AK5ac","H2BK20ac","H4K5ac","H4K12ac","H2A.Z","H3K23me2","H2AK9ac","H3T11ph")

f_metagenePlotsForComparison=function(chrommark,twopointElements, TSSElements, TESElements, savePlotPath=getwd())
{
	psc <- 1; # pseudocount # required to avoid log2 of 0
	if (chrommark %in% Hlist)
	{
		for (endung in c("TWO","TES","TSS"))
		{
			print(endung)
			if (endung=="TWO")
			{
				iframe=log2(do.call(rbind,twopointElements$input)+psc)
				cframe=log2(do.call(rbind,twopointElements$chip)+psc)				
			}else if(endung=="TES"){				
				iframe=log2(do.call(rbind,TESElements$input)+psc)
				cframe=log2(do.call(rbind,TESElements$chip)+psc)	
			}else{
				iframe=log2(do.call(rbind,TSSElements$input)+psc)
				cframe=log2(do.call(rbind,TSSElements$chip)+psc)	
			}

			frame=NULL
			#load average dataframe
			name=paste(chrommark,"_",endung,"norm.RData", sep="")
			a= tryCatch({load(file.path(profilePath,name))},warning=function(z){print("error")})
			normMin=min(frame$mean-frame$sderr)
			normMax= max(frame$mean+frame$sderr)
			n_mean=frame

			frame=NULL
			#load average dataframe
			name=paste(chrommark,"_",endung,"chip.RData", sep="")
			a= tryCatch({load(file.path(profilePath,name))},warning=function(z){print("error")})
			absoluteMin=min(frame$mean-frame$sderr)
			absoluteMax= max(frame$mean+frame$sderr)
			c_mean=frame

			frame=NULL
			#load average dataframe
			name= paste(chrommark,"_",endung,"input.RData", sep="")
			a= tryCatch({load(file.path(profilePath,name))},warning=function(z){print("error")})
			if ((min(frame$mean-frame$sderr))<absoluteMin){absoluteMin=min(frame$mean-frame$sderr)}
			if ((max(frame$mean+frame$sderr))>absoluteMax){absoluteMax= max(frame$mean+frame$sderr)}
			i_mean=frame

			
			##prepare individual data
			#nframe<-colMeans(t(t(iframe)-t(cframe)),na.rm=T)
			nframe<-colMeans(t(t(cframe)-t(iframe)),na.rm=T)
			iframe=colMeans(iframe,na.rm=T)
			cframe= colMeans(cframe,na.rm=T)
			
			finalframe=NULL
			finalframe=cbind(c_mean$x,iframe)
			rownames(finalframe)=NULL
			finalframe=as.data.frame(finalframe)	
			colnames(finalframe)=c("x","mean") ##q for: 0% 25% 50% 75% 100%
			iframe=finalframe

			finalframe=NULL
			finalframe=cbind(c_mean$x,cframe)
			rownames(finalframe)=NULL
			finalframe=as.data.frame(finalframe)	
			colnames(finalframe)=c("x","mean") ##q for: 0% 25% 50% 75% 100%
			cframe=finalframe

			newMin=min(cframe$mean,absoluteMin,iframe$mean)
			newMax= max(cframe$mean,absoluteMax,iframe$mean)


			finalframe=NULL
			finalframe=cbind(c_mean$x,nframe)
			rownames(finalframe)=NULL
			finalframe=as.data.frame(finalframe)	
			colnames(finalframe)=c("x","mean") ##q for: 0% 25% 50% 75% 100%
			nframe=finalframe
			
			normMin=min(nframe$mean,normMin)
			normMax=max(nframe$mean,normMax)

			f_plotProfiles(i_mean, iframe, endung, c(newMin-0.001,newMax+0.001),maintitel=paste(chrommark,endung,"Input",sep="_"),savePlotPath)
			f_plotProfiles(c_mean, cframe, endung, c(newMin-0.001,newMax+0.001),maintitel=paste(chrommark,endung,"Chip",sep="_"),savePlotPath)
			f_plotProfiles(n_mean, nframe, endung, c(normMin-0.001,normMax+0.001),maintitel=paste(chrommark,endung,"norm",sep="_"),ylab="mean log2 enrichment (signal/input)")
			print(paste("pdf saved under",savePlotPath,sep=" "))
		}
	}else{
		print("Chromatin marks has to be one of the following:")
		print(Hlist)
	}
}


#'@title Function to create reference distribution plot for comparison
#'
#' @description
#' Creates a density plot (in pdf) for the sample against the reference distribution (density plots) of the compendium values stratified by chromatin marks.
#' 
#'
#' f_metagenePlotsForComparison
#' @param chrommark
#' @param 
#' @param valueToBePlotted Default="RSC"
#' @param currentValue
#' @param savePlotPath Default=NULL, when set saves the density plot as pdf under the gieven path
#'
#' @examples
#'\{dontrun
#' WHATEVER
#'}


f_plotReferenceDistribution=function(chrommark,valueToBePlotted="RSC",currentValue,savePlotPath=NULL)
{
	if (chrommark %in% Hlist)
	{
	
		filename=file.path(dataDirectory,"numbersHistone_AllGenes_ENCODE.txt")
		d1=read.table(filename,stringsAsFactors=TRUE,header=TRUE)
		filename=file.path(dataDirectory,"numbersHistone_AllGenes_unconsolidated.txt")
		d2=read.table(filename,stringsAsFactors=TRUE,header=TRUE)
		values=rbind(d1,d2)
		#get values for chrommark
		fake=paste("CC",valueToBePlotted,sep="_")
		subset= values[which(values$CC_TF==chrommark),fake]
		
		f_plotValueDistribution(subset,valueToBePlotted,currentValue,savePlotPath)
	}
}
