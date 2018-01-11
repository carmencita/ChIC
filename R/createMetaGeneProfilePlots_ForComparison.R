
# f_plot(i_mean, iframe, endung, c(newMin-0.001,newMax+0.001),maintitel=paste(ll,endung,"Input",sep="_"))


# meanFrame=i_mean
# currentFrame =iframe
# endung="TWO"
# absoluteMinMax=c(newMin-0.001,newMax+0.001)
# maintitel=paste(ll,endung,"Input",sep="_")
# ylab="mean of log2 tag density"

f_plot <- function(meanFrame, currentFrame , endung="TWO", absoluteMinMax, maintitel="title",ylab="mean of log2 tag density",savePlotPath=getwd())#, color="green") 
{
	filename=paste(maintitel,".pdf", sep="")
	pdf(file=file.path(savePlotPath,filename),width=10, height=7)
	#absoluteMin=min(frame$mean,r_o_frame$mean)
	#absoluteMax= max(frame$mean,r_o_frame$mean)
	
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

f_metagenePlotsForComparison=function(chrommark,twopointElements, TSSElements, TESElements,plotName="NA",profilePath=getwd())
{
	psc <- 1; # pseudocount # required to avoid log2 of 0

	Hlist=c("H3K36me3","POLR2A","H3K4me3","POLR3G","H3K79me2","H4K20me1","H2AFZ","H3K27me3","H3K9me3","H3K27ac","POLR2AphosphoS5","H3K9ac","H3K4me2",
"H3K9me1","H3K4me1","POLR2AphosphoS2","H3K79me1","H3K4ac","H3K14ac","H2BK5ac","H2BK120ac","H2BK15ac","H4K91ac","H4K8ac","H3K18ac","H2BK12ac","H3K56ac",
"H3K23ac","H2AK5ac","H2BK20ac","H4K5ac","H4K12ac","H2A.Z","H3K23me2","H2AK9ac","H3T11ph")

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

		f_plot(i_mean, iframe, endung, c(newMin-0.001,newMax+0.001),maintitel=paste(plotName,chrommark,endung,"Input",sep="_"),savePlotPath=getwd())
		f_plot(c_mean, cframe, endung, c(newMin-0.001,newMax+0.001),maintitel=paste(plotName,chrommark,endung,"Chip",sep="_"),savePlotPath=getwd())
		f_plot(n_mean, nframe, endung, c(normMin-0.001,normMax+0.001),maintitel=paste(plotName,chrommark,endung,"norm",sep="_"),ylab="mean log2 enrichment (signal/input)",savePlotPath=getwd())
	}
}