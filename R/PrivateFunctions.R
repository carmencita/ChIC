#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for narrow binding PROFILES 	###########################
#######														###########################
#######################################################################################

load("Settings.RData")
##format conversion
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

##format conversion
f_converNarrowPeakFormat =function(bd, margin = bd$whs) 
{
	if (is.null(margin)) 
	{
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

##get strand shift
f_getCustomStrandShift= function(x,y){
	#x=binding.characteristics$cross.correlation$x
	#y=binding.characteristics_cross.correlation_y_smoothed
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

##read bam file of tagalign file 
f_readFile=function(f_filename,f_path=getwd())
{
	reads.aligner.type="bam"
	read_tags_current_function<-get(paste("read", reads.aligner.type , "tags", sep="."))
	#if (f_reads.aligner.type=="bam")
	#{
	    data<-read_tags_current_function(file.path(f_path,paste(f_filename,".bam",sep="")))
	#}
	#if (f_reads.aligner.type=="tagalign")
	#{
	#    data<-read_tags_current_function(file.path(f_path,paste(f_filename,".tagAlign",sep="")))
	#}
	readCount=sum(sapply(data$tags, length))
	print(readCount)


	return(data)
}

##caluclate QC tag
f_qflag=function(RSC)
{
	qflag=NA
	if ( (RSC >= 0) & (RSC < 0.25) ){
		qflag <- -2
	} else if ( (RSC >= 0.25) & (RSC < 0.5) ) {
		qflag <- -1
	} else if ( (RSC >= 0.5) & (RSC < 1) ) {
		qflag <- 0
	} else if ( (RSC >= 1) & (RSC < 1.5) ) {
		qflag <- 1
	} else if ( (RSC >= 1.5) ) {
		qflag <- 2
	}
	return(qflag)
}

# #My version for binding detection
# get.binding.characteristicsMy = function(data,srange=c(50,500),bin=5,cluster=NULL,debug=F,min.tag.count=1e3,acceptance.z.score=3,remove.tag.anomalies=T,anomalies.z=5,accept.all.tags=F) 
# {
#   if(remove.tag.anomalies) {
#     data <- remove.tag.anomalies(data,z=anomalies.z);
#   }
  
#   # take highest quality tag bin
#   if(!is.null(data$quality) && !accept.all.tags) {
#     min.bin <- min(unlist(lapply(data$quality,min)))
#     chrl <- names(data$tags); names(chrl) <- chrl;
#     otl <- lapply(chrl,function(chr) data$tags[[chr]][data$quality[[chr]]==min.bin]);
#   } else {
#     otl <- data$tags;
#   }
#   # remove empty chromosomes
#   otl <- otl[unlist(lapply(otl,length))!=0];


#   # calculate strand scc
#   if(!is.null(cluster)) {
#     cc <- clusterApplyLB(cluster,otl,tag.scc,srange=srange,bin=bin);
#     names(cc) <- names(otl); 
#   } else {
#     cc <- lapply(otl,tag.scc,srange=srange,bin=bin);
#   }
#   ccl<-list(sample=cc);
#   ccl.av <- lapply(names(ccl),t.plotavcc,type='l',ccl=ccl,return.ac=T,ttl=list(sample=otl),plot=F)[[1]]
#   ccl.av <- data.frame(x=as.numeric(names(ccl.av)),y=as.numeric(ccl.av));
  
#   # find peak
#   pi <- which.max(ccl.av$y);
  
#   # determine width at third-height
#   th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/3+ccl.av$y[length(ccl.av$y)]
#   whs <- max(ccl.av$x[ccl.av$y>=th]);

# #  if (! is.integer(whs)) { # Anshul: added this to avoid situations where whs ends up being -Inf
#   if (!is.finite(whs)) { # fixed to avoid a TRUE with numeric values
#     whs <- ccl.av$x[ min(c(2*pi,length(ccl.av$y))) ]
#   }

#   # determine acceptance of different quality bins
#   # calculates tag scc for the best tags, and combinations of best tag category with every other category
#   # for subsequent selection of acceptable categories
#   scc.acceptance.calc <- function() {

#     qr <- range(unlist(lapply(data$quality,range)))

#     # start with best tags

#     # determine half-width for scc calculations
#     pi <- which.max(ccl.av$y);

#     # determine width at half-height
#     th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/2+ccl.av$y[length(ccl.av$y)]
#     lwhs <- max(ccl.av$x[ccl.av$y>=th])-ccl.av$x[pi];
#     lwhs <- max(c(20,bin*10,lwhs));
#     srange <- ccl.av$x[pi]+c(-lwhs,lwhs)

#     # calculate chromosome-average scc
#     t.scc <- function(tags) {
#       if(is.null(cluster)) {
#         cc <- lapply(tags,tag.scc,srange=srange,bin=bin);
#       } else {
#         cc <- clusterApplyLB(cluster,tags,tag.scc,srange=srange,bin=bin); names(cc) <- names(tags);
#       }
#       return(t.plotavcc(1,type='l',ccl=list(cc),ttl=list(tags),plot=F,return.ac=T))
#     }


#     # returns info list for a given tag length (lv), mismatch count (nv)
#     t.cat <- function(qual) {
#       # construct tag set
#       if(qual==qr[1]) {
#         ts <- otl;
#       } else {
#         nts <- names(otl); names(nts) <- nts;
#         # select tags
#         at <- lapply(nts,function(chr) data$tags[[chr]][data$quality[[chr]]==qual]);
#         ntags <- sum(unlist(lapply(at,length)));
#         if(ntags<min.tag.count) { return(NULL); }

#         # append to otl
#         ts <- lapply(nts,function(nam) c(otl[[nam]],at[[nam]]));
#       }
#       return(t.scc(ts));
#     }

#     # calculate cross-correlation values for each quality bin
#     ql <- sort(unique(unlist(lapply(data$quality,unique)))); names(ql) <- ql;

#     qccl <- lapply(ql,t.cat);

#     # acceptance tests
#     ac <- c(T,unlist(lapply(qccl[-1],function(d) if(is.null(d)) { return(F) } else { t.test(d-qccl[[as.character(min.bin)]],alternative="greater")$p.value<pnorm(acceptance.z.score,lower.tail=F) }))); names(ac) <- names(qccl);
#     return(list(informative.bins=ac,quality.cc=qccl))
#   }

#   if(accept.all.tags | is.null(data$quality)) {
#     return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs))    
#   } else {
#     acc <- scc.acceptance.calc();
#     return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs,quality.bin.acceptance=acc));
#   }

# }

##selet informative tags
f_selectInformativeTag=function(chip,input,chip_b.characteristics,input_b.characteristics,select.informative.tags=FALSE,remove.local.tag.anomalies=TRUE)
{
	print("Filter tags")

	if (select.informative.tags) {
	      print("select.informative.tags filter")
	     #load(paste("sppdata", "binding", chip.data.samplename, "RData", sep="."))
	      chip.dataSelected <- select.informative.tags(chip, chip_b.characteristics)
	      #load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
	      input.dataSelected <- select.informative.tags(input.data, input_b.characteristics)
	} else {
	      print("SKIP select.informative.tags filter")
	      chip.dataSelected<-chip$tags
	      input.dataSelected<-input$tags
	}
	if (remove.local.tag.anomalies) {
		print("remove.local.tag.anomalies filter")
	      # restrict or remove singular positions with very high tag counts
		chip.dataSelected <- remove.local.tag.anomalies(chip.dataSelected)
		input.dataSelected <- remove.local.tag.anomalies(input.dataSelected)
	} else {
		print("SKIP remove.local.tag.anomalies filter")
		input.dataSelected<-input$tags
		chip.dataSelected=chip$tags
	}

	finalList=list("input.dataSelected"=input.dataSelected,
		"chip.dataSelected"=chip.dataSelected)
	return(finalList)
}

##takes dataSelected as input, parallel is the number of CPUs used for parallelization
f_tagDensity=function(data,tag.shift,rngl,mc=1)
{
	##Smoothing parameters
	smoothingStep<-20	##is size of sw ##before it was 10
	smoothingBandwidth<-50
	## density distribution for data
	print("Smooth tag density")
	ts <- sum(unlist(lapply(data,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
	##parallelisation
	chromosomes_list<-names(data)
	##creates a list of lists
	data<-lapply(chromosomes_list, FUN=function(x) {
		return(data[x])
	})
	if (mc>1)
	{
		smoothed.density<-mclapply(data, FUN=function(current_chr_list)
		{
		    current_chr<-names(current_chr_list)
		    str(current_chr_list)
		    if (length(current_chr) != 1) 
		    {
		        stop("unexpected dataSelected structure")
		    }
		    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
		}, mc.preschedule = FALSE,mc.cores=mc)
	}else{
		smoothed.density<-lapply(data, FUN=function(current_chr_list)
		{
		    current_chr<-names(current_chr_list)
		    str(current_chr_list)
		    if (length(current_chr) != 1) 
		    {
		        stop("unexpected dataSelected structure")
		    }
		    get.smoothed.tag.density(current_chr_list, bandwidth=smoothingBandwidth, step=smoothingStep,tag.shift=tag.shift, rngl=rngl[current_chr])
		})
	}
	smoothed.density=(unlist(smoothed.density,recursive=FALSE))
	#normalizing smoothed tag density by library size
	smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })
	return(smoothed.density)
}



# # -------- ROUTINES FOR WRITING OUT TAG DENSITY AND ENRICHMENT PROFILES  ------------
# # calculate smoothed tag density, optionally subtracting the background
# get.smoothed.tag.densityMy <- function(signal.tags,control.tags=NULL,bandwidth=150,bg.weight=NULL,tag.shift=146/2,step=round(bandwidth/3),background.density.scaling=T,rngl=NULL,scale.by.dataset.size=F) {
#   chrl <- names(signal.tags);
#   #if(!is.null(rngl)) { chrl <- names(rngl); }
#   names(chrl) <- chrl;
  

#   if(!is.null(control.tags) && is.null(bg.weight)) {
#     bg.weight <- dataset.density.ratio(signal.tags,control.tags,background.density.scaling=background.density.scaling);
#   }

#   if(scale.by.dataset.size) {
#     den.scaling <- 1/(dataset.density.size(signal.tags,background.density.scaling=background.density.scaling)/1e6);
#   } else {
#     den.scaling <- 1;
#   }
  
#   lapply(chrl,function(chr) {
#     ad <- abs(signal.tags[[chr]]+tag.shift);
#     rng <- NULL;
#     if(!is.null(rngl)) {
#       rng <- rngl[[chr]];
#     }
#     if(is.null(rng)) {
#       rng <- range(ad);
#     }

#     ds <- densum(ad,bw=bandwidth,from=rng[1],to=rng[2],return.x=T,step=step);
#     if(!is.null(control.tags)) {
#       if(!is.null(control.tags[[chr]])) {
#         bsd <- densum(abs(control.tags[[chr]]+tag.shift),bw=bandwidth,from=rng[1],to=rng[2],return.x=F,step=step);
#         ds$y <- ds$y-bsd*bg.weight;
#       }
#     }
#     return(data.frame(x=seq(ds$x[1],ds$x[2],by=step),y=den.scaling*ds$y))
#   })
# }


#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for global read distribution 	###########################
#######														###########################
#######################################################################################


f_shortenFrame=function(smoothedDensity)
{
	##shorten frame with cumulative distribution
	newSmoothedDensity=NULL
	chrl=names(smoothedDensity)
	for (i in chrl)
	{
		#print(i)
		##appending tags and quality elements of all inputs to newControl
		newSmoothedDensity[[i]]$x=smoothedDensity[[i]]$x[seq.int(1, length(smoothedDensity[[i]]$x), 6)]
		newSmoothedDensity[[i]]$y=smoothedDensity[[i]]$y[seq.int(1, length(smoothedDensity[[i]]$y), 6)]
	}
	return(newSmoothedDensity)
}

##sorts and bins the dataframe
f_sortAndBinning=function(shortframe)
{
	shortframeSorted=sort(shortframe)
	BINS_NUMBER<-1e4

	cumSum<-cumsum(shortframeSorted)
	cumSumBins<-quantile(cumSum, probs=seq(0,1,(1/BINS_NUMBER)))
	pj<-(cumSumBins/cumSum[length(cumSum)])
	normalizedCumSum=data.frame(x=seq(0,1,(1/BINS_NUMBER)),pj=pj)
	rownames(normalizedCumSum)<-NULL
	return(normalizedCumSum)
}

#The plot function creates the "Fingerprint plot" i.e. the cumulative distribution of the read counts across genomic bins for Input and ChIP.
#  cumChip The cumulative distribution of the read counts for the ChIP
#  cumInput The cumulative distribution of the read counts for the Input
# savePlotPath 


f_fingerPrintPlot=function(cumChip,cumInput,savePlotPath=NULL)
{	
	
	if (!is.null(savePlotPath))
	{
		filename=file.path(savePlotPath,"FingerPrintPlot.pdf")
		pdf(file=filename,width=7, height=7)
	}
	plot(cumChip,type="l",col="blue",lwd=2,xlab="Percentage of bins",ylab="Percentage of tags",main="Fingerprint: global read distribution")
	lines(cumInput,col="red",lwd=2)
	arrowx=cumChip[which.max(abs(cumInput$pj-cumChip$pj)),]$x
	abline(v =arrowx, col = "green",lty=2,lwd=2)
	#abline(h=schneidePunktY,col='cyan',lty=2,lwd=2)
	#abline(v=schneidePunktX,col='cyan',lty=2,lwd=2)
	legend("topleft", legend = c("Input","ChIP"), fill = c("red","blue"))
	if (!is.null(savePlotPath))
	{
		dev.off()
    	print(paste("pdf saved under",filename,sep=" "))

	}
}



#######################################################################################
#######													 	###########################
####### FUNCTIONS QC-metrics for local ENRICHMENT 			###########################
#######														###########################
#######################################################################################



###############################################################
#### BEGINF OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
###############################################################
#### - f_sfeature.bin.averages
#### - two.point.scaling
#### - one.point.scaling
###############################################################



# calculates bin average values for a list of one- or two-point features
# dat - $x/$y data frame
# feat - data frame with either $x/{$strand} or $s/$e/{$strand}
# return.scaling causes function to just calculate, cleanup and return scaling mapping for further use
f_feature.bin.averages <- function(dat,feat,nu.feat.omit=F,nu.point.omit=T, scaling=NULL, return.scaling=F, trim=0, min.feature.size=NULL, ... ) {
  if(is.null(feat)) { return(NULL) };
  if(dim(feat)[1]<1) { return(NULL) };

  if(is.null(scaling)) { # calculate scaling
    if(!is.null(feat$e)) { # two-point
      scaling <- f_two.point.scaling(dat$x,feat, ... )
      #scaling <- two.point.scaling(dat$x,feat,om=2e3,im=1e3,bs=2e3,nbins=100)
    } else { # one-point
      scaling <- f_one.point.scaling(dat$x,feat$x,feat$strand, ... );
      #scaling <- one.point.scaling(dat$x,feat$x,feat$strand, m=2e3,nbins=50 );
    }
    # clean up
    if(nu.feat.omit) {
      scaling <- scaling[!scaling$si %in% unique(scaling$si[scaling$nu>1]),]
    } else {
      if(nu.point.omit) {
        scaling <- scaling[scaling$nu==1,]
      }
    }
    if(return.scaling) {
      return(scaling);
    }
  }
  
  if(!is.null(min.feature.size)) {
    if(!is.null(feat$e)) {
      scaling <- scaling[scaling$si %in% which(feat$e-feat$s>=min.feature.size),]
    }
  }
  
  # determine and return gene bin average table
  return(tapply(dat$y[scaling$i],list(factor(scaling$si,levels=c(1:dim(feat)[1])),scaling$bin),function(x) mean(na.omit(x))))
  #return(tapply(dat$y[scaling$i],list(scaling$si,scaling$bin),mean,trim=trim,na.rm=T))
}



# identifies points located within margins of given segments
# and returns rescaled relative positions
# m - margin; om - outer margin; im - inner margin; rom - right outer margin, etc.
# x - x positions being mapped/rescaled
# seg - $s/$e/{$strand} data frame
# return $i/$rx/$nu data frame corresponding to indecies of original x values ($i) and new relative $rx positions, and whether a point maps to several segments
# $r5x/$r3x give distances relative to the 5' and 3' ends, $si gives segment index
# $bin index factor is added if nbins is supplied
#two.point.scaling <- function(x,seg,m=1e3,bs=gene_body,om=m,im=m,rom=om,lom=om,rim=im,lim=im,nbins=predefnum_bins/2) {
f_two.point.scaling <- function(x,seg,bs=gene_body,om=m,im=m,rom=om,lom=om,rim=im,lim=im,nbins=predefnum_bins) {
  # map points to the segments defined by outer margins
  nseg <- length(seg$s);
  if(!is.null(seg$strand)) {
    nsi <- which(seg$strand=="-");
    ml <- rep(lom,nseg); ml[nsi] <- rom;
    mr <- rep(rom,nseg); mr[nsi] <- lom;
    sg5 <- seg$s; sg5[nsi] <- seg$e[nsi];
    sg3 <- seg$e; sg3[nsi] <- seg$s[nsi];
    sstrand <- ifelse(seg$strand=="-",-1,1); sstrand[is.na(sstrand)] <- 1;
  } else {
    ml <- rep(lom,nseg);
    mr <- rep(rom,nseg);
    sg5 <- seg$s;
    sg3 <- seg$e;
    sstrand <- rep(1,length(sg5));
  }
  spi <- points.within(x,seg$s-ml,seg$e+mr,return.list=T);
  lspi <- unlist(lapply(spi,length))
  df <- data.frame(i=rep(1:length(x),lspi),si=unlist(spi),nu=rep(lspi,lspi))
  df$r5x <- (x[df$i]-sg5[df$si])*sstrand[df$si];
  df$r3x <- (x[df$i]-sg3[df$si])*sstrand[df$si];

  # between g3-rim and g3+rom - unscaled
  df$rx <- rep(NA,length(df$i));
  vi <- which(df$r3x >= -rim);
  df$rx[vi] <- df$r3x[vi]+lim+rim+bs;
  
  #df$rx <- df$r3x+lim+rim+bs;
  
  # body scaling
  vi <- which(df$r3x <  -rim & df$r5x > lim);
  df$rx[vi] <- ((df$r5x[vi]-lim)/((seg$e-seg$s)[df$si[vi]]-lim-rim))*bs+lim;
  # between g5-lom and g5+lim
  vi <- which(df$r5x <= lim);
  df$rx[vi] <- df$r5x[vi];

  if(!is.null(nbins)) {
    breaks <- seq(-lom, lim+bs+rim+rom, by=(lom+lim+bs+rim+rom)/nbins);
    rn <- breaks[-1]-diff(breaks)/2;
    breaks[1] <- breaks[1]-1; breaks[length(breaks)] <- breaks[length(breaks)]+1;
    df$bin <- cut(df$rx, breaks, labels=rn)
  }

  return(df);
}



# much simpler, one-point scaling that returns relative positions in the format similar to the two-point scaling ($i/$rx/$nu/$si)
#one.point.scaling <- function(x, pos, strand=NULL, m=1e3, lm=m, rm=m, nbins=predefnum_bins/2) {  m=2010
f_one.point.scaling <- function(x, pos, strand=NULL,m=up_downStream/2, lm=m, rm=m, nbins=predefnum_bins/2) {
 # print(nbins)
 # print (up_downStream)
#  print(m)
  if(is.null(pos)) { return(NULL) }
  nseg <- length(pos);
  if(nseg<1) { return(NULL) }
  ml <- rep(lm,nseg); 
  mr <- rep(rm,nseg); 
  if(!is.null(strand)) {
    nsi <- which(strand=="-");
    ml[nsi] <- rm;
    mr[nsi] <- lm;
    sstrand <- ifelse(strand=="-",-1,1); sstrand[is.na(sstrand)] <- 1;
  } else {
    sstrand <- rep(1,length(pos));
  }
  spi <- points.within(x,pos-ml,pos+mr,return.list=T);
  lspi <- unlist(lapply(spi,length))
  df <- data.frame(i=rep(1:length(x),lspi),si=unlist(spi),nu=rep(lspi,lspi))
  df$rx <- (x[df$i]-pos[df$si])*sstrand[df$si];
  if(!is.null(nbins)) {
    if(is.null(strand)) {
      breaks <- seq(-lm,rm,by=(lm+rm)/nbins);
    } else {
      breaks <- seq(-max(lm,rm),max(lm,rm),by=(lm+rm)/nbins);
    }
    rn <- breaks[-1]-diff(breaks)/2;
    breaks[1] <- breaks[1]-1; breaks[length(breaks)] <- breaks[length(breaks)]+1;
    df$bin <- cut(df$rx,breaks,labels=rn)
  }
  return(df);
}

###############################################################
#### END OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
###############################################################

##MODIFIED VERSION OF t.get.gene.av.density for TWO.POINT
# a function for getting bin-averaged profiles for individual genes TWO.POINT
#t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=500,lom=5e3, rom=1e3, bs=2e3, nbins=400,separate.strands=F) {
	#min.feature.size=min.feature.sizeMy=3000
f_t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=inner_margin,lom=left_outer_margin, rom=right_outer_margin, bs=gene_body, nbins=predefnum_bins,separate.strands=F, min.feature.size=3000,mc=1) {
  chrl <- names(gdl);
  names(chrl) <- chrl;
  #lapply(chrl[chrl %in% names(chipTags_current$td)],function(chr) {
  mclapply(chrl[chrl %in% names(chipTags_current$td)],  mc.preschedule = FALSE,mc.cores=mc, FUN=function(chr) {
    print(chr)
    nsi <- gdl[[chr]]$strand=="-";
    print(length(nsi))
    current_gene_names<-gdl[[chr]]$geneSymbol
    if ((sum(!nsi)>0)) {
    #if ((sum(!nsi)>1)) {
      px <- f_feature.bin.averages(chipTags_current$td[[chr]],data.frame(s=gdl[[chr]]$txStart[!nsi],e=gdl[[chr]]$txEnd[!nsi]),lom=lom,rom=rom,im=im,bs=bs, nbins=nbins, min.feature.size=min.feature.size, nu.point.omit=FALSE)
      rownames(px)<-current_gene_names[!nsi]
      } else { 
        px<-NULL
      }
    if ((sum(nsi)>0)) {
    #if ((sum(nsi)>1)) {
        nd <- chipTags_current$td[[chr]]; nd$x <- -1*nd$x;
        nx <- f_feature.bin.averages(nd,data.frame(s=-1*gdl[[chr]]$txEnd[nsi],e=-1*gdl[[chr]]$txStart[nsi]), lom=lom,rom=rom,im=im,bs=bs, nbins=nbins,min.feature.size=min.feature.size, nu.point.omit=FALSE)
        rownames(nx)<-current_gene_names[nsi]
           } else { 
        nx<-NULL
    }

    if(separate.strands) {
      return(p=px,n=nx);
          } else {
      return(rbind(px,nx));
          }
    })
  }


##MODIFIED VERSION OF t.get.gene.av.density FOR TSS
# a function for getting bin-averaged profiles for individual genes TSS ONE.POINT
f_t.get.gene.av.density_TSS <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F,mc=1) { ##binning of the frame
    chrl <- names(gdl);
    names(chrl) <- chrl;
    #lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
    mclapply(chrl[chrl %in% names(tl_current$td)],  mc.preschedule = FALSE,mc.cores=mc, FUN=function(chr) {
      nsi <- gdl[[chr]]$strand=="-";
      current_gene_names<-gdl[[chr]]$geneSymbol

       # if ((sum(!nsi)>0)) {
  if ((sum(!nsi)>0)) {
          px <- f_feature.bin.averages(tl_current$td[[chr]],data.frame(x=gdl[[chr]]$txStart[!nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
    #str(px)
          rownames(px)<-current_gene_names[!nsi]
        } else { 
          px<-NULL
        }
  
  if ((sum(nsi)>0)) {
       #if ((sum(nsi)>0)) {
          nd <- tl_current$td[[chr]]; nd$x <- -1*nd$x;
          nx <- f_feature.bin.averages(nd,data.frame(x=-1*gdl[[chr]]$txEnd[nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
          rownames(nx)<-current_gene_names[nsi]
       } else { 
          nx<-NULL
       }

      if(separate.strands) {
        return(p=px,n=nx);
      } else {
        return(rbind(px,nx));
      }
    })
  }

##MODIFIED VERSION OF t.get.gene.av.density FOR TES
# a function for getting bin-averaged profiles for individual genes TES ONE.POINT 
f_t.get.gene.av.density_TES <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F,mc=1) {
    chrl <- names(gdl);
    names(chrl) <- chrl;
    #lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
    mclapply(chrl[chrl %in% names(tl_current$td)],  mc.preschedule = FALSE,mc.cores=mc, FUN=function(chr) {
      nsi <- gdl[[chr]]$strand=="-";
      current_gene_names<-gdl[[chr]]$geneSymbol

        #if ((sum(!nsi)>0)) {
        if ((sum(!nsi)>0)) {
          px <- f_feature.bin.averages(tl_current$td[[chr]],data.frame(x=gdl[[chr]]$txEnd[!nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
          rownames(px)<-current_gene_names[!nsi]
        } else { 
          px<-NULL
        }

       #if ((sum(nsi)>1)) {
       if ((sum(nsi)>0)) {
          nd <- tl_current$td[[chr]]; nd$x <- -1*nd$x;
          nx <- f_feature.bin.averages(nd,data.frame(x=-1*gdl[[chr]]$txStart[nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
          rownames(nx)<-current_gene_names[nsi]
       } else { 
          nx<-NULL
       }

      if(separate.strands) {
        return(p=px,n=nx);
      } else {
        return(rbind(px,nx));
      }
    })
  }


##helper function for global settings
f_metaGeneDefinition=function(selection="break_points_2P")
{
	predefnum_bins=301   # 151

	##TWO point scaling
	inner_margin=500   # 500
	right_outer_margin=1010   # 1020
	left_outer_margin=2010    # 2020
	gene_body=2000   # 2000


	totalGeneLength=gene_body+2*inner_margin
	inner_margin2=totalGeneLength-inner_margin

	break_points_2P=c(-2000,0,inner_margin,gene_body+inner_margin,totalGeneLength,totalGeneLength+1000)

	total_output_gene_length<-left_outer_margin+inner_margin+gene_body+inner_margin+right_outer_margin
	estimated_bin_size_2P<-total_output_gene_length/predefnum_bins


	##ONE point scaling
	up_downStream=4020
	predefnum_bins_1P=201 ##binsize=20
	break_points=c(-2000,-1000,-500,0,500,1000,2000)
	### %%% need to use this for one point scaling estimated bin size....
	estimated_bin_size_1P<-up_downStream/predefnum_bins_1P

	if (selection=="break_points_2P")
	{return(break_points_2P)}
	if (selection=="break_points")
	{return(break_points)}

}


f_spotfunction=function(dframe,breaks, estimatedBinSize,tag="twopoints")
{
  hotSpots=NULL ##takes values at different predefined points
  for (i  in breaks){
    print(i)
    if (length(which(as.integer(row.names(dframe))==i))<1){
      x_bin= which(as.integer(row.names(dframe))==i+estimatedBinSize/2)}else{
      x_bin= which(as.integer(row.names(dframe))==i)
    }
    hotSpots=rbind(hotSpots,c(i, dframe[x_bin,][1], dframe[x_bin,][2]))
  }
  ##save hotSpots values  ChiP and input
  #write.table(cbind("chip",hotSpots[,1],hotSpots[,2],0,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  hotSpotsChip=cbind(paste("chip","hotSpots",tag,hotSpots[,1],sep="_"),round(hotSpots[,2],3)) # chip_hotSpots_twopoints_-2000 chip_hotSpots_twopoints_0 
  #write.table(cbind("input",hotSpots[,1],hotSpots[,3],0,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  hotSpotsInput=cbind(paste("input","hotSpots",tag,hotSpots[,1],sep="_"),round(hotSpots[,3],3)) #input_hotSpots_twopoints_-2000  input_hotSpots_twopoints_0 
  result=data.frame(rbind(hotSpotsChip,hotSpotsInput))
  colnames(result)=c("Feature","Value")
  return(result)
}


f_spotfunctionNorm=function(dframe,breaks, estimatedBinSize,tag="norm")
{

  hotSpotsNorm=NULL ##takes values at different predefined points
  for (i  in breaks){
    if (length(which(as.integer(names(dframe))==i))<1){
      x_bin= which(as.integer(names(dframe))==i+estimatedBinSize/2)}else{
      x_bin= which(as.integer(names(dframe))==i)
    }
    hotSpotsNorm=rbind(hotSpotsNorm,c(i, dframe[x_bin][1]))
  }
  ##save hotSpots values for normalized profile
  #write.table(cbind("norm",hotSpotsNorm[,1],hotSpotsNorm[,2],1,"hotSpots twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  hotSpotsN=data.frame(cbind(paste("norm","hotSpots",tag,hotSpotsNorm[,1],sep="_"),round(hotSpotsNorm[,2],3))) #norm_hotSpots_twopoints_-2000 
  colnames(hotSpotsN)=c("Feature","Value")
  return(hotSpotsN)
}


f_maximaAucfunction=function(dframe,breaks, estimatedBinSize,tag="twopoint")
{
  ##local maxima and area in all the predefined regions
  localMaxima_auc=NULL
  for (i in seq(length(breaks)-1))
  {
    if (length(which(as.integer(row.names(dframe))==breaks[i]))<1){
      xi= which(as.integer(row.names(dframe))==breaks[i]+estimatedBinSize/2)}else{
      xi= which(as.integer(row.names(dframe))==breaks[i])
    }
    if (length(which(as.integer(row.names(dframe))==breaks[i+1]))<1){
      xj= which(as.integer(row.names(dframe))==breaks[i+1]+estimatedBinSize/2)}else{
      xj= which(as.integer(row.names(dframe))==breaks[i+1])
    }

    sub_set=dframe[xi:xj,]
    l1=max(sub_set[,1])
    l2=max(sub_set[,2])
    print(breaks[i])
    localMaxima_auc=rbind(localMaxima_auc,c(rownames(sub_set)[which(sub_set[,1]==l1)], l1, sum((sub_set[,1])*estimatedBinSize),rownames(sub_set)[which(sub_set[,2]==l2)],l2,sum((sub_set[,2])*estimated_bin_size_2P)))
  }

  ##save x and y coordinate for local Maxima ChiP
  # write.table(cbind("chip",localMaxima_auc[,1],localMaxima_auc[,2],0,"localMax twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  xvaluesChip=cbind(paste("chip","localMax",tag,rownames(data.frame(localMaxima_auc)),"x",sep="_"),round(as.numeric(localMaxima_auc[,1]),3)) #chip_localMax_twopoint_1_x 
  yvaluesChip=cbind(paste("chip","localMax",tag,rownames(data.frame(localMaxima_auc)),"y",sep="_"),round(as.numeric(localMaxima_auc[,2]),3)) #chip_localMax_twopoint_1_y 
  ##save x and y coordinate for local Maxima Input
  # write.table(cbind("input",localMaxima_auc[,4],localMaxima_auc[,5],0,"localMax twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  xvaluesInput=cbind(paste("input","localMax",tag,rownames(data.frame(localMaxima_auc)),"x",sep="_"),round(as.numeric(localMaxima_auc[,4]),3)) #input_localMax_twopoint_1_x 
  yvaluesInput=cbind(paste("input","localMax",tag,rownames(data.frame(localMaxima_auc)),"y",sep="_"),round(as.numeric(localMaxima_auc[,5]),3)) #input_localMax_twopoint_1_y  

  ##save auc for Chip and Input
  #write.table(cbind("chip",localMaxima_auc[,3],0,"auc twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  aucChip=cbind(paste("chip","auc",tag,rownames(data.frame(localMaxima_auc)),sep="_"),round(as.numeric(localMaxima_auc[,3]),3)) #chip_auc_twopoint_1 
  #write.table(cbind("input",localMaxima_auc[,6],0,"auc twopoint"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  aucInput=cbind(paste("input","auc",tag,rownames(data.frame(localMaxima_auc)),sep="_"),round(as.numeric(localMaxima_auc[,6]),3)) #ipnut_auc_twopoint_1 

  resultsnoNorm=NULL
  resultsnoNorm= data.frame(rbind(xvaluesChip,yvaluesChip,xvaluesInput,yvaluesInput,aucChip,aucInput))
  colnames(resultsnoNorm)=c("Feature","Value")
  return(resultsnoNorm)
}

f_maximaAucfunctionNorm=function(dframe,breaks, estimatedBinSize,tag="norm")
{
  ##local maxima and area in all the predefined regions
  localMaxima_aucNorm=NULL
  for (i in seq(length(breaks)-1)){
    if (length(which(as.integer(names(dframe))==breaks[i]))<1){
      xi= which(as.integer(names(dframe))==breaks[i]+estimatedBinSize/2)}else{
      xi= which(as.integer(names(dframe))==breaks[i])
    }
    if (length(which(as.integer(names(dframe))==breaks[i+1]))<1){
      xj= which(as.integer(names(dframe))==breaks[i+1]+estimatedBinSize/2)}else{
      xj= which(as.integer(names(dframe))==breaks[i+1])
    }

    sub_set= data.frame(dframe[xi:xj])
    l1=max(sub_set)
    localMaxima_aucNorm=rbind(localMaxima_aucNorm,c(rownames(sub_set)[which(sub_set[,1]==l1)], l1, sum((sub_set)*estimatedBinSize)))
  }
  ##save x and y coordinate for local Maxima Normalized profile
  #write.table(cbind("norm",localMaxima_aucNorm[,1],localMaxima_aucNorm[,2],1,"localMax twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  xvaluesNorm=cbind(paste("norm","localMax",tag,rownames(data.frame(localMaxima_aucNorm)),"x",sep="_"),round(as.numeric(localMaxima_aucNorm[,1]),3)) #norm_localMax_twopoints_1_x
  yvaluesNorm=cbind(paste("norm","localMax",tag,rownames(data.frame(localMaxima_aucNorm)),"y",sep="_"),round(as.numeric(localMaxima_aucNorm[,2]),3)) #norm_localMax_twopoints_1_y
  ##save auc for normalized profile
  #write.table(cbind("norm",localMaxima_aucNorm[,3],1,"auc twopoints"),file=outname,row.names = FALSE,col.names=FALSE,append=TRUE, quote = FALSE)
  aucNorm=cbind(paste("norm","auc",tag,rownames(data.frame(localMaxima_aucNorm)),sep="_"),round(as.numeric(localMaxima_aucNorm[,3]),3)) #norm_auc_twopoints_1 
  
  resultsNorm=NULL
  resultsNorm= data.frame(rbind(xvaluesNorm,yvaluesNorm,aucNorm))
  colnames(resultsNorm)=c("Feature","Value")
  return(resultsNorm)
}





#######################################################################################
#######													 	###########################
####### FUNCTIONS comparison with compendium	 			###########################
#######														###########################
#######################################################################################

#helper functions to incorporate profile chromatins
f_Hlist=function(giveBack=TRUE)
{

	##GLOBAL VARIABLES
	Hlist<-c("H3K36me3","POLR2A","H3K4me3","H3K79me2","H4K20me1","H2AFZ","H3K27me3","H3K9me3","H3K27ac","POLR2AphosphoS5","H3K9ac","H3K4me2",
	"H3K9me1","H3K4me1","POLR2AphosphoS2","H3K79me1","H3K4ac","H3K14ac","H2BK5ac","H2BK120ac","H2BK15ac","H4K91ac","H4K8ac","H3K18ac","H2BK12ac","H3K56ac",
	"H3K23ac","H2AK5ac","H2BK20ac","H4K5ac","H4K12ac","H2A.Z","H3K23me2","H2AK9ac","H3T11ph")
	#"POLR3G"
	return(Hlist)
}

#helper functions to define chromating
f_chromatinMarkClasses=function(giveBack=TRUE)
{
	allSharp<-c("H3K27ac","H3K9ac","H3K14ac","H2BK5ac","H4K91ac","H3K18ac","H3K23ac","H2AK9ac","H3K4me3","H3K4me2","H3K79me1","H2AFZ","H2A.Z","H4K12ac"
		,"H4K8ac","H3K4ac","H2BK12ac","H4K5ac","H2BK20ac","H2BK120ac","H2AK5ac","H2BK15ac")
	allBroad<-c("H3K23me2","H3K9me2","H3K9me3","H3K27me3","H4K20me1","H3K36me3","H3K56ac","H3K9me1","H3K79me2","H3K4me1","H3T11ph")
	##USE ONLY POLR2A for Pol2 class
	RNAPol2<-"POLR2A"
	final=list("allSharp"=allSharp,
		"allBroad"=allBroad,
		"RNAPol2"=RNAPol2)
	return(final)
}

##plot profiles compendium versus current dataset
f_plotProfiles <- function(meanFrame, currentFrame , endung="TWO", absoluteMinMax, maintitel="title",ylab="mean of log2 tag density",savePlotPath=NULL)#, color="green") 
{
	if (!is.null(savePlotPath))
    {
		filename=paste(maintitel,".pdf", sep="")
		pdf(file=file.path(savePlotPath,filename),width=10, height=7)
	}
	break_points_2P=f_metaGeneDefinition("break_points_2P")
	break_points=f_metaGeneDefinition("break_points")
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
		currBreak_points=break_points_2P[c(-2,-5)] #c(-2000,500,2500,4000)
		abline(v=c( break_points_2P[c(2,5)]),lty=2,col="darkgrey", lwd=3)
	   	abline(v=currBreak_points,lty=3,col="darkgrey", lwd=2)
	   	axis(side = 1, at = break_points_2P, labels = c("-2KB","TSS","TSS+500","TES-500","TES","+1KB"))
	}else{
	   	TSSbreak_points=break_points[-c(4)]# c(-2000,-1000,-500,500,1000,2000)
	   	if (endung=="TSS")
		{
			axis(side = 1, at = break_points, labels = c("-2KB","-1KB","-500","TSS","+500","+1KB","+2KB"))#,cex.axis=1.3)
		}else{
			axis(side = 1, at = break_points, labels = c("-2KB","-1KB","-500","TES","+500","+1KB","+2KB"))#,cex.axis=1.3)
		}
		abline(v=0,lty=2,col="darkgrey", lwd=2)
	    abline(v=TSSbreak_points,lty=3,col="darkgrey", lwd=2)
	}
	legend("topleft",legend=c("mean","2*stdErr"),fill=c("black","lightblue"),bg="white")
	if (!is.null(savePlotPath))
    {
    	dev.off()
    	print(paste("pdf saved under",filename,sep=" "))
	}
}

##density plot for QC-value distribution versus a single value
f_plotValueDistribution = function(compendium,title,coordinateLine,savePlotPath=NULL)
{   
	if (!is.null(savePlotPath))
    {
    	filename=file.path(savePlotPath,"PlotValueDistribution.pdf")
		pdf(file=filename,width=10, height=7)
	}
    mycolors= c("lightsteelblue1","lightsteelblue3")
	maxi=max(compendium)
    mini=min(compendium)

    ##get density
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

    abline(v=coordinateLine,,col="red",lwd=3,lty=2)

    if (!is.null(savePlotPath))
    {
    	dev.off()
    	print(paste("pdf saved under",filename,sep=" "))
	}    
}
