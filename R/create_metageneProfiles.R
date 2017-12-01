###############################################################
###############################################################
####
#### BEGINF OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
####
###############################################################
###############################################################
####
#### - feature.bin.averages
#### - two.point.scaling
#### - one.point.scaling
####
###############################################################
###############################################################

feature.bin.averages <- function(dat,feat,nu.feat.omit=F,nu.point.omit=T, scaling=NULL, return.scaling=F, trim=0, min.feature.size=NULL, ... ) {
  if(is.null(feat)) { return(NULL) };
  if(dim(feat)[1]<1) { return(NULL) };

  if(is.null(scaling)) { # calculate scaling
    if(!is.null(feat$e)) { # two-point
      scaling <- two.point.scaling(dat$x,feat, ... )
      #scaling <- two.point.scaling(dat$x,feat,om=2e3,im=1e3,bs=2e3,nbins=100)
    } else { # one-point
      scaling <- one.point.scaling(dat$x,feat$x,feat$strand, ... );
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

two.point.scaling <- function(x,seg,bs=gene_body,om=m,im=m,rom=om,lom=om,rim=im,lim=im,nbins=predefnum_bins) {
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

one.point.scaling <- function(x, pos, strand=NULL,m=up_downStream/2, lm=m, rm=m, nbins=predefnum_bins/2) {
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

t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=inner_margin,lom=left_outer_margin, rom=right_outer_margin, bs=gene_body, nbins=predefnum_bins,separate.strands=F, min.feature.size=min.feature.size_carmen) {
	chrl <- names(gdl);
	names(chrl) <- chrl;
	lapply(chrl[chrl %in% names(chipTags_current$td)],function(chr) {
		print(chr)
		nsi <- gdl[[chr]]$strand=="-";
		print(length(nsi))
		current_gene_names<-gdl[[chr]]$geneSymbol
		if ((sum(!nsi)>0)) {
		#if ((sum(!nsi)>1)) {
			px <- feature.bin.averages(chipTags_current$td[[chr]],data.frame(s=gdl[[chr]]$txStart[!nsi],e=gdl[[chr]]$txEnd[!nsi]),lom=lom,rom=rom,im=im,bs=bs, nbins=nbins, min.feature.size=min.feature.size, nu.point.omit=FALSE)
			rownames(px)<-current_gene_names[!nsi]
			} else { 
			  px<-NULL
			}
		if ((sum(nsi)>0)) {
		#if ((sum(nsi)>1)) {
			  nd <- chipTags_current$td[[chr]]; nd$x <- -1*nd$x;
			  nx <- feature.bin.averages(nd,data.frame(s=-1*gdl[[chr]]$txEnd[nsi],e=-1*gdl[[chr]]$txStart[nsi]), lom=lom,rom=rom,im=im,bs=bs, nbins=nbins,min.feature.size=min.feature.size, nu.point.omit=FALSE)
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

t.get.gene.av.density_TSS <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F){ 
    chrl <- names(gdl);
    names(chrl) <- chrl;
    lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
      nsi <- gdl[[chr]]$strand=="-";
      current_gene_names<-gdl[[chr]]$geneSymbol
	if ((sum(!nsi)>0)) {
          px <- feature.bin.averages(tl_current$td[[chr]],data.frame(x=gdl[[chr]]$txStart[!nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
	  #str(px)
          rownames(px)<-current_gene_names[!nsi]
  }else { 
          px<-NULL
    }
	if ((sum(nsi)>0)) {
       #if ((sum(nsi)>0)) {
          nd <- tl_current$td[[chr]]; nd$x <- -1*nd$x;
          nx <- feature.bin.averages(nd,data.frame(x=-1*gdl[[chr]]$txEnd[nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
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


# a function for getting bin-averaged profiles for individual genes TES ONE.POINT 
t.get.gene.av.density_TES <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F) {
    chrl <- names(gdl);
    names(chrl) <- chrl;
    lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
      nsi <- gdl[[chr]]$strand=="-";
      current_gene_names<-gdl[[chr]]$geneSymbol

        #if ((sum(!nsi)>0)) {
        if ((sum(!nsi)>0)) {
          px <- feature.bin.averages(tl_current$td[[chr]],data.frame(x=gdl[[chr]]$txEnd[!nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
          rownames(px)<-current_gene_names[!nsi]
        } else { 
          px<-NULL
        }

       #if ((sum(nsi)>1)) {
       if ((sum(nsi)>0)) {
          nd <- tl_current$td[[chr]]; nd$x <- -1*nd$x;
          nx <- feature.bin.averages(nd,data.frame(x=-1*gdl[[chr]]$txStart[nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
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



###############################################################
###############################################################
####
#### END OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
####
###############################################################
###############################################################
CreateMetageneProfile = function(path=NULL,datafilename=NULL){

  cat ("[1] Getting strandshift...")
  ssfile=read.delim(paste(path,"CrossCorrelation/out/",datafilename,".results",sep=""))
  ssindex=grep("Substitution of StrandShift from",ssfile[,1]) 
  strandShift=as.integer(strsplit(as.character(ssfile[ssindex,])," ")[[1]][10])
  print(strandShift)

  source(paste(path,"GlobalParameters.R",sep=""))
  sampleinfo_file<-paste(path,"matchlist.txt",sep="")
  sampleinfo<-read.table(sampleinfo_file,  header=TRUE, quote="", stringsAsFactors=FALSE)

  t0<-proc.time()[3]


  workingdir<-paste(path,"/MetaGene",sep="")
  outputdir<-paste(workingdir,"/out",sep="")
  timedir<-paste(workingdir,"/time_stamps",sep="")
  plotsdir<-paste(workingdir,"/plots",sep="")

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

  chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
  rownames(chrominfo)<-chrominfo$chrom

  rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(c(1,x))})

  require(spp)
  cluster<-NULL

  tag.shift= round(strandShift/2)


  ## load gene annotations
  require(girafe)

  print("Load geneannotation")
  load(geneAnnotations_file) #RefSeqGenesAll_object
  current_annotations_type<-gsub(pattern=".RData", replacement="", fixed=TRUE, x=basename(geneAnnotations_file))
  #current_annotations_object_name<-paste(current_annotations_type, "_genomeIntervals_object", sep="")
  #current_annotations_object<-get(current_annotations_object_name)
  current_annotations_object=Genes_genomeIntervals_object

  # format annotations as chromosome lists
  current_annotations_object<-data.frame(current_annotations_object@.Data, annotation(current_annotations_object), stringsAsFactors=FALSE)
  current_annotations_object$interval_starts<-as.integer(current_annotations_object$interval_starts)
  current_annotations_object$interval_ends<-as.integer(current_annotations_object$interval_ends)
  current_annotations_object$seq_name<-as.character(current_annotations_object$seq_name)

  annotatedGenesPerChr <-split(current_annotations_object, f=current_annotations_object$seq_name)

  print("Load TagDensisty input")
  load(paste(path,"TagDensity_",datafilename,"_",inputname,".RData",sep=""))
  #input.smoothed.density=smoothed.density
  sd.gen=smoothed.density
  sd.gen=list(td=sd.gen)
  #names(sd.gen)<-inputname
  rm(smoothed.density)
  gd.gen <- t.get.gene.av.density(sd.gen)

  print("Load TagDensisty chip")
  load(paste(path,"TagDensity_",datafilename,".RData",sep=""))
  sd.r1=smoothed.density
  sd.r1 <-list(td=sd.r1) ## need a list of list structure (nested list)
  rm(smoothed.density)
  #names(sd.r1)="td"
  #sd.r1=list(sd.r1)
  #names(sd.r1)<-datafilename

  ##two.point.scaling\
  print("Calculate two point scaling")
  #gd.r1 <- lapply(sd.r1, t.get.gene.av.density) ##frame with chromosome:genes:startsites:density
  gd.r1 <- t.get.gene.av.density(sd.r1) ##frame with chromosome:genes:startsites:density
  gd.gfp<-gd.r1#[[1]] 	##%%%### ONLY ONE SAMPLE ### the first... there is only one sample btw
  #nam<-datafilename
  file.out<-paste(outputdir,"/two.point.scaling_",datafilename, ".RData", sep="")
  save(gd.gfp, gd.gen, file=file.out)
  #rm(gd.gen)

  print("Calculate one point scaling")
  ##one.point.scaling o mi calcolo in numero di bin cosi: up_downStream /10
  #gd.r1_TSS <- lapply(sd.r1, t.get.gene.av.density_TSS)
  gd.r1_TSS <- t.get.gene.av.density_TSS(sd.r1)
  #gd.gen_TSS <- t.get.gene.av.density_TSS(list(td=sd.gen))
  gd.gen_TSS <- t.get.gene.av.density_TSS(sd.gen)
  gd.gfp<-gd.r1_TSS#[[1]] 
  file.out<-paste(outputdir,"/one.point.scalingTSS_",datafilename, ".RData", sep="")
  print(file.out)
  
  #rm(gd.gen_TSS)

  gd.r1_TES <- t.get.gene.av.density_TES(sd.r1)
  #gd.gen_TES <- t.get.gene.av.density_TES(list(td=sd.gen))
  gd.gen_TES <- t.get.gene.av.density_TES(sd.gen)
  gd.gfp<-gd.r1_TES#[[1]] 
  file.out<-paste(outputdir,"/one.point.scalingTES_",datafilename, ".RData", sep="")
  print(file.out)
  save(gd.gfp, gd.gen_TES, file=file.out)
  rm(gd.gen_TES)

  t1<-proc.time()[3]
  deltat=t1-t0
  write(paste("prepare_profiles",datafilename,deltat,sep=" "),file=paste(timedir,"/timing",datafilename,"_",TFname,".txt",sep=""),append=TRUE)

  file.remove((paste(path,"TagDensity_",datafilename,"_",inputname,".RData",sep="")))
  file.remove((paste(path,"TagDensity_",datafilename,"_",inputname,".RData",sep="")))
  return(list(gd.gfp, gd.gen_TSS))
}





