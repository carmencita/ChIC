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


require(spp)

## load gene annotations
require(girafe)
require(snow)
require(parallel)

# calculates bin average values for a list of one- or two-point features
# dat - $x/$y data frame
# feat - data frame with either $x/{$strand} or $s/$e/{$strand}
# return.scaling causes function to just calculate, cleanup and return scaling mapping for further use
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



# identifies points located within margins of given segments
# and returns rescaled relative positions
# m - margin; om - outer margin; im - inner margin; rom - right outer margin, etc.
# x - x positions being mapped/rescaled
# seg - $s/$e/{$strand} data frame
# return $i/$rx/$nu data frame corresponding to indecies of original x values ($i) and new relative $rx positions, and whether a point maps to several segments
# $r5x/$r3x give distances relative to the 5' and 3' ends, $si gives segment index
# $bin index factor is added if nbins is supplied
#two.point.scaling <- function(x,seg,m=1e3,bs=gene_body,om=m,im=m,rom=om,lom=om,rim=im,lim=im,nbins=predefnum_bins/2) {
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



# much simpler, one-point scaling that returns relative positions in the format similar to the two-point scaling ($i/$rx/$nu/$si)
#one.point.scaling <- function(x, pos, strand=NULL, m=1e3, lm=m, rm=m, nbins=predefnum_bins/2) {  m=2010
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




# a function for getting bin-averaged profiles for individual genes TWO.POINT
#t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=500,lom=5e3, rom=1e3, bs=2e3, nbins=400,separate.strands=F) {

# a function for getting bin-averaged profiles for individual genes TWO.POINT
#t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=500,lom=5e3, rom=1e3, bs=2e3, nbins=400,separate.strands=F) {
t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=inner_margin,lom=left_outer_margin, rom=right_outer_margin, bs=gene_body, nbins=predefnum_bins,separate.strands=F, min.feature.size=min.feature.size_carmen) {
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




#gd.r1_TSS <- lapply(sd.r1, t.get.gene.av.density_TSS)
#tl_current=sd.r1$ENCFF001GVI
# a function for getting bin-averaged profiles for individual genes TSS ONE.POINT
t.get.gene.av.density_TSS <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F) { ##binning of the frame
    chrl <- names(gdl);
    names(chrl) <- chrl;
    #lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
    mclapply(chrl[chrl %in% names(tl_current$td)],  mc.preschedule = FALSE,mc.cores=mc, FUN=function(chr) {
      nsi <- gdl[[chr]]$strand=="-";
      current_gene_names<-gdl[[chr]]$geneSymbol

       # if ((sum(!nsi)>0)) {
	if ((sum(!nsi)>0)) {
          px <- feature.bin.averages(tl_current$td[[chr]],data.frame(x=gdl[[chr]]$txStart[!nsi]),m=m,nbins=nbins,nu.point.omit=FALSE)
	  #str(px)
          rownames(px)<-current_gene_names[!nsi]
        } else { 
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
    #lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
    mclapply(chrl[chrl %in% names(tl_current$td)],  mc.preschedule = FALSE,mc.cores=mc, FUN=function(chr) {
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

#path="/gpfs/work/IscrC_CONCEPt/QC_pipeline/"
#datafilename="ENCFF000OYO"
#read_length=50


###############################################################
###############################################################
####
#### END OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
####
###############################################################
###############################################################


#################################
#################################

#setwd(workingdir)
#sampleinfo<-read.table(sampleinfo_file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
#rownames(sampleinfo)<-sampleinfo[,"Samplename"]


arg <- commandArgs()
print(arg)

#sampleIndex=as.integer(arg[7])
path=arg[6]
datafilename=arg[7]
print(datafilename)

print ("get strandshift")
#strandShift=as.integer(arg[8])
ssfile=read.delim(paste(path,"CrossCorrelation/",datafilename,".results",sep=""))
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


#################################
#################################
## scaling over total chromsome size (fixed)

chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
rownames(chrominfo)<-chrominfo$chrom

rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(c(1,x))})

#################################
#################################
###
### Analysis LOOP
###
#################################
#################################

#strandShift.matrix= read.table("/data/FF/Carmen/Pipeline/Summary.Results",header=T)
#current_tagShift=strandShift.matrix[which(strandShift.matrix$Name==datafilename),]$StrandShift
#tag.shift <- round(current_tagShift/2)
tag.shift= round(strandShift/2)




print("Load geneannotation")
load(geneAnnotations_file) #RefSeqGenesAll_object
current_annotations_type<-gsub(pattern=".RData", replacement="", fixed=TRUE, x=basename(geneAnnotations_file))
#current_annotations_object_name<-paste(current_annotations_type, "_genomeIntervals_object", sep="")
#current_annotations_object<-get(current_annotations_object_name)
#current_annotations_object=Genes_genomeIntervals_object
current_annotations_object=RefSeqGenes_annotated_filteredByOverlap_geneLength

##FILTER genes, otherwise it takes too much time
#filterGenes=read.delim(paste(path,"HouseKeepingGenes_wikicell.org.txt",sep=""),header=TRUE) #X.Gene and geneSymbol
#filtered_current_annotations_object=current_annotations_object[which(current_annotations_object$geneSymbol %in% filterGenes$X.Gene==TRUE),]


# format annotations as chromosome lists
current_annotations_object<-data.frame(current_annotations_object@.Data, annotation(current_annotations_object), stringsAsFactors=FALSE)
current_annotations_object$interval_starts<-as.integer(current_annotations_object$interval_starts)
current_annotations_object$interval_ends<-as.integer(current_annotations_object$interval_ends)
current_annotations_object$seq_name<-as.character(current_annotations_object$seq_name)

annotatedGenesPerChr <-split(current_annotations_object, f=current_annotations_object$seq_name)



#mc <- getOption("mc.cores",mpi.universe.size() )
mc=1

print("Load TagDensisty input")
load(paste(storagedir,"TagDensity_",datafilename,"_",inputname,".RData",sep=""))
#input.smoothed.density=smoothed.density
sd.gen=smoothed.density
sd.gen=list(td=sd.gen)
#names(sd.gen)<-inputname
rm(smoothed.density)
gd.gen <- t.get.gene.av.density(sd.gen)


#save(smoothed.density,file=paste(path,"TagDensity_",datafilename,".RData",sep=""))
print("Load TagDensisty chip")
load(paste(storagedir,"TagDensity_",datafilename,".RData",sep=""))
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
rm(gd.gen)




print("Calculate one point scaling")
##one.point.scaling o mi calcolo in numero di bin cosi: up_downStream /10
#gd.r1_TSS <- lapply(sd.r1, t.get.gene.av.density_TSS)
gd.r1_TSS <- t.get.gene.av.density_TSS(sd.r1)
#gd.gen_TSS <- t.get.gene.av.density_TSS(list(td=sd.gen))
gd.gen_TSS <- t.get.gene.av.density_TSS(sd.gen)
gd.gfp<-gd.r1_TSS#[[1]] 
file.out<-paste(outputdir,"/one.point.scalingTSS_",datafilename, ".RData", sep="")
print(file.out)
save(gd.gfp, gd.gen_TSS, file=file.out)
rm(gd.gen_TSS)

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





#file.remove((paste(storagedir,"TagDensity_",datafilename,"_",inputname,".RData",sep="")))
#file.remove((paste(storagedir,"TagDensity_",datafilename,".RData",sep="")))


