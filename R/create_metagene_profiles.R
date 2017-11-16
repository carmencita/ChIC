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




#gd.r1_TSS <- lapply(sd.r1, t.get.gene.av.density_TSS)
#tl_current=sd.r1$ENCFF001GVI
# a function for getting bin-averaged profiles for individual genes TSS ONE.POINT
t.get.gene.av.density_TSS <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F) { ##binning of the frame
    chrl <- names(gdl);
    names(chrl) <- chrl;
    lapply(chrl[chrl %in% names(tl_current$td)],function(chr) {
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

require(spp)
#require(snow) # first load snow
#require(Rmpi) # thel load Rmpi (the reverse sequence won't work)
cluster<-NULL


#strandShift.matrix= read.table("/data/FF/Carmen/Pipeline/Summary.Results",header=T)
#current_tagShift=strandShift.matrix[which(strandShift.matrix$Name==datafilename),]$StrandShift
#tag.shift <- round(current_tagShift/2)
tag.shift= round(strandShift/2)

#get readcount using specific aligner 
#read.tags.current_function<-get(paste("read", reads.aligner.type , "tags", sep="."))
#chip.data<-read.tags.current_function(file.path(path,paste(sampleinfo$Filename[sampleIndex],"bam",sep=".")))
###step 1: tag distribution 
#get binding characteristics 
#chip.data_tagdistribution<-sapply(chip.data$tags, length)
#binding.characteristics<-get.binding.characteristics(chip.data, srange=estimating_fragment_length_range, bin=estimating_fragment_length_bin, cluster=cluster)


## format as tag lists
#chipTags <- list(chip.dataSelected )
#names(chipTags)<-chip.data.samplename
#names(chipTags)<-datafilename
#inputTags<-input.dataSelected 
#rm(chip.data,input.data)
#gc()



## load gene annotations
require(girafe)

print("Load geneannotation")
load(geneAnnotations_file) #RefSeqGenesAll_object
current_annotations_type<-gsub(pattern=".RData", replacement="", fixed=TRUE, x=basename(geneAnnotations_file))
#current_annotations_object_name<-paste(current_annotations_type, "_genomeIntervals_object", sep="")
#current_annotations_object<-get(current_annotations_object_name)
current_annotations_object=Genes_genomeIntervals_object

##FILTER genes, otherwise it takes too much time
#filterGenes=read.delim(paste(path,"HouseKeepingGenes_wikicell.org.txt",sep=""),header=TRUE) #X.Gene and geneSymbol
#filtered_current_annotations_object=current_annotations_object[which(current_annotations_object$geneSymbol %in% filterGenes$X.Gene==TRUE),]


# format annotations as chromosome lists
current_annotations_object<-data.frame(current_annotations_object@.Data, annotation(current_annotations_object), stringsAsFactors=FALSE)
current_annotations_object$interval_starts<-as.integer(current_annotations_object$interval_starts)
current_annotations_object$interval_ends<-as.integer(current_annotations_object$interval_ends)
current_annotations_object$seq_name<-as.character(current_annotations_object$seq_name)

annotatedGenesPerChr <-split(current_annotations_object, f=current_annotations_object$seq_name)


#MetaGeneDir="/data/FF/Carmen/Pipeline/MetaGenome"

# genomic density, density (x and y coordinates per chromosome)
#sd.gen <- get.smoothed.tag.density(inputTags,bandwidth=bw,step=step,tag.shift=tag.shift, rngl=rngl)
#ts <- sum(unlist(lapply(inputTags,length)))/1e6
#sd.gen <- lapply(sd.gen,function(d) { d$y <- d$y/ts; return(d); }) 

#ts <- sum(unlist(lapply(input.dataSelected,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
#smoothed.density <- get.smoothed.tag.density(input.dataSelected, bandwidth=wig.bw, step=wig.step,tag.shift=tag.shift)
#smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })
#save(smoothed.density,file=paste(path,"TagDensity_",inputname,".RData",sep=""))

print("Load TagDensisty input")
load(paste(path,"TagDensity_",datafilename,"_",inputname,".RData",sep=""))
#input.smoothed.density=smoothed.density
sd.gen=smoothed.density
sd.gen=list(td=sd.gen)
#names(sd.gen)<-inputname
rm(smoothed.density)
gd.gen <- t.get.gene.av.density(sd.gen)


#gd.gen <- t.get.gene.av.density(list(td=sd.gen)) 

# calculate ChIP samples densities
#bw = MetaprofilesSmoothingBandwidth  ##is slidingwindow shift
#step <- MetaprofilesSmoothingStep ##is slidingwindow size

#sd.r1 <- lapply(chipTags,function(chipTags_intern, current_tagshift=tag.shift) { ## splits the tag distribution per chromosmoem
#	ts <- sum(unlist(lapply(chipTags_intern,length)))/1e6;
#	td <- get.smoothed.tag.density(chipTags_intern, bandwidth=bw,step=step,tag.shift=current_tagshift, rngl=rngl); # tagsfhift is half strandshift#
#	td <- lapply(td,function(d) { d$y <- d$y/ts; return(d); })
#    	return(list(td=td));
#})

#ts <- sum(unlist(lapply(chip.dataSelected,length)))/1e6 ##tag smoothing, (sum of tags in all chr)/1e6
#smoothed.density <- get.smoothed.tag.density(chip.dataSelected, bandwidth=wig.bw, step=wig.step,tag.shift=tag.shift)
#smoothed.density<-lapply(smoothed.density,function(d) { d$y <- d$y/ts; return(d); })

#save(smoothed.density,file=paste(path,"TagDensity_",datafilename,".RData",sep=""))
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
save(gd.gfp, gd.gen_TSS, file=file.out)
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
#print("Load control file...")
#newControl=NULL
#if (length(unlist(strsplit(inputname,"\\.")))>1)
#{
#	print("Concat multiple inputs...")
#	for (element in unlist(strsplit(inputname,"\\.")))
#	{	
#		if (element!="")
#		{
#			if (is.null(newControl))
#			{	
#				##if newControl is empty, fill it with the first input-data
#				print(element)
#				load(file.path(datadir, paste(element,".RData",sep="")))
#				newControl<-input.data
#				print("original size")
#				for (i in seq(length(newControl$tags)))
#				{
#					print(length(newControl$tags[[i]]))
#				}
#				#newControl.tags=(sapply(newControl$tags,c))
#				#newControl.q=(sapply(newControl$quality,c))
#				#print(length(newControl.tags[[25]]))
#				#print(length(newControl.q[[25]]))
#			}else{
#				##otherwise add the input-data and sort
#				print(element)
#				load(file.path(datadir, paste(element,".RData",sep="")))
#				temp<-input.data
#				#temp.tags=(sapply(temp$tags,c))
#				#temp.q=(sapply(temp$quality,c))
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

#	##sort tags and quality of the final NewControl to get an increasing list independently from the sign:  -57 -95 -112 151 159 166 169 -217...
#	print("Sort tags and quality flags...")
#	for (i in seq(length(newControl$tags)))
#	{
#		index=sort.list(abs(newControl$tags[[i]]))
#		newControl$tags[[i]]=newControl$tags[[i]][index]	
#		newControl$quality[[i]]=newControl$quality[[i]][index]	
#	}
#
#	input.data=newControl
#}else{
#
#	#inputIndex=as.integer(arg[8])
#	print("Single input...")
#	load(file.path(datadir, paste(inputname,".RData",sep="")))
#	input.data<-input.data
#}



##   print("Filter tags")
##   # select informative tags based on the binding characteristics
##   chip.data <- select.informative.tags(chip.data,binding.characteristics)
##  input.data <- select.informative.tags(input.data,binding.characteristics)
## 
## here we skip this step because it requires using of binding characteristics that we have manually modified to force the selection of tag shift within user specified window (100-200)
#
#     print("Filter tags")
#    if (select.informative.tags_filter_metaprofiles) {
#      print("select.informative.tags filter")
#      load(paste("sppdata", "binding", chip.data.samplename, "RData", sep="."))
#      chip.data <- select.informative.tags(chip.data, binding.characteristics)
#       load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
#       input.data <- select.informative.tags(input.data, binding.characteristics)
#    } else {
#      print("SKIP select.informative.tags filter")
#      chip.data<-chip.data$tags
#       input.data<-input.data$tags
#    }
#    if (remove.local.tag.anomalies_filter_metaprofiles) {
#      print("remove.local.tag.anomalies filter")
#      # restrict or remove singular positions with very high tag counts
#      chip.data <- remove.local.tag.anomalies(chip.data)
#      input.data <- remove.local.tag.anomalies(input.data)
#    } else {
#      print("SKIP remove.local.tag.anomalies filter")
#    }


#print("Filter tags")
#if (select.informative.tags_filter_metaprofiles) {
#      print("select.informative.tags filter")
#     #load(paste("sppdata", "binding", chip.data.samplename, "RData", sep="."))
#      chip.dataSelected <- select.informative.tags(chip.data, binding.characteristics)
#      #load(paste("sppdata", "binding", input.data.samplename, "RData", sep="."))
#      input.dataSelected <- select.informative.tags(input.data, binding.characteristics)
#} else {
#      print("SKIP select.informative.tags filter")
#      chip.dataSelected<-chip.data$tags
#      input.dataSelected<-input.data$tags
#}
#
#
#if (remove.local.tag.anomalies_filter_metaprofiles) {
#	print("remove.local.tag.anomalies filter")
#      # restrict or remove singular positions with very high tag counts
#	chip.dataSelected <- remove.local.tag.anomalies(chip.dataSelected)
#	input.dataSelected <- remove.local.tag.anomalies(input.dataSelected)
#} else {
#	print("SKIP remove.local.tag.anomalies filter")
#}





