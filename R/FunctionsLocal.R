## load gene annotations
require(girafe)
require(snow)
require(parallel)

##################
#private variables
##################
min.feature.sizeMy=3000
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

###############################################################
###############################################################
####
#### BEGINF OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
####
###############################################################
###############################################################
####
#### - f_sfeature.bin.averages
#### - two.point.scaling
#### - one.point.scaling
####
###############################################################
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
###############################################################
####
#### END OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
####
###############################################################
###############################################################


##MODIFIED VERSION OF t.get.gene.av.density for TWO.POINT
# a function for getting bin-averaged profiles for individual genes TWO.POINT
#t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=500,lom=5e3, rom=1e3, bs=2e3, nbins=400,separate.strands=F) {
f_t.get.gene.av.density <- function(chipTags_current,gdl=annotatedGenesPerChr,im=inner_margin,lom=left_outer_margin, rom=right_outer_margin, bs=gene_body, nbins=predefnum_bins,separate.strands=F, min.feature.size=min.feature.sizeMy) {
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
f_t.get.gene.av.density_TSS <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F) { ##binning of the frame
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
f_t.get.gene.av.density_TES <- function(tl_current,gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F) {
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
