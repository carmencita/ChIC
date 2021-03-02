#####################################################################  
#####FUNCTIONS QC-metrics for narrow binding PROFILES ###########
##################################################################### 

#' @keywords internal 
## format conversion (spp; Koustav)
f_convertFormatBroadPeak <- function(given.clusters) {
    chrl <- names(given.clusters)
    names(chrl) <- chrl
    chrl <- chrl[unlist(lapply(given.clusters, function(d) length(d$s))) > 0]
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- given.clusters[[chr]]
        data.frame(chr = chr, start = df$s, end = df$e, name = ".", 
            score = "0", 
            strand = ".", rv = df$rv, A1 = -1, A2 = -1)
    }))
    md <- md[order(as.numeric(md[, 7]), decreasing = TRUE), ]
    md <- data.frame(md)
    md.minimal <- md[, c(1, 2, 3)]
    return(md.minimal)
}


#' @keywords internal 
## format conversion (spp;Koustav)
f_converNarrowPeakFormat <- function(bd, margin = bd$whs) {
    if (is.null(margin)) {
        margin <- 50
    }
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        x <- df$x
        rs <- df$rs
        if (is.null(rs)) {
            rs <- rep(NA, length(x))
        }
        re <- df$re
        if (is.null(re)) {
            re <- rep(NA, length(x))
        }
        ivi <- which(is.na(rs))
        if (any(ivi)) {
            rs[ivi] <- x[ivi] - margin
        }
        ivi <- which(is.na(re))
        if (any(ivi)) {
            re[ivi] <- x[ivi] + margin
        }
        data.frame(chr = chr, start = rs, end = re, name = ".", 
            score = "0", strand = ".", 
            y = df$y, A1 = -1, 
            fdr = format(df$fdr, scientific = TRUE, digits = 3), 
            dist = x - rs)
    }))
    md <- md[order(as.numeric(md[, 7]), decreasing = TRUE), ]
    md <- data.frame(md)
    md.minimal <- md[, c(1, 2, 3)]
    return(md.minimal)
}





#' @keywords internal 
## writing wig files 
## writewig function of spp with different track type header
f_writewig <- function (dat, fname, feature, threshold = 5, zip = FALSE) 
{
    chrl <- names(dat)
    names(chrl) <- chrl
    invisible(lapply(chrl, function(chr) {
        bdiff <- dat[[chr]]
        ind <- seq(1, length(bdiff$x))
        ind <- ind[!is.na(bdiff$y[ind])]
        header <- chr == chrl[1]
        f_write.probe.wig(chr, bdiff$x[ind], bdiff$y[ind], fname, 
            append = !header, feature = feature, header = header)
    }))
    if (zip) {
    zf <- paste(fname, "zip", sep = ".")
        system(paste("zip \"", zf, "\" \"", fname, "\"", sep = ""))
        system(paste("rm \"", fname, "\"", sep = ""))
        return(zf)
    } else {
        return(fname)
    }
}


#' @keywords internal 
## function used by writewig function 
##Here I am using a slightly changed version of write.probe.wig from spp.
##This version changes the track type in the header to "bedGraph"
##ORIGINAL Version can be found in spp package!!!
f_write.probe.wig <- function(chr, pos, val, fname, append=FALSE, feature="M", 
    probe.length=35, header=TRUE) 
{
    min.dist <- min(diff(pos));
    if(probe.length>=min.dist) {
        probe.length <- min.dist-1;
        cat("warning: adjusted down wig segment length to", probe.length,"\n");
    }
    mdat <- data.frame(chr, as.integer(pos), as.integer(pos+probe.length), val)

    if(header) {
        write(paste("track type=bedGraph name=\"Bed Format\" description=\"", 
            feature,"\" visibility=dense color=200,100,0 altColor=0,100,200 
            priority=20",sep=""),file=fname,append=append)
        write.table(mdat, file=fname, col.names=FALSE, row.names=FALSE, 
            quote=FALSE, sep=" ", append=TRUE);
    } else {
        write.table(mdat, file=fname, col.names=FALSE, row.names=FALSE, 
        quote=FALSE, sep=" ", append=append);
    }
}


#' @keywords internal 
## reads bam file or tagalign file
f_readFile <- function(filename, reads.aligner.type) {
    if (reads.aligner.type == "bam") {
        #calling our internal f_read.bam.tags instead than spp default one to handle paired end reads BAM files
        currentFormat <- get(paste("f_read", reads.aligner.type, "tags", sep = ".")) #calling our internal f_read.bam.tags
        data <- currentFormat(file.path(paste(filename, ".bam", sep = "")))
    } else if (reads.aligner.type == "tagAlign") {
        currentFormat <- get(paste("read", "tagalign", "tags", sep = "."))
        data <- currentFormat(file.path(paste(filename,".tagAlign", sep = "")))
    } else {
        stop("Only BAM or tagalign formats are currently supported")
    }
    
    ## readCount=sum(sapply(data$tags, length))
    ## readCount <- sum(unlist(lapply(data$tags, length)))
    readCount <- sum(lengths(data$tags))
    message(readCount," reads")
    ##double check data structure to make sure the structure contains
    ##two lists called $tags and $quality
    helper=data
    data=NULL
    data$tags=helper$tags
    data$quality=helper$quality
    return(data)
}

                               
#' @keywords internal 
## this function reads BAM files into a taglist object of SPP 
## this function code is derived from the original "read.bam.tags" from spp package by Peter Kharchenko 
## the version included here has been revised to hanlde BAM files containing paired end reads  
f_read.bam.tags <- function(filename,read.tag.names=FALSE,fix.chromosome.names=F) {
  #require(Rsamtools)
  if(!is.element("Rsamtools", installed.packages()[, 1])) {
    stop("Rsamtools Bioconductor package is now required for BAM file support. Please install")
  }

## this is setting the list of fileds to be extracted from the BAM file  
## note that "pos" is the mapping position on the reference sequence (chromosome name stored in "rname")
## whereas the "qname" field contains the ID of each sequencing reads (or read pairs for paired end reads in the bam file)
## check if paired end reads
  checkIfPAired<-any(bitwAnd((Rsamtools::scanBam(filename,param=Rsamtools::ScanBamParam(what="flag",flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]])$flag,0x1));
  
  ## which fileds will be read from the BAM file    
  ww <- c("flag","rname","pos","isize","strand","mapq","qwidth");
  if(read.tag.names || checkIfPAired) { ww <- c(ww,"qname") }; ## need this also to handle paired ends
  bam <- Rsamtools::scanBam(filename,param=Rsamtools::ScanBamParam(what=ww,flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]];
  ## removed unused levels (this is needed otehrwise by default we will have an element for each chromosome listed in the BAM header even if there are no reads for that chromosome)
  bam$rname <- droplevels(bam$rname)
   
## this is returning an empty tagglist object if the BAM file contains no valid alignment
  if(is.null(bam$pos) || length(bam$pos)==0) { return(list(tags=c(),quality=c())) }

## this is just creating a 1/0 vector for postiive/negative strand mapped reads
  strm <- as.integer(bam$strand=="+")

## this is checking if the BAM file is actually containing paired end reads
  if(checkIfPAired) { 
    # paired-end data
    ## for paired end data, we can select one (random) read out of the pair, so as to have equally represented both the positive and the negative strand mapped reads
    ## we must design the code so as to take 1 random read for each read ID (qname) so that we get 1 even if we have only one read mapped in the pair
    oneSelectedInPair<-unlist(tapply(X=1:length(bam$pos), INDEX=bam$qname, FUN=function(ii)  {
        if (length(ii)>1) {
            return(sample(ii, size=1))
         } else if (length(ii)==1) {
            return(ii)
         } else {
            stop("unexpected BAM file content format")
         }
    }))  # return only one for each pair (or one for each group if more than 2 alignments are present)
    ## the selection of indexes (oneSelectedInPair) is performed on the full vectors, thus we can use these indexes to perfrm subselections on the full BAM vectors
    rl <- list(tags=tapply(X=oneSelectedInPair, INDEX=bam$rname[oneSelectedInPair],function(ii) bam$pos[ii]*strm[ii]  - (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii])))
    rl <- c(rl,list(quality=tapply(X=oneSelectedInPair,INDEX=bam$rname[oneSelectedInPair],function(ii) bam$mapq[ii])))
    
    ## return also the read IDS (query name = "qname") if required
    if (read.tag.names) {
        rl <- c(rl,list(names=tapply(X=oneSelectedInPair,INDEX=bam$rname[oneSelectedInPair],function(ii) bam$qname[ii])))
    }
      
  } else {
    ## this is the "standard" workflow in case the BAM file contains only single end reads
    ## most of the ChIP-seq peaks calling alogorithms expect single end reads, disributed on both positive and negative strand
    ## this line of code is traversing the BAM file content (1:length(bam$pos)), chromosome by chromosome (bam$rname)
    ## and keeping the annotated position as 5'end of positive strand mapped reads (bam$pos[ii]*strm[ii] )
    ## or the "-" 3'end (i.e. the 5'-end of the reads mapped on hte negative strand (- (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii])))
    ## this "ifelse" condition for positive and negative strand reads is actually managed by the 1/0 vectors for strand,
    ## as it will change to "zero" either the first or the second element in the subtraction below
    rl <- list(tags=tapply(1:length(bam$pos),bam$rname,function(ii) bam$pos[ii]*strm[ii]  - (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii])))
    rl <- c(rl,list(quality=tapply(1:length(bam$pos),bam$rname,function(ii) bam$mapq[ii])))
    ## return also the read IDS (query name = "qname") if required
    if(read.tag.names) {
        rl <- c(rl,list(names=tapply(1:length(bam$pos),bam$rname,function(ii) bam$qname[ii])))
    }
  }

  if(fix.chromosome.names) {
    # remove ".fa"
    names(rl) <- gsub("\\.fa","",names(rl))
  }
  return(rl)
}

                                 
                                 
                                 
#' @keywords internal 
## filters canonical chromosomes 
f_clearChromStructure <- function(structure, annotationID) {
    ##delete empty chromosomes
    message("load chrom_info")

    chInf <- f_chromInfoLoad(annotationID)
    dCheck <- NULL
    ##check on regulare chromosomes
    ss=structure
    if (is.null(ss$tags)){
        dCheck <- structure[ which( names(ss) %in%  names(chInf) )]
        #if (sum(sapply(dCheck, is.null ))>0){
        if (sum ( unlist (lapply(dCheck, is.null )))>0){
            #dCheck<- dCheck[ -which( sapply(dCheck, is.null ))]
            dCheck<- dCheck[ -which( unlist (lapply(dCheck, is.null )))]
        }
    }else{
        dCheck$tags <- ss$tags[ which( names(ss$tags) %in%  names(chInf) )]
        dCheck$quality <- 
            ss$quality[ which( names(ss$quality) %in%  names(chInf) )]
        #if (sum(sapply(dCheck$tags, is.null ))>0) {
        if (sum (unlist (lapply (dCheck$tags, is.null )))>0) {
            #dCheck$tags <- 
            #dCheck$tags[ -which( sapply(dCheck$tags, is.null ))]
            dCheck$tags <- 
                dCheck$tags[ -which( unlist(lapply(dCheck$tags, is.null )))]
            dCheck$quality <- 
                dCheck$quality[ -which(  
                    unlist(lapply(dCheck$quality, is.null )))]
                #dCheck$quality[ -which( sapply(dCheck$quality, is.null ))]
        }
    }
    return(dCheck)
}


#' @keywords internal 
## caluclate QC tag
f_qflag <- function(RSC) {
    qflag <- NA
    if ((RSC >= 0) & (RSC < 0.25)) {
        qflag <- -2
    } else if ((RSC >= 0.25) & (RSC < 0.5)) {
        qflag <- -1
    } else if ((RSC >= 0.5) & (RSC < 1)) {
        qflag <- 0
    } else if ((RSC >= 1) & (RSC < 1.5)) {
        qflag <- 1
    } else if ((RSC >= 1.5)) {
        qflag <- 2
    }
    return(qflag)
}


#' @keywords internal 
## selet informative tags
f_removeLocalTagAnomalies <- function(chip, input, chip_b.characteristics, 
    #input_b.characteristics,
    remove.local.tag.anomalies, select.informative.tags) 
{
    message("Filter tags")
    
    if (select.informative.tags) {
        message("select.informative.tags filter")
        chipSelected <- select.informative.tags(chip, 
            chip_b.characteristics)
        #inputSelected <- select.informative.tags(input, 
        #    input_b.characteristics)
        inputSelected <- select.informative.tags(input, 
            chip_b.characteristics)

    } else {
        message("SKIP select.informative.tags filter")
        chipSelected <- chip$tags
        inputSelected <- input$tags
    }
    if (remove.local.tag.anomalies) {
        message("remove.local.tag.anomalies filter")
        ## restrict or remove singular positions with very high tag counts
        chipSelected <- remove.local.tag.anomalies(chipSelected)
        inputSelected <- remove.local.tag.anomalies(inputSelected)
    } else {
        message("SKIP remove.local.tag.anomalies filter")
        inputSelected <- input$tags
        chipSelected <- chip$tags
    }
    
    finalList <- list(input.dataSelected = inputSelected, 
        chip.dataSelected = chipSelected)
    return(finalList)
}

#' @keywords internal 
## takes dataSelected as input, parallel is the number of CPUs used for
## parallelization
f_tagDensity <- function(data, tag.shift, chromDef, mc = 1) {
    ## Smoothing parameters
    smoothingStep <- 20  ###is size of sw###before it was 10
    smoothingBandwidth <- 50
    
    ## density distribution for data
    message("Smooth tag density")
    ## tag smoothing, (sum of tags in all chr)/1e6
    ## ts <- sum(unlist(lapply(data, length)))/1e+06
    ts <- sum(lengths(data)) / 1e6
    ### parallelisation
    chromosomes_list <- names(data)
    ### creates a list of lists
    data <- lapply(chromosomes_list, FUN = function(x) {
        return(data[x])
    })
    
    #smoothed.density <- BiocParallel::bplapply(data, 
    #    BPPARAM = BiocParallel::MulticoreParam(workers = mc), 
    smoothed.density<-mclapply(data, 
        mc.preschedule = FALSE, 
        mc.cores  = mc,   
        FUN = function(current_chr_list) {
            current_chr <- names(current_chr_list)
            str(current_chr_list)
            if (length(current_chr) != 1) {
                stop("unexpected dataSelected structure")
            }
            spp::get.smoothed.tag.density(current_chr_list, 
                bandwidth = smoothingBandwidth, 
                step = smoothingStep, tag.shift = tag.shift, 
                rngl = chromDef[current_chr])
        })

    smoothed.density <- (unlist(smoothed.density, recursive = FALSE))
    ## normalizing smoothed tag density by library size
    smoothed.density <- lapply(smoothed.density, function(d) {
        d$y <- d$y/ts
        return(d)
    })
    return(smoothed.density)
}

#' @keywords internal 
# Author: Koustav Pal, Contact: koustav.pal@ifom.eu, Affiliation: IFOM - FIRC
# Institute of Molecular Oncology Why use this function and not the 
# GenomicRanges
# function makeGRangesFromDataFrame?  1. This function creates names out of
# coordinates and does not inherit names from a data.frame as 
# the default package
# function does. If it must be inherited, it can be specified with Names=Args.
# This provides better control on what the names are for each row and in some
# cases ensures that users are proactive when using the %in% construct. Cause
# the default behaviour is to construct new names, when users use an %in%
# construct they will encounter an error when older names are not recycled,
# forcing them to use the Names=Args parameter.  2. Using cbind > rbind >
# as.data.frame creates a scenario where first a matrix is created and then 
# cast
# as a data.frame. In a matrix all points must be of the same type. Therefore,
# when a single row/col contains a character type value all values are cast as
# factors.  makeGRangesFromDataFrame sticks to R principles and and produces no
# errors. Rather, as it expects to encounter a numeric type in start and end
# cols, all factors are cast back to numeric using the as.character >as.numeric
# construct. While this behaviour is helpful, it needlessly uses two additional
# steps to circumvent a probable bug.  MakeGRangesObject does no such thing and
# will return an error as IRanges does not accept factor values.
f_makeGRangesObject <- function(Chrom = NULL, Start = NULL, End = NULL, 
    Strand = NULL, Names = NULL, Sep = ":") 
{
    # require(GenomicRanges)
    if (is.null(Names)) {
        Names <- paste(Chrom, as.integer(Start), as.integer(End), sep = Sep)
    }
    if (is.null(Strand)) {
        Strand <- rep("*", length(Chrom))
    }
    Object <- GenomicRanges::GRanges(seqnames = Rle(Chrom), 
        ranges = IRanges(Start, end = End, names = Names), 
        strand = Rle(strand(Strand)))
    return(Object)
}
# Author: Koustav Pal, Contact: koustav.pal@ifom.eu, Affiliation: IFOM - FIRC
# Institute of Molecular Oncology Creates non-overlapping regions from
#' @keywords internal 
f_reduceOverlappingRegions <- function(ranges = NULL) 
{
    # Try it out.
    rangesByChrom <- split(ranges,seqnames(ranges))
    ##loop over chromosomes
    nonOverlappingRanges <- lapply(rangesByChrom,function(chromFrame){
        myRange <- chromFrame[order(start(chromFrame))]
        Chrom <- unique(as.vector(seqnames(myRange)))
        #print(Chrom)
        Q.name <- "Queries"
        S.name <- "Subject"
        ##getting overlap object with the hits
        HitsObject <- findOverlaps(query = myRange, subject = myRange)
        PairList.hits <- HitsObject[
            queryHits(HitsObject) != subjectHits(HitsObject)]
        PairList <- data.frame(queryHits(PairList.hits), 
            subjectHits(PairList.hits))
        colnames(PairList) <- c(Q.name,S.name)
        #3number of non overlapping reagions
        No.overlaps <- myRange[!(seq_along(myRange) %in% PairList[,Q.name])]
        if (nrow(PairList) == 0) {
            return(No.overlaps)
        }
        PairList.gt <- PairList[PairList[,S.name] > PairList[,Q.name],]
        PairList.lt <- PairList[PairList[,S.name] < PairList[,Q.name],]
        colnames(PairList.lt) <- c(S.name,Q.name)
        PairList.lt <- PairList.lt[,c(Q.name,S.name)]
        PairList <- rbind(PairList.lt,PairList.gt)
        PairList <- PairList[order(PairList[,Q.name]),]

        unique.queries <- unique(PairList[,Q.name])
        unique.subjects <- unique(PairList[,S.name])

        Which.q <- unique.queries[which(
            !(unique.queries %in% unique.subjects))]
        Which.s <- unique.subjects[which(
            !(unique.subjects %in% unique.queries))]

        #if(length(Which.q)!=length(Which.s)){
        #    stop("Cannot resolve overlaps. Contact the writer of the function
        #        to deconvolute the logic!\n")
        #}
        if(length(Which.q)!=length(Which.s)){
            message("Skipping for overlapping analysis ",Chrom)
            return(No.overlaps)
        }else{
            Starts <- start(myRange[Which.q])
            Ends <- sapply(seq_along(Which.s),function(x){
                Index <- Which.q[x]:Which.s[x]
                max(end(myRange[Index]))
            })

            NewRanges <- f_makeGRangesObject(Chrom=rep(Chrom,length(Starts)), 
                Start= Starts, End= Ends)
            return(c(No.overlaps,NewRanges))
        }
    })
    nonOverlappingRangesFinal <- do.call(c, 
        unlist(nonOverlappingRanges,use.names = FALSE))
    return(nonOverlappingRangesFinal)
}


###Orig function
#ReduceOverlappingRegions <- function(Ranges = NULL) 
#{
#    # Try it out.
#    Ranges.split <- split(Ranges,seqnames(Ranges))
#    Non.overlapping.ranges.list <- lapply(Ranges.split,function(Range){
#        My.Range <- Range[order(start(Range))]
#        Chrom <- unique(as.vector(seqnames(My.Range)))
#        Q.name <- "Queries"
#        S.name <- "Subject"
#        HitsObject <- findOverlaps(query = My.Range, subject = My.Range)
#        PairList.hits <- 
#        HitsObject[queryHits(HitsObject) != subjectHits(HitsObject)]
#        PairList <- data.frame(queryHits(PairList.hits), 
#            subjectHits(PairList.hits))
#        colnames(PairList) <- c(Q.name,S.name)
#        No.overlaps <- My.Range[!(seq_along(My.Range) 
#       %in% PairList[,Q.name])]
#        if (nrow(PairList) == 0) {
#            return(No.overlaps)
#        }
#        PairList.gt <- PairList[PairList[,S.name] > PairList[,Q.name],]
#        PairList.lt <- PairList[PairList[,S.name] < PairList[,Q.name],]
#        colnames(PairList.lt) <- c(S.name,Q.name)
#        PairList.lt <- PairList.lt[,c(Q.name,S.name)]
#        PairList <- rbind(PairList.lt,PairList.gt)
#        PairList <- PairList[order(PairList[,Q.name]),]
#
#       unique.queries <- unique(PairList[,Q.name])
#        unique.subjects <- unique(PairList[,S.name])
#
#        Which.q <- unique.queries[
#               which(!(unique.queries %in% unique.subjects))]
#        Which.s <- unique.subjects[
#               which(!(unique.subjects %in% unique.queries))]
#
#        if(length(Which.q)!=length(Which.s)){
#            stop("Cannot resolve overlaps. Contact the writer of the function
#                to deconvolute the logic!\n")
#        }
#
#        Starts <- start(My.Range[Which.q])
#        Ends <- sapply(seq_along(Which.s),function(x){
#            Index <- Which.q[x]:Which.s[x]
#            max(end(My.Range[Index]))
#        })
#
#        NewRanges <- MakeGRangesObject(Chrom=
#               rep(Chrom,length(Starts)),Start= Starts, End= Ends)
#        return(c(No.overlaps,NewRanges))
#    })
#    Non.overlapping.ranges <- do.call(c, 
#       unlist(Non.overlapping.ranges.list,use.names = FALSE))
#    return(Non.overlapping.ranges)
#}


##################################################################### 
######### FUNCTIONS QC-metrics for global read distribution ######### 


#' @keywords internal 
f_shortenFrame <- function(smothDens) {
    ## shorten frame with cumulative distribution 
    new <- lapply(smothDens, function(chromelement) {
        tt <- NULL
        tt$x <- chromelement$x[seq.int(1, length(chromelement$x), 6)]
        tt$y <- chromelement$y[seq.int(1, length(chromelement$y), 6)]
        return(tt)
    })
    
    # for (i in chrl) { ##print(i) 
    ###appending tags and quality elements of all
    # inputs to newControl new[[i]]$x=smothDens[[i]]$x[seq.int(1,
    # length(smothDens[[i]]$x), 6)] new[[i]]$y=smothDens[[i]]$y[seq.int(1,
    # length(smothDens[[i]]$y), 6)] }
    return(new)
}

#' @keywords internal 
## sorts and bins the dataframe
f_sortAndBinning <- function(shortframe) {
    shortframeSorted <- sort(shortframe)
    BINS_NUMBER <- 10000
    
    cumSum <- cumsum(shortframeSorted)
    cumSumBins <- quantile(cumSum, probs = seq(0, 1, (1/BINS_NUMBER)))
    pj <- (cumSumBins/cumSum[length(cumSum)])
    normalizedCumSum <- data.frame(x = seq(0, 1, (1/BINS_NUMBER)), pj = pj)
    rownames(normalizedCumSum) <- NULL
    return(normalizedCumSum)
}

## The plot function creates the 'Fingerprint plot' i.e. the cumulative
## distribution of the read counts across genomic bins for Input and ChIP. 
## cumChip
## The cumulative distribution of the read counts for the ChIP cumInput The
## cumulative distribution of the read counts for the Input savePlotPath

#' @keywords internal 
f_fingerPrintPlot <- function(cumChip, cumInput, savePlotPath = NULL) {
    
    if (!is.null(savePlotPath)) {
        filename <- file.path(savePlotPath, "FingerPrintPlot.pdf")
        pdf(file = filename, width = 7, height = 7)
    }
    plot(cumChip, type = "l", col = "blue", lwd = 2, 
        xlab = "Fraction of bins", 
        ylab = "Fr. of total reads coverage", 
        main = "Fingerprint: global read distribution")
    
    lines(cumInput, col = "red", lwd = 2)
    arrowx <- cumChip[which.max(abs(cumInput$pj - cumChip$pj)), ]$x
    abline(v = arrowx, col = "green", lty = 2, lwd = 2)
    ## abline(h=schneidePunktY,col='cyan',lty=2,lwd=2)
    ## abline(v=schneidePunktX,col='cyan',lty=2,lwd=2)
    legend("topleft", legend = c("Input", "ChIP"), fill = c("red", "blue"))
    if (!is.null(savePlotPath)) {
        dev.off()        
    }
}


##################################################################### 
########## FUNCTIONS QC-metrics for local ENRICHMENT 



############################################################
# BEGINF OF PETER Kharchenko's 
# FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES
############################################################
# - f_feature.bin.averages - f_two.point.scaling - f_one.point.scaling


#' @keywords internal 
## calculates bin average values for a list of one- or two-point 
## features dat -
## $x/$y data frame feat - data frame with either 
## $x/{$strand} or $s/$e/{$strand}
## return.scaling causes function to just calculate, cleanup and return scaling
## mapping for further use
f_feature.bin.averages <- function(dat, feat, nu.feat.omit = FALSE, 
    nu.point.omit = TRUE, scaling = NULL, return.scaling = FALSE, 
    trim = 0, min.feature.size = NULL, m = 4020/2, ...) 
{
    if (is.null(feat)) {
        return(NULL)
    }
    if (dim(feat)[1] < 1) {
        return(NULL)
    }
    
    if (is.null(scaling)) {
        ## calculate scaling
        if (!is.null(feat$e)) {
            ## two-point
            scaling <- f_two.point.scaling(dat$x, feat, ...)
            
        } else {
            scaling <- f_one.point.scaling(dat$x, feat$x, feat$strand, 
                m = m, ...)
        }
        ## clean up
        if (nu.feat.omit) {
            scaling <- scaling[!scaling$si %in% 
                unique(scaling$si[scaling$nu > 1]), ]
        } else {
            if (nu.point.omit) {
                scaling <- scaling[scaling$nu == 1, ]
            }
        }
        if (return.scaling) {
            return(scaling)
        }
    }
    
    if (!is.null(min.feature.size)) {
        if (!is.null(feat$e)) {
            iindex <- which(feat$e - feat$s >= min.feature.size)
            scaling <- scaling[scaling$si %in% iindex, ]
        }
    }
    ## determine and return gene bin average table
    return(tapply(dat$y[scaling$i], 
        list(factor(scaling$si, levels = c(1:dim(feat)[1])), scaling$bin),
        function(x) mean(na.omit(x))))
    
    ## return(tapply(dat$y[scaling$i],list(scaling$si,scaling$bin),
    ## mean,trim=trim,na.rm=T))
}

#' @keywords internal 
## identifies points located within margins of given segments and 
## returns rescaled
## relative positions m - margin; om - outer margin; 
## im - inner margin; rom - right outer margin, etc. 
## x - x positions being mapped/rescaled 
## seg - $s/$e/{$strand} data frame return $i/$rx/$nu data 
## frame corresponding to indecies of original x values ($i) 
## and new relative $rx positions, and whether
## a point maps to several segments $r5x/$r3x give distances relative to the 5'
## and 3' ends, $si gives segment index $bin index factor is added if nbins is
## supplied two.point.scaling <- function(x,seg,m=1e3,bs=gene_body,om=m,im=m,
## rom=om,lom=om, rim=im,lim=im,nbins=predefnum_bins/2) {
f_two.point.scaling <- function(x, seg, bs = 2000, im, rom, lom, nbins = 301) 
{
    rim <- im
    lim <- im
    ## map points to the segments defined by outer margins
    nseg <- length(seg$s)
    if (!is.null(seg$strand)) {
        nsi <- which(seg$strand == "-")
        ml <- rep(lom, nseg)
        ml[nsi] <- rom
        mr <- rep(rom, nseg)
        mr[nsi] <- lom
        sg5 <- seg$s
        sg5[nsi] <- seg$e[nsi]
        sg3 <- seg$e
        sg3[nsi] <- seg$s[nsi]
        sstrand <- ifelse(seg$strand == "-", -1, 1)
        sstrand[is.na(sstrand)] <- 1
    } else {
        ml <- rep(lom, nseg)
        mr <- rep(rom, nseg)
        sg5 <- seg$s
        sg3 <- seg$e
        sstrand <- rep(1, length(sg5))
    }
    
    
    spi <- spp::points_within(x, seg$s - ml, seg$e + mr, 
        return.list = TRUE)

    #spi <- spp::points_withinFunction(x, seg$s - ml, seg$e + mr, 
        #return.list = TRUE)
    
    lspi <- unlist(lapply(spi, length))
    
    #df <- data.frame(i = rep(1:length(x), lspi), si = unlist(spi), 
    #    nu = rep(lspi, lspi))
    df <- data.frame(i = rep(seq_along(x), lspi), si = unlist(spi), 
        nu = rep(lspi, lspi))

    df$r5x <- (x[df$i] - sg5[df$si]) * sstrand[df$si]
    df$r3x <- (x[df$i] - sg3[df$si]) * sstrand[df$si]
    
    ## between g3-rim and g3+rom - unscaled
    df$rx <- rep(NA, length(df$i))
    vi <- which(df$r3x >= -rim)
    df$rx[vi] <- df$r3x[vi] + lim + rim + bs
    
    ## df$rx <- df$r3x+lim+rim+bs; body scaling
    vi <- which(df$r3x < -rim & df$r5x > lim)
    df$rx[vi] <- ((df$r5x[vi] - lim)/
        ((seg$e - seg$s)[df$si[vi]] - lim - rim)) * bs + lim
    ## between g5-lom and g5+lim
    vi <- which(df$r5x <= lim)
    df$rx[vi] <- df$r5x[vi]
    
    if (!is.null(nbins)) {
        breaks <- seq(-lom, lim + bs + rim + rom, 
            by = (lom + lim + bs + rim + rom)/nbins)
        rn <- breaks[-1] - diff(breaks)/2
        breaks[1] <- breaks[1] - 1
        breaks[length(breaks)] <- breaks[length(breaks)] + 1
        df$bin <- cut(df$rx, breaks, labels = rn)
    }
    return(df)
}


#' @keywords internal 
## much simpler one-point scaling that returns relative positions in the format
## similar to the two-point scaling ($i/$rx/$nu/$si) one.point.scaling <-
## function(x, pos, strand=NULL, m=1e3, lm=m, rm=m, nbins=predefnum_bins/2) {
## m=2010 
## f_one.point.scaling <- function(x, pos, strand=NULL,m=up_downStream/2,
## lm=m, rm=m, nbins=predefnum_bins/2) 
## { f_one.point.scaling <- function(x, pos,
## strand=NULL,m=4020/2, lm=m, rm=m, nbins=301/2) {

f_one.point.scaling <- function(x, pos, strand = NULL,  nbins = 201, m = 4020/2)
{
    lm <- m
    rm <- m
    
    if (is.null(pos)) {
        return(NULL)
    }
    nseg <- length(pos)
    if (nseg < 1) {
        return(NULL)
    }
    ml <- rep(lm, nseg)
    mr <- rep(rm, nseg)
    if (!is.null(strand)) {
        nsi <- which(strand == "-")
        ml[nsi] <- rm
        mr[nsi] <- lm
        sstrand <- ifelse(strand == "-", -1, 1)
        sstrand[is.na(sstrand)] <- 1
    } else {
        sstrand <- rep(1, length(pos))
    }

    
    spi <- spp::points_within(x, pos - ml, pos + mr, 
        return.list = TRUE)


    #spi <- spp::points_withinFunction(x, pos - ml, pos + mr, 
        #return.list = TRUE)

    lspi <- unlist(lapply(spi, length))
    #df <- data.frame(i = rep(1:length(x), lspi), si = unlist(spi), 
    #    nu = rep(lspi, lspi))
    df <- data.frame(i = rep(seq_along(x), lspi), si = unlist(spi), 
        nu = rep(lspi, lspi))
    df$rx <- (x[df$i] - pos[df$si]) * sstrand[df$si]
    if (!is.null(nbins)) {
        if (is.null(strand)) {
            breaks <- seq(-lm, rm, by = (lm + rm)/nbins)
        } else {
            breaks <- seq(-max(lm, rm), max(lm, rm), by = (lm + rm)/nbins)
        }
        rn <- breaks[-1] - diff(breaks)/2
        breaks[1] <- breaks[1] - 1
        breaks[length(breaks)] <- breaks[length(breaks)] + 1
        df$bin <- cut(df$rx, breaks, labels = rn)
    }
    return(df)
}

########################################################################
# END OF PETER K.'s FUNCTIONS FOR BIN AVERAGES + SCALING OF METAGENES ##

#' @keywords internal 
## helper function for global settings
f_metaGeneDefinition <- function(selection = "Settings") 
{
    predefnum_bins <- 301  ## 151
    
    ## TWO point scaling
    inner_margin <- 500  ## 500
    right_outer_margin <- 1010  ## 1020
    left_outer_margin <- 2010  ## 2020
    gene_body <- 2000  ## 2000
    
    
    totalGeneLength <- gene_body + 2 * inner_margin
    inner_margin2 <- totalGeneLength - inner_margin
    
    break_points_2P <- c(-2000, 0, inner_margin, 
        gene_body + inner_margin, totalGeneLength, 
        totalGeneLength + 1000)
    
    total_output_gene_length <- sum(left_outer_margin, 
        inner_margin, gene_body, inner_margin, 
        right_outer_margin)
    
    estimated_bin_size_2P <- total_output_gene_length/predefnum_bins
    
    
    ## ONE point scaling
    up_downStream <- 4020
    predefnum_bins_1P <- 201  ###binsize=20
    break_points <- c(-2000, -1000, -500, 0, 500, 1000, 2000)
    ### %%% need to use this for one point scaling estimated bin size....
    estimated_bin_size_1P <- up_downStream/predefnum_bins_1P
    
    if (selection == "Settings") {
        settings <- NULL
        settings$break_points_2P <- break_points_2P
        settings$estimated_bin_size_2P <- estimated_bin_size_2P
        settings$break_points <- break_points
        settings$estimated_bin_size_1P <- estimated_bin_size_1P
        settings$predefnum_bins <- predefnum_bins
        settings$inner_margin <- inner_margin
        settings$right_outer_margin <- right_outer_margin
        settings$left_outer_margin <- left_outer_margin
        settings$gene_body <- gene_body
        settings$up_downStream <- up_downStream
        settings$predefnum_bins_1P <- predefnum_bins_1P
        return(settings)
    }
    data(classesDefList, package = "ChIC.data", envir = environment())

    if (selection == "Hlist") {
        ## GLOBAL VARIABLES
        #Hlist <- c("H3K36me3", "POLR2A", "H3K4me3", "H3K79me2", "H4K20me1",
        #    "H2AFZ", "H3K27me3", "H3K9me3", "H3K27ac", "POLR2AphosphoS5", 
        #    "H3K9ac", "H3K4me2", "H3K9me1", "H3K4me1", "POLR2AphosphoS2", 
        #    "H3K79me1", "H3K4ac", "H3K14ac", "H2BK5ac", "H2BK120ac", 
        #    "H2BK15ac", "H4K91ac", "H4K8ac", "H3K18ac", "H2BK12ac", 
        #    "H3K56ac", "H3K23ac", "H2AK5ac", "H2BK20ac", "H4K5ac", 
        #    "H4K12ac", "H2A.Z", "H3K23me2", "H2AK9ac", "H3T11ph")
        ##'POLR3G'
        Hlist <- classesDefList$Hlist
        return(Hlist)
    }

    if (selection == "TFlist") {
        ## GLOBAL VARIABLES
        #Hlist <- c("H3K36me3", "POLR2A", "H3K4me3", "H3K79me2", "H4K20me1",
        #    "H2AFZ", "H3K27me3", "H3K9me3", "H3K27ac", "POLR2AphosphoS5", 
        #    "H3K9ac", "H3K4me2", "H3K9me1", "H3K4me1", "POLR2AphosphoS2", 
        #    "H3K79me1", "H3K4ac", "H3K14ac", "H2BK5ac", "H2BK120ac", 
        #    "H2BK15ac", "H4K91ac", "H4K8ac", "H3K18ac", "H2BK12ac", 
        #    "H3K56ac", "H3K23ac", "H2AK5ac", "H2BK20ac", "H4K5ac", 
        #    "H4K12ac", "H2A.Z", "H3K23me2", "H2AK9ac", "H3T11ph")
        ##'POLR3G'
        TFlist <- classesDefList$TF
        return(TFlist)
    }

    if (selection == "Classes") {
        ## helper functions to define chromating
        allSharp <- classesDefList$allSharp
        # c("H3K27ac", "H3K9ac", "H3K14ac", "H2BK5ac", "H4K91ac", 
        #    "H3K18ac", "H3K23ac", "H2AK9ac", "H3K4me3", "H3K4me2", "H3K79me1", 
        #    "H2AFZ", "H2A.Z", "H4K12ac", "H4K8ac", "H3K4ac", "H2BK12ac", 
        #    "H4K5ac", "H2BK20ac", "H2BK120ac", "H2AK5ac", "H2BK15ac")
        allBroad <- classesDefList$allBroad
        #c("H3K23me2", "H3K9me2", "H3K9me3", "H3K27me3", "H4K20me1",
        #    "H3K36me3", "H3K56ac", "H3K9me1", "H3K79me2", "H3K4me1", "H3T11ph")
        ## USE ONLY POLR2A for Pol2 class
        RNAPol2 <- classesDefList$RNAPol2
        final <- list(allSharp = allSharp, 
            allBroad = allBroad, 
            RNAPol2 = RNAPol2)
        return(final)
    }
}

#' @keywords internal 
## helper function to check if annotationID is valid
f_annotationCheck <- function(annotationID)
{
    checkMe <- ((annotationID == "hg19") | 
        (annotationID == "mm9")| (annotationID == "dm3")
        | (annotationID == "mm10")| (annotationID == "hg38"))
    if (is.character(annotationID) & checkMe)
    {
            message("\n",annotationID, " valid annotation...")
    }else{
        warning("annotationID not valid. Setting it back to default value 
            (hg19). Currently supported annotations are hg19, 
            hg38, mm9 and mm10.")
        annotationID <- "hg19"
    }
    return(annotationID)        
}

#' @keywords internal 
## helper function to check if annotationID is valid
f_annotationLoad <- function(annotationID)
{
    message("Load gene annotation")
    ## require(ChIC.data)
    if (annotationID == "hg19") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("hg19_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- hg19_refseq_genes_filtered_granges
    } else if (annotationID == "hg38") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("hg38_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- hg38_refseq_genes_filtered_granges
    } else if (annotationID == "mm9") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("mm9_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- mm9_refseq_genes_filtered_granges
    } else if (annotationID == "mm10") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("mm10_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- mm10_refseq_genes_filtered_granges
    } else if (annotationID == "dm3") {
        # hg19_refseq_genes_filtered_granges=NULL
        data("dm3_refseq_genes_filtered_granges", 
            package = "ChIC.data", envir = environment())
        annotObject <- dm3_refseq_genes_filtered_granges
    } else {
        stop(paste("Annotations for", annotationID, "currently not supported"))
    }

    return(annotObject)        
}


#' @keywords internal 
## helper function to check if annotationID is valid
f_chromInfoLoad <- function(annotationID)
{
    ## load chrom_info
    ##message("load chrom_info")
    if (annotationID == "hg19") {
        # hg19_chrom_info=NULL
        data("hg19_chrom_info", package = "ChIC.data", envir = environment())
        # this is keeping only chromsomes 1-22 (excluding X/Y and other non nuclear (Mitocondrial genome) or not in the  main assembly
        chromInfo <- hg19_chrom_info[paste("chr", c(seq_len(22)), sep = "")]
    } else if (annotationID == "hg38") {
        # hg19_chrom_info=NULL
        data("hg38_chrom_info", package = "ChIC.data", envir = environment())
        # this is keeping only chromsomes 1-22 (excluding X/Y and other non nuclear (Mitocondrial genome) or not in the  main assembly
        chromInfo <- hg38_chrom_info[paste("chr", c(seq_len(22)), sep = "")]
    } else if (annotationID == "mm9") {
        # hg19_chrom_info=NULL
        data("mm9_chrom_info", package = "ChIC.data", envir = environment())
        # this is keeping only chromsomes 1-19 (excluding X/Y and other non nuclear (Mitocondrial genome) or not in the  main assembly
        chromInfo <- mm9_chrom_info[paste("chr", c(seq_len(19)), sep = "")]
    } else if (annotationID == "mm10") {
        # hg19_chrom_info=NULL
        data("mm10_chrom_info", package = "ChIC.data", envir = environment())
        # this is keeping only chromsomes 1-19 (excluding X/Y and other non nuclear (Mitocondrial genome) or not in the  main assembly
        chromInfo <- mm10_chrom_info[paste("chr", c(seq_len(19)), sep = "")]
    } else if (annotationID == "dm3") {
        # hg19_chrom_info=NULL
        data("dm3_chrom_info", package = "ChIC.data", envir = environment())
        chromInfo <- dm3_chrom_info[c("chr2L","chr2R","chr3L","chr3R","chr4")]
    } else {
        stop(paste("Annotations for", annotationID, "currently not supported"))
    }

    return(chromInfo)        
}


#' @keywords internal 
##helper function to check if chromosome names contain "chr"
f_checkOfChrNames = function( data )
{
    checkOfChr <- (grep( "chr", names( data$tags )))
    if ( length(checkOfChr) < 1L)
    {
        newnames <- paste("chr", names(data$tags), sep="")
        names(data$tags) <- newnames
        names(data$quality) <- newnames
    }
    return(data)
}

#' @keywords internal 
## masked_t.get.gene.av.density(smoothed.densityInput,
## gdl=annotatedGenesPerChr,mc=mc)
masked_t.get.gene.av.density <- function(chipTags_current, gdl, mc = 1) 
{
    settings <- NULL
    settings <- f_metaGeneDefinition(selection = "Settings")
    message("\nloading metaGene settings")
    ## print(str(settings))
    result <- f_t.get.gene.av.density(chipTags_current, 
        gdl = gdl, im = settings$inner_margin, 
        lom = settings$left_outer_margin, 
        rom = settings$right_outer_margin, bs = settings$gene_body, 
        nbins = settings$predefnum_bins, mc = mc)
    return(result)
}

## binnedInput_TSS <- masked_getGeneAvDensity_TES_TSS
## (smoothed.densityInput,gdl=annotatedGenesPerChr,mc=mc,tag='TSS')
#' @keywords internal 
masked_getGeneAvDensity_TES_TSS <- function(smoothed.density, gdl, 
    tag = "TSS", mc = 1) 
{
    settings <- NULL
    settings <- f_metaGeneDefinition(selection = "Settings")
    message("\nloading metaGene settings")
    ## print(str(settings))
    if (tag == "TSS") {
        result <- f_t.get.gene.av.density_TSS(tl_current = smoothed.density, 
            gdl = gdl, 
            m = (settings$up_downStream / 2), 
            nbins = settings$predefnum_bins_1P, 
            separate.strands = FALSE, 
            mc = mc)
        message("TSS")
    } else {
        result <- f_t.get.gene.av.density_TES(tl_current = smoothed.density, 
            gdl = gdl, 
            m = (settings$up_downStream / 2), 
            nbins = settings$predefnum_bins_1P, 
            separate.strands = FALSE, 
            mc = mc)
        message("TES")
    }
    return(result)
}


#' @keywords internal 
## MODIFIED VERSION OF 
## t.get.gene.av.density for TWO.POINT a function for getting
## bin-averaged profiles for individual genes 
## TWO.POINT t.get.gene.av.density <-
## function(chipTags_current,gdl=annotatedGenesPerChr, im=500,lom=5e3, rom=1e3,
## bs=2e3, nbins=400,separate.strands=F) {
## min.feature.size=min.feature.sizeMy=3000 f_t.get.gene.av.density <-
## function(chipTags_current, gdl=annotatedGenesPerChr, im=inner_margin,
## lom=left_outer_margin, rom=right_outer_margin, bs=gene_body,
## nbins=predefnum_bins,separate.strands=F, min.feature.size=3000,mc=1) {
## f_t.get.gene.av.density <- function(chipTags_current,gdl, im=500, lom=2010,
## rom=1010, bs=2000, nbins=301,separate.strands=F, min.feature.size=3000,mc=1)
f_t.get.gene.av.density <- function(chipTags_current, gdl, im, lom, 
    rom, bs, nbins, separate.strands = FALSE, min.feature.size = 3000, mc = 1) 
{
    chrl <- names(gdl)
    names(chrl) <- chrl
    ## lapply(chrl[chrl %in% names(chipTags_current$td)],function(chr) {
    ## BiocParallel::bplapply(chrl[chrl %in% names(chipTags_current$td)], 
    ## BPPARAM = BiocParallel::MulticoreParam(workers = mc), 
    #mclapply(chrl[chrl %in% names(chipTags_current$td)], 
    mclapply(chrl[chrl %in% names(chipTags_current)], 
        mc.preschedule = FALSE, 
        mc.cores = mc, 
        FUN = function(chr) {
            # print(chr)
            nsi <- gdl[[chr]]$strand == "-"
            # print(length(nsi))
            current_gene_names <- gdl[[chr]]$geneSymbol
            if ((sum(!nsi) > 0)) {
                ## if ((sum(!nsi)>1)) {
                #px <- f_feature.bin.averages(chipTags_current$td[[chr]], 
                px <- f_feature.bin.averages(chipTags_current[[chr]], 
                    data.frame(s = gdl[[chr]]$txStart[!nsi], 
                        e = gdl[[chr]]$txEnd[!nsi]), 
                    lom = lom, rom = rom, im = im, bs = bs, 
                    nbins = nbins, min.feature.size = min.feature.size, 
                    nu.point.omit = FALSE)
                rownames(px) <- current_gene_names[!nsi]
            } else {
                px <- NULL
            }
            if ((sum(nsi) > 0)) {
                ## if ((sum(nsi)>1)) {
                #nd <- chipTags_current$td[[chr]]
                nd <- chipTags_current[[chr]]
                nd$x <- -1 * nd$x
                nx <- f_feature.bin.averages(nd, 
                    data.frame(s = -1 * gdl[[chr]]$txEnd[nsi], 
                        e = -1 * gdl[[chr]]$txStart[nsi]), 
                    lom = lom, rom = rom, im = im, 
                    bs = bs, nbins = nbins, 
                    min.feature.size = min.feature.size, 
                    nu.point.omit = FALSE)
                rownames(nx) <- current_gene_names[nsi]
            } else {
                nx <- NULL
            }
            
            if (separate.strands) {
                return(p = px, n = nx)
            } else {
                return(rbind(px, nx))
            }
        })
}

#' @keywords internal 
## MODIFIED VERSION OF t.get.gene.av.density FOR TSS a function for getting
## bin-averaged profiles for individual genes TSS ONE.POINT
## f_t.get.gene.av.density_TSS <- function(tl_current,
## gdl=annotatedGenesPerChr,m=up_downStream, nbins=predefnum_bins_1P,
## separate.strands=F,mc=1) {###binning of the frame
f_t.get.gene.av.density_TSS <- function(tl_current, gdl, m = 4020, 
    nbins = 201, separate.strands = FALSE, mc = 1) 
{
    ## binning of the frame
    chrl <- names(gdl)
    names(chrl) <- chrl
    ## lapply(chrl[chrl %in% names(tl_current$td)],function(chr) 
    ## BiocParallel::bplapply(chrl[chrl %in% names(tl_current$td)], 
    ## BPPARAM = BiocParallel::MulticoreParam(workers = mc), 
    
    mclapply(chrl[chrl %in% names(tl_current)], 
        mc.preschedule = FALSE, 
        mc.cores = mc,     
        FUN = function(chr) {
            nsi <- gdl[[chr]]$strand == "-"
            current_gene_names <- gdl[[chr]]$geneSymbol
            
            ## if ((sum(!nsi)>0)) {
            if ((sum(!nsi) > 0)) {
                ## px <- f_feature.bin.averages(tl_current$td[[chr]],
                ## data.frame(x=gdl[[chr]]$txStart[!nsi]),m=m, 
                ## nbins=nbins,nu.point.omit=FALSE)
                px <- f_feature.bin.averages(tl_current[[chr]], 
                    data.frame(x = gdl[[chr]]$txStart[!nsi]), 
                m = m,
                min.feature.size=NULL,
                nbins = nbins, nu.point.omit = FALSE)
                rownames(px) <- current_gene_names[!nsi]
            } else {
                px <- NULL
            }
            
            if ((sum(nsi) > 0)) {
                ## if ((sum(nsi)>0)) {
                nd <- tl_current[[chr]]
                nd$x <- -1 * nd$x
                ## nx <- f_feature.bin.averages(nd,data.frame
                ## (x=-1*gdl[[chr]]$txEnd[nsi]),m=m,nbins=nbins,
                ## nu.point.omit=FALSE)
                nx <- f_feature.bin.averages(nd, 
                    data.frame(x = -1 * gdl[[chr]]$txEnd[nsi]), 
                m = m, 
                min.feature.size=NULL,
                nbins = nbins, nu.point.omit = FALSE)
                rownames(nx) <- current_gene_names[nsi]
            } else {
                nx <- NULL
            }
            
            if (separate.strands) {
                return(p = px, n = nx)
            } else {
                return(rbind(px, nx))
            }
        })
}

#' @keywords internal 
## MODIFIED VERSION OF t.get.gene.av.density FOR TES a function for getting
## bin-averaged profiles for individual genes TES ONE.POINT
## f_t.get.gene.av.density_TES <- function(tl_current,gdl=annotatedGenesPerChr,
## m=up_downStream, nbins=predefnum_bins_1P,separate.strands=F,mc=1) {
f_t.get.gene.av.density_TES <- function(tl_current, gdl, m = 4020, 
    nbins = 201, separate.strands = FALSE, mc = 1) 
{
    
    chrl <- names(gdl)
    names(chrl) <- chrl
    ## lapply(chrl[chrl %in% names(tl_current$td)],function(chr) 
    ## BiocParallel::bplapply(chrl[chrl %in% names(tl_current$td)], 
    ## BPPARAM = BiocParallel::MulticoreParam(workers = mc), 
    
    mclapply(chrl[chrl %in% names(tl_current)], 
        mc.preschedule = FALSE, 
        mc.cores = mc,
        FUN = function(chr) {
            nsi <- gdl[[chr]]$strand == "-"
            current_gene_names <- gdl[[chr]]$geneSymbol
            
            ## if ((sum(!nsi)>0)) {
            if ((sum(!nsi) > 0)) {
                ## f_feature.bin.averages <- 
                ##function(dat,feat,nu.feat.omit=F, nu.point.omit=T,
                ## scaling=NULL, return.scaling=F, trim=0,
                ## min.feature.size=NULL, ... ) {
                px <- f_feature.bin.averages(tl_current[[chr]], 
                    data.frame(x = gdl[[chr]]$txEnd[!nsi]), 
                m = m, 
                min.feature.size=NULL,
                nbins = nbins, nu.point.omit = FALSE)
                rownames(px) <- current_gene_names[!nsi]
            } else {
                px <- NULL
            }
            ## if ((sum(nsi)>1)) {
            if ((sum(nsi) > 0)) {
                nd <- tl_current[[chr]]
                nd$x <- -1 * nd$x
                nx <- f_feature.bin.averages(nd, 
                    data.frame(x = -1 * gdl[[chr]]$txStart[nsi]), 
                m = m, 
                min.feature.size=NULL,
                nbins = nbins, nu.point.omit = FALSE)
                rownames(nx) <- current_gene_names[nsi]
            } else {
                nx <- NULL
            }
            
            if (separate.strands) {
                return(p = px, n = nx)
            } else {
                return(rbind(px, nx))
            }
        })
}


#' @keywords internal 
f_spotfunction <- function(dframe, breaks, tag) 
{
    # hotSpots=NULL###takes values at different predefined points
    fframe <- dframe[rownames(dframe) %in% breaks, ]
    rownames(fframe) <- c(paste("hotSpots", tag, breaks, sep = "_"))
    return(fframe)
}

#' @keywords internal 
f_spotfunctionNorm <- function(dframe, breaks, tag) 
{
    newframe <- data.frame(dframe)
    Norm <- newframe[rownames(newframe) %in% breaks, ]
    fframe <- data.frame(breaks, Norm)
    rownames(fframe) <- c(paste("hotSpots", tag, breaks, sep = "_"))
    fframe$breaks <- NULL
    
    return(fframe)
}


#' @keywords internal 
f_maximaAucfunction <- function(dframe, breaks, estBinSize, tag) 
{
    ## local maxima and area in all the predefined regions
    ## settings=f_metaGeneDefinition(selection='Settings')
    ## estimated_bin_size_2P=settings$estimated_bin_size_2P
    dframe <- as.data.frame(dframe)
    
    localMaxima_auc <- lapply(seq(length(breaks) - 1), FUN = function(i) {
        xi <- breaks[i]
        xj <- breaks[i + 1]
        sub_set <- dframe[which((as.numeric(rownames(dframe)) <= xj) & 
            (as.numeric(rownames(dframe)) >= xi)), ]
        l1 <- max(sub_set$Chip)
        l2 <- max(sub_set$Input)
        chip_frame <- cbind(i, "chip", 
            rownames(sub_set[which(sub_set$Chip == l1), ]), l1, 
            sum((sub_set$Chip) * estBinSize))
        input_frame <- cbind("input", 
            rownames(sub_set[which(sub_set$Input == l2), ]), l2, 
            sum((sub_set$Input) * estBinSize))
        tt <- (cbind(chip_frame, input_frame))
        
        return(tt)
    })
    
    if (tag == "geneBody") {
        fframe <- data.frame(matrix(unlist(localMaxima_auc), 
            nrow = 5, byrow = TRUE))
    } else {
        fframe <- data.frame(matrix(unlist(localMaxima_auc), 
            nrow = 6, byrow = TRUE))
    }
    colnames(fframe) <- c("i", "chip", "chipX", "chipY", "chipAUC", 
        "input", "inputX", "inputY", "inputAUC")
    
    ## save x-coordiante for max value for Chip and Input
    xFrame <- data.frame(as.numeric(as.character(fframe$chipX)), 
        as.numeric(as.character(fframe$inputX)))
    colnames(xFrame) <- c("Chip", "Input")
    rownames(xFrame) <- paste("localMax", tag, fframe$i, "x", sep = "_")
    
    ## save y-coordiante for max value for Chip and Input
    yFrame <- data.frame(as.numeric(as.character(fframe$chipY)), 
        as.numeric(as.character(fframe$inputY)))
    colnames(yFrame) <- c("Chip", "Input")
    rownames(yFrame) <- paste("localMax", tag, fframe$i, "y", sep = "_")
    
    ## save auc for Chip and Input
    aucFrame <- data.frame(as.numeric(as.character(fframe$chipAUC)), 
        as.numeric(as.character(fframe$inputAUC)))
    colnames(aucFrame) <- c("Chip", "Input")
    rownames(aucFrame) <- paste("auc", tag, fframe$i, sep = "_")
    
    
    finalReturn <- NULL
    finalReturn <- data.frame(rbind(xFrame, yFrame, aucFrame))
    
    return(finalReturn)
}


#' @keywords internal 
f_maximaAucfunctionNorm <- function(dframe, breaks, estBinSize, tag) 
{
    
    newframe <- as.data.frame(dframe)
    newframe$x <- rownames(newframe)
    colnames(newframe) <- c("Norm", "break")
    
    localMaxima_auc <- lapply(seq(length(breaks) - 1), FUN = function(i) {
        # print(i)
        xi <- breaks[i]
        xj <- breaks[i + 1]
        sub_set <- newframe[which((as.numeric(rownames(newframe)) <= xj) & 
            (as.numeric(rownames(newframe)) >= xi)), ]
        l1 <- max(sub_set$Norm)
        
        norm_frame <- cbind(i, "norm", 
            rownames(sub_set[which(sub_set$Norm == l1), ]), 
            l1, sum((sub_set$Norm) * estBinSize))
        return(norm_frame)
    })
    
    if (tag == "geneBody") {
        fframe <- data.frame(matrix(unlist(localMaxima_auc), 
            nrow = 5, byrow = TRUE))
    } else {
        fframe <- data.frame(matrix(unlist(localMaxima_auc), 
            nrow = 6, byrow = TRUE))
    }
    colnames(fframe) <- c("i", "norm", "X", "Y", "AUC")
    
    ## save x-coordiante for max value for Chip and Input
    xFrame <- data.frame(as.numeric(as.character(fframe$X)))
    colnames(xFrame) <- c("Norm")
    rownames(xFrame) <- paste("localMax", tag, fframe$i, "x", sep = "_")
    
    ## save y-coordiante for max value for Chip and Input
    yFrame <- data.frame(as.numeric(as.character(fframe$Y)))
    colnames(yFrame) <- c("Norm")
    rownames(yFrame) <- paste("localMax", tag, fframe$i, "y", sep = "_")
    
    ## save auc for Chip and Input
    aucFrame <- data.frame(as.numeric(as.character(fframe$AUC)))
    colnames(aucFrame) <- c("Norm")
    rownames(aucFrame) <- paste("auc", tag, fframe$i, sep = "_")
    
    
    finalReturn <- NULL
    finalReturn <- data.frame(rbind(xFrame, yFrame, aucFrame))
    return(finalReturn)
}

#' @keywords internal
## calculating variance, sd and qartiles of the value ditribution in different
## intervals
f_variabilityValues <- function(dframe, breaks, tag) {
    #variabilityValues <- lapply(breaks[1:3], FUN = function(start) {
    variabilityValues <- lapply(breaks[seq_len(3)], FUN = function(start) {
        #print(start)
        end <- abs(start)
        sub_set <- dframe[which((as.numeric(rownames(dframe)) <= end) & 
            (as.numeric(rownames(dframe)) >= start)), ]
        sub_set <- data.frame(sub_set)
        name <- paste("dispersion", tag, start, "variance", sep = "_")
        varI <- var(sub_set$Input)
        varC <- var(sub_set$Chip)
        name <- c(name, paste("dispersion", tag, start, "sd", sep = "_"))
        sdI <- sd(sub_set$Input)
        sdC <- sd(sub_set$Chip)
        #valueFrameI <- data.frame(quantile(sub_set$Input)[1:4])
        valueFrameI <- data.frame(quantile(sub_set$Input)[seq_len(4)])
        colnames(valueFrameI) <- c("value")
        #valueFrameC <- data.frame(quantile(sub_set$Chip)[1:4])
        valueFrameC <- data.frame(quantile(sub_set$Chip)[seq_len(4)])
        colnames(valueFrameC) <- c("value")
        name <- c(name, paste("dispersion", tag, start, 
            rownames(valueFrameI), sep = "_"))
        back <- data.frame(cbind(c(varI, sdI, valueFrameI$value), 
            c(varC, sdC, valueFrameC$value)))
        
        colnames(back) <- c("Input", "Chip")
        rownames(back) <- c(name)
        return(back)
    })
    return(variabilityValues)
}


#' @keywords internal 
# calculating variance, sd and qartiles of the value ditribution in different
# intervals
f_variabilityValuesNorm <- function(dframe, breaks, tag) {
    
    newframe <- as.data.frame(dframe)
    newframe$x <- rownames(newframe)
    colnames(newframe) <- c("Norm", "break")
    #variabilityValues <- lapply(breaks[1:3], FUN = function(start) {
    variabilityValues <- lapply(breaks[seq_len(3)], FUN = function(start) {
        #print(start)
        end <- abs(start)
        sub_set <- newframe[which((as.numeric(rownames(newframe)) <= end) & 
            (as.numeric(rownames(newframe)) >= start)), ]
        sub_set <- data.frame(sub_set)
        name <- paste("dispersion", tag, start, "variance", sep = "_")
        var <- var(sub_set$Norm)
        name <- c(name, paste("dispersion", tag, start, "sd", sep = "_"))
        sd <- sd(sub_set$Norm)
        
        #valueFrame <- data.frame(quantile(sub_set$Norm)[1:4])
        valueFrame <- data.frame(quantile(sub_set$Norm)[seq_len(4)])
        colnames(valueFrame) <- c("value")
        name <- c(name, paste("dispersion", tag, start, 
            rownames(valueFrame), sep = "_"))
        
        back <- data.frame(cbind(c(var, sd, valueFrame$value)))
        
        colnames(back) <- c("Norm")
        rownames(back) <- c(name)
        return(back)
    })
    return(variabilityValues)
}


##################################################################### 
######################### FUNCTIONS comparison with compendium
#####################################################################

#' @keywords internal 
## helper function to load profiles from ChIC.data
f_loadDataCompendium <- function(endung, target, tag) 
{
    # compendium_profiles=ChIC.data::compendium_profiles
    if (tag == "geneBody") {
        name <- paste(target, "_", "TWO", endung, sep = "")
    } else {
        name <- paste(target, "_", tag, endung, sep = "")
    }
    #load profiles
    if (target %in% f_metaGeneDefinition("Hlist")){
        data("compendium_profiles", 
            package = "ChIC.data", 
            envir = environment())
        frame=compendium_profiles[[name]]
    }else{
        data("compendium_profiles_TF", 
            package = "ChIC.data", 
            envir = environment())
        frame=compendium_profiles_TF[[name]]
    }
    return(frame)
}

#' @keywords internal 
## helper function to prepare dataframe
f_prepareData <- function(fmean, frame) 
{
    finalframe <- cbind(fmean$x, frame)
    rownames(finalframe) <- NULL
    finalframe <- as.data.frame(finalframe)
    colnames(finalframe) <- c("x", "mean")
    return(finalframe)
}



#' @keywords internal 
## plot profiles compendium versus current dataset
f_plotProfiles <- function(meanFrame, currentFrame, endung = "geneBody", 
    absoluteMinMax, maintitel = "title", ylab = "mean of log2 read density", 
    savePlotPath = NULL) 
{
    message("Load settings")
    settings <- f_metaGeneDefinition(selection = "Settings")
    
    if (!is.null(savePlotPath)) {
        filename <- paste(maintitel, ".pdf", sep = "")
        pdf(file = file.path(savePlotPath, filename), width = 10, height = 7)
    }
    break_points_2P <- settings$break_points_2P
    break_points <- settings$break_points
    ## The standard error of the mean (SEM) is the standard deviation of the
    ## sample-mean's estimate of a population mean.  
    ## (It can also be seen as the standard deviation of the error in the 
    ## sample mean with respect to the true
    ## mean, since the sample mean is an unbiased estimator.)  SEM is usually
    ## estimated by the sample estimate of the population 
    ## standard deviation (sample standard deviation) divided by the 
    ## square root of the sample size (assuming
    ## statistical independence of the values in the sample)
    plot(x = c(min(meanFrame$x), max(meanFrame$x)), 
        y = c(absoluteMinMax[1], absoluteMinMax[2]), 
        type = "n", xlab = "metagene coordinates", ylab = ylab, 
        main = maintitel, 
        xaxt = "n")
    polygon(x = c(meanFrame$x, rev(meanFrame$x)), 
        y = c(meanFrame$mean + 2 * meanFrame$sderr, 
        rev(meanFrame$mean)), col = "lightblue", border = NA)
    polygon(x = c(meanFrame$x, rev(meanFrame$x)), 
        y = c(meanFrame$mean - 2 * meanFrame$sderr, 
        rev(meanFrame$mean)), col = "lightblue", border = NA)
    lines(x = meanFrame$x, y = meanFrame$mean, col = "black", lwd = 2)
    lines(x = currentFrame$x, y = currentFrame$mean, col = "red", lwd = 2)
    if (endung == "geneBody") {
        currBreak_points <- break_points_2P[c(-2, -5)]  
        ##c(-2000,500,2500,4000)
        abline(v = c(break_points_2P[c(2, 5)]), lty = 2, 
            col = "darkgrey", lwd = 3)
        abline(v = currBreak_points, lty = 3, 
            col = "darkgrey", lwd = 2)
        axis(side = 1, at = break_points_2P, 
            labels = c("-2KB", "TSS", "TSS+500", "TES-500", "TES", "+1KB"))
    } else {
        TSSbreak_points <- break_points[-c(4)]  
        # c(-2000,-1000,-500,500,1000,2000)
        if (endung == "TSS") {
            axis(side = 1, at = break_points, 
                labels=c("-2KB","-1KB", "-500", "TSS", "+500", "+1KB", "+2KB"))
        } else {
            axis(side = 1, at = break_points, 
                labels=c("-2KB","-1KB", "-500", "TES", "+500", "+1KB", "+2KB"))
        }
        abline(v = 0, lty = 2, col = "darkgrey", lwd = 2)
        abline(v = TSSbreak_points, lty = 3, col = "darkgrey", lwd = 2)
    }
    legend("topleft", legend = c("mean", "2*stdErr"), 
        fill = c("black", "lightblue"), 
        bg = "white")
    if (!is.null(savePlotPath)) {
        dev.off()
    }
}

#' @keywords internal 
## helper function to get binding class
f_getBindingClass <- function(target) {
    allChrom <- f_metaGeneDefinition("Classes")
    if (target %in% allChrom$allSharp) {
        profileSet <- allChrom$allSharp
        tag <- "(Sharp class)"
    }
    if (target %in% allChrom$allBroad) {
        profileSet <- allChrom$allBroad
        tag <- "(Broad class)"
    }
    if (target %in% allChrom$RNAPol2) {
        profileSet <- allChrom$RNAPol2
        tag <- "(RNAPol2 class)"
    }

    return(list(profileSet = profileSet, tag = tag))
}


#' @keywords internal 
## density plot for QC-value distribution versus a single value
f_plotValueDistribution <- function(compendium, title, coordinateLine, 
    savePlotPath = NULL) 
{
    if (!is.null(savePlotPath)) {
        filename <- file.path(savePlotPath, "PlotValueDistribution.pdf")
        pdf(file = filename, width = 10, height = 7)
    }
    mycolors <- c("lightsteelblue1", "lightsteelblue3")
    maxi <- max(compendium)
    mini <- min(compendium)
    
    ## get density
    d <- density(compendium)
    plot(d, main = title, xlim = c(mini, maxi))
    
    median <- median(compendium)
    qq <- quantile(compendium, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    
    i <- d$x >= (qq[1]) & d$x <= (qq[2])
    lines(d$x, d$y)
    polygon(c(qq[1], d$x[i], qq[2]), c(0, d$y[i], 0), col = mycolors[1])
    
    i <- d$x >= (qq[2]) & d$x <= (qq[3])
    lines(d$x, d$y)
    polygon(c(qq[2], d$x[i], qq[3]), c(0, d$y[i], 0), col = mycolors[2])
    
    i <- d$x >= (qq[3]) & d$x <= (qq[4])
    lines(d$x, d$y)
    polygon(c(qq[3], d$x[i], qq[4]), c(0, d$y[i], 0), col = mycolors[2])
    
    i <- d$x >= (qq[4]) & d$x <= (qq[5])
    lines(d$x, d$y)
    polygon(c(qq[4], d$x[i], qq[5]), c(0, d$y[i], 0), col = mycolors[1])
    
    abline(v = coordinateLine, , col = "red", lwd = 3, lty = 2)
    
    if (!is.null(savePlotPath)) {
        dev.off()
    }
}



#' @keywords internal 
## helper function to select the random forest model for the respective
## chromatinmark or TF
f_getPredictionModel <- function(id) {
    # library(randomForest)
    allChrom <- f_metaGeneDefinition("Classes")
    data("rf_models", package = "ChIC.data", envir = environment())
    
    if (id %in% f_metaGeneDefinition("Hlist")) {
        message("Load chromatinmark model")
        if (id %in% allChrom$allSharp) {
            model <- rf_models[["sharpEncode"]]
        }
        
        if (id %in% allChrom$allBroad) {
            model <- rf_models[["broadEncode"]]
        }
        
        if (id %in% allChrom$RNAPol2) {
            model <- rf_models[["RNAPol2Encode"]]
        }
        
        if (id == "H3K9me3") {
            model <- rf_models[["H3K9Encode"]]
        }
        
        if (id == "H3K27me3") {
            model <- rf_models[["H3K27Encode"]]
        }
        
        if (id == "H3K36me3") {
            model <- rf_models[["H3K36Encode"]]
        }
    } else if ((id %in% f_metaGeneDefinition("TFlist")) | (id== "TF"))
    {
        message("Load TF model")
        model <- rf_models$TFmodel
    } else {
        message(id, "not found")
        model=NULL
    }
    return(model)
}


#' @keywords internal 
## helper function that converts frame with chip and normalized values to one
## column frame to further be procecced.
f_convertframe <- function(oldframe) {
    values <- c(oldframe$Chip, oldframe$Norm)
    newframe <- data.frame(values)
    nn <- c(paste("chip", rownames(oldframe), sep = "_"), 
        paste("norm", rownames(oldframe), 
        sep = "_"))
    rownames(newframe) <- nn
    return(newframe)
}



