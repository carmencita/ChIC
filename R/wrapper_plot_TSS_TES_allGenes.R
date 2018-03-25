#'@title Wrapper function that plots non-scaled profiles for 
#' TSS or TES and to collects the QC-metrics
#'
#'@description 
#' The non-scaled profile is constructed around the TSS/TES,
#' with 2KB up- and downstream regions respectively. Different values 
#' are taken at the TSS/TES and surroundings with +/-2KB, +/-1KB 
#' and +/-500 sizes. For all the genomic 
#' positions, we kept the values for the ChIP and the normalized profile,
#' as the normalization already contains information from the input. 
#' Additionally, we calculated for all of the intervals between 
#' the predefined positions the area under the profile, 
#' the local maxima (x, y coordinates), the variance, 
#' the standard deviation and the quantiles at 0%, 25%, 
#' 50% and 75%. In total the function returns 43 QC-metrics.
#'
#' qualityScores_LM
#'
#' @param data metagene-list for input and chip sample for
#' TSS or TES returned by createMetageneProfile()
#' @param tag String that can be "TSS" or "TES",indicating if the TSS or 
#' the TES profile should be calcualted (Default="TSS")
#' @param savePlotPath if set the plot will be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout. 
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @export
#'
#' @return result Dataframe with QC-values for chip, input and 
#' normalized metagene profile
#'
#'@examples
#' data(TSSProfile)
#' TSS_Plot=qualityScores_LM(data=TSSProfile, tag="TSS")
#'


qualityScores_LM<-function(data, tag, savePlotPath=NULL, debug=FALSE)
{    
    stopifnot(tag %in% c("TES","TSS"))
    stopifnot(length(data) == 2L)

    binnedChip=data$chip
    binnedInput=data$chip
    message("load metagene setting")
    #load("Settings.RData")
    settings=f_metaGeneDefinition(selection="Settings")
    break_points=settings$break_points
    estimated_bin_size_1P=settings$estimated_bin_size_1P

    psc <- 1;## pseudocount## required to avoid log2 of 0
    chip <-log2(do.call(rbind,binnedChip)+psc)
    input<-log2(do.call(rbind,binnedInput)+psc)


    input.noNorm<-colMeans(input,na.rm=TRUE)
    chip.noNorm <- colMeans(chip,na.rm=TRUE)
    all.noNorm=NULL
    all.noNorm<-cbind(chip.noNorm, input.noNorm)
    colnames(all.noNorm)<-c("Chip","Input")
    

    ##values at specific predefined points
    hotSpotsValues=f_spotfunction(all.noNorm, break_points, 
        tag=tag)
    ##local maxima and area in all the predefined regions
    maxAucValues=f_maximaAucfunction(all.noNorm, break_points, 
        estimated_bin_size_1P, tag=tag)

    ## chip_dispersion_TES_-1000_0% 
    ## chip_dispersion_TES_-1000_25% 
    ## chip_dispersion_TES_-1000_50% 
    ## chip_dispersion_TES_-1000_75% 
    ## chip_dispersion_TES_-1000_sd 
    ## chip_dispersion_TES_-1000_variance 

    
    variabilityValues=f_variabilityValues(all.noNorm, break_points, 
        tag=tag)

    ##make plots
    colori<-c(rev(rainbow(ncol(all.noNorm)-1)), "black")
    if (!is.null(savePlotPath))
    {
        filename=file.path(savePlotPath,paste("ChIP_Input_",tag,".pdf",sep=""))
        pdf(file=filename,width=10, height=7)
    }    
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    matplot(x=as.numeric(rownames(all.noNorm)), 
        y=all.noNorm, type="l", lwd=2, lty=1,
        col=colori,xlab="metagene coordinates",
        ylab="mean of log2 tag density",main=tag,xaxt='n')


    newbreak_points=break_points[-c(4)]# c(-2000,-1000,-500,500,1000,2000)
    abline(v=0,lty=2,col="darkgrey", lwd=2)
    abline(v=newbreak_points,lty=3,col="darkgrey", lwd=2)
    ##abline(v=0,lty=2,col="darkgrey", lwd=2)###plot TSS        
    ##plotPoints=c(-2000,-1000,-500,500,1000,2000)###plot remaining be
    ##abline(v=plotPoints,lty=3,col="darkgrey", lwd=2)
    axis(side = 1, at = break_points, 
        labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
    legend(x="topleft", fill=colori, legend=colnames(all.noNorm), 
        bg="white",cex=0.8)

    if (!is.null(savePlotPath))
    {
        dev.off()
        message("pdf saved under ",filename)
    }


    ##normalized plot and values
    common_genes<-rownames(input)[rownames(input) %in% rownames(chip)]

    all.Norm<-colMeans(t
        (t(chip[common_genes,])-t(input[common_genes,])),na.rm=TRUE)

    hotSpotsValuesNorm=f_spotfunctionNorm(all.Norm, break_points, 
        tag=tag)

    maxAucValuesNorm=f_maximaAucfunctionNorm(all.Norm, break_points, 
        estimated_bin_size_1P, tag=tag)

    variabilityValuesNorm=f_variabilityValuesNorm(all.noNorm, break_points, 
        tag=tag)


    if (!is.null(savePlotPath))
    {
        filename=file.path(savePlotPath,paste("Normalized_",tag,".pdf",sep=""))
        pdf(file=filename,width=10, height=7)
    }
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    plot(x=as.numeric(names(all.Norm)),y=all.Norm, type="l", lwd=2, lty=1, 
        col="orange",xlab="metagene coordinates",
        ylab="mean log2 enrichment (signal/input)",
        main=paste("normalized",tag,sep=" "), xaxt='n')

    abline(v=0,lty=2,col="darkgrey", lwd=2)
    abline(v=newbreak_points,lty=3,col="darkgrey", lwd=2)
    axis(side = 1, at = break_points, 
        labels = c("-2KB","-1KB","-500",tag,"500","+1KB","+2KB"))
    legend(x="topleft", fill=colori, legend=colnames(all.noNorm), 
        bg="white", cex=0.8)
    
    if (!is.null(savePlotPath))
    {
        dev.off()
        message("pdf saved under ",filename)
    }

    #convert values and features to one frame
    p1=rbind(hotSpotsValues,maxAucValues)
    p2=rbind(variabilityValues[[1]], 
        variabilityValues[[2]],
        variabilityValues[[3]])
    p3=rbind(hotSpotsValuesNorm, maxAucValuesNorm)
    p4=rbind(variabilityValuesNorm[[1]], 
        variabilityValuesNorm[[2]],
        variabilityValuesNorm[[3]])
    result=data.frame(cbind(rbind(p1,p2),rbind(p3,p4)))

    if (debug)
    {
        message("Debugging mode ON")
        outname=file.path(getwd(), paste(tag,"onepoints.result",sep="_"))
        file.remove(outname)
        write.table(result,file=outname,row.names = FALSE, 
            col.names=FALSE,append=TRUE, quote = FALSE)
    }
    return(result)
}