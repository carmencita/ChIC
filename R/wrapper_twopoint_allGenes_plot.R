#'@title Wrapper function to plot the scaled metagene- profile  
#' and to collect the QC-metrics
#'
#'@description The scaled metagene profile that includes the gene body, 
#' the signal is captured on a real scale from the TSS and an upstream 
#' region of 2KB. From the TSS, the gene body is constructed with 0.5KB 
#' in real scale at the gene start (TSS + 0.5KB) and the gene end 
#' (TES - 0.5KB), whereas the remaining gene body is 
#' scaled to a virtual length of 2000. Considering the length 
#' of these regions, the minimum gene length required is 3KB and 
#' shorter genes are filtered out. From the profile, we take enrichment 
#' values at different coordinates: at 
#' -2KB, at the TSS, inner margin (0.5KB), gene body 
#' (2KB + 2 * inner margin), 
#' gene body+1KB. We collect in total 42 QC-metrics from the ChIP and 
#' normalized profile. 
#'
#' qualityScores_LMgenebody
#'
#' @param data metagene-list for input and chip sample 
#' of the genebody profile returned by createMetageneProfile()
#' @param savePlotPath if set the plot will be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout. 
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @export
#'
#' @return returnList
#'
#'@examples
#' print("Example")
#' data(geneBodyProfile)
#' geneBody_Plot=qualityScores_LMgenebody(geneBodyProfile)
#'

qualityScores_LMgenebody<-function(data, savePlotPath=NULL, debug=FALSE)
{
    stopifnot(length(data) == 2L)

    binnedChip=data$chip
    binnedInput=data$input
    message("load metagene setting")
    settings=f_metaGeneDefinition(selection="Settings")
    psc <- 1;## pseudocount## required to avoid log2 of 0
    break_points_2P=settings$break_points_2P
    estimated_bin_size_2P=settings$estimated_bin_size_2P
    chip <- log2(do.call(rbind,binnedChip)+psc)
    input<-log2(do.call(rbind,binnedInput)+psc)    

    input.noNorm<-colMeans(input,na.rm=TRUE)
    chip.noNorm <- colMeans(chip,na.rm=TRUE)
    all.noNorm<-cbind(chip.noNorm, input.noNorm)
    colnames(all.noNorm)<-c("Chip","Input")


    ##values at specific predefined points
    hotSpotsValues=f_spotfunction(all.noNorm, break_points_2P, 
        tag="geneBody")
    ##local maxima and area in all the predefined regions
    maxAucValues=f_maximaAucfunction(all.noNorm, breaks=break_points_2P, 
        estBinSize=estimated_bin_size_2P, tag="geneBody")

    ##make plots
    colori<-c(rev(rainbow(ncol(all.noNorm)-1)), "black")
    if (!is.null(savePlotPath))
    {
        filename=file.path(savePlotPath,"ScaledMetaGene_ChIP_Input.pdf")
        pdf(file=filename,width=10, height=7)
    }    
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    matplot(x=as.numeric(rownames(all.noNorm)),y=all.noNorm, type="l", 
        lwd=2, lty=1, col=colori, xlab="metagene coordinates", 
        ylab="mean of log2 tag density", main="metagene", xaxt='n')

    currBreak_points=break_points_2P[c(-2,-5)]
    abline(v=c( break_points_2P[c(2,5)]),lty=2,col="darkgrey", lwd=3)
    abline(v=currBreak_points,lty=3,col="darkgrey", lwd=2)
    axis(side = 1, at = break_points_2P, 
        labels = c("-2KB","TSS","TSS+500","TES-500","TES","+1KB"))

    legend(x="topleft", fill=colori, 
        legend=colnames(all.noNorm),bg="white",cex=0.8)

    if (!is.null(savePlotPath))
    {
        dev.off()
        message("pdf saved under ",filename)
    }

    ####################
    ###normalized plot and values
    ####################

    common_genes<-rownames(input)[rownames(input) %in% rownames(chip)]
    frameNormalized<-colMeans(
        t(t(chip[common_genes,])-t(input[common_genes,])), 
        na.rm=TRUE)

    hotSpotsValuesNorm<-f_spotfunctionNorm(frameNormalized, 
        breaks=break_points_2P, tag="geneBody")
    maxAucValuesNorm<-f_maximaAucfunctionNorm(frameNormalized,
        breaks=break_points_2P, 
        estBinSize=estimated_bin_size_2P, tag="geneBody")
    if (!is.null(savePlotPath))
    {
        filename=file.path(savePlotPath,"ScaledMetaGene_normalized.pdf")
        pdf(filename,width=10, height=7)
    }
    par(mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1);
    plot(x=as.numeric(names(frameNormalized)),y=frameNormalized, type="l", 
        lwd=2, lty=1, col="orange",xlab="metagene coordinates", 
        ylab="mean log2 enrichment (signal/input)",
        main="normalized metagene", xaxt='n')#,cex.axis=1.3,cex.lab=1.3)
    
    currBreak_points=break_points_2P[c(-2,-5)]##c(-2000,500,2500,4000)
    abline(v=c( break_points_2P[c(2,5)]),lty=2,col="darkgrey", lwd=3)
    abline(v=currBreak_points,lty=3,col="darkgrey", lwd=2)
    axis(side = 1, at = break_points_2P, 
        labels = c("-2KB","TSS","TSS+500","TES-500","TES","+1KB"))

    legend(x="topleft", fill=colori, legend=colnames(all.noNorm), 
        bg="white", cex=0.8)

    if (!is.null(savePlotPath))
    {
        dev.off()
        message("pdf saved under ",filename)
    }

    p1=rbind(hotSpotsValues,maxAucValues)
    p3=rbind(hotSpotsValuesNorm, maxAucValuesNorm)
    result=data.frame(cbind(p1,p3)) 

    if (debug)
    {
        message("Debugging mode ON")
        outname=file.path(getwd(), "geneBody.result")
        file.remove(outname)
        write.table(result,file=outname,
            row.names = TRUE, col.names=TRUE, quote = FALSE)
    }

    return(result)   
}
