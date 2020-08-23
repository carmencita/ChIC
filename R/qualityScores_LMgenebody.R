#' @title Wrapper function to plot the scaled metagene- profile and to 
#' collect the QC-metrics
#'
#' @description The scaled metagene profile that includes the gene body, the 
#' signal is captured on a real scale from the TSS and an upstream region of 
#' 2KB. From the TSS, the gene body is constructed with 0.5KB in real scale 
#' at the gene start (TSS + 0.5KB) and the gene end (TES - 0.5KB), whereas 
#' the remaining gene body is scaled to a virtual length of 2000. Considering 
#' the length of these regions, the minimum gene length required is 3KB and 
#' shorter genes are filtered out. From the profile, we take enrichment values 
#' at different coordinates: at -2KB, at the TSS, inner margin (0.5KB), gene 
#' body (2KB + 2 * inner margin), gene body+1KB. We collect in total 42 
#' QC-metrics from the ChIP and normalized profile. 
#'
#' qualityScores_LMgenebody
#'
#' @param data metagene-list for input and chip sample of the genebody profile
#' returned by createMetageneProfile()
#' @param savePlotPath if set the plot will be saved under 'savePlotPath'. 
#' Default=NULL and plot will be forwarded to stdout. 
#'
#' @export
#'
#' @return returnList
#'
#' @examples
#'
#' ## This command is time intensive to run
#' ##
#' ## To run the example code the user must provide two bam files for the ChIP
#' ## and the input and read them with the readBamFile() function. To make it
#' ## easier for the user to run the example code we provide tow bam examples 
#' ## (chip and input) in our ChIC.data package that have already been loaded 
#' ## with the readBamFile() function.
#'
#' mc=4
#' finalTagShift=82
#' \dontrun{
#'
#' filepath=tempdir()
#' setwd(filepath)
#'
#' data("chipBam", package = "ChIC.data", envir = environment())
#' data("inputBam", package = "ChIC.data", envir = environment())
#'
#' ## calculate binding characteristics 
#' chip_binding.characteristics<-spp::get.binding.characteristics(
#'     chipBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
#'
#' input_binding.characteristics<-spp::get.binding.characteristics(
#'     inputBam, srange=c(0,500), bin=5,accept.all.tags=TRUE)
#'
#' ##get chromosome information and order chip and input by it
#' chrl_final=intersect(names(chipBam$tags), names(inputBam$tags))
#' chipBam$tags=chipBam$tags[chrl_final]
#' chipBam$quality=chipBam$quality[chrl_final]
#' inputBam$tags=inputBam$tags[chrl_final]
#' inputBam$quality=inputBam$quality[chrl_final]
#'
#' ##remove sigular positions with extremely high read counts with 
#' ##respect to the neighbourhood
#' selectedTags=removeLocalTagAnomalies(chipBam, inputBam, 
#' chip_binding.characteristics, input_binding.characteristics)
#'
#' inputBamSelected=selectedTags$input.dataSelected
#' chipBamSelected=selectedTags$chip.dataSelected
#'
#' ##smooth input and chip tags
#' smoothedChip=tagDensity(chipBamSelected, tag.shift=finalTagShift)
#' smoothedInput=tagDensity(inputBamSelected, tag.shift=finalTagShift)
#' 
#' ##calculate metagene profiles
#' Meta_Result=createMetageneProfile(smoothed.densityChip=smoothedChip, 
#'     smoothedInput,tag.shift=finalTagShift, mc=mc)
#'
#' #create scaled metagene profile
#' geneBody_Scores=qualityScores_LMgenebody(Meta_Result$geneBody,
#' savePlotPath=filepath)
#'}


qualityScores_LMgenebody <- function(data, savePlotPath = NULL)
{
    stopifnot(length(data) == 3L)
    
    binnedChip <- data$chip
    binnedInput <- data$input
    binnedNorm <- data$norm
    message("load metagene setting")
    settings <- f_metaGeneDefinition(selection = "Settings")
    ## pseudocount## required to avoid log2 of 0
    psc <- 1
    break_points_2P <- settings$break_points_2P
    estimated_bin_size_2P <- settings$estimated_bin_size_2P
    chip <- log2(do.call(rbind, binnedChip) + psc)
    input <- log2(do.call(rbind, binnedInput) + psc)
    
    input.noNorm <- colMeans(input, na.rm = TRUE)
    chip.noNorm <- colMeans(chip, na.rm = TRUE)
    all.noNorm <- cbind(chip.noNorm, input.noNorm)
    colnames(all.noNorm) <- c("Chip", "Input")
    
    
    ## values at specific predefined points
    hotSpotsValues <- f_spotfunction(all.noNorm, 
        break_points_2P, tag = "geneBody")
    ## local maxima and area in all the predefined regions
    maxAucValues <- f_maximaAucfunction(all.noNorm, 
        breaks = break_points_2P, 
        estBinSize = estimated_bin_size_2P, 
        tag = "geneBody")
    
    ## make plots
    colori <- c(rev(rainbow(ncol(all.noNorm) - 1)), "black")
    if (!is.null(savePlotPath)) {
        filename <- file.path(savePlotPath, "ScaledMetaGene_ChIP_Input.pdf")
        pdf(file = filename, width = 10, height = 7)
    }
    par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.65, 0), cex = 1)
    matplot(x = as.numeric(rownames(all.noNorm)), 
        y = all.noNorm, type = "l", lwd = 2, 
        lty = 1, col = colori, xlab = "metagene coordinates", 
        ylab = "mean of log2 read density", 
        main = "scaled metagene profile", xaxt = "n")
    
    currBreak_points <- break_points_2P[c(-2, -5)]
    abline(v = c(break_points_2P[c(2, 5)]), lty = 2, 
        col = "darkgrey", lwd = 3)
    abline(v = currBreak_points, lty = 3, 
        col = "darkgrey", lwd = 2)
    axis(side = 1, at = break_points_2P, 
        labels = c("-2KB", "TSS", "TSS+500", "TES-500", "TES", "+1KB"))
    
    legend(x = "topleft", fill = colori, legend = colnames(all.noNorm), 
        bg = "white", 
        cex = 0.8)
    
    if (!is.null(savePlotPath)) {
        dev.off()
        message("pdf saved under ", filename)
    }
    
    #################### normalized plot and values
    
    # common_genes <- rownames(input)[rownames(input) %in% rownames(chip)]
    # frameNormalized <- colMeans(
    #     t(t(chip[common_genes, ]) - t(input[common_genes, ])), 
    #     na.rm = TRUE)
    norm <- do.call(rbind, binnedNorm)
    normFinal <- data.frame(norm=colMeans(norm, na.rm = TRUE))
    
    hotSpotsValuesNorm <- f_spotfunctionNorm(normFinal, 
        breaks = break_points_2P, 
        tag = "geneBody")
    maxAucValuesNorm <- f_maximaAucfunctionNorm(normFinal, 
        breaks = break_points_2P, 
        estBinSize = estimated_bin_size_2P, tag = "geneBody")

    if (!is.null(savePlotPath)) {
        filename <- file.path(savePlotPath, "ScaledMetaGene_normalized.pdf")
        pdf(filename, width = 10, height = 7)
    }
    par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.65, 0), cex = 1)
    plot(x = as.numeric(rownames(normFinal)), 
        y = normFinal$norm, 
        type = "l", 
        lwd = 2, lty = 1, col = "orange", 
        xlab = "metagene coordinates", 
        ylab = "mean of log2 enrichment (signal/input)", 
        main = "normalized scaled metagene profile", xaxt = "n")  
            #,cex.axis=1.3,cex.lab=1.3)
    
    currBreak_points <- break_points_2P[c(-2, -5)]  ##c(-2000,500,2500,4000)
    abline(v = c(break_points_2P[c(2, 5)]), lty = 2, 
        col = "darkgrey", lwd = 3)
    abline(v = currBreak_points, lty = 3, 
        col = "darkgrey", lwd = 2)
    axis(side = 1, at = break_points_2P, 
        labels = c("-2KB", "TSS", "TSS+500", "TES-500", "TES", "+1KB"))
    
    legend(x = "topleft", fill = "orange", 
        legend = colnames(normFinal), bg = "white", 
        cex = 0.8)
    
    if (!is.null(savePlotPath)) {
        dev.off()
        message("pdf saved under ", filename)
    }
    
    result=data.frame(rbind(
        cbind(round(hotSpotsValues,3),round(hotSpotsValuesNorm,3)),
        cbind(round(maxAucValues,3),round(maxAucValuesNorm,3))))

    message("Calculation of LM for scaled profile done!")
    return(result)
}


