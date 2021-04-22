
#' @title Wrapper function that plots non-scaled profiles for TSS of TES and 
#' to collect the QC-metrics
#'
#' @description 
#' The non-scaled profile is constructed around the TSS/TES, with 2KB up- and 
#' downstream regions respectively. Different values are taken at the TSS/TES 
#' and surroundings with +/-2KB, +/-1KB and +/-500 sizes. For all the genomic 
#' positions, we kept the values for the ChIP and the normalized profile, as 
#' the normalization already contains information from the input. Additionally,
#' we calculated for all of the intervals between the predefined positions the
#' area under the profile, the local maxima (x, y coordinates), the variance, 
#' the standard deviation and the quantiles at 0%, 25%, 50% and 75%. In total 
#' the function returns 43 QC-metrics.
#'
#' qualityScores_LM
#'
#' @param data metagene-list for input and chip samples as returned 
#' by createMetageneProfile()
#' @param tag String that can be 'TSS' or 'TES',indicating if the TSS or the 
#' TES profile should be calcualted (Default='TSS')
#' @param savePlotPath if set the plot will be saved under 'savePlotPath'. 
#' Default=NULL and plot will be forwarded to stdout. 
#' @param plot character, possible values "norm", "split", "all" (default) to plot metaprofiles for 
#' normalized ChIP/input enrichment, or the two samples reads densities separated, or both plots (respectively)
#'
#' @export
#'
#' @return result Dataframe with QC-values for chip, input and normalized 
#' metagene profile
#'
#' @examples
#'
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
#' ##estract QC-values and plot metageneprofile for TSS
#' TSS_Scores=qualityScores_LM(data=Meta_Result$TSS, tag="TSS",
#' savePlotPath=filepath))
#'}

qualityScores_LM <- function(data, tag, savePlotPath = NULL, plot="all") 
{
    if (!(tag %in% c("geneBody", "TES", "TSS"))) {
        stop("tag not valid! Please use: geneBody, TES or TSS") }
    stopifnot(length(data[[tag]]) == 3L)

    if (!(plot %in% c("all","norm","split")) ) {
        stop("wrong plot parameter value")
    }

    if (tag=="geneBody") {
           result<- qualityScores_LMgenebody(data=data, tag=tag, savePlotPath = savePlotPath, plot=plot)
        } else {
    
            binnedChip <- data[[tag]]$chip
            binnedInput <- data[[tag]]$input
            binnedNorm <- data[[tag]]$norm
            message("load metagene setting")
            # load('Settings.RData')
            settings <- f_metaGeneDefinition(selection = "Settings")
            break_points <- settings$break_points
            estimated_bin_size_1P <- settings$estimated_bin_size_1P

            ## pseudocount## required to avoid log2 of 0
            psc <- 1
            chip <- log2(do.call(rbind, binnedChip) + psc)
            input <- log2(do.call(rbind, binnedInput) + psc)


            input.noNorm <- colMeans(input, na.rm = TRUE)
            chip.noNorm <- colMeans(chip, na.rm = TRUE)
            all.noNorm <- NULL
            all.noNorm <- cbind(chip.noNorm, input.noNorm)
            colnames(all.noNorm) <- c("Chip", "Input")


            ## values at specific predefined points
            hotSpotsValues <- f_spotfunction(all.noNorm, break_points, tag = tag)
            ## local maxima and area in all the predefined regions
            maxAucValues <- f_maximaAucfunction(all.noNorm, break_points, 
                estimated_bin_size_1P, 
                tag = tag)

            ## chip_dispersion_TES_-1000_0% chip_dispersion_TES_-1000_25%
            ## chip_dispersion_TES_-1000_50% chip_dispersion_TES_-1000_75%
            ## chip_dispersion_TES_-1000_sd chip_dispersion_TES_-1000_variance


            variabilityValues <- f_variabilityValues(all.noNorm, 
                break_points, 
                tag = tag)

            ## make plots
            colori <- c(rev(rainbow(ncol(all.noNorm) - 1)), "black")
            if (plot %in% c("all", "split")) {
                if (!is.null(savePlotPath)) {
                    filename <- file.path(savePlotPath, 
                        paste("ChIP_Input_", tag, ".pdf", sep = ""))
                    pdf(file = filename, width = 10, height = 7)
                }
                par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.65, 0), cex = 1)
                matplot(x = as.numeric(rownames(all.noNorm)),
                    y=all.noNorm,
                    type = "l", lwd = 2, 
                    lty = 1, col = colori, xlab = "metagene coordinates", 
                    ylab = "mean of log2 read density", 
                    main = tag, xaxt = "n")


                newbreak_points <- break_points[-c(4)]  
                # c(-2000,-1000,-500,500,1000,2000)
                abline(v = 0, lty = 2, col = "darkgrey", lwd = 2)
                abline(v = newbreak_points, lty = 3, col = "darkgrey", lwd = 2)
                ## abline(v=0,lty=2,col='darkgrey', lwd=2)###plot TSS
                ## plotPoints=c(-2000,-1000,-500,500,1000,2000)###plot remaining be
                ## abline(v=plotPoints,lty=3,col='darkgrey', lwd=2)
                axis(side = 1, at = break_points, 
                    labels = c("-2KB", "-1KB", "-500", tag, "500", "+1KB", "+2KB"))
                legend(x = "topleft", fill = colori, legend = colnames(all.noNorm), 
                    bg = "white", 
                    cex = 0.8)

                if (!is.null(savePlotPath)) {
                    dev.off()
                    message("pdf saved under ", filename)
                }
            }

            ## normalized plot and values
            #common_genes <- rownames(input)[rownames(input) %in% rownames(chip)]

            # all.Norm <- colMeans(
            #     t(t(chip[common_genes, ]) - t(input[common_genes, ])), 
            #     na.rm = TRUE)

            norm <- do.call(rbind, binnedNorm)
            normFinal <- data.frame(norm=colMeans(norm, na.rm = TRUE))


            hotSpotsValuesNorm <- f_spotfunctionNorm(normFinal, 
                break_points, tag = tag)

            maxAucValuesNorm <- f_maximaAucfunctionNorm(normFinal, 
                break_points, estimated_bin_size_1P, 
                tag = tag)

            variabilityValuesNorm <- f_variabilityValuesNorm(normFinal, 
                break_points, tag = tag)

            if (plot %in% c("all", "norm")) {
                if (!is.null(savePlotPath)) {
                    filename <- file.path(savePlotPath, 
                        paste("Normalized_", tag, ".pdf", sep = ""))
                    pdf(file = filename, width = 10, height = 7)
                }
                par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.65, 0), cex = 1)
                plot(x = as.numeric(rownames(normFinal)),
                    y=normFinal$norm, 
                    type = "l", lwd = 2, lty = 1, 
                    col = "orange", xlab = "metagene coordinates", 
                    ylab = "mean of log2 enrichment (signal/input)", 
                    main = paste("normalized", tag, sep = " "), xaxt = "n")

                abline(v = 0, lty = 2, col = "darkgrey", lwd = 2)
                abline(v = newbreak_points, lty = 3, col = "darkgrey", lwd = 2)
                axis(side = 1, at = break_points, 
                    labels = c("-2KB", "-1KB", "-500", tag, "500", "+1KB", "+2KB"))
                legend(x = "topleft", fill = "orange", 
                    legend = "Normalized", bg = "white", 
                    cex = 0.8)

                if (!is.null(savePlotPath)) {
                    dev.off()
                    message("pdf saved under ", filename)
                }
            }

            # convert values and features to one frame
            result=data.frame(rbind(
                cbind(round(hotSpotsValues,3),round(hotSpotsValuesNorm,3)),
                cbind(round(variabilityValues[[1]],3),round(variabilityValuesNorm[[1]],3)),
                cbind(round(variabilityValues[[2]],3),round(variabilityValuesNorm[[2]],3)),
                cbind(round(variabilityValues[[3]],3),round(variabilityValuesNorm[[3]],3)),
                cbind(round(maxAucValues,3),round(maxAucValuesNorm,3))))

            message("Calculation of LM for non-scaled profiles done!")
        
        } ## closing the upper if for switching between qualityScores_LM and qualityScores_LMgenebody
    return(result)
}

