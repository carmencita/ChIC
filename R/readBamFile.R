
#'@title Read bam file
#'
#' @description
#' Reading bam file format
#' 
#' readBamFile
#'
#' @param filename, name/path of the bam file to be read (without extension)
#' @param readAlignerType, format of the input file. Currently only 'bam' 
#' (default) and 'tagAlign' are supported
#'
#' @return result list of lists, every list corresponds to a chromosome and 
#' contains a vector of coordinates of the 5' ends of the aligned tags.
#'
#' @export
#'
#' @examples
#'
#' ## To run this example code the user MUST provide a bam file: The user can 
#' ## download a ChIP-seq bam file for example from ENCODE:
#' ## https://www.encodeproject.org/files/ENCFF000BFX/
#' ## and save it in the working directory 
#'
#' bamID="ENCFF000BFX"
#' \dontrun{
#' filepath=tempdir()
#' setwd(filepath)
#' system("wget 
#' https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")
#'
#' bamName=file.path(filepath,bamID)
#' chipBam=readBamFile(bamName)
#' }

readBamFile <- function(filename, readAlignerType = "bam") {
    ########## check if input format is ok
    if (!is.character(filename)) 
        stop("Invalid filename (String required)")
    if (!(readAlignerType %in% c("bam","tagAlign"))) 
        stop("Invalid filename (String required)")
    ######### 
    fileContent <- f_readFile(filename = filename, reads.aligner.type = readAlignerType)
    result=f_checkOfChrNames(fileContent)
    return(result)
}

