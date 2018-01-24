#'@title Wrapper function to create scaled and non-scaled metageneprofiles 
#'
#' @description
#' Metagene plots show the signal enrichment around a region of interest like the TSS or over a predefined set of genes.
#' The tag density of the immunoprecipitation is taken over all RefSeg annotated human genes, averaged and log2 transformed. The same is done for the 
#' input. The normalized profile is calculated as the signal enrichment (immunoprecipitation over the input). 
#' We created two types of metagene profiles: a non-scaled profile for the TSS  and TES, and a scaled profile 
#' for the entire gene, including the gene body. 
#' The non-scaled profile is constructed around the TSS/TES, with 2KB up- and downstream regions respectively. 
#' 
#' f_CreateMetageneProfile
#'
#' @param smoothed.densityChip DESCRIBE
#' @param smoothed.densityInput DESCRIBE
#' @param tag.shift Integer, tag shift returned by f_CrossCorrelation()
#' @param annotationID String, indicating the genome assembly (Default="hg19")
#' @param debug Boolean to enter in debugging mode (default= FALSE)
#'
#' @return list with 3 objects: scaled profile ("twopoint"), non-scaled profile for TSS (TSS) and TES (TES). each object is made of two lists 
#' the chip and the input profile
#'
#' @examples
#'\{dontrun
#' Meta_Result=f_CreateMetageneProfile(smoothedDensityChip,smoothedDensityInput,tag.shift,annotationID="hg19",debug=FALSE)
#'}

annotationPath="/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/"

f_CreateMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,annotationID,debug=FALSE)
{

  annotation=file.path(annotationID,"RefSeqAllGenesFiltered.RData")
  geneAnnotations_file<-file.path(annotationPath,annotation)
  #loading genome annotation
  source("FunctionsLocal.R")
  print("Load geneannotation")
  load(geneAnnotations_file) #RefSeqGenesAll_object
  current_annotations_type<-gsub(pattern=".RData", replacement="", fixed=TRUE, x=basename(geneAnnotations_file))
  current_annotations_object=RefSeqGenes_annotated_filteredByOverlap_geneLength
  # format annotations as chromosome lists
  current_annotations_object<-data.frame(current_annotations_object@.Data, annotation(current_annotations_object), stringsAsFactors=FALSE)
  current_annotations_object$interval_starts<-as.integer(current_annotations_object$interval_starts)
  current_annotations_object$interval_ends<-as.integer(current_annotations_object$interval_ends)
  current_annotations_object$seq_name<-as.character(current_annotations_object$seq_name)
  annotatedGenesPerChr <-split(current_annotations_object, f=current_annotations_object$seq_name)

  ##two.point.scaling
  #create scaled metageneprofile
  #input
  print("Calculating scaled metageneprofile ...")
  smoothed.densityInput=list(td=smoothed.densityInput)
  print("process input")
  binned_Input = f_t.get.gene.av.density(smoothed.densityInput,gdl=annotatedGenesPerChr)
  #Chip
  smoothed.densityChip=list(td=smoothed.densityChip)
  print("process ChIP")
  binned_Chip= f_t.get.gene.av.density(smoothed.densityChip,gdl=annotatedGenesPerChr)

  twopoint=list(chip=binned_Chip,input= binned_Input)

  if (debug)
  {
    save(binned_Chip, binned_Input,file=file.path(getwd(), paste(chipName,inputName,"Twopoint.RData",sep="_")))
  }

   ##one.point.scaling
  #create non-scaled metageneprofile for TSS
  print("Creating scaled metageneprofiles...")

  print("...TSS")
  binnedInput_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityInput,gdl=annotatedGenesPerChr)
  binnedChip_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityChip,gdl=annotatedGenesPerChr)

  onepointTSS=list(chip=binnedChip_TSS,input= binnedInput_TSS)
  if (debug)
  {
    save(binnedChip_TSS, binnedInput_TSS,file=file.path(getwd(), paste(chipName,inputName,"OnePointTSS.RData",sep="_")))
  }

  ##one.point.scaling
  #create non-scaled metageneprofile for TES

  print("...TES")
  binnedInput_TES <- f_t.get.gene.av.density_TES(smoothed.densityInput,gdl=annotatedGenesPerChr)
  binnedChip_TES <- f_t.get.gene.av.density_TES(smoothed.densityChip,gdl=annotatedGenesPerChr)

  onepointTES=list(chip=binnedChip_TES,input= binnedInput_TES)
  if (debug)
  {
    save(binnedChip_TES, binnedInput_TES,file=file.path(getwd(), paste(chipName,inputName,"OnePointTES.RData",sep="_")))
  }


  return(list("twopoint"=twopoint,"TSS"=onepointTSS,"TES"=onepointTES))
}

