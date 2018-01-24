# #########################
# ##### FOR DEVEL ONLY

# ## load gene annotations
# require(girafe)
# require(snow)
# require(parallel)


# source("FunctionsLocal.R")
# path=getwd()
# dataPath<-"/lustre/data/FF/Carmen/BitBucket/chic/data/"
# source(file.path(path,"GlobalParameters.R"))


# chrominfo<-read.table(chrominfo_file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
# rownames(chrominfo)<-chrominfo$chrom
# rngl<-lapply(split(x=chrominfo$size, f=chrominfo$chrom), FUN=function(x) {return(c(1,x))})


# chipName="ENCFF000BBB"
# inputName="ENCFF000BAF"
# debug=TRUE
# cluster=NULL
# dataPath="/lustre//data/FF/Carmen/BitBucket/chic/data"



# tag.shift= 95
# mc=1

# #load Tagdensity Input smoothed.densityInput
# load(file.path(path,paste(chipName,inputName,"TagDensityInput.RData",sep="_")))

# #load Tagdensity ChIP smoo
# load(file.path(path,paste(chipName,inputName,"TagDensityChip.RData",sep="_")))


# #########################
# ##### FOR DEVEL ONLY END


f_CreateMetageneProfile = function(smoothed.densityChip,smoothed.densityInput,tag.shift,geneAnnotations_file,debug=FALSE)
{

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
  #Input smoothed.densityInput
  print("Calculate two point scaling")
  smoothed.densityInput=list(td=smoothed.densityInput)
  print("process input")
  binned_Input = f_t.get.gene.av.density(smoothed.densityInput,gdl=annotatedGenesPerChr)
  #Chip smoothed.densityChip
  smoothed.densityChip=list(td=smoothed.densityChip)
  print("process ChIP")
  binned_Chip= f_t.get.gene.av.density(smoothed.densityChip,gdl=annotatedGenesPerChr)

  twopoint=list(chip=binned_Chip,input= binned_Input)
  if (debug)
  {
    save(binned_Chip, binned_Input,file=file.path(getwd(), paste(chipName,inputName,"Twopoint.RData",sep="_")))
  }

   ##one.point.scaling
  print("Calculate one point scaling...")

  print("...TSS")
  binnedInput_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityInput,gdl=annotatedGenesPerChr)
  binnedChip_TSS <- f_t.get.gene.av.density_TSS(smoothed.densityChip,gdl=annotatedGenesPerChr)

  onepointTSS=list(chip=binnedChip_TSS,input= binnedInput_TSS)
  if (debug)
  {
    save(binnedChip_TSS, binnedInput_TSS,file=file.path(getwd(), paste(chipName,inputName,"OnePointTSS.RData",sep="_")))
  }

   ##one.point.scaling
  print("...TSS")
  binnedInput_TES <- f_t.get.gene.av.density_TES(smoothed.densityInput,gdl=annotatedGenesPerChr)
  binnedChip_TES <- f_t.get.gene.av.density_TES(smoothed.densityChip,gdl=annotatedGenesPerChr)

  onepointTES=list(chip=binnedChip_TES,input= binnedInput_TES)
  if (debug)
  {
    save(binnedChip_TES, binnedInput_TES,file=file.path(getwd(), paste(chipName,inputName,"OnePointTES.RData",sep="_")))
  }


  return(list("twopoint"=twopoint,"TSS"=onepointTSS,"TES"=onepointTES))
}

