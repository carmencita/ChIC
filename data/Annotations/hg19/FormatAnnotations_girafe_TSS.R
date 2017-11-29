
workingdir<-"."
outdir<-workingdir
RefSeqGenes_file<-"refGene.txt"#"hg19.refGene.txt"#"mm9.RefSeqGenes.txt"
RefSeqGenes_colClasses_custom<-c("NULL", "character", "character", "character", "integer", "integer", "integer", "integer", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "character", "NULL")

refLink_file<-"refLink.txt"#"hg19.refLink.txt"#"mm9.refLink.txt"
refLink_colClasses_custom<-c("character", "character", "character", "character", "NULL", "NULL", "character", "NULL")



require(girafe)
require(genomeIntervals)
setwd(workingdir)






#####################################################
#####################################################
###
### FORMAT REFSEQ BASED GENE INTERVALS ANNOTATIONS
###
#####################################################
#####################################################

RefSeqGenes<-read.table(file=RefSeqGenes_file, colClasses=RefSeqGenes_colClasses_custom, comment.char = "", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

refLink<-read.table(file=refLink_file, colClasses=refLink_colClasses_custom, comment.char = "", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")


RefSeqGenes_annotated<-merge(x=RefSeqGenes, y=refLink, 
           by.x = "name", by.y ="mrnaAcc", all = FALSE, all.x = TRUE, all.y = FALSE,
           sort = TRUE)

# For NON protein coding gene the field product and protAcc are empy (==""). Moreover the fields cdsStartStat and cdsEndStat are =="unk"
# we can use these fields to filter out non protein coding genes
proteinCodingSelection<- (!((RefSeqGenes_annotated$protAcc=="") & (RefSeqGenes_annotated$product == "") & (RefSeqGenes_annotated$cdsStartStat == "unk") & (RefSeqGenes_annotated$cdsEndStat == "unk")))
RefSeqGenes_annotated<-RefSeqGenes_annotated[proteinCodingSelection, ]





## TSS +/- 2Kb intervals
print("TSS +/- 2Kb intervals")


stradSpecific_TSS<-ifelse(RefSeqGenes_annotated$strand=="+", yes=RefSeqGenes_annotated$txStart, no=RefSeqGenes_annotated$txEnd)

RefSeqGenesTSS2kbwindow_genomeIntervals_object<-new("Genome_intervals_stranded",
              .Data=(cbind((stradSpecific_TSS-2000), (stradSpecific_TSS+2000)) ),
              closed=TRUE,
              annotation = data.frame(
                    seq_name = factor(RefSeqGenes_annotated$chrom),
                    inter_base=FALSE,
                    strand = factor(RefSeqGenes_annotated$strand, levels=c("-","+")),
                    txStart = RefSeqGenes_annotated$txStart,
                    txEnd = RefSeqGenes_annotated$txEnd,
                    RefSeq = RefSeqGenes_annotated$name,
                    geneSymbol = RefSeqGenes_annotated$X.name,
                    description = RefSeqGenes_annotated$product,
                    entrezID= RefSeqGenes_annotated$locusLinkId
              )
)
# PLEASE NOTICE THAT ORIGINAL REFSEQ txStart AND txEnd ARE INCORPORATED AS PART OF ANNOTATION SLOT INTO Genome_intervals_stranded OBJECT
save(RefSeqGenes_annotated, RefSeqGenesTSS2kbwindow_genomeIntervals_object, file="RefSeqGenesTSS2kbwindow.RData")









