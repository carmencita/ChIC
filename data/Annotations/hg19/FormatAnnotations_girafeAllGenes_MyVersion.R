RefSeqGenes_file<-"refGene.txt"#"mm9.RefSeqGenes.txt"
RefSeqGenes_colClasses_custom<-c("NULL", "character", "character", "character", "integer", "integer", "integer", "integer", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "character", "NULL")

refLink_file<-"refLink.txt"#"mm9.refLink.txt"
refLink_colClasses_custom<-c("character", "character", "character", "character", "NULL", "NULL", "character", "NULL")



require(girafe)
require(genomeIntervals)

source("/data/FF/Carmen/Pipeline_allEncodeBeta/GlobalParameters.R")
#####################################################
#####################################################
###
### FORMAT REFSEQ BASED GENE INTERVALS ANNOTATIONS
###
#####################################################
#####################################################

##I am taking ALL the genes. No filtering for protein coding or whatever


RefSeqGenes<-read.table(file=RefSeqGenes_file, colClasses=RefSeqGenes_colClasses_custom, comment.char = "", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
refLink<-read.table(file=refLink_file, colClasses=refLink_colClasses_custom, comment.char = "", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
RefSeqGenes_annotated<-merge(x=RefSeqGenes, y=refLink, by.x = "name", by.y ="mrnaAcc", all = FALSE, all.x = TRUE, all.y = FALSE,sort = TRUE)

## RefSeq Genes intervals
print("RefSeq Genes intervals")
RefSeqAllGenes_object<-new("Genome_intervals_stranded",
              .Data=(cbind(RefSeqGenes_annotated$txStart, RefSeqGenes_annotated$txEnd) ),
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
save(RefSeqAllGenes_object, file="RefSeqAllGenes.RData")
length(RefSeqAllGenes_object@annotation$geneSymbol)


#########################################################
#########################################################
#########################################################


### STRAND SPECIFIC START AND END SPECIFIC INTERVAL
###enlarges the gene length by half_window_width_upstream and half_window_width_downstream to have the regions around the TSS and TES
interval_5end<-ifelse(RefSeqGenes_annotated$strand=="+", yes=(RefSeqGenes_annotated$txStart - half_window_width_upstream), no=(RefSeqGenes_annotated$txEnd + half_window_width_upstream))
interval_3end<-ifelse(RefSeqGenes_annotated$strand=="+", yes=(RefSeqGenes_annotated$txEnd + half_window_width_downstream), no=(RefSeqGenes_annotated$txStart - half_window_width_downstream))


interval_starts<-ifelse(RefSeqGenes_annotated$strand=="+", yes=interval_5end, no=interval_3end)
interval_ends<-ifelse(RefSeqGenes_annotated$strand=="+", yes=interval_3end, no=interval_5end)

RefSeqAllGenesFiltered_genomeIntervals_object<-new("Genome_intervals_stranded",
              .Data=cbind(interval_starts, interval_ends),
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
                    entrezID= RefSeqGenes_annotated$locusLinkId,
                    stringsAsFactors=FALSE
              )
)


##Filter for overlapping genes


#AnnotationsOverlap<-interval_overlap(RefSeq25kbGenes2kbUp1kbDown_genomeIntervals_object, RefSeq25kbGenes2kbUp1kbDown_genomeIntervals_object)
AnnotationsOverlap<-interval_overlap(RefSeqAllGenesFiltered_genomeIntervals_object, RefSeqAllGenesFiltered_genomeIntervals_object)
count_overlaps<-(sapply(AnnotationsOverlap, length))
select_non_overlappings<-(count_overlaps == 1)


RefSeqGenes_annotated_filteredByOverlap<-RefSeqAllGenesFiltered_genomeIntervals_object[select_non_overlappings,]

##filtering by minimal gene lengh min_gene_length
#min_length_selection<-((RefSeqGenes_annotated$txEnd-RefSeqGenes_annotated$txStart)>=min_gene_length)
#RefSeqGenes_annotated_filtered<-RefSeqGenes_annotated[min_length_selection,]
min_length_selection<-((RefSeqGenes_annotated_filteredByOverlap$txEnd-RefSeqGenes_annotated_filteredByOverlap$txStart)>=min_gene_length)
RefSeqGenes_annotated_filteredByOverlap_geneLength<-RefSeqGenes_annotated_filteredByOverlap[min_length_selection,]



length(RefSeqGenes_annotated_filteredByOverlap_geneLength@annotation$geneSymbol)

save(RefSeqGenes_annotated_filteredByOverlap_geneLength, file="RefSeqAllGenesFiltered.RData")



