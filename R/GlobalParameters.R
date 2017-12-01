
###
### GLOBAL PARAMETERS
###


# workingdir<-"/gpfs/work/IscrC_CONCEPt/QC_pipeline/"
# bamdir="/pico/scratch/userexternal/clivi000/bam/"
# storagedir="/pico/scratch/userexternal/clivi000/Storage/"
#sampleinfo_file<-paste(workingdir,"/SampleinfoAlpha.txt",sep="")

chrominfo_file<-"/lustre/data/FF/Carmen/BitBucket/chic/data/Annotations/hg19/hg19.chromInfo.txt"

###
###CLUSTER settings
cluster_ON_OFF=TRUE
#if TRUE sets cluster=makeCluster(mc, type="MPI")
###



###
### savedata parameters
###
#reads.aligner.type<-"bowtie"
#reads.aligner.type<-"bam"
#reads.aligner.type<-"tagalign"


###
### wig parameters
###
#wig.bw <- 100;
#wig.step <- 50;



###
### cross_correlation parameters
###
custom_chrorder<-paste("chr", c(1:19, "X","Y"), sep="")
#custom_chrorder<-paste("chr", c(1:22, "X","Y"), sep="")
estimating_fragment_length_range<-c(0,500)
estimating_fragment_length_bin<-5



###
### Phantom peaks
###

PhantomPeak_range<-c(-500, 1500)
PhantomPeak_bin<-5
PhantomPeak_exclude.min<-10


###
### cross_correlation_customShift_withSmoothing parameters
###
cross_correlation_range_subset<-100:500
cross_correlation_smoothing_k<-10
cross_correlation_smoothing_color<-"blue"

##Smoothing parameters
smoothingBandwidth<-50 ##is slidingwindow shift 
smoothingStep<-20	##is size of sw ##before it was 10


###
### broadRegionsOfEnrichment_loop parameters
###
window_sizes_list<-c(1000, 500, 200)[1]
remove.local.tag.anomalies_filter<-TRUE
select.informative.tags_filter<-FALSE




###
###metagene annotation settings
###
### example of hashtable-like data structure to handle genome specific annotation files
#Annotation_hashTable<-list(
#"hg19"=c(
#	"chrominfo"="/data/FF/Carmen/Annotations/hg19/hg19.chromInfo.txt",
#	"geneAnnotations"= "/data/FF/Carmen/Annotations/hg19/RefSeqAllGenesFiltered.RData"
#	),

#"hg18"=c(
#	"chrominfo"="/data/FF/Carmen/Annotations/hg18/hg18.chromInfo.txt",
#	"geneAnnotations"= "/data/FF/Carmen/Annotations/hg18/RefSeqAllGenesFiltered.RData"
#	),

#"dm3"=c(
#	"chrominfo"="/data/FF/Carmen/Annotations/dm3/dm3.chromInfo.txt",
#	"geneAnnotations"= "/data/FF/Carmen/Annotations/dm3/RefSeqAllGenesFiltered.RData"
#	))

interval_starts<-NULL
interval_ends<-NULL

half_window_width_upstream<-2000
half_window_width_downstream<-2000
min_gene_length<-2000 ##select only genes with at least min_gene_length 


###
### MetageneProfiles/TSSprofiles parameters
###PLOT for TWO POINT SCALING
###
psc <- 1; # pseudocount # required to avoid log2 of 0
geneAnnotations_file<-"/gpfs/work/IscrC_CONCEPt/Annotations/hg19/RefSeqAllGenesFiltered.RData"


#### this is the less "painful" way to reduce the smoothing step in the output without changing all parameters
#### however this will result in dropping a few more genes (if possible double check how many genes are lost)
min.feature.size_carmen=3000 # min gene length filtering applied inside function feature.bin.averages called by t.get.gene.av.density 
### REMEMBER this might result in some discrepancies in the list of genes used for the two point scaling (where this filter is applied) and the one point scaling (where this filter is not applied?? double check if this is the case)
#min.feature.size_carmen=2000 # min gene length filtering applied inside function feature.bin.averages called by t.get.gene.av.density 


##TWO point scaling
predefnum_bins=301   # 151
inner_margin=500   # 500
right_outer_margin=1010   # 1020
left_outer_margin=2010    # 2020
gene_body=2000   # 2000

totalGeneLength=gene_body+2*inner_margin
inner_margin2=totalGeneLength-inner_margin

break_points_2P=c(-2000,0,inner_margin,gene_body+inner_margin,totalGeneLength,totalGeneLength+1000)

total_output_gene_length<-left_outer_margin+inner_margin+gene_body+inner_margin+right_outer_margin
estimated_bin_size_2P<-total_output_gene_length/predefnum_bins

##ONE point scaling
up_downStream=4020
predefnum_bins_1P=201 ##binsize=20
break_points=c(-2000,-1000,-500,0,500,1000,2000)
### %%% need to use this for one point scaling estimated bin size....
estimated_bin_size_1P<-up_downStream/predefnum_bins_1P


### %%% if you want one point and two point scaling parameters to be connected you can do something like this
# up_downStream=left_outer_margin*2
# predefnum_bins_onePS=up_downStream/20 ##binsize=20
# estimated_bin_size_onePoint<-up_downStream/predefnum_bins_onePS



