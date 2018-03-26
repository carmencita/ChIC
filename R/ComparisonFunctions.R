#'@title Function to create metage plots for comparison
#'
#' @description
#' QC-metrics of newly analysed ChIP-seq samples can be compared with 
#' the reference values of the compendium and enrichment profiles 
#' can be plotted against pre-computed profiles of published 
#' datasets. The metagene profiles show the problematic samples signal 
#' (red line) for the ChIP, for the input and their relative 
#' enrichment when compared to the compendium’s mean signal 
#' (black line) and its 2 x standard error (blue shadow). 
#' Additionally the function plots the desired QC-metric
#' as a red dashed line for the sample plotted against the 
#' reference distribution (density plots) of the 
#' compendium values stratified by chromatin marks.
#'
#' metagenePlotsForComparison
#' @param chrommark String, chromatin mark to be analysed. 
#' Has to be one of the following: H3K36me3, H3K4me3, 
#' H3K79me2, H4K20me1,H2AFZ,H3K27me3, H3K9me3,H3K27ac,
#' POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1, H3K4me1,
#' H3K79me1, H3K4ac, H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac,
#' H4K91ac, H4K8ac, H3K18ac, H2BK12ac, H3K56ac, 
#' H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, 
#' H2A.Z, H3K23me2, H2AK9ac, H3T11ph. 
#' For RNAPOL2 different variants are available: POLR2A 
#' (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param data metagene-object of metagene profile by 
#' createMetageneProfile() containing 
#' input and chip profile
#' @param tag indicating the kind of profile to plot. 
#' Can be either: "geneBody", "TES" or "TSS".
#' @param savePlotPath if set the plot will be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout. 
#'
#' @return nothing, creates a figure under 'savePlotPath'
#'
#' @export
#'
#' @examples
#' data(TSSProfile)
#' metagenePlotsForComparison(data=TSSProfile,chrommark="H3K36me3",tag="TSS")



metagenePlotsForComparison<-function(data, chrommark, tag, savePlotPath=NULL)
{
    Hlist=f_metaGeneDefinition("Hlist")
    psc <- 1; # pseudocount, required to avoid log2 of 0

    stopifnot(chrommark %in% Hlist)
    stopifnot(tag %in% c("geneBody","TES","TSS"))
    stopifnot(length(data) == 2L)

    iframe=log2(do.call(rbind,data$input)+psc)
    cframe=log2(do.call(rbind,data$chip)+psc)

    #load average dataframe normalized
    n_mean=f_loadDataCompendium(endung="norm",chrommark=chrommark,tag=tag)
    normMin=min(n_mean$mean-n_mean$sderr)
    normMax= max(n_mean$mean+n_mean$sderr)
    ##load average dataframe chip
    c_mean=f_loadDataCompendium(endung="chip",chrommark=chrommark,tag=tag)
    absoluteMin=min(c_mean$mean-c_mean$sderr)
    absoluteMax= max(c_mean$mean+c_mean$sderr)
    ##load average dataframe input
    i_mean=f_loadDataCompendium(endung="input",chrommark=chrommark,tag=tag)
    ##get the range for x-y axis fro the final plot
    if ((min(i_mean$mean-i_mean$sderr))<absoluteMin)
    {absoluteMin=min(i_mean$mean-i_mean$sderr)}
    if ((max(i_mean$mean+i_mean$sderr))>absoluteMax)
    {absoluteMax=max(i_mean$mean+i_mean$sderr)}
    
    nframe<-colMeans(t(t(cframe)-t(iframe)),na.rm=TRUE)
    iframe=colMeans(iframe,na.rm=TRUE)
    cframe= colMeans(cframe,na.rm=TRUE)

    iframe<-f_prepareData(c_mean,iframe)
    cframe<-f_prepareData(c_mean,cframe)
    nframe<-f_prepareData(c_mean,nframe)
    
    ##get max and min for same y-axis values for chip and input
    newMin=min(cframe$mean,absoluteMin,iframe$mean)
    newMax=max(cframe$mean,absoluteMax,iframe$mean)

    ##get max and min for y-axis values for input
    normMin=min(nframe$mean,normMin)
    normMax=max(nframe$mean,normMax)

    #create comparison plots
    f_plotProfiles(i_mean, iframe, tag,c(newMin-0.001,newMax+0.001),
    maintitel=paste(chrommark,tag,"Input",sep="_"),savePlotPath=savePlotPath)

    f_plotProfiles(c_mean, cframe, tag, c(newMin-0.001,newMax+0.001), 
    maintitel=paste(chrommark,tag,"Chip",sep="_"),savePlotPath=savePlotPath)

    f_plotProfiles(n_mean, nframe, tag, c(normMin-0.001,normMax+0.001), 
    maintitel=paste(chrommark,tag,"norm",sep="_"),
    ylab="mean log2 enrichment (signal/input)", savePlotPath=savePlotPath)
}


#'@title Function to create reference distribution plot for comparison
#'
#' @description
#' Creates a density plot (in pdf) for the sample against the reference 
#' distribution (density plots) of the compendium values stratified by 
#' chromatin marks.
#'
#' plotReferenceDistribution
#'
#' @param chrommark String, chromatin mark to be analysed. Has to be one 
#' of the following: H3K36me3, H3K4me3, H3K79me2, H4K20me1,H2AFZ,H3K27me3, 
#' H3K9me3,H3K27ac, 
#' POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1, H3K4me1, H3K79me1, H3K4ac, 
#' H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac, H4K91ac, H4K8ac, H3K18ac, 
#' H2BK12ac, H3K56ac, H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, H2A.Z, 
#' H3K23me2,
#' H2AK9ac, H3T11ph. For RNAPOL2 different variants are available: POLR2A 
#' (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param metricToBePlotted The metric to be plotted (Default="RSC")
#' @param currentValue The value of the current sample
#' @param savePlotPath if set the plot will be saved under 
#' "savePlotPath". Default=NULL and plot will be forwarded to stdout. 
#'
#' @export
#'
#' @return nothing, creates a figure under 'savePlotPath'
#'@examples
#' print ("Plot distribution of RSC")
#' plotReferenceDistribution(chrommark="H3K4me1",metricToBePlotted="RSC",
#' currentValue=0.8)
#'

plotReferenceDistribution<-function(chrommark,metricToBePlotted="RSC",
    currentValue,savePlotPath=NULL)
{
    Hlist=f_metaGeneDefinition("Hlist")
    
    stopifnot(chrommark %in% Hlist)
    stopifnot(is.numeric(currentValue))
    
    allChrom=f_metaGeneDefinition("Classes")
    ##reading compendium
    #compendium_db=NULL
    #data(compendium_db,package="ChIC.data",envir = environment())
    compendium_db=ChIC.data::compendium_db
    ##select the class for the respective chromatin mark
    #load("data/compendium_db.rda")
    profileInfo=f_getBindingClass(chrommark)

    message(profileInfo$tag)
    #get values for chrommark
    #alias=paste("CC",metricToBePlotted,sep="_")
    alias=NULL
    if (paste("CC",metricToBePlotted,sep="_") %in% colnames(compendium_db))
    { alias=paste("CC",metricToBePlotted,sep="_")
    }else{ 
        helpi= paste("Ch",metricToBePlotted,sep="_")
        if (helpi %in% colnames(compendium_db))
        { alias=paste("Ch",metricToBePlotted,sep="_")}
    }
    ##get the values of respective set
    subset= compendium_db[which(compendium_db$CC_TF %in% 
        profileInfo$profileSet),alias]
    ##plot distribution
    f_plotValueDistribution(subset,
        title=paste(metricToBePlotted,"\n", chrommark,profileInfo$tag,set=" "),
        currentValue,savePlotPath)
}



#'@title Predict score
#'
#' @description
#' BLABLA 
#' plotPredictionScore
#'
#' @param chrommark String, chromatin mark to be analysed. Has to be 
#' one of the following: H3K36me3, H3K4me3, H3K79me2, H4K20me1,H2AFZ,
#' H3K27me3, H3K9me3,H3K27ac, POLR2AphosphoS5, H3K9ac, H3K4me2, H3K9me1,
#' H3K4me1, H3K79me1, H3K4ac, 
#' H3K14ac, H2BK5ac, H2BK120ac, H2BK15ac, H4K91ac, H4K8ac, H3K18ac, 
#' H2BK12ac, H3K56ac, H3K23ac, H2AK5ac, H2BK20ac, H4K5ac, H4K12ac, 
#' H2A.Z, H3K23me2,
#' H2AK9ac, H3T11ph. For RNAPOL2 different variants are available: 
#' POLR2A (for RNAPol2), POLR3G and POLR2AphosphoS2
#' @param features_cc list, with QC-metrics returned from 
#' qualityScores_EM()
#' @param features_global list, list with QC-metrics returned from 
#' qualityScores_GM()
#' @param features_local_TSS list, list with QC-metrics returned from 
#' qualityScores_LM() with option TSS
#' @param features_local_TES list,  list with QC-metrics returned from 
#' qualityScores_LM() with option TES
#' @param features_local_scaled list, list with QC-metrics returned from 
#' qualityScores_LMgenebody()
#' @param savePlotPath Default=NULL, when set saves the density plot 
#' (pdf format) under the given path.
#'
#' @export
#'
#' @return something something
#'@examples
#' print("Prediction")
#'\dontrun{
#' something somethign
#'}
plotPredictionScore<-function(chrommark="H3K36me3", features_cc,
    features_global,features_local_TSS, features_local_TES, 
    features_local_scaled, savePlotPath=NULL)
{
    #library(randomForest)
    message("load profile classes...")
    
    allChrom=f_metaGeneDefinition("Classes")
    #rf_models=NULL
    #data(rf_models,package="ChIC.data",envir = environment())
    rf_models=ChIC.data::rf_models

    #load("data/rf_models.rda")

    Hlist=f_metaGeneDefinition("Hlist")    
    model=NULL
    if (chrommark %in% Hlist){

        if (chrommark %in% allChrom$allSharp)
        {model=rf_models[["sharpEncode"]]}
        
        if (chrommark%in%allChrom$allBroad)
        {model=rf_models[["broadEncode"]]}

        if (chrommark%in%allChrom$RNAPol2)
        {model=rf_models[["RNAPol2Encode"]]}
        
        if (chrommark=="H3K9me3")
        {model=rf_models[["H3K9Encode"]]}
        if (chrommark=="H3K27me3")
        {model=rf_models[["H3K27Encode"]]}
        if (chrommark=="H3K36me3")
        {model=rf_models[["H3K36Encode"]]}

    }else{
        model=rf_models$TFmodel
    }


    part1=data.frame(unlist(c(features_cc,features_global)))
    colnames(part1)="Value"
    part1$Feature=rownames(part1)
    frame=rbind(part1,features_local_TSS,features_local_TES,
        features_local_scaled)
    rownames(frame)=NULL


    selectedFeatures=c(colnames(model$trainingData))
    selectedFeatures=selectedFeatures[-(which(selectedFeatures==".outcome"))]

    featureVector=frame
    a=lapply(as.list(featureVector$Feature),function(element){
        word=strsplit(element,"-")[[1]]
        if (length(word)>1)
        {        
            new=paste(word[1],word[2],sep=".")
            word=new
        }
                
        if (length(grep("%",word))>0)
        {
            word=strsplit(word,"%")[[1]]
            new=paste(word,".",sep="")            
            word=new
        }
        return(word)
    })

    featureVector$FeaturesNew=c(unlist(a))

    helpi=selectedFeatures[!(selectedFeatures %in% featureVector$FeaturesNew)]
    if (length(helpi) != 0)
    {
        print ("ERROR in features")
    }else{
        inn=featureVector[featureVector$FeaturesNew %in% selectedFeatures,]
        test=data.frame(t(inn$Value))
        colnames(test)=c(inn$FeaturesNew)
        test =test[order(colnames(test))]
        prediction=predict(model, newdata=test,type="prob")
    }    
}
