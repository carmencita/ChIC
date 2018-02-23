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
#' @param twopointElements  DESCRIBE!!!
#' @param TSSElements DESCRIBE!!!
#' @param TESElements DESCRIBE!!!
#' @param savePlotPath Path in which the profiles (in pdf) 
#' should be saved (Default=working directory)
#'
#' @return returnList
#'@examples
#' print("possible usage: ")
#'\dontrun{
#' metagenePlotsForComparison(chrommark="H3K4me1", 
#' Meta_Result$twopoint, 
#' Meta_Result$TSS, 
#' Result_Meta$TES,plotName=chipName,
#' savePlotPath=getwd())
#' }

metagenePlotsForComparison<-function(chrommark, twopointElements, TSSElements, 
TESElements, savePlotPath=NULL)
{
    compendium_profiles=NULL
    data(compendium_profiles, package="chic.data", envir=environment())
    #compendium_profiles=chic.data::compendium_profiles

    Hlist=f_metaGeneDefinition("Hlist")

    psc <- 1; # pseudocount, required to avoid log2 of 0
    if (chrommark %in% Hlist)
    {
        for (endung in c("TWO","TES","TSS"))
        {
            print(endung)
            #loading current profiles 
            if (endung == "TWO")
            {
                iframe=log2(do.call(rbind,twopointElements$input)+psc)
                cframe=log2(do.call(rbind,twopointElements$chip)+psc)
            }else if(endung == "TES"){                
                iframe=log2(do.call(rbind,TESElements$input)+psc)
                cframe=log2(do.call(rbind,TESElements$chip)+psc)    
            }else{
                iframe=log2(do.call(rbind,TSSElements$input)+psc)
                cframe=log2(do.call(rbind,TSSElements$chip)+psc)
            }

            #load data from compendium 
            frame=NULL
            #load average dataframe normalized
            name=paste(chrommark,"_",endung,"norm", sep="")
            frame=compendium_profiles[[name]]
            normMin=min(frame$mean-frame$sderr)
            normMax= max(frame$mean+frame$sderr)
            n_mean=frame

            frame=NULL
            #load average dataframe chip
            name=paste(chrommark,"_",endung,"chip", sep="")
            frame=compendium_profiles[[name]]
            absoluteMin=min(frame$mean-frame$sderr)
            absoluteMax= max(frame$mean+frame$sderr)
            c_mean=frame

            frame=NULL
            #load average dataframe input
            name= paste(chrommark,"_",endung,"input", sep="")
            frame=compendium_profiles[[name]]
            if ((min(frame$mean-frame$sderr))<absoluteMin)
            {
                absoluteMin=min(frame$mean-frame$sderr)
            }
            if ((max(frame$mean+frame$sderr))>absoluteMax)
            {
                absoluteMax=max(frame$mean+frame$sderr)
            }
            i_mean=frame

            
            ##prepare individual data normalized, input and chip
            nframe<-colMeans(t(t(cframe)-t(iframe)),na.rm=TRUE)
            iframe=colMeans(iframe,na.rm=TRUE)
            cframe= colMeans(cframe,na.rm=TRUE)
            
            finalframe=NULL
            finalframe=cbind(c_mean$x,iframe)
            rownames(finalframe)=NULL
            finalframe=as.data.frame(finalframe)    
            colnames(finalframe)=c("x","mean") 
            iframe=finalframe

            finalframe=NULL
            finalframe=cbind(c_mean$x,cframe)
            rownames(finalframe)=NULL
            finalframe=as.data.frame(finalframe)    
            colnames(finalframe)=c("x","mean") 
            cframe=finalframe

            ##get max and min for same y-axis values for chip and input
            newMin=min(cframe$mean,absoluteMin,iframe$mean)
            newMax=max(cframe$mean,absoluteMax,iframe$mean)

            finalframe=NULL
            finalframe=cbind(c_mean$x,nframe)
            rownames(finalframe)=NULL
            finalframe=as.data.frame(finalframe)    
            colnames(finalframe)=c("x","mean") 
            nframe=finalframe
            
            ##get max and min for y-axis values for input
            normMin=min(nframe$mean,normMin)
            normMax=max(nframe$mean,normMax)

            #create comparison plots
            f_plotProfiles(i_mean, iframe, endung, 
                c(newMin-0.001,newMax+0.001), 
                maintitel=paste(chrommark,endung,"Input",sep="_"),
                savePlotPath=savePlotPath)
            f_plotProfiles(c_mean, cframe, endung, 
                c(newMin-0.001,newMax+0.001), 
                maintitel=paste(chrommark,endung,"Chip",sep="_"),
                savePlotPath=savePlotPath)
            f_plotProfiles(n_mean, nframe, endung, 
                c(normMin-0.001,normMax+0.001), 
                maintitel=paste(chrommark,endung,"norm",sep="_"),
                ylab="mean log2 enrichment (signal/input)",
                savePlotPath=savePlotPath)
            if(!is.null(savePlotPath))
            { 
                print(paste("pdf saved under",savePlotPath,sep=" "))
            }
        }
    }else{
        print("Chromatin marks has to be one of the following:")
        print(Hlist)
    }
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
#' @param savePlotPath Default=NULL, when set saves the density plot 
#' (pdf format) under the given path.
#'
#'
#' @return nothing, creates a figure under 'savePlotPath'
#'@examples
#'\dontrun{
#' plotReferenceDistribution(chrommark="H3K4me1",metricToBePlotted="RSC",
#' currentValue=crossvalues_Chip$RSC,savePlotPath=getwd())
#'}

plotReferenceDistribution<-function(chrommark,metricToBePlotted="RSC",
currentValue,savePlotPath=NULL)
{
    Hlist=f_metaGeneDefinition("Hlist")
    allChrom=f_metaGeneDefinition("Classes")
    if (chrommark %in% Hlist)
    {
        ##reading compendium
        compendium_db=NULL
        data(compendium_db,package="chic.data",envir = environment())
        #compendium_db=chic.data::compendium_db
        #load("data/compendium_db.rda")
        ##select the class for the respective chromatin mark
        if (chrommark %in% allChrom$allSharp)
        {
            profileSet=allChrom$allSharp
            tag="(Sharp class)"
        }
        if (chrommark %in% allChrom$allBroad)
        {
            profileSet=allChrom$allBroad
            tag="(Broad class)"
        }
        if (chrommark %in% allChrom$RNAPol2)
        {
            profileSet=allChrom$RNAPol2
            tag="(RNAPol2 class)"
        }
        print(tag)
        #get values for chrommark
        #alias=paste("CC",metricToBePlotted,sep="_")
        alias=NULL
        if (paste("CC",metricToBePlotted,sep="_") %in% colnames(compendium_db))
        {
            alias=paste("CC",metricToBePlotted,sep="_")
        }else {
            helpi= paste("Ch",metricToBePlotted,sep="_")
            if (helpi %in% colnames(compendium_db))
            {    
                alias=paste("Ch",metricToBePlotted,sep="_")
            }
        }        
        
        ##get the values of respective set
        subset= compendium_db[which(compendium_db$CC_TF %in% profileSet),alias]
        ##plot distribution
        f_plotValueDistribution(subset,
            title=paste(metricToBePlotted,"\n", chrommark,tag,set=" "),
            currentValue,savePlotPath)
    }else{
        print("Chromatin marks has to be one of the following:")
        print(Hlist)
    }
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
#' crossCorrelation()
#' @param features_global list, list with QC-metrics returned from 
#' QCscores_global()
#' @param features_local_TSS list, list with QC-metrics returned from 
#' nonScaledMetageneProfile() with option TSS
#' @param features_local_TES list,  list with QC-metrics returned from 
#' nonScaledMetageneProfile() with option TES
#' @param features_local_scaled list, list with QC-metrics returned from 
#' scaledMetageneProfile()
#' @param savePlotPath Default=NULL, when set saves the density plot 
#' (pdf format) under the given path.
#'
#' @return something something
#'@examples
#'\dontrun{
#' something somethign
#'}
plotPredictionScore<-function(chrommark="H3K36me3", features_cc,
features_global,features_local_TSS, features_local_TES, 
features_local_scaled, savePlotPath=NULL)
{
    #library(randomForest)
    print("load profile classes...")
    
    allChrom=f_metaGeneDefinition("Classes")
    rf_models=NULL
    data(rf_models,package="chic.data",envir = environment())
    #rf_models=chic.data::rf_models

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
