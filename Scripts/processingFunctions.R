#####functions for data processing prior to Linear Modeling#####

## combining replicates (adds values from featureCounts)
combineReplicates<-function(metaData, countData) {
  meta_sub<-metaData[rownames(metaData) %in% colnames(countData),]
  combinedMetadata<-as.data.frame(matrix(ncol=ncol(meta_sub), nrow=0))
  colnames(combinedMetadata)<-colnames(meta_sub)
  combinedCounts<-data.frame(row.names=rownames(countData))
  for (sampleID in unique(meta_sub$ID)) {
    metaToCombine<-meta_sub[meta_sub$ID == sampleID,]
    countsToCombine<-countData[,colnames(countData) %in% rownames(metaToCombine)]
    if (nrow(metaToCombine) == 1) {
      combinedMetadata<-rbind(combinedMetadata, metaToCombine)
      old_col_names<-colnames(combinedCounts)
      combinedCounts<-cbind(combinedCounts, countsToCombine)
      colnames(combinedCounts)<-c(old_col_names, rownames(metaToCombine))
    }
    else {
      countsEntry<-rowSums(countsToCombine)
      old_col_names<-colnames(combinedCounts)
      combinedCounts<-data.frame(cbind(combinedCounts, countsEntry), check.names = FALSE)
      colnames(combinedCounts)<-c(old_col_names, paste0(sampleID,"_combined"))
      metaEntry<-metaToCombine[1,]
      rownames(metaEntry)<-paste0(sampleID,"_combined")
      metaEntry$AutosomalReads<-sum(metaToCombine$AutosomalReads)
      if ("ecoliCounts" %in% colnames(metaData)) {
        metaEntry$ecoliCounts<-sum(metaToCombine$ecoliCounts)
      }
      newLibPrep<-''
      newSeqBatch<-''
      for (sample in rownames(metaToCombine)) {
        newLibPrep<-paste0(newLibPrep,'_', metaToCombine[sample,]$LibPrepBatch)
        newSeqBatch<-paste0(newSeqBatch, '_', metaToCombine[sample,]$SequencingBatch)
      }
      metaEntry$LibPrepBatch<-newLibPrep
      metaEntry$SequencingBatch<-newSeqBatch
      combinedMetadata<-rbind(combinedMetadata, metaEntry)
    }
  }
  combinedMetadata<-combinedMetadata[colnames(combinedCounts),]
  return(list(combinedMetadata, combinedCounts))
} 

## running deseq2

# implement fallback fit types
run_deseq_with_fallback <- function(dds) {
  
  try_fit <- function(fitType) {
    message("Trying fitType = ", fitType)
    
    result <- withCallingHandlers(
      {
        DESeq(dds, fitType = fitType)
      },
      warning = function(w) {
        message("Warning with fitType = ", fitType, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
        stop("Treating warning as failure")
      }
    )
    
    return(result)
  }
  
  for (fit in c("parametric", "local", "mean")) {
    res <- tryCatch(
      try_fit(fit),
      error = function(e) {
        message("Failed with fitType = ", fit)
        return(NULL)
      }
    )
    
    if (!is.null(res)) {
      message("Success with fitType = ", fit)
      return(list(dds = res, fitType = fit))
    }
  }
  
  stop("All fitTypes failed.")
}


# deseq2 without sequencing batch correction
deseqRun_noseqBatch<-function(countFile, metaFile, prefixFileName=""){
  suppressPackageStartupMessages(library(DESeq2))
  #load the counts
  expdata<-countFile
  
  #load the metadata
  expmeta<-metaFile
  
  #deseq2 analysis-X
  dds_exp<-DESeqDataSetFromMatrix(countData=expdata, colData=expmeta, design= ~y_count + xi_count)
  sizeFactors(dds_exp)<-expmeta$AutosomalReads/10000000
  tryDeseq<-run_deseq_with_fallback(dds_exp)
  dds_exp<-tryDeseq$dds
  norm_counts_exp<-counts(dds_exp, normalized=TRUE)
  

  if (exists("myPath")==FALSE) {
    myPath<-getwd()
  }
  
  if (file.exists(paste0(myPath,"/deseq2"))==FALSE) {
    dir.create(paste0(myPath,"/deseq2"), recursive = TRUE)
  }
  
  #write file of normalized counts to table
  expNormCounts<-counts(dds_exp, normalized=TRUE)
  write.table(expNormCounts,file=paste0(myPath,'/deseq2/', as.character(prefixFileName),"_normalizedcounts.txt"),quote=FALSE,sep="\t",row.names=TRUE)
  
 
  return(expNormCounts)
  
}


# deseq2 with sequencing batch correction
deseqRun_seqBatch<-function(countFile, metaFile, prefixFileName=""){
  suppressPackageStartupMessages(library(DESeq2))
  #load the counts
  expdata<-countFile
  
  #load the metadata
  expmeta<-metaFile
  
  #deseq2 analysis-X
  dds_exp<-DESeqDataSetFromMatrix(countData=expdata, colData=expmeta, design= ~SequencingBatch+y_count + xi_count)
  sizeFactors(dds_exp)<-expmeta$AutosomalReads/10000000
  tryDeseq<-run_deseq_with_fallback(dds_exp)
  dds_exp<-tryDeseq$dds
  norm_counts_exp<-counts(dds_exp, normalized=TRUE)
  
  
  #check if deseq2 directory exists and create it if not
  
  if (exists("myPath")==FALSE) {
    myPath<-getwd()
  }
  
  if (file.exists(paste0(myPath,"/deseq2"))==FALSE) {
    dir.create(paste0(myPath,"/deseq2"), recursive = TRUE)
  }

  #write file of normalized counts to table
  expNormCounts<-counts(dds_exp, normalized=TRUE)
  write.table(expNormCounts,file=paste0(myPath,'/deseq2/', as.character(prefixFileName),"_normalizedcounts.txt"),quote=FALSE,sep="\t",row.names=TRUE)
  

  return(expNormCounts)
  
}

# decide which deseq model to use and run it
to_run_deseq<-function(countFile1, metaFile1, prefixFileName1=""){
  if ("SequencingBatch" %in% colnames(metaFile1) & length(unique(metaFile1$SequencingBatch)) > 1) {
    message('deseq run - correction for sequencingbatch')
    a<-deseqRun_seqBatch(countFile1, metaFile1, prefixFileName1)
  }
  else {
    message('deseq run - no batch correction')
    a<-deseqRun_noseqBatch(countFile1, metaFile1, prefixFileName1)
  }
}



