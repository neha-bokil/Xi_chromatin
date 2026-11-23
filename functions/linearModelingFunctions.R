#####functions for Linear Modeling#####
suppressPackageStartupMessages(library(broom))


#linear modeling - no sequencing batch corrections
lmanalysis_noseqBatch<-function(countTable, resTable, metaTable){
  region_res<-character()
  gene_res<-character()
  for (gene in rownames(countTable)){
    options(scipen=0, digits=5)
    
    geneCounts_experimental<-t(countTable[gene,])
    metaTable$xi_count<-metaTable$x_count -1
    
    exp_Data<-cbind(metaTable, geneCounts_experimental)
    
    colnames(exp_Data)[ncol(exp_Data)]<-'Counts'
    
    #linear regression 
    
    #formulas
    formula_lm<-formula(Counts ~ xi_count + y_count)
    formula_int<-formula(Counts ~ xi_count + y_count)
    
    #linear modeling 
    linmodelExp<-lm(formula_lm, data=exp_Data)
    linmodelExp_int<-lm(formula_int, data=exp_Data)
    
    #find intercept, beta values, and pvalues

    intercept_his<-summary(linmodelExp_int)$coefficients['(Intercept)', 1]
    interr_his<-summary(linmodelExp_int)$coefficients['(Intercept)', 2]
    xbeta_his<-summary(linmodelExp)$coefficients['xi_count',1]
    xpval_his<-summary(linmodelExp)$coefficients['xi_count',4] 
    xerr_his<-summary(linmodelExp)$coefficients['xi_count',2]
    ybeta_his<-summary(linmodelExp)$coefficients['y_count', 1]
    ypval_his<-summary(linmodelExp)$coefficients['y_count',4]
    yerr_his<-summary(linmodelExp)$coefficients['y_count',2]
    adj_rsq_his<-glance(linmodelExp)$adj.r.squared
    pval_his<-glance(linmodelExp)$p.value
    
    genename<-sub("_.*", "", gene)
    # label genes with regions
    if (genename %in% NPX_genes){
      region<-'NPX'
    }
    if (genename %in% PAR1_genes) {
      region<-'PAR1'
    }
    if (genename %in% PAR2_genes) {
      region<-'PAR2'
    }
    if (genename %in% NPY_genes) {
      region<-'NPY'
    }
    if (genename %in% sexChrom_genes == FALSE) {
      region<-NA
    }
    #build results into table 
    
    resTable[gene,]<-c(intercept_his, interr_his, xbeta_his, xpval_his, xerr_his, ybeta_his, ypval_his, yerr_his, adj_rsq_his, pval_his)
    region_res<-c(region_res, region)
    gene_res<-c(gene_res, genename)
  }
  
  resTable<-cbind(resTable, region_res, gene_res)
  
  
  #correct for multiple hypothesis testing
  bh_adj_x_his<-p.adjust(resTable$xpval_his, method='BH')
  bh_adj_y_his<-p.adjust(resTable$ypval_his, method="BH")
  bh_adj_model_his<-p.adjust(resTable$pval_his, method='BH')
  
  res_final<-resTable
  colnames(res_final)<-c('intercept_his', 'intercept_his_err', 'x_his_coeff','x_his_pval', 'x_his_coeff_err','y_his_coeff', 'y_his_pval', 'y_his_coeff_err', 'adj_rsq_his', 'his_regression_pval','region', 'Gene')
  
  res_final$gene<-rownames(resTable)
  
  res_final$x_his_pval_adj<-bh_adj_x_his
  res_final$y_his_pval_adj<-bh_adj_y_his
  
  res_final$his_regression_pval_adj<-bh_adj_model_his
  
  col_order<-c('intercept_his', 'intercept_his_err', 'x_his_coeff', 'x_his_coeff_err','x_his_pval', 'x_his_pval_adj',
               'y_his_coeff', 'y_his_coeff_err','y_his_pval', 'y_his_pval_adj', 'adj_rsq_his', 'his_regression_pval', 'his_regression_pval_adj','region', 'Gene')
  res_final<-res_final[,col_order]
  return(res_final)
}


#linear modeling - sequencing batch correction
lmanalysis_seqBatch<-function(countTable, resTable, metaTable){
  region_res<-character()
  gene_res<-character()
  for (gene in rownames(countTable)){
    options(scipen=0, digits=5)
    
    geneCounts_experimental<-t(countTable[gene,])
    metaTable$xi_count<-metaTable$x_count -1
    
    exp_Data<-cbind(metaTable, geneCounts_experimental)
    
    colnames(exp_Data)[ncol(exp_Data)]<-'Counts'
    #linear regression 
    
    #formulas
    formula_lm<-formula(Counts ~ xi_count + y_count + SequencingBatch)
    formula_int<-formula(Counts ~ xi_count + y_count)
    
    #linear modeling 
    linmodelExp<-lm(formula_lm, data=exp_Data)
    linmodelExp_int<-lm(formula_int, data=exp_Data)
    
    #find intercept, beta values, and pvalues

    intercept_his<-summary(linmodelExp_int)$coefficients['(Intercept)', 1]
    interr_his<-summary(linmodelExp_int)$coefficients['(Intercept)', 2]
    xbeta_his<-summary(linmodelExp)$coefficients['xi_count',1]
    xpval_his<-summary(linmodelExp)$coefficients['xi_count',4] 
    xerr_his<-summary(linmodelExp)$coefficients['xi_count',2]
    ybeta_his<-summary(linmodelExp)$coefficients['y_count', 1]
    ypval_his<-summary(linmodelExp)$coefficients['y_count',4]
    yerr_his<-summary(linmodelExp)$coefficients['y_count',2]
    adj_rsq_his<-glance(linmodelExp)$adj.r.squared
    pval_his<-glance(linmodelExp)$p.value
    
    genename<-sub("_.*", "", gene)
    #label genes with regions
    if (genename %in% NPX_genes){
      region<-'NPX'
    }
    if (genename %in% PAR1_genes) {
      region<-'PAR1'
    }
    if (genename %in% PAR2_genes) {
      region<-'PAR2'
    }
    if (genename %in% NPY_genes) {
      region<-'NPY'
    }
    if (genename %in% sexChrom_genes == FALSE) {
      region<-NA
    }
    #build results into table 
    
    resTable[gene,]<-c(intercept_his, interr_his, xbeta_his, xpval_his, xerr_his, ybeta_his, ypval_his, yerr_his, adj_rsq_his, pval_his)
    region_res<-c(region_res, region)
    gene_res<-c(gene_res, genename)
  }
  
  resTable<-cbind(resTable, region_res, gene_res)
  
  
  #correct for multiple hypothesis testing
  bh_adj_x_his<-p.adjust(resTable$xpval_his, method='BH')
  bh_adj_y_his<-p.adjust(resTable$ypval_his, method="BH")
  bh_adj_model_his<-p.adjust(resTable$pval_his, method='BH')
  
  res_final<-resTable
  colnames(res_final)<-c('intercept_his', 'intercept_his_err', 'x_his_coeff','x_his_pval', 'x_his_coeff_err','y_his_coeff', 'y_his_pval', 'y_his_coeff_err', 'adj_rsq_his', 'his_regression_pval','region', 'Gene')
  
  res_final$gene<-rownames(resTable)
  
  res_final$x_his_pval_adj<-bh_adj_x_his
  res_final$y_his_pval_adj<-bh_adj_y_his
  
  res_final$his_regression_pval_adj<-bh_adj_model_his
  
  col_order<-c('intercept_his', 'intercept_his_err', 'x_his_coeff', 'x_his_coeff_err','x_his_pval', 'x_his_pval_adj',
               'y_his_coeff', 'y_his_coeff_err','y_his_pval', 'y_his_pval_adj', 'adj_rsq_his', 'his_regression_pval', 'his_regression_pval_adj','region', 'Gene')
  res_final<-res_final[,col_order]
  return(res_final)
}


#use this to actually run the linear modeling--it will decide whether you can correct for seqBatch or not
to_run_lmanalysis<-function(countTable1, metaTable1, filePrefix, selectXY=FALSE){
  # ensure samples are in the correct order
  commonSamps<-intersect(rownames(metaTable1), colnames(countTable1))
  metaTable1<-metaTable1[commonSamps,]
  countTable1<-countTable1[,commonSamps]
  resTable1<-data.frame(intercept_his=numeric(), interr_his=numeric(), xbeta_his=numeric(), xpval_his=numeric(), xerr_his=numeric(), 
                        ybeta_his=numeric(), ypval_his=numeric(), yerr_his=numeric(), adj_rsq_his=numeric(), pval_his=numeric())
  
  #extract only X and Y regions
  if (selectXY==TRUE) {
    countTable1<-countTable1[gsub('_.*', '', rownames(countTable1)) %in% sexChrom_genes,]
  }
  if ("SequencingBatch" %in% colnames(metaTable1) & length(unique(metaTable1$SequencingBatch))>1) {
    a<-lmanalysis_seqBatch(countTable1, resTable1, metaTable1)
    message('lm analysis run - correction for sequencingbatch')
  }
  else {
    a<-lmanalysis_noseqBatch(countTable1, resTable1, metaTable1)
    message('lm analysis run - only one sequencing batch')
  }

  return(a)
}

#compute the scaling factor to downsample background for further analyses
compute_scaling_factor<-function(signal, background, noOverlapPeaks,bgThreshold=0.05) {
  
  
  #use only fragments that have positive values for the signal and above the threshold for background
  fragments_signalPos<-rownames(signal[signal$intercept_his>0,])
  fragments_bgPos<-rownames(background[background$intercept_his>bgThreshold,])
  fragments_bothPos<-fragments_signalPos[fragments_signalPos %in% fragments_bgPos]
  
  #keep only those fragments that do not overlap peaks
  keepFrags<-noOverlapPeaks$V4
  fragments_bothPos<-fragments_bothPos[fragments_bothPos %in% keepFrags]
  
  signal_bothPos<-signal[rownames(signal) %in% fragments_bothPos,]
  bg_bothPos<-background[rownames(background) %in% fragments_bothPos,]
  
  ordrow<-match(rownames(signal_bothPos), rownames(bg_bothPos))
  
  signal_bothPos$b0Ratio<-signal_bothPos$intercept_his/bg_bothPos[ordrow,]$intercept_his
  
  Q1 <- quantile(signal_bothPos$b0Ratio, probs = 0.25)
  
  return(list(signal_df=signal_bothPos, scalingFactor=Q1))
  
}
