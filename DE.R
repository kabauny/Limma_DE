library(limma)
parent = "/Users/zhongningchen/Data/miRNA"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")
ReadmiRNA = function(cancerType, sample, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index
  
  f = paste(parent, cancerType, sample, tumorStage, sep = "/")
  
  manif = paste(f, "MANIFEST.txt", sep = "/")
  manifest = read.table(file = manif, sep = "\t", header = TRUE)
  fileName = paste(f, toString(manifest$filename[1]), sep = "/")
  dft = read.table(file = fileName, sep = "\t", header = TRUE)
  df = data.frame(dft$read_count)
  colnames(df)[1] = paste(cancerType, tumorStage, "1", sep = "_")
  
  rownames(df) = dft$miRNA_ID
  
  for(i in 2:length(manifest$filename)){
    fileName = paste(f, toString(manifest$filename[i]), sep = "/")
    dft = read.table(file = fileName, sep = "\t", header = TRUE)
    tempSampleName = paste(cancerType, tumorStage, toString(i), sep = "_")
    df[tempSampleName] = dft$read_count
  }
  
  return(df)
}
add_Label = function(DF, n){
  #transpose dataframe to n x p
  #adds label column 
  DF = as.data.frame(t(DF))
  DF$label = rep(n, nrow(DF))
  return(DF)
}
near_zero_var = function(X, i){
  if(all(X[,i] == 0)){
    return(FALSE)
  }else{
    n = length(X[,i])
    u = as.double(length(unique(X[,i]))/n)
    
    return(((u > .1)))
  }
}
DE_voom = function(ct, a,b, p.threshold = 0.05){

  LUAD.1 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[a])
  LUAD.2 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[b])
  LUAD.2 = add_Label(LUAD.2, tumorStageList[a])
  LUAD.1 = add_Label(LUAD.1, tumorStageList[b])
  LUAD = rbind(LUAD.1, LUAD.2)
  LUAD.12 = rbind(LUAD.1, LUAD.2)
  
  X = LUAD.12[,names(LUAD.12) != 'label']
  
  NZV_gene = data.frame(matrix(TRUE, 1882, 1))
  for(i in 1:1881){
    NZV_gene[i,] = near_zero_var(X, i)
  }
  NZV = colnames(LUAD.12)[NZV_gene == TRUE]
  LUAD.12 = LUAD.12[,colnames(LUAD.12) %in% NZV]
  
  # Apply voom transformation
  LUAD.design <- model.matrix(~ LUAD.12$label)
  LUAD.12$label = NULL
  LUAD.12 = t(LUAD.12)
  LUAD.nf <- calcNormFactors(LUAD.12)
  LUAD.v <- voom(LUAD.12, LUAD.design, lib.size=colSums(LUAD.12)*LUAD.nf, 
            normalize.method="quantile", plot=TRUE)
  
  # Usual limma pipeline
  LUAD.fit.voom <- lmFit(LUAD.v, LUAD.design)
  LUAD.fit.voom <- eBayes(LUAD.fit.voom)
  
  LUAD.voom_results <- topTable(LUAD.fit.voom, coef=2,  adjust="BH", number = nrow(LUAD.12))
  LUAD.voom_results$threshold <- as.logical(LUAD.voom_results$adj.P.Val < p.threshold)
  LUAD.genes.voom <- row.names(LUAD.voom_results)[which(LUAD.voom_results$threshold)]
  
  return(LUAD.genes.voom)
}


geneList.LUAD.1.2 = DE_voom("LUAD", 1, 2)
geneList.LUAD.1.3 = DE_voom("LUAD", 1, 3)
geneList.LUAD.1.4 = DE_voom("LUAD", 1, 4)
geneList.LUAD.1.5 = DE_voom("LUAD", 1, 5)
geneList.LUAD.1.6 = DE_voom("LUAD", 1, 6)
geneList.LUAD.2.3 = DE_voom("LUAD", 2, 3)
geneList.LUAD.2.4 = DE_voom("LUAD", 2, 4)
geneList.LUAD.2.5 = DE_voom("LUAD", 2, 5)
geneList.LUAD.2.6 = DE_voom("LUAD", 2, 6)
geneList.LUAD.3.4 = DE_voom("LUAD", 3, 4)
geneList.LUAD.3.5 = DE_voom("LUAD", 3, 5)
geneList.LUAD.3.6 = DE_voom("LUAD", 3, 6)
geneList.LUAD.4.5 = DE_voom("LUAD", 4, 5)
geneList.LUAD.4.6 = DE_voom("LUAD", 4, 6)
geneList.LUAD.5.6 = DE_voom("LUAD", 5, 6)

geneList.LUSC.1.2 = DE_voom("LUSC", 1, 2, p.threshold = .1)
geneList.LUSC.1.3 = DE_voom("LUSC", 1, 3, p.threshold = .1)
geneList.LUSC.1.4 = DE_voom("LUSC", 1, 4, p.threshold = .1)
geneList.LUSC.1.5 = DE_voom("LUSC", 1, 5, p.threshold = .1)
geneList.LUSC.1.6 = DE_voom("LUSC", 1, 6, p.threshold = .1)
geneList.LUSC.2.3 = DE_voom("LUSC", 2, 3, p.threshold = .1)
geneList.LUSC.2.4 = DE_voom("LUSC", 2, 4, p.threshold = .1)
geneList.LUSC.2.5 = DE_voom("LUSC", 2, 5, p.threshold = .1)
geneList.LUSC.2.6 = DE_voom("LUSC", 2, 6, p.threshold = .1)
geneList.LUSC.3.4 = DE_voom("LUSC", 3, 4, p.threshold = .1)
geneList.LUSC.3.5 = DE_voom("LUSC", 3, 5, p.threshold = .1)
geneList.LUSC.3.6 = DE_voom("LUSC", 3, 6, p.threshold = .1)
geneList.LUSC.4.5 = DE_voom("LUSC", 4, 5, p.threshold = .1)
geneList.LUSC.4.6 = DE_voom("LUSC", 4, 6, p.threshold = .1)
geneList.LUSC.5.6 = DE_voom("LUSC", 5, 6, p.threshold = .1)
