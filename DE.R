library(limma)
library(edgeR)
parent = "C:/Users/kabau/Desktop/Data/miRNA"
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
getDF = function(ct, mi){
  df1 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[1])[mi,]
  df2 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[2])[mi,]
  df3 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[3])[mi,]
  df4 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[4])[mi,]
  df5 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[5])[mi,]
  df6 = ReadmiRNA(ct, "PrimaryTumor", tumorStageList[6])[mi,]
  
  return(list(df.1 = df1, df.2 = df2, df.3 = df3, df.4 = df4, df.5= df5, df.6= df6))
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
{
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

geneList.LUSC.1.2 = DE_voom("LUSC", 1, 2)
geneList.LUSC.1.3 = DE_voom("LUSC", 1, 3)
geneList.LUSC.1.4 = DE_voom("LUSC", 1, 4)
geneList.LUSC.1.5 = DE_voom("LUSC", 1, 5)
geneList.LUSC.1.6 = DE_voom("LUSC", 1, 6)
geneList.LUSC.2.3 = DE_voom("LUSC", 2, 3)
geneList.LUSC.2.4 = DE_voom("LUSC", 2, 4)
geneList.LUSC.2.5 = DE_voom("LUSC", 2, 5)
geneList.LUSC.2.6 = DE_voom("LUSC", 2, 6)
geneList.LUSC.3.4 = DE_voom("LUSC", 3, 4)
geneList.LUSC.3.5 = DE_voom("LUSC", 3, 5)
geneList.LUSC.3.6 = DE_voom("LUSC", 3, 6)
geneList.LUSC.4.5 = DE_voom("LUSC", 4, 5)
geneList.LUSC.4.6 = DE_voom("LUSC", 4, 6)
geneList.LUSC.5.6 = DE_voom("LUSC", 5, 6)
}
{
  parent = "C:/Users/kabau/Desktop/Data/miRNA/Limma_to_ENSG/LUAD"
  library(biomaRt)
  mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  desired_attributes = c("mirbase_id", "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")
  ensg = getBM(attributes= desired_attributes, filters = 'mirbase_id', 
               values = geneList.LUAD.4.6, mart = mart)
  write.table(ensg, paste(parent, "LUAD_46.txt", sep = "/"))
}
{
parent = "C:/Users/kabau/Desktop/Data/miRNA/Limma_results"

write.table(geneList.LUAD.1.2, paste(parent, "LUAD", "LUAD_12.txt", sep = "/"))
write.table(geneList.LUAD.1.3, paste(parent, "LUAD", "LUAD_13.txt", sep = "/"))
write.table(geneList.LUAD.1.4, paste(parent, "LUAD", "LUAD_14.txt", sep = "/"))
write.table(geneList.LUAD.1.5, paste(parent, "LUAD", "LUAD_15.txt", sep = "/"))
write.table(geneList.LUAD.1.6, paste(parent, "LUAD", "LUAD_16.txt", sep = "/"))

write.table(geneList.LUAD.2.3, paste(parent, "LUAD", "LUAD_23.txt", sep = "/"))
write.table(geneList.LUAD.2.4, paste(parent, "LUAD", "LUAD_24.txt", sep = "/"))
write.table(geneList.LUAD.2.5, paste(parent, "LUAD", "LUAD_25.txt", sep = "/"))
write.table(geneList.LUAD.2.6, paste(parent, "LUAD", "LUAD_26.txt", sep = "/"))

write.table(geneList.LUAD.3.4, paste(parent, "LUAD", "LUAD_34.txt", sep = "/"))
write.table(geneList.LUAD.3.5, paste(parent, "LUAD", "LUAD_35.txt", sep = "/"))
write.table(geneList.LUAD.3.6, paste(parent, "LUAD", "LUAD_36.txt", sep = "/"))

write.table(geneList.LUAD.4.5, paste(parent, "LUAD", "LUAD_45.txt", sep = "/"))
write.table(geneList.LUAD.4.6, paste(parent, "LUAD", "LUAD_46.txt", sep = "/"))

write.table(geneList.LUAD.5.6, paste(parent, "LUAD", "LUAD_56.txt", sep = "/"))

write.table(geneList.LUSC.1.2, paste(parent, "LUSC", "LUSC_12.txt", sep = "/"))
write.table(geneList.LUSC.1.3, paste(parent, "LUSC", "LUSC_13.txt", sep = "/"))
write.table(geneList.LUSC.1.4, paste(parent, "LUSC", "LUSC_14.txt", sep = "/"))
write.table(geneList.LUSC.1.5, paste(parent, "LUSC", "LUSC_15.txt", sep = "/"))
write.table(geneList.LUSC.1.6, paste(parent, "LUSC", "LUSC_16.txt", sep = "/"))

write.table(geneList.LUSC.2.3, paste(parent, "LUSC", "LUSC_23.txt", sep = "/"))
write.table(geneList.LUSC.2.4, paste(parent, "LUSC", "LUSC_24.txt", sep = "/"))
write.table(geneList.LUSC.2.5, paste(parent, "LUSC", "LUSC_25.txt", sep = "/"))
write.table(geneList.LUSC.2.6, paste(parent, "LUSC", "LUSC_26.txt", sep = "/"))

write.table(geneList.LUSC.3.4, paste(parent, "LUSC", "LUSC_34.txt", sep = "/"))
write.table(geneList.LUSC.3.5, paste(parent, "LUSC", "LUSC_35.txt", sep = "/"))
write.table(geneList.LUSC.3.6, paste(parent, "LUSC", "LUSC_36.txt", sep = "/"))

write.table(geneList.LUSC.4.5, paste(parent, "LUSC", "LUSC_45.txt", sep = "/"))
write.table(geneList.LUSC.4.6, paste(parent, "LUSC", "LUSC_46.txt", sep = "/"))

write.table(geneList.LUSC.5.6, paste(parent, "LUSC", "LUSC_56.txt", sep = "/"))

}
{
common.12.16 = Reduce(intersect, list(geneList.LUAD.1.2, geneList.LUAD.1.4, geneList.LUAD.1.5, geneList.LUAD.1.6))
common.12.15 = Reduce(intersect, list(geneList.LUAD.1.2, geneList.LUAD.1.4, geneList.LUAD.1.5))
common.12.14 = Reduce(intersect, list(geneList.LUAD.1.2, geneList.LUAD.1.4))
common.14.16 = Reduce(intersect, list(geneList.LUAD.1.4, geneList.LUAD.1.5, geneList.LUAD.1.6))
common.14.15 = Reduce(intersect, list(geneList.LUAD.1.4, geneList.LUAD.1.5))
common.15.16 = Reduce(intersect, list(geneList.LUAD.1.5, geneList.LUAD.1.6))
}

miRNA = getDF(ct = "LUAD", mi = common.12.16)
{
  m = c(apply(miRNA$df.1, 1, mean), 
        apply(miRNA$df.2, 1, mean),
        apply(miRNA$df.3, 1, mean),
        apply(miRNA$df.4, 1, mean),
        apply(miRNA$df.5, 1, mean),
        apply(miRNA$df.6, 1, mean))
  
}

