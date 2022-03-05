rm(list = ls());gc();rm(list = ls())
Num = "S003"

library(data.table)
library(magrittr)
library(stringr)
library(dplyr)
library(DESeq2)

wd = "/home/qians/MamDC/Data/Seqdata/NicholasTony2020NatureRNAseq/"
a <- fread(file.path(wd,"GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv")) %>% as.data.frame()
info <- read.csv(file.path(wd,"GSE132040_MACA_Bulk_metadata.csv"))
info$characteristics..age %<>% as.numeric()

info$source.name %<>% lapply(.,function(x)strsplit(x,split = "_",fixed = T)[[1]][1]) %>% unlist() %>%
  gsub("Limb","Limb Muscle",.) %>% gsub("Small","Small Intestine",.)
info %<>% dplyr::filter(.,!grepl("NA",source.name))
info$sample = paste(info$source.name,info$characteristics..sex,info$characteristics..age,sep = "_") # tissue + sex + stage
row.names(a) = a$gene
colnames(a) %<>% gsub(".gencode.vM19","",.)
b = read.csv("~/MamDC/Data/Seqdata/NicholasTony2020NatureRNAseq/Mouse.gencode.gene.bed",
             header = FALSE,sep = "\t",stringsAsFactors = F)
for (Tissue in unique(info$source.name)) {
  focal.tissue = dplyr::filter(info,source.name == Tissue & characteristics..age == "18")
  focal.express = a[,focal.tissue$Sample.name]
  dds = DESeq2::DESeqDataSetFromMatrix(focal.express,
                                       colData = focal.tissue,
                                       design = ~characteristics..sex)
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  
  res = results(dds,contrast = c("characteristics..sex","m","f")) %>% as.data.frame() %>%
    #dplyr::filter(.,abs(log2FoldChange) >=0.5 & padj <= 0.05) %>% 
    dplyr::select(.,c(log2FoldChange,padj))
  res[,c("type","chr")] = b[match(row.names(res),b$V6),c(5,1)]
  res %<>% na.omit #%>% 
  res %<>% dplyr::filter(., grepl("pseudogene|protein_coding",type) & !chr %in% c("chrM","chrY"))
  res[res$type!="protein_coding","type"] = "Pseudogene"
  res[res$chr!="chrX","chr"] = "chrA"
  res$bias = "Unbiased"
  res[res$log2FoldChange>=0.5 & res$padj <= 0.05, "bias"] = "Male-biased"
  res[res$log2FoldChange<= -0.5 & res$padj <= 0.05, "bias"] = "Female-biased"
  
  table(res$bias,paste(res$type,res$chr))
  
  
wd = "/home/qians/Pseudo/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
Species = "Mouse" #"Human"
a <- fread(file.path(wd,"counts.tsv")) %>% as.data.frame()
coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
coldata$name %<>% gsub("KidneyTestis","Kidney",.)
coldata$time <- lapply(coldata$stage2,function(x){substr(x,1,str_length(x)-3)}) %>% 
  unlist() %>% as.numeric()

if (all(row.names(coldata) %in% colnames(a)[-1])) {
  a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  row.names(a) = lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  a = a[,-1]
  a = a[,row.names(coldata)]
}
b = read.csv(paste0("/home/qians/Pseudo/Data/Ref/",Species,"/","gene",Species,".bed"),
             header = F,sep = "\t",stringsAsFactors = FALSE)
b %<>% dplyr::filter(., grepl("pseudogene|protein_coding",V5) & V1 %in% c(1:100,"X"))
b[b$V1!="X","V1"]="A"
b$V5 %<>% gsub(" ","",.)
b[b$V5!="protein_coding","V5"] ="Pseudogene"
b[b$V5=="protein_coding","V5"] ="Coding"

coldata$condition %<>% gsub("Testis","Gonad",.) %>% gsub("Ovary","Gonad",.)
df = c()
for ( Tissue in unique(coldata$condition) ) {
  coldata2 = dplyr::filter(coldata, condition == Tissue )
  #coldata2$stage4 %<>% factor(.,levels = deve[Tissue,][deve[Tissue,]!=" "] )
  onetissue = a[,row.names(coldata2)] #one tissue matrix
  if (all(row.names(coldata2) == colnames(onetissue))) {
    dds = DESeq2::DESeqDataSetFromMatrix(onetissue,
                                         colData = coldata2,
                                         design = ~sex)
  }
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  
  res = results(dds,contrast = c("sex","Male","Female")) %>% as.data.frame() %>%
    dplyr::filter(.,abs(log2FoldChange) >=0.5 & padj <= 0.05) %>% 
    dplyr::select(.,c(log2FoldChange,padj))
  res[,c("type","chr")] = b[match(row.names(res),b$V4),c(5,1)]
  res %<>% na.omit #%>% dplyr::filter(., grepl("pseudogene|protein_coding",type) & chr %in% c(1:100,"X")) 
  sexbias = table(paste(res$type,res$chr,sep = "_"),res$log2FoldChange>0) %>% as.data.frame()
  sexbias$Var2 %<>% gsub("FALSE","Female-biased",.) %>% gsub("TRUE","Male-biased",.) 
  all = table(paste(b$V5,b$V1,sep="_")) %>% as.data.frame()
  sexbias$all = all[match(sexbias$Var1,all$Var1),2]
  sexbias$not = sexbias$all - sexbias$Freq
  c = data.frame(Type=unique(paste(sexbias$Var1,sexbias$Var2,sep = "_") %>% gsub("_A","",.) %>% gsub("_X","",.)))
  for (i in 1:4) {
    p = fisher.test(sexbias[c(2*i-1,2*i),c(3,5)])$p.value
    odd = fisher.test(sexbias[c(2*i-1,2*i),c(3,5)])$estimate %>% as.numeric()
    c[i,2] = p 
    c[i,3] = odd 
  }
  c$tissue = Tissue
  df = rbind(df,c)
}
df$type = lapply(df$Type,function(x)strsplit(x,split = "_")[[1]][1]) %>% unlist()
df$sex = lapply(df$Type,function(x)strsplit(x,split = "_")[[1]][2]) %>% unlist()
ggplot(df,aes())

