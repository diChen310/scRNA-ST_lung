####Downstream analysis for the inferCNV results,2021-10-11####
####Install: manually install JAGS-4.3.0.exe, 
####BiocManager::install('rjags') ##do not compile for the pop-up window
####BiocManager::install('infercnv')
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library("RColorBrewer")
setwd('H:/project/single cell/MPSC&ST/')
rm(list = ls())
epi.data = readRDS(file = './variables/epi.data.V2.Rds')
####Manually change the patients, TIM_R_P1, TIM_R_P2, TIM_R_P3, TIS_L_P4
infercnv_obj = readRDS("./documents/inferCNV/res3/TIS_L_P4_refT/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$`T&NK cell`
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc.ti <- test_loc$TI_L_P4
test_loc.tm = c(test_loc$TS1_L_P4,test_loc$TS2_L_P4)
test_loc = c(test_loc.ti,test_loc.tm)
anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc))),
  from=c(rep("T&NK",length(normal_loc)),rep('TI',length(test_loc.ti)),
         rep('TM',length(test_loc.tm)))
)
head(anno.df)

####Clustering
set.seed(20210418)
epi.data.subtype = data.frame(sample = colnames(epi.data),subtype = epi.data$subType,
                              row.names = colnames(epi.data))
kmeans.result <- kmeans(t(expr), 7)#7 clusters, 
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #merge
kmeans_df_s=arrange(kmeans_df,kmeans_class) #sorting
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #clusters
kmeans_df_s$subtype = as.character(epi.data.subtype[rownames(kmeans_df_s),'subtype'])
kmeans_df_s[is.na(kmeans_df_s$subtype),'subtype']='Not'
head(kmeans_df_s)
table(kmeans_df_s$kmeans_class,kmeans_df_s$class)
####Define heatmap and colors
subtype.cols = pal_npg('nrc')(length(unique(epi.data$subType))+1)
names(subtype.cols)=c('AT2','CXCL14+','Ciliated',
                       'Club','Basal','AT1','Not')

top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #number of clusters
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,
                           col=list(class=c("test"="red","normal" = "blue"),
                                    kmeans_class=color_v,
                                    from=c('T&NK'="blue","TI"=brewer.pal(4,'Set2')[1],
                                           'TM'=brewer.pal(4,'Set2')[2]),
                                    subtype = subtype.cols))

#plot
pdf("./documents/inferCNV/res3/TIS_L_P4_refT/try1-v2.pdf",width = 16,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(infercnv_obj@gene_order$chr, paste("chr",1:22,sep = "")), #
             column_gap = unit(2, "mm"),
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()
##Manually revised the class ids for different patients
kmeans_df_s$cellMN = ifelse(kmeans_df_s$kmeans_class %in% c(3),'Normal','Malignant' )
kmeans_df_s[kmeans_df_s$class == 'normal','cellMN'] = 'Normal'

write.table(kmeans_df_s, file = "./documents/inferCNV/res3/TIS_L_P4_refT/kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

kmeans_df_s = read.delim(file = "./documents/inferCNV/res3/TIS_L_P4_refT/kmeans_df_s.txt",row.names = 1)
####Only display the cells of tumor tissue sample 
pdf("./documents/inferCNV/res3/TIS_L_P4_refT/try2-v2.pdf",width = 15,height = 10)
kmeans_df_s_malignant = kmeans_df_s[kmeans_df_s$class=='test',]
set.seed(20210418)
kmeans.result2 <- kmeans(t(expr)[rownames(kmeans_df_s_malignant),], 4)
kmeans_df_s_malignant$km_class = kmeans.result2$cluster
kmeans_df_s_malignant$from = substr(rownames(kmeans_df_s_malignant),1,3)
kmeans_df_s_malignant = arrange(kmeans_df_s_malignant,from,cellMN,km_class,subtype)
left_anno <- rowAnnotation(df = kmeans_df_s_malignant[,c('from','cellMN','km_class','subtype')],
                           col=list(cellMN=c("Malignant"="red","Normal" = "blue"),
                                    km_class=color_v,
                                    from=c("TI-"=brewer.pal(4,'Set2')[1],
                                           'TS1'=brewer.pal(4,'Set2')[2],
                                           "TS2"=brewer.pal(4,'Set2')[3]),
                                    subtype = subtype.cols))

ht = Heatmap(t(expr)[rownames(kmeans_df_s_malignant),], 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(infercnv_obj@gene_order$chr, paste("chr",1:22,sep = "")), #
             column_gap = unit(2, "mm"),
             
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()
write.table(kmeans_df_s_malignant, file = "./documents/inferCNV/res3/TIS_L_P4_refT/kmeans_df_s_malignant.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

####Check on the marker genes and cnv genes.2021-10-20####
library(dplyr)
cosmic.info <- read.csv(file='E:/data/COSMIC/cancer_gene_census-v85 - noFusion.csv')
cosmic.info <- select(cosmic.info,Gene.Symbol,Role.in.Cancer,Tumour.Types.Somatic.)
cosmic.gene = as.character(cosmic.info$Gene.Symbol)
patients.markers = readRDS(file='variables/patients.diff.markers.Rds')

####P1
res.markers.p1=patients.markers[[1]]
res.markers.p1 = res.markers.p1[res.markers.p1$p_val_adj<0.05 & res.markers.p1$cellType == 'Epithelial cell',]
p1.infercnv = read.delim(file='documents/inferCNV/res3/TIM_R_P1_refT/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_genes.dat')
p1.infercnv.sig = p1.infercnv[p1.infercnv$chr=='chr3' & p1.infercnv$state == 4,]
p1.infercnv.ti = p1.infercnv.sig[grepl('TI',p1.infercnv.sig$cell_group_name),]
p1.infercnv.tm = p1.infercnv.sig[grepl('TM',p1.infercnv.sig$cell_group_name),]
intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))

ti.spc.cnas = as.character(p1.infercnv.ti$gene)[as.character(p1.infercnv.ti$gene) %in% as.character(p1.infercnv.tm$gene) == F]
markers.p1.cnas = res.markers.p1[res.markers.p1$Gene %in% ti.spc.cnas & res.markers.p1$Compare=='TIVSNI',]
p1.infercnv.ti$Tag = ifelse(as.character(p1.infercnv.ti$gene) %in% markers.p1.cnas$Gene,'Tag','N')
hist(markers.p1.cnas$avg_logFC,breaks = 50)
ti.spc.paths = enricher(ti.spc.cnas,pvalueCutoff = 1.2,
                        TERM2GENE = path2Gene[,c(1,2)])@result

write.csv(ti.spc.paths,file = 'variables/P1 TI spc cnas pathways.csv')
write.csv(ti.spc.cnas,file = 'variables/P1 TI spc cnas.csv')

cosmic.p1 = cosmic.info[cosmic.info$Gene.Symbol %in% ti.spc.cnas,]
table(cosmic.p1$Role.in.Cancer)
#"RPN1"   "CNBP"   "STAG1"  "PIK3CB" "ATR"    "WWTR1"  "GMPS"   "MLF1"  
####P2
res.markers.p2=patients.markers[[2]]
res.markers.p2 = res.markers.p2[res.markers.p2$p_val_adj<0.05 & res.markers.p2$cellType == 'Epithelial cell',]
p2.infercnv = read.delim(file='documents/inferCNV/res3/TIM_R_P2_refT/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_genes.dat')

p2.infercnv$tag = paste(p2.infercnv$chr,p2.infercnv$state,p2.infercnv$gene)
dup.p2 = p2.infercnv[duplicated(p2.infercnv$tag),]
uni.p2= p2.infercnv[p2.infercnv$tag %in% as.character(dup.p2$tag)==F,]
table(uni.p2$cell_group_name, uni.p2$chr)

p2.infercnv.sig = p2.infercnv[p2.infercnv$chr=='chr1' & p2.infercnv$state >= 4,]
p2.infercnv.ti = p2.infercnv.sig[grepl('TI',p2.infercnv.sig$cell_group_name),]
p2.infercnv.tm = p2.infercnv.sig[grepl('TM',p2.infercnv.sig$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p2.infercnv.ti = p2.infercnv.ti[as.character(p2.infercnv.ti$gene) %in% as.character(p2.infercnv.tm$gene) == F, ]
markers.p2.cnas = res.markers.p2[res.markers.p2$Gene %in% as.character(p2.infercnv.ti$gene) & res.markers.p2$Compare=='TIVSNI',]
markers.p2.cnas.tm = res.markers.p2[res.markers.p2$Gene %in% as.character(p2.infercnv.ti$gene) & res.markers.p2$Compare=='TMVSNM',]

hist(markers.p2.cnas$avg_logFC,breaks = 50)
hist(markers.p2.cnas.tm$avg_logFC,breaks = 50)
boxplot(as.double(markers.p2.cnas$avg_logFC),as.double(markers.p2.cnas.tm$avg_logFC))

p2.infercnv.sig2 = p2.infercnv[p2.infercnv$chr=='chr1' & p2.infercnv$state <= 2,]
p2.infercnv.ti2 = p2.infercnv.sig2[grepl('TI',p2.infercnv.sig2$cell_group_name),]
p2.infercnv.tm2 = p2.infercnv.sig2[grepl('TM',p2.infercnv.sig2$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p2.infercnv.ti2 = p2.infercnv.ti2[as.character(p2.infercnv.ti2$gene) %in% as.character(p2.infercnv.tm2$gene) == F, ]
markers.p2.cnds = res.markers.p2[res.markers.p2$Gene %in% as.character(p2.infercnv.ti2$gene) & res.markers.p2$Compare=='TIVSNI',]
hist(markers.p2.cnds$avg_logFC,breaks = 50)

cosmic.p2.a = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p2.infercnv.ti$gene),]
table(cosmic.p2.a$Role.in.Cancer)
cosmic.p2.a$Gene.Symbol
# AKT3     ARNT     CAMTA1   CASP9    CSF3R    ELF3     FH       H3F3A    JAK1     LCK      MTOR     PBX1     PRCC    
# [14] PRDM16   PRDM2    SDHB     SDHC     SFPQ     SKI      SPEN     THRAP3   TNFRSF14

cosmic.p2.d = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p2.infercnv.ti2$gene),]
table(cosmic.p2.d$Role.in.Cancer)
cosmic.p2.d$Gene.Symbol
cosmic.uni.p2 = uni.p2[uni.p2$gene %in% cosmic.info$Gene.Symbol,]
#ATP1A1 FAM46C NOTCH2 NRAS   RBM15  TRIM33

ti2.spc.paths = enricher(as.character(p2.infercnv.ti$gene),pvalueCutoff = 1.2,
                        TERM2GENE = path2Gene[,c(1,2)])@result
ti2.spc.paths.d = enricher(as.character(p2.infercnv.ti2$gene),pvalueCutoff = 1.2,
                         TERM2GENE = path2Gene[,c(1,2)])@result

####P3
res.markers.p3=patients.markers[[3]]
res.markers.p3 = res.markers.p3[res.markers.p3$p_val_adj<0.05 & res.markers.p3$cellType == 'Epithelial cell',]
p3.infercnv = read.delim(file='documents/inferCNV/res3/TIM_R_P3_refT/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_genes.dat')
p3.infercnv.sig = p3.infercnv[p3.infercnv$chr %in% c('chr5','chr17') & p3.infercnv$state == 4,]
p3.infercnv.ti = p3.infercnv.sig[grepl('TI1',p3.infercnv.sig$cell_group_name),]
p3.infercnv.to = p3.infercnv.sig[!grepl('TI1',p3.infercnv.sig$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p3.infercnv.ti = p3.infercnv.ti[as.character(p3.infercnv.ti$gene) %in% as.character(p3.infercnv.to$gene) == F, ]
markers.p3.cnas = res.markers.p3[res.markers.p3$Gene %in% as.character(p3.infercnv.ti$gene) & res.markers.p3$Compare=='TIVSNI',]
markers.p3.cnas.to = res.markers.p3[res.markers.p3$Gene %in% as.character(p3.infercnv.ti$gene) & res.markers.p3$Compare=='TMVSNM',]

hist(markers.p3.cnas$avg_logFC,breaks = 50)
hist(markers.p3.cnas.to$avg_logFC,breaks = 50)
boxplot(as.double(markers.p3.cnas$avg_logFC),as.double(markers.p3.cnas.to$avg_logFC))

p3.infercnv.sig2 = p3.infercnv[p3.infercnv$chr=='chr7' & p3.infercnv$state == 2,]
p3.infercnv.ti2 = p3.infercnv.sig2[grepl('TI1',p3.infercnv.sig2$cell_group_name),]
p3.infercnv.tm2 = p3.infercnv.sig2[!grepl('TI1',p3.infercnv.sig2$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p3.infercnv.ti2 = p3.infercnv.ti2[as.character(p3.infercnv.ti2$gene) %in% as.character(p3.infercnv.tm2$gene) == F, ]
markers.p3.cnds = res.markers.p3[res.markers.p3$Gene %in% as.character(p3.infercnv.ti2$gene) & res.markers.p3$Compare=='TIVSNI',]
hist(markers.p3.cnds$avg_logFC,breaks = 50)

cosmic.p3.a = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p3.infercnv.ti$gene),]
table(cosmic.p3.a$Role.in.Cancer)
cosmic.p3.a$Gene.Symbol
#ARHGAP26 ASPSCR1  CANT1    CD74     CLTC     HLF      MSI2     RNF213   RNF43    SRSF2 

cosmic.p3.d = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p3.infercnv.ti2$gene),]
table(cosmic.p3.d$Role.in.Cancer)
cosmic.p3.d$Gene.Symbol
#BRAF    CREB3L2 EZH2    KMT2C   TRIM24 

ti3.spc.paths = enricher(as.character(p3.infercnv.ti$gene),pvalueCutoff = 1.2,
                         TERM2GENE = path2Gene[,c(1,2)])@result
ti3.spc.paths.d = enricher(as.character(p3.infercnv.ti2$gene),pvalueCutoff = 1.2,
                           TERM2GENE = path2Gene[,c(1,2)])@result

####P4
res.markers.p4=patients.markers[[4]]
res.markers.p4 = res.markers.p4[res.markers.p4$p_val_adj<0.05 & res.markers.p4$cellType == 'Epithelial cell',]
p4.infercnv = read.delim(file='documents/inferCNV/res3/TIS_L_P4_refT/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_genes.dat')
p4.infercnv$tag = paste(p4.infercnv$chr,p4.infercnv$state,p4.infercnv$gene)
dup.p4 = p4.infercnv[duplicated(p4.infercnv$tag),]
uni.p4= p4.infercnv[p4.infercnv$tag %in% as.character(dup.p4$tag)==F,]
table(uni.p4$cell_group_name, uni.p4$chr)
p4.infercnv.sig = p4.infercnv[p4.infercnv$chr %in% c('chr1','chr1','chr16','chr17') & p4.infercnv$state == 4,]
p4.infercnv.ti = p4.infercnv.sig[grepl('TI',p4.infercnv.sig$cell_group_name),]
p4.infercnv.to = p4.infercnv.sig[!grepl('TI',p4.infercnv.sig$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p4.infercnv.ti = p4.infercnv.ti[as.character(p4.infercnv.ti$gene) %in% as.character(p4.infercnv.to$gene) == F, ]
markers.p4.cnas = res.markers.p4[res.markers.p4$Gene %in% as.character(p4.infercnv.ti$gene) & res.markers.p4$Compare=='TIVSNI',]
markers.p4.cnas.to = res.markers.p4[res.markers.p4$Gene %in% as.character(p4.infercnv.ti$gene) & res.markers.p4$Compare=='TSVSNS',]

hist(markers.p4.cnas$avg_logFC,breaks = 50)
hist(markers.p4.cnas.to$avg_logFC,breaks = 50)
boxplot(as.double(markers.p3.cnas$avg_logFC),as.double(markers.p4.cnas.to$avg_logFC))

p4.infercnv.sig2 = p4.infercnv[p4.infercnv$chr %in% c('chr5','chr12') & p4.infercnv$state == 2,]
p4.infercnv.ti2 = p4.infercnv.sig2[grepl('TI',p4.infercnv.sig2$cell_group_name),]
p4.infercnv.ts = p4.infercnv.sig2[!grepl('TI',p4.infercnv.sig2$cell_group_name),]
#intersect(as.character(p1.infercnv.ti$gene),as.character(p1.infercnv.tm$gene))
p4.infercnv.ts = p4.infercnv.ts[as.character(p4.infercnv.ts$gene) %in% as.character(p4.infercnv.ti2$gene) == F, ]
markers.p4.cnds = res.markers.p4[res.markers.p4$Gene %in% as.character(p4.infercnv.ts$gene) & res.markers.p4$Compare=='TSVSNS',]
hist(markers.p4.cnds$avg_logFC,breaks = 50)

cosmic.p4.a = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p4.infercnv.ti$gene),]
table(cosmic.p4.a$Role.in.Cancer)
cosmic.p4.a$Gene.Symbol
#CIITA  CLTC   COL1A1 CREBBP ERCC4  H3F3A  HLF    KAT7   LMNA   MSI2   MUC1   RNF43  SOCS1

cosmic.p4.d = cosmic.info[cosmic.info$Gene.Symbol %in% as.character(p4.infercnv.ts$gene),]
table(cosmic.p4.d$Role.in.Cancer)
cosmic.p4.d$Gene.Symbol
#ALDH2  BCL7A  BTG1   CCND2  CDKN1B CHD4   CHST11 CLIP1  ERC1   ETV6   KDM5A  NCOR2  PTPN11 PTPN6  USP44  ZCCHC8
#[17] ZNF384

ti4.spc.paths = enricher(as.character(p4.infercnv.ti$gene),pvalueCutoff = 1.2,
                         TERM2GENE = path2Gene[,c(1,2)])@result
ts4.spc.paths.d = enricher(as.character(p4.infercnv.ts$gene),pvalueCutoff = 1.2,
                           TERM2GENE = path2Gene[,c(1,2)])@result
save.image('inferCNV&markerGene.RData')
