####Sub-clustering of the Epithelial cell, 2021-09-30####
setwd('H:/Project/single cell/MPSC&ST/')
library(Seurat)
library(cowplot)
library(ggplot2)

alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
fib.data = alldata[,alldata$CellTypeManully %in% c('Fibroblast cell')]
dim(fib.data)#33261  2675

fib.data = NormalizeData(fib.data)
fib.data = FindVariableFeatures(fib.data, selection.method = 'mvp',verbose = T)
fib.data = ScaleData(fib.data, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                      verbose = T)

fib.data = RunPCA(fib.data, verbose = T, npcs = 50)
fib.data = RunHarmony(fib.data,dims.use = 1:50,sigma = 0.30,
                       group.by.vars = c( "orig.ident" )
)

fib.data <- RunUMAP(fib.data,reduction = "harmony", dims = 1:30)
fib.data <- FindNeighbors(fib.data, reduction = "harmony",dims = 1:30, verbose = FALSE)
names(fib.data@graphs)

for (res in c(0.1,0.3,0.5,1)) {
  fib.data <- FindClusters(fib.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, 
          DimPlot(fib.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), 
          DimPlot(fib.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), 
          DimPlot(fib.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))



sel.clust = "RNA_snn_res.0.3"
fib.data <- SetIdent(fib.data, value = sel.clust)
table(fib.data@active.ident)
DimPlot(fib.data, reduction = "umap", label = T)
VlnPlot(fib.data, features = c("percent_ribo",'percent_mito','nFeature_RNA'), 
        ncol = 2, pt.size = 0)

fib.data$orig.ident = factor(fib.data$orig.ident,
                              levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                         'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                         'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                         'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))


#fib.data@assays$RNA = alldata@assays$RNA[,colnames(fib.data)]

markers_genes <- FindAllMarkers(fib.data, logfc.threshold = 0.1, 
                                test.use = "wilcox", 
                                min.pct = 0.2, 
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_fibrCellSubClustering.csv')

markers_genes = split(markers_genes,markers_genes$cluster)

markers_genes = read.csv(file = './variables/marker_genes_fibrcellSubClustering.csv',row.names = 1)


fib.data<-RenameIdents(fib.data, '0'='COL6A3+','1'='notD','2'='MF',
                        '3'='Tcell','4'='SFRP1+','5'='Pericyte','6'='Lipofibroblast',
                       '7'='SMC','8'='EC','9'='Mesothelia')
fib.data <- subset(fib.data, idents = c('Tcell','EC'), invert = T)
# fib.data$subType = Idents(fib.data)
# table(fib.data$CellTypeManully,fib.data$subType)
fib.data$subType = Idents(fib.data)
table(fib.data$orig.ident,fib.data$subType)

marker.genes = c('COL6A3',"CDH11",
  'TAGLN','ACTA2','ACTG2','MYH11','MYLK',
  'SYNPO2','CRYAB','CNN1','DES',
  'UPK3B','MSLN','CALB2','WT1',
  'FABP4','FABP5',
  'RGS5','CSPG4','ABCC9','KCNJ8',
  'SFRP1'
)
DotPlot(fib.data, features = marker.genes, assay = "RNA",cluster.idents = T) + 
  coord_flip()+ RotatedAxis() +
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')+
  theme(text = element_text(size = 8))
DimPlot(fib.data,group.by = 'subType',reduction = 'umap',label = T,cols = pal_npg('nrc')(10))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Fibr/Dimplot cellType labeled.pdf
plot_grid(ncol = 3, 
          DimPlot(fib.data, group.by = 'subType',label = T) + NoAxes(), 
          DimPlot(fib.data, group.by = "orig.ident") + NoAxes(), 
          DimPlot(fib.data, group.by = "type") + NoAxes())
fib.data$orig.ident = factor(fib.data$orig.ident,
                              levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                         'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                         'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                         'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))
DimPlot(fib.data,group.by = 'subType',split.by = 'orig.ident',ncol=4)
####

annotation = fib.data@meta.data

annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = c('TI','TM','TI1','TI2','TM1','TM2','TS1','TS2','NI','NM','NS'))
subType.counts = table(annotation$subType)
annotation$subType = factor(annotation$subType,levels =names(subType.counts[order(subType.counts,decreasing = T)]) )
library(ggplot2)
library(RColorBrewer)
ggplot(data = annotation,aes(x=Patient,fill=type))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = brewer.pal(3,'Paired')[-1])
ggplot(data = annotation,aes(x=Patient,fill=From))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))+theme( strip.background = element_blank())
ggplot(data = annotation,aes(x=Patient))+geom_bar(fill='blue',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+theme( strip.background = element_blank())

saveRDS(fib.data,file = './variables/fib.data.Rds')
fib.data = readRDS(file = './variables/fib.data.Rds')
####
library(pheatmap)
com2orig = table(fib.data$orig.ident,fib.data$subType)
type.anno = data.frame(Patient= unlist(lapply(rownames(com2orig),function(a){
  return(strsplit(a,'_')[[1]][3])
})),
type = unlist(lapply(rownames(com2orig),function(a){
  return(substr(a,1,1))
})))
rownames(type.anno)=rownames(com2orig)
#pheatmap(cor(t(com2orig)),annotation_row = type.anno)

####PCA analysis about the celltype composition
com2orig = com2orig/apply(com2orig, 1, sum)
mat = t(com2orig)
mat = mat[apply(mat, 1, sum)>0,]
mat = scale(t(mat))
mat = mat[apply(mat, 1, sd)!=0,]
mat[is.na(mat)]
pheatmap(cor(t(mat),method = 'spearman'),clustering_distance_rows = "correlation",
         ,annotation_row = type.anno)

pcr.int = prcomp(mat)
biplot(pcr.int)
library(RColorBrewer)
pcr.cluster = data.frame(From = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][3]
})),
name = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][1]
})),
row.names  = rownames(mat),
PC1 = pcr.int$x[,1],
PC2 = pcr.int$x[,2])
pcr.cluster[is.na(pcr.cluster$From),'From']= paste('GEO', unlist(lapply(rownames(pcr.cluster[is.na(pcr.cluster$From),]), function(a){
  substr( strsplit(a,'_')[[1]][2],2,3)
})))
pcr.cluster[pcr.cluster$name == 'LUNG','name']=  unlist(lapply(rownames(pcr.cluster[pcr.cluster$name == 'LUNG',]), function(a){
  strsplit(a,'_')[[1]][2]
}))
pcr.cluster$Type = ifelse(grepl('T',rownames(pcr.cluster)),'Tumor','Normal')
library(ggrepel)
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=From,shape=Type),size = 4)+
  scale_color_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(8,'Set3')))+
  theme_bw()+geom_text_repel(aes(label = name),size=4,nudge_y = 0.1,nudge_x = 0.1)


####Cell type compostion differences####
com2orig = table(fib.data$orig.ident,fib.data$subType)
com2orig = com2orig/apply(com2orig, 1, sum)

diff.ratios.1 = data.frame(Diff = c(as.double(abs(com2orig['NI_R_P1',]-com2orig['NM_R_P1',])),
                                    as.double(abs(com2orig['NI_R_P2',]-com2orig['NM_R_P2',])),
                                    as.double(abs(com2orig['NI_R_P3',]-com2orig['NM_R_P3',])),
                                    as.double(abs((com2orig['NS_L_P4',]-com2orig['NI_L_P4',])))),
                           CellType = rep(colnames(com2orig),4),
                           Patient = rep(c('P1','P2','P3','P4'),each = ncol(com2orig)))
diff.ratios.1$Type = 'Normal'

diff.ratios = data.frame(Diff = c(as.double(abs(com2orig['TI_R_P1',]-com2orig['TM_R_P1',])),
                                  as.double(abs(com2orig['TI_R_P2',]-com2orig['TM_R_P2',])),
                                  as.double(abs((com2orig['TI1_R_P3',]+com2orig['TI2_R_P3',])*0.5-
                                                  (com2orig['TM1_R_P3',]+com2orig['TM2_R_P3',]*0.5))),
                                  as.double(abs((com2orig['TS1_L_P4',]+com2orig['TS2_L_P4',])*0.5-
                                                  com2orig['TI_L_P4',]))),
                         CellType = rep(colnames(com2orig),4),
                         Patient = rep(c('P1','P2','P3','P4'),each = ncol(com2orig)))
diff.ratios$Type = 'Tumor'
diff.ratios = rbind(diff.ratios,diff.ratios.1)
diff.ratios$CellType = factor(diff.ratios$CellType,
                              levels = names(sort(table(fib.data$subType))))
library(ggpubr)
library(dplyr)
ptest <- diff.ratios %>% group_by(CellType) %>% summarize(p.value = kruskal.test(Diff ~Type)$p.value)
#write.csv(ptest,file = './variables/ptest for immu subtype diff ration.csv')
ggbarplot(diff.ratios,
          x='CellType',y = 'Diff',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  geom_text(data = ptest,mapping = aes(x=CellType,y=0.15,label = ifelse(signif(p.value,2)<0.05,'*','')),inherit.aes = F,size=3)+
  coord_flip()+
  theme_cowplot() + 
  ylab("Diff Ratio") + 
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Fibr/subtype diffRatio barplot.pdf


####
library(reshape2)
ratio.mat = melt(com2orig,measure.vars = colnames(com2orig))
ratio.mat$Type = substr(ratio.mat$Var1,1,1)
ptest2 <- ratio.mat %>% group_by(Var2) %>% summarize(p.value = kruskal.test(value ~Type)$p.value)
ratio.mat$Var2 = factor(ratio.mat$Var2,levels = rev(levels(ratio.mat$Var2)))

ggbarplot(ratio.mat,
          x='Var2',y = 'value',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  geom_text(data = ptest2,mapping = aes(x=Var2,y=0.6,label = ifelse(signif(p.value,2)<0.05,'*','')),inherit.aes = F,size=3)+
  coord_flip()+
  theme_cowplot() + 
  ylab("Proportion") + 
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Fibr/Subtype T&N compare barplot.pdf
