####Sub-clustering of the Epithelial cell, 2021-09-27####
setwd('H:/Project/single cell/MPSC&ST/')
library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony) 
library(ggsci)
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
Bcell.data = alldata[,alldata$CellTypeManully %in% c('B cell')]
dim(Bcell.data)#33261  4857
Bcell.data = NormalizeData(Bcell.data)
Bcell.data = FindVariableFeatures(Bcell.data, verbose = T)
Bcell.data = ScaleData(Bcell.data, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                      verbose = T)

Bcell.data = RunPCA(Bcell.data, verbose = T, npcs = 50)
Bcell.data = RunUMAP(Bcell.data, dims = 1:50, verbose = T)
Bcell.data = RunHarmony(Bcell.data,
                        c( "orig.ident" ))
Bcell.data <- RunUMAP(Bcell.data,reduction = "harmony", dims = 1:20)
Bcell.data <- FindNeighbors(Bcell.data, reduction = "harmony",dims = 1:20, verbose = FALSE)
names(Bcell.data@graphs)

for (res in c(0.1,0.3,0.5,1)) {
  Bcell.data <- FindClusters(Bcell.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, DimPlot(Bcell.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), DimPlot(Bcell.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), DimPlot(Bcell.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))


sel.clust = "RNA_snn_res.0.1"
Bcell.data <- SetIdent(Bcell.data, value = sel.clust)
table(Bcell.data@active.ident)
DimPlot(Bcell.data, reduction = "umap", label = T)
VlnPlot(Bcell.data, features = c("percent_ribo",'percent_mito','nFeature_RNA'), 
        ncol = 2, pt.size = 0)

Bcell.data$orig.ident = factor(Bcell.data$orig.ident,
                              levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                         'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                         'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                         'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))

DimPlot(Bcell.data,group.by = 'RNA_snn_res.0.1',split.by = 'orig.ident',ncol = 4)

#Bcell.data@assays$RNA = alldata@assays$RNA[,colnames(Bcell.data)]

markers_genes <- FindAllMarkers(Bcell.data, logfc.threshold = 0.1, 
                                test.use = "wilcox", 
                                min.pct = 0.2, 
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_BCellSubClustering.csv')
markers_genes = read.csv(file = './variables/marker_genes_BCellSubClustering.csv')
markers_genes = split(markers_genes,markers_genes$cluster)
#


#GPR183: memory B
#MZB1: Plasma B
#CD37: Mature B
FeaturePlot(Bcell.data,features = c('CD79A','CD79B','MZB1','CD37','GPR183','SFTPB'),slot = 'counts')

Bcell.data<-RenameIdents(Bcell.data, '0'='Mature B','1'='Plasma B','2'='CD2+ B',
                         '3'='Memory B','4'='Epi','5'='Plasmablast')
Bcell.data <- subset(Bcell.data, idents = c('Epi'), invert = T)
Bcell.data$subType = Idents(Bcell.data)
# table(Bcell.data$CellTypeManully,Bcell.data$subType)
marker.genes = c('GPR183','MZB1','CD37','CD2','CD19','MKI67')
DotPlot(Bcell.data, features = marker.genes, assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() + 
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')
# file:///H:/project/single cell/MPSC&ST/figure results/four samples/bcell/Dotplot markers for subtype.pdf
# 
DimPlot(Bcell.data,group.by = 'subType',reduction = 'umap',label = T,cols = pal_npg('nrc')(10))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/bcell/Dimplot subtype labeled.pdf

# plot_grid(ncol = 3, 
#           DimPlot(Bcell.data, group.by = 'subType',label = T) + NoAxes(), 
#           DimPlot(Bcell.data, group.by = "orig.ident") + NoAxes(), 
#           DimPlot(Bcell.data, group.by = "type") + NoAxes())
# DimPlot(Bcell.data,group.by = 'subType',reduction = 'umap',label = T)
# Bcell.data$orig.ident = factor(Bcell.data$orig.ident,
#                                levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
#                                           'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
#                                           'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
#                                           'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))
DimPlot(Bcell.data,group.by = 'subType',split.by = 'orig.ident',ncol=4)

####

# annotation = Bcell.data@meta.data
# 
# annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
#   return(strsplit(a,'_')[[1]][1])
# }))
# annotation$From = factor(annotation$From ,levels = c('TI','TM','TI1','TI2','TM1','TM2','TS1','TS2','NI','NM','NS'))
# library(ggplot2)
# library(RColorBrewer)
# ggplot(data = annotation,aes(x=Patient,fill=type))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
#   facet_grid( rows= vars(subType))+theme_bw()+theme( strip.background = element_blank())+
#   scale_fill_manual(values = brewer.pal(3,'Paired')[-1])
# ggplot(data = annotation,aes(x=Patient,fill=From))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
#   facet_grid( rows= vars(subType))+theme_bw()+
#   scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))+theme( strip.background = element_blank())
# ggplot(data = annotation,aes(x=Patient))+geom_bar(fill='blue',width = 0.7)+coord_flip()+
#   facet_grid( rows= vars(subType))+theme_bw()+theme( strip.background = element_blank())
# 
# alldata.p1 =  Bcell.data[, Bcell.data$orig.ident %in% c('TI_P1','NI_P1','NM_P1','TM_P1')]
# alldata.p2 =  Bcell.data[, Bcell.data$orig.ident %in% c('TI_P2','NI_P2','NM_P2','TM_P2')]
# 
# par(mfrow=c(1,2), mar = c(4, 6, 3, 1))
# 
# com2orig.1 = table(alldata.p1$orig.ident,alldata.p1$subType)[c('TI_P1','NI_P1','NM_P1','TM_P1'),]
# cor(com2orig.1)
# cor(t(com2orig.1))
# library(corrplot)
# com2orig.1=com2orig.1/apply(com2orig.1, 1, sum)
# corrplot(cor(t(com2orig.1),method='pearson'), method="pie")
# 
# com2orig.2 = table(alldata.p2$orig.ident,alldata.p2$subType)[c('TI_P2','NI_P2','NM_P2','TM_P2'),]
# com2orig.2=com2orig.2/apply(com2orig.2, 1, sum)
# corrplot(cor(t(com2orig.2),method = 'pearson'), method="pie")
# 
# ggplot(annotation,aes(x=subType,fill=From))+geom_bar(width = 0.7,position = 'dodge')+
#   facet_grid(rows = vars(Patient))+coord_flip()+theme_bw()+theme( strip.background = element_blank())+
#   scale_fill_manual(values = brewer.pal(4,'Set2'))
saveRDS(Bcell.data,file = './variables/Bcell.data.Rds')
Bcell.data = readRDS(file = './variables/Bcell.data.Rds')
####
library(pheatmap)
Bcell.data = readRDS(file = './variables/Bcell.data.Rds')
com2orig = table(Bcell.data$orig.ident,Bcell.data$subType)
type.anno = data.frame(Patient= unlist(lapply(rownames(com2orig),function(a){
  return(strsplit(a,'_')[[1]][3])
})),
type = unlist(lapply(rownames(com2orig),function(a){
  return(substr(a,1,1))
})))
rownames(type.anno)=rownames(com2orig)
pheatmap(cor(t(com2orig)),annotation_row = type.anno)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
f.data = t(f.data)

pheatmap(cor(f.data,method = 'spearman'),clustering_distance_rows = "correlation",annotation_row = type.anno)
pheatmap(cor(f.data,method = 'spearman'),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_row = type.anno)
###

annotation =Bcell.data@meta.data

annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = c('TI','TM','TI1','TI2','TM1','TM2','TS1','TS2','NI','NM','NS'))
subType.counts = table(annotation$subType)
annotation$subType = factor(annotation$subType,levels =names(subType.counts[order(subType.counts,decreasing = T)]) )
library(ggplot2)
library(RColorBrewer)


ggplot(data = annotation,aes(x=From,fill=subType))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free')+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = pal_npg('nrc')(10))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/bcell/cell type summary part4.pdf
