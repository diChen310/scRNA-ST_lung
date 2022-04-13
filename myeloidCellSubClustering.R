####Sub-clustering of the myeloid cell, 2021-09-29####
setwd('H:/Project/single cell/MPSC&ST/')
library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony)
library(ggsci)
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
mcell.data = alldata[,alldata$CellTypeManully %in% c('Myeloid cell')]

dim(mcell.data)#33261  38174
rm(alldata)
mcell.data = NormalizeData(mcell.data)
mcell.data = FindVariableFeatures(mcell.data, selection.method = 'mvp',verbose = T)
mcell.data = ScaleData(mcell.data, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                       verbose = T)

mcell.data = RunPCA(mcell.data, verbose = T, npcs = 50)
mcell.data = RunHarmony(mcell.data,dims.use = 1:50,sigma = 0.30,
                        group.by.vars = c( "orig.ident" )
)

mcell.data <- RunUMAP(mcell.data,reduction = "harmony", dims = 1:30)
mcell.data <- FindNeighbors(mcell.data, reduction = "harmony",dims = 1:30, verbose = FALSE)
names(mcell.data@graphs)

for (res in c(0.1,0.3,0.5,1)) {
  mcell.data <- FindClusters(mcell.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, 
          DimPlot(mcell.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), 
          DimPlot(mcell.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), 
          DimPlot(mcell.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))


sel.clust = "RNA_snn_res.0.3"
mcell.data <- SetIdent(mcell.data, value = sel.clust)
table(mcell.data@active.ident)
DimPlot(mcell.data, reduction = "umap", label = T)
DimPlot(mcell.data,group.by = 'Patient')
VlnPlot(mcell.data, 
        features = c("percent_ribo",'percent_mito','nFeature_RNA'), 
        ncol = 2, pt.size = 0)
#cluster 10: low nFeature, probably null cells
mcell.data <- subset(mcell.data, idents = c('10'), invert = T)

#mcell.data@assays$RNA = alldata@assays$RNA[,colnames(mcell.data)]
marker.genes = c('MARCO',"FABP4","C1QB","C1QA","MRC1","CD68",
                 'CCL5',
                 "LYZ",'VCAN','CD14',"FCN1",
                 'FCGR3B',
                 "SPP1","CTSB",'CTSL',"CTSZ",
                 "CD1C", "FCER1A", "CLEC10A", "IRF8",
                 'S100A14',
                 "MT1G",'MT1M','MT1X','MT1E',
                 "LAMP3"
                 )
DotPlot(mcell.data, features = marker.genes,assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis()

markers_genes <- FindAllMarkers(mcell.data, logfc.threshold = 0.25, 
                                test.use = "wilcox", 
                                min.pct = 0.2, 
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_mcellSubClustering.csv')
markers_genes = read.csv(file = './variables/marker_genes_mcellSubClustering.csv',row.names = 1)
markers_genes = split(markers_genes,markers_genes$cluster)

table(mcell.data@active.ident)

mcell.data<-RenameIdents(mcell.data, '0'='classical macrophage','1'='CCL5+IMC','2'='Monocytes',
                         '3'='Neutrophils','4'='SPP1+macrophage','5'='cDC','6'='notD',
                         '7'='MT+macrophage','8'='classical macrophage','9'='mDC')
mcell.data <- subset(mcell.data, idents = c('notD'), invert = T)
mcell.data$subType = Idents(mcell.data)
table(mcell.data$orig.ident,mcell.data$subType)

DotPlot(mcell.data, features = marker.genes, assay = "RNA",cluster.idents = T) + 
  coord_flip()+ 
  RotatedAxis()+
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')+
  theme(text = element_text(size = 8))
DoHeatmap(mcell.data, features = marker.genes, assay = "RNA",slot = 'counts',
          size = 2)+scale_fill_gradient2(low='blue',mid='white',high = 'red')


DimPlot(mcell.data,group.by = 'subType',reduction = 'umap',label = T,cols = pal_npg('nrc')(10))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/myeloid/Dimplot cellType labeled.pdf
plot_grid(ncol = 3, 
          DimPlot(mcell.data, group.by = 'subType',label = T) + NoAxes(), 
          DimPlot(mcell.data, group.by = "orig.ident") + NoAxes(), 
          DimPlot(mcell.data, group.by = "type") + NoAxes())
mcell.data$orig.ident = factor(mcell.data$orig.ident,
                               levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                          'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                          'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                          'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))
DimPlot(mcell.data,group.by = 'subType',split.by = 'orig.ident',ncol=4)
####

annotation = mcell.data@meta.data

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


ggplot(data = annotation,aes(x=From,fill=subType))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free')+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = pal_npg('nrc')(10))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/myeloid/cell type summary part4.pdf


saveRDS(mcell.data,file = './variables/mcell.data.Rds')
mcell.data = readRDS(file = './variables/mcell.data.Rds')
####
library(pheatmap)
mcell.data = readRDS(file = './variables/mcell.data.Rds')
com2orig = table(mcell.data$orig.ident,mcell.data$subType)
type.anno = data.frame(Patient= unlist(lapply(rownames(com2orig),function(a){
  return(strsplit(a,'_')[[1]][3])
})),
type = unlist(lapply(rownames(com2orig),function(a){
  return(substr(a,1,1))
})))
rownames(type.anno)=rownames(com2orig)
#pheatmap(cor(t(com2orig)),annotation_row = type.anno)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
f.data = t(f.data)

pheatmap(cor(f.data,method = 'spearman'),clustering_distance_rows = "correlation",annotation_row = type.anno)
pheatmap(cor(f.data,method = 'spearman'),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_row = type.anno)
pheatmap(scale(t(f.data)),clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_row = type.anno)
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

