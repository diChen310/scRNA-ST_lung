####Sub-clustering of the Epithelial cell, 2021-09-22,2021-09-28####
setwd('H:/Project/single cell/MPSC&ST/')
library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony)
library(ggsci)
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')

tcell.data = alldata[,alldata$CellTypeManully %in% c('T&NK cell')]

dim(tcell.data)#33261  45369
rm(alldata)
tcell.data = NormalizeData(tcell.data)
tcell.data = FindVariableFeatures(tcell.data, selection.method = 'mvp',verbose = T)
tcell.data = ScaleData(tcell.data, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                       verbose = T)

tcell.data = RunPCA(tcell.data, verbose = T, npcs = 50)
#tcell.data = RunUMAP(tcell.data, dims = 1:50, verbose = T)
tcell.data = RunHarmony(tcell.data,dims.use = 1:50,sigma = 0.30,
                        group.by.vars = c( "orig.ident" ),
                        max.iter.harmony = 30,
                        max.iter.cluster = 50
                        )

tcell.data <- RunUMAP(tcell.data,reduction = "harmony", dims = 1:30)
tcell.data <- FindNeighbors(tcell.data, reduction = "harmony",dims = 1:30, verbose = FALSE)
names(tcell.data@graphs)

for (res in c(0.1,0.3,0.5,1)) {
  tcell.data <- FindClusters(tcell.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, 
          DimPlot(tcell.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), 
          DimPlot(tcell.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), 
          DimPlot(tcell.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))


sel.clust = "RNA_snn_res.0.5"
tcell.data <- SetIdent(tcell.data, value = sel.clust)
table(tcell.data@active.ident)
DimPlot(tcell.data, reduction = "umap", label = T)
DimPlot(tcell.data,group.by = 'Patient')
VlnPlot(tcell.data, 
        features = c("percent_ribo",'percent_mito','nFeature_RNA'), 
        ncol = 2, pt.size = 0)


#tcell.data@assays$RNA = alldata@assays$RNA[,colnames(tcell.data)]
marker.genes = c("CD3D","CD3E","CD3G","XCL1","FCGR3A",'KLRD1',"KLRF1",
                 "CD8A","CD8B","CD4","IL7R",
                 "TCF7",'SELL','LEF1','CCR7',
                 'LAG3','TIGIT','PDCD1','HAVCR2','CTLA4',
                 'IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7',
                 'IL2RA','FOXP3','IKZF2','TGFB1','TGFB3','TGFBI','TGFBR1',
                 'MAF','CXCL13','CXCR5',
                 'IRF4','CREM','NR4A2','STAT4','IL12RB2',
                 'GATA3','STAT6','IL4','TRDC','TRGC2','TRGC1',
                 'CD27','ZNF683',
                 'CD45RA', 'CD45RO', 'CD69', 'CXCR3')
DotPlot(tcell.data, features = marker.genes,assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis()

markers_genes <- FindAllMarkers(tcell.data, logfc.threshold = 0.25, 
                                test.use = "wilcox", 
                                min.pct = 0.2, 
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_tcellSubClustering.V2.csv')
markers_genes = read.csv(file = './variables/marker_genes_tcellSubClustering.V2.csv',row.names = 1)
markers_genes = split(markers_genes,markers_genes$cluster)

markers <- read.csv("documents/T cell markers.csv")
celltype_list <- lapply(unique(markers$CellType), function(x) {
  x <- markers$Gene[markers$CellType == x]
})
names(celltype_list) <- unique(markers$CellType)

library(fgsea)

# run fgsea for each of the clusters in the list
res <- lapply(markers_genes, function(x) {
  if(nrow(x)>1){
    gene_rank <- setNames(x$avg_logFC, x$gene)
    fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank,nperm = 500)
  }else{
    fgseaRes = data.frame()
  }
  return(fgseaRes)
})

table(tcell.data@active.ident)

tcell.data<-RenameIdents(tcell.data, '0'='CD4+Trm','1'='CD4+Tn','2'='Treg',
                       '3'='NK1','4'='CD8+Tem','5'='CD8+CTL','6'='CD8+Trm',
                       '7'='CD8+Teff','8'='NK2','9'='Tfh','10'='NK3',
                       '11'='Th1-like','12'='CD8+CD161+','13'='NK2','14'='CD4+Tn','15'='CD4+Trm')
#tcell.data = RenameIdents(tcell.data,"KLRG1+CD8+Teff"='CD8+Teff')
#tcell.data <- subset(tcell.data, idents = c('Myeloid'), invert = T)

tcell.data$subType = Idents(tcell.data)
table(tcell.data$orig.ident,tcell.data$subType)
 marker.genes = c("XCL1","XCL2","FCGR3A",'KLRD1',"KLRF1",'KLRC1','KLRB1',
                   "CD8A","CD8B","CD4","IL7R", "CX3CR1", "FGFBP2",'PTGDS','SPON2',
                    "TCF7",'SELL','LEF1','CCR7',
                    'TIGIT','PDCD1','HAVCR2','CTLA4',
                    'GZMA','GNLY','PRF1','GZMB','GZMK','NKG7',
                     'IL2RA','FOXP3','NR4A1',"NR4A2","CXCL13",
                  "CD200", "ICA1",  "TOX2", "BTLA",'KLRG1')

DotPlot(tcell.data, features = marker.genes, 
        assay = "RNA",cluster.idents = T,
        dot.scale = 4) + 
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')+
  coord_flip()+ 
  RotatedAxis()+
  theme(text = element_text(size = 8))

DoHeatmap(tcell.data, features = marker.genes, assay = "RNA",slot = 'counts',
          size = 2)+scale_fill_gradient2(low='blue',mid='white',high = 'red')


DimPlot(tcell.data,group.by = 'subType',reduction = 'umap',label = T,cols = c(pal_npg('nrc')(10),"#E6AB02" ,"#A6761D", "#666666"))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/tcell/V2/Dimplot cell subtype.pdf
plot_grid(ncol = 3, 
          DimPlot(tcell.data, group.by = 'subType',label = T) + NoAxes(), 
          DimPlot(tcell.data, group.by = "orig.ident") + NoAxes(), 
          DimPlot(tcell.data, group.by = "type") + NoAxes())
DimPlot(tcell.data,group.by = 'subType',reduction = 'umap',label = T)
tcell.data$orig.ident = factor(tcell.data$orig.ident,
                               levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                          'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                          'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                          'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))
DimPlot(tcell.data,group.by = 'subType',split.by = 'orig.ident',ncol=4)
epi.data.tumor = tcell.data[,tcell.data$type=='TUMOR']
p1 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P1'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p2 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P2'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p3 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P3'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p4 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P4'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p1/p2/p3/p4
####

annotation = tcell.data@meta.data

annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = c('TI','TM','TI1','TI2','TM1','TM2','TS1','TS2','NI','NM','NS'))
subType.counts = table(annotation$subType)
annotation$subType = factor(annotation$subType,levels =names(subType.counts[order(subType.counts,decreasing = T)]) )
library(ggplot2)
library(RColorBrewer)
ggplot(data = annotation,aes(x=Patient,fill=type))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+
  scale_fill_manual(values = brewer.pal(3,'Paired')[-1])+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))
ggplot(data = annotation,aes(x=Patient,fill=From))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))
ggplot(data = annotation,aes(x=Patient))+geom_bar(fill='blue',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(subType))+theme_bw()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))


ggplot(data = annotation,aes(x=From,fill=subType))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free')+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = c(pal_npg('nrc')(10),"#E6AB02" ,"#A6761D", "#666666"))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/tcell/V2/cell type summary part4.pdf

ggplot(data = annotation,aes(x=From,fill=From))+geom_bar(width = 0.7)+coord_flip()+
  facet_grid( Patient~subType,scales = 'free',space = 'free')+theme_cowplot()+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,angle = 90),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  scale_y_log10()


saveRDS(tcell.data,file = './variables/tcell.data.Rds')

tcell.data = readRDS(file = './variables/tcell.data.Rds')
####
library(pheatmap)
tcell.data = readRDS(file = './variables/tcell.data.Rds')
com2orig = table(tcell.data$orig.ident,tcell.data$subType)
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
pcr.cluster$Source = ifelse(grepl('LUNG',rownames(pcr.cluster)),'GEO','Our')
library(ggrepel)
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=From,shape=Source),size = 4)+
  scale_color_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(8,'Set3')))+
  theme_bw()+geom_text_repel(aes(label = name),size=4,nudge_y = 0.1,nudge_x = 0.1)
####top genes check
int.genes = c("ALOX15B", "CD79B" ,  "ITSN2"  , "LCN2"   , "PRSS12" , "TNFRSF18" ,"WIF1"    )
VlnPlot(tcell.data,group.by = 'subType',
        cols = c(pal_npg('nrc')(10),"#E6AB02" ,"#A6761D", "#666666"),
        features = 'TNFRSF18',sort = 'counts',
        ncol = 1,
        pt.size = 0)
VlnPlot(tcell.data,group.by = 'orig.ident',
        features = c('TNFRSF18'),sort = 'counts',
        ncol = 1,
        pt.size = 0)
VlnPlot(alldata,group.by = 'CellTypeManully',
        cols = c(pal_npg('nrc')(10),"#E6AB02" ,"#A6761D", "#666666"),
        features = 'TNFRSF18',sort = 'counts',
        ncol = 1,
        pt.size = 0)
