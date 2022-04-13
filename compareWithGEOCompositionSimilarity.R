####Compare the sample-sample similarity in terms of cell compositions, 2021-08-19, Di Chen
setwd('H:/project/single cell/MPSC&ST/')
options(stringsAsFactors = F)
library(Seurat)
library(ggsci)
library(ggplot2)
library(cowplot)
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
com2orig = table(alldata$orig.ident,alldata$CellTypeManully)
rm(alldata)

ref.geo = read.delim(file='../documents/GEO/GSE131907_Lung_Cancer_cell_annotation.txt/GSE131907_Lung_Cancer_cell_annotation.txt')

ref.geo.sub.tummors = ref.geo[ref.geo$Sample_Origin %in% c('nLung','tLung'),]
ref.geo.sub.tummors$Cell_type[ref.geo.sub.tummors$Cell_type %in% c('T lymphocytes','NK cells')]='T&NK cell'

ref.geo.sub.tummors$Patient = substr(ref.geo.sub.tummors$Sample,7,8)

cellType.cols = pal_npg('nrc')(length(unique(ref.geo.sub.tummors$Cell_type)))
names(cellType.cols)=c('Epithelial cells','Endothelial cells','Fibroblasts',
                       'T&NK cell','B lymphocytes','Myeloid cells','MAST cells')
ref.geo.sub.tummors$Cell_type = factor(ref.geo.sub.tummors$Cell_type,levels = 
                                         c('Epithelial cells','Endothelial cells','Fibroblasts',
                                           'T&NK cell','B lymphocytes','Myeloid cells','MAST cells'))
ggplot(data = ref.geo.sub.tummors,aes(x=Sample,fill=Cell_type))+geom_bar(position = 'stack',width = 0.7)+
  coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free',space = 'free')+
  theme_cowplot()+
  scale_fill_manual(values = cellType.cols)+
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_blank())
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/compareWithGEO/T&NK/GEO cell type sum part4.pdf

com2orig.2 = table(as.character(ref.geo.sub.tummors$Sample),ref.geo.sub.tummors$Cell_type)

colnames(com2orig.2)=c('B cell',"Endothelial cell","Epithelial cell",
                       "Fibroblast cell","Mast cell", "Myeloid cell",
                       "T&NK cell")
both.cols = intersect(colnames(com2orig),colnames(com2orig.2))

com2orig.both = rbind(com2orig[,both.cols],com2orig.2[,both.cols])
library(corrplot)
corrplot(cor(t(com2orig.both),method = 'spearman'))
library(pheatmap)
type.anno = data.frame(From=c(rep('Our',19),rep('GEO',nrow(com2orig.2))),
                       Type=c(rep('Normal',8),rep('Tumor',11),rep('Normal',11),rep('Tumor',11)),
                       row.names = rownames(com2orig.both))
pheatmap(cor(t(com2orig.both)),annotation_row = type.anno)
 
f.data<-as.matrix(com2orig.both)
f.data = f.data/apply(f.data, 1, sum)
f.data = t(f.data)

pheatmap(cor(f.data,method = 'spearman'),clustering_distance_rows = "correlation",annotation_row = type.anno)
pheatmap(cor(t(com2orig/apply(com2orig, 1, sum)),method = 'spearman'),clustering_distance_rows = "correlation")
#pheatmap(t(f.data),clustering_distance_rows = "correlation",annotation_row = type.anno)


####PCA analysis####
#com2orig = table(alldata$orig.ident,alldata$RNA_snn_res.1)
com2orig = com2orig.both
com2orig = com2orig/apply(com2orig, 1, sum)
mat = t(com2orig)
mat = mat[apply(mat, 1, sum)>0,]
mat = scale(t(mat))
mat = mat[apply(mat, 1, sd)!=0,]
mat[is.na(mat)]
pheatmap(cor(t(mat),method = 'spearman'),clustering_distance_rows = "correlation",
         ,annotation_row = type.anno)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/compareWithGEO/T&NK/pheatmap pca utilized mat.pdf

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
pcr.cluster$Type = substr(pcr.cluster$name,1,1)
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=From,shape=Type),size = 3)+
  scale_color_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(4,'Set3'),pal_npg('nrc')(4)))+
  theme_cowplot()+geom_text_repel(aes(label = name),size=2,nudge_y = 0.1,nudge_x = 0.1)+
theme(axis.text = element_text(size=8),
      text = element_text(size=8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/compareWithGEO/T&NK/PCA reduction about com2orig for both GEO&Our data.pdf
####Read GEOdata martirx####
library(harmony)
geo.data = readRDS(file='GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds')
geo.data = CreateSeuratObject(counts = geo.data)
rownames(ref.geo)=ref.geo$Index
geo.data@meta.data$orig.ident = ref.geo[rownames(geo.data@meta.data),'Sample']
geo.data@meta.data = cbind(geo.data@meta.data,ref.geo[rownames(geo.data@meta.data),c(-1,-3)])

geo.data = geo.data[,geo.data$Sample_Origin %in% c('nLung','tLung')]
geo.data = FindVariableFeatures(geo.data, verbose = T)
geo.data = ScaleData(geo.data, verbose = T)
geo.data = RunPCA(geo.data, verbose = T, npcs = 50)
geo.data = RunUMAP(geo.data, dims = 1:50)
DimPlot(geo.data, reduction = "umap",group.by = 'orig.ident')
DimPlot(geo.data, reduction = "umap",group.by = 'Sample_Origin')

geo.data = RunHarmony(geo.data,group.by.vars = 'orig.ident')
geo.data = RunUMAP(geo.data,reduction = "harmony", dims = 1:30)
DimPlot(geo.data, reduction = "umap",group.by = 'orig.ident')
DimPlot(geo.data, reduction = "umap",group.by = 'Sample_Origin')
DimPlot(geo.data, reduction = "umap",group.by = 'Cell_type')
cxcl14.markers = c('CXCL14','CLDN2','CEACAM5','CEACAM6','MDK')
FeaturePlot(geo.data, reduction = "umap", 
            features = c("EPCAM","KRT19","KRT18","CDH1"), 
            order = T, slot = "counts", combine = T)

FeaturePlot(geo.data[,geo.data$Cell_type == 'Epithelial cells'], reduction = "umap", 
            features = cxcl14.markers, 
            order = T, slot = "counts", combine = T)


epi.geo.data = geo.data[,geo.data$Cell_type == 'Epithelial cells']
#epi.data = NormalizeData(epi.geo.data)
rm(epi.geo.data)
epi.data <- epi.data[!grepl("MALAT1", rownames(epi.data)), ]
# Filter Mitocondrial
epi.data <- epi.data[!grepl("^MT-", rownames(epi.data)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) epi.data
epi.data <- epi.data[ ! grepl('^RP[SL]', rownames(epi.data)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
epi.data <- epi.data[!grepl("^HB[^(P)]", rownames(epi.data)), ]

epi.data = FindVariableFeatures(epi.data)
epi.data = ScaleData(epi.data)
epi.data = RunPCA(epi.data, verbose = T, npcs = 50)
epi.data = RunHarmony(epi.data,c( "orig.ident" ))
epi.data <- RunUMAP(epi.data,reduction = "harmony", dims = 1:30)

DimPlot(epi.data, reduction = "umap",group.by = 'Cell_subtype')

FeaturePlot(epi.data, reduction = "umap", 
            features = cxcl14.markers, 
            order = T, slot = "counts", combine = T)

epi.data <- FindNeighbors(epi.data, reduction = "harmony",dims = 1:30, verbose = FALSE)
names(epi.data@graphs)

for (res in c(0.1,0.2,0.3,0.5)) {
  epi.data <- FindClusters(epi.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))
sel.clust = "RNA_snn_res.0.5"
epi.data <- SetIdent(epi.data, value = sel.clust)
table(epi.data@active.ident)
DimPlot(epi.data, reduction = "umap", label = T)
VlnPlot(epi.data, features = cxcl14.markers, 
        group.by = "RNA_snn_res.0.5", 
        ncol = 2, pt.size = 0)
table(epi.data@active.ident,epi.data$Cell_subtype)

markers_genes <- FindAllMarkers(epi.data, logfc.threshold = 0.25, 
                                test.use = "wilcox", 
                                min.pct = 0.1, 
                                only.pos = TRUE, 
                                assay = "RNA")
saveRDS(geo.data,file = './variables/GSE131907_Lung.Rds')

top25 <- markers_genes %>% group_by(cluster) %>% top_n(-15, p_val_adj)%>% top_n(15, avg_log2FC) 

marker.our = read.csv(file = './variables/marker_genes_epiSubClustering-v2.csv')
marker.cxcl14 = marker.our[marker.our$cluster == '1',]
top25.cxcl14 = marker.cxcl14 %>% top_n(-15, p_val_adj)%>% top_n(15, avg_log2FC) 
top25.cxcl14$cluster='CXCL14+'
top25 = rbind(top25,top25.cxcl14)
top25.genes.cxcl14 = as.character(top25.cxcl14$gene)

library(clusterProfiler)

top25.enrich = enricher(top25.genes.cxcl14,pvalueCutoff = 1.2,TERM2GENE = top25[,c('cluster','gene')],
                        minGSSize = 0,maxGSSize = 30)
top25.enrich.res = top25.enrich@result[-1,]

top25.enrich.res$Description = factor(top25.enrich.res$Description,levels = rev(as.character(top25.enrich.res$Description)))
ggplot(top25.enrich.res,mapping = aes(x=-log10(pvalue),y=Description))+
  geom_bar(stat = 'Identity',width=0.75,fill='grey')+
  geom_text(mapping = aes(x=1,y=ID,label=geneID),color='blue')+theme_cowplot()+
  geom_vline(xintercept = 2,lty=2,color='red')


top.mark = top25[top25$cluster == '2',][1:15,]
top.mark = top.mark[order(top.mark$avg_log2FC,decreasing = T),]
top.mark$gene = factor(top.mark$gene,levels = as.character(top.mark$gene[15:1]))
ggplot(top.mark,mapping = aes(x=avg_log2FC,y=gene))+geom_point(mapping = aes(size = pct.1,color=pct.1),fill='grey')+
  theme_cowplot()+scale_color_gradient2(low='white',high = 'red')

