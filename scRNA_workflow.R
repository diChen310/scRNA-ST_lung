####Workflow for the four samples 2021-07-29 to 2021-08#######
#BiocManager::install("SydneyBioX/scdney")
library(Seurat)
#library(DoubletFinder)
library(Matrix)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)
#library(scdney)
setwd('H:/project/single cell/MPSC&ST/')

##################################################Data input##############################################
lung.ti.p1 = Read10X(data.dir = '../ourData/scRNA/P1/P1_T_R_I/filtered_feature_bc_matrix')
lung.tm.p1 = Read10X(data.dir = '../ourData/scRNA/P1/P1_T_R_M/filtered_feature_bc_matrix')
lung.ni.p1 = Read10X(data.dir = '../ourData/scRNA/P1/P1_N_R_I/filtered_feature_bc_matrix')
lung.nm.p1 = Read10X(data.dir = '../ourData/scRNA/P1/P1_N_R_M/filtered_feature_bc_matrix')
lung.ti.p2 = Read10X(data.dir = '../ourData/scRNA/P2/P2_T_R_I/filtered_feature_bc_matrix')
lung.tm.p2 = Read10X(data.dir = '../ourData/scRNA/P2/P2_T_R_M/filtered_feature_bc_matrix')
lung.ni.p2 = Read10X(data.dir = '../ourData/scRNA/P2/P2_N_R_I/filtered_feature_bc_matrix')
lung.nm.p2 = Read10X(data.dir = '../ourData/scRNA/P2/P2_N_R_M/filtered_feature_bc_matrix')
lung.ti1.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_T1_R_I/filtered_feature_bc_matrix')
lung.tm1.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_T1_R_M/filtered_feature_bc_matrix')
lung.ni.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_N_R_I/filtered_feature_bc_matrix')
lung.nm.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_N_R_M/filtered_feature_bc_matrix')
lung.ti2.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_T2_R_I/filtered_feature_bc_matrix')
lung.tm2.P3 = Read10X(data.dir = '../ourData/scRNA/P8/P8_T2_R_M/filtered_feature_bc_matrix')
#P4-T-L-I??P4-1T????P4-N-L-I??P4-1N????P4-T1-L-S??P4-2T1????P4-T2-L-S??P4-2T2????P4-N-L-S??P4-2N??
lung.ti.p4 = read.delim(file = 'H:/project/single cell/ourData/scRNA/V2/P4-1T_matrix.tsv')
lung.ni.p4 = read.delim(file = 'H:/project/single cell/ourData/scRNA/V2/P4-1N_matrix.tsv')
lung.ts1.p4 = read.delim(file = 'H:/project/single cell/ourData/scRNA/V2/P4-2T1_matrix.tsv')
lung.ts2.p4 = read.delim(file = 'H:/project/single cell/ourData/scRNA/V2/P4-2T2_matrix.tsv')
lung.ns.p4 = read.delim(file = 'H:/project/single cell/ourData/scRNA/V2/P4-2N_matrix.tsv')

sdata.ti.p1 <- CreateSeuratObject(lung.ti.p1,min.cells=3, project = "TI_R_P1")
sdata.tm.p1 <- CreateSeuratObject(lung.tm.p1,min.cells=3, project = "TM_R_P1")
sdata.ni.p1 <- CreateSeuratObject(lung.ni.p1,min.cells=3, project = "NI_R_P1")
sdata.nm.p1 <- CreateSeuratObject(lung.nm.p1,min.cells=3, project = "NM_R_P1")
sdata.ti.p2 <- CreateSeuratObject(lung.ti.p2,min.cells=3, project = "TI_R_P2")
sdata.tm.p2 <- CreateSeuratObject(lung.tm.p2,min.cells=3, project = "TM_R_P2")
sdata.ni.p2 <- CreateSeuratObject(lung.ni.p2,min.cells=3, project = "NI_R_P2")
sdata.nm.p2 <- CreateSeuratObject(lung.nm.p2,min.cells=3, project = "NM_R_P2")
sdata.ti1.P3 <- CreateSeuratObject(lung.ti1.P3,min.cells=3, project = "TI1_R_P3")
sdata.tm1.P3 <- CreateSeuratObject(lung.tm1.P3,min.cells=3, project = "TM1_R_P3")
sdata.ni.P3 <- CreateSeuratObject(lung.ni.P3,min.cells=3, project = "NI_R_P3")
sdata.nm.P3 <- CreateSeuratObject(lung.nm.P3,min.cells=3, project = "NM_R_P3")
sdata.ti2.P3 <- CreateSeuratObject(lung.ti2.P3,min.cells=3, project = "TI2_R_P3")
sdata.tm2.P3 <- CreateSeuratObject(lung.tm2.P3,min.cells=3, project = "TM2_R_P3")
sdata.ti.p4 <- CreateSeuratObject(counts = lung.ti.p4,min.cells=3, project = "TI_L_P4")
sdata.ni.p4 <- CreateSeuratObject(counts = lung.ni.p4,min.cells=3, project = "NI_L_P4")
sdata.ts1.p4<- CreateSeuratObject(counts = lung.ts1.p4,min.cells=3, project = "TS1_L_P4")
sdata.ts2.p4  <- CreateSeuratObject(counts = lung.ts2.p4,min.cells=3, project = "TS2_L_P4")
sdata.ns.p4 <- CreateSeuratObject(counts = lung.ns.p4,min.cells=3, project = "NS_L_P4")

sdata.ti.p1$type = "TUMOR"
sdata.tm.p1$type = "TUMOR"
sdata.ni.p1$type = "NORMAL"
sdata.nm.p1$type = "NORMAL"
sdata.ti.p2$type = "TUMOR"
sdata.tm.p2$type = "TUMOR"
sdata.ni.p2$type = "NORMAL"
sdata.nm.p2$type = "NORMAL"
sdata.ti1.P3$type = "TUMOR"
sdata.tm1.P3$type = "TUMOR"
sdata.ni.P3$type = "NORMAL"
sdata.nm.P3$type = "NORMAL"
sdata.ti2.P3$type = "TUMOR"
sdata.tm2.P3$type = "TUMOR"
sdata.ti.p4$type = "TUMOR"
sdata.ts1.p4$type = "TUMOR"
sdata.ts2.p4$type = "TUMOR"
sdata.ni.p4$type = "NORMAL"
sdata.ns.p4$type = "NORMAL"


alldata <- merge(sdata.ti.p1, c(sdata.tm.p1, sdata.ni.p1, sdata.nm.p1,
                                sdata.ti.p2, sdata.tm.p2, sdata.ni.p2, sdata.nm.p2,
                                sdata.ti1.P3,sdata.tm1.P3,sdata.ni.P3,sdata.nm.P3,sdata.ti2.P3,sdata.tm2.P3,
                                sdata.ti.p4,sdata.ts1.p4,sdata.ts2.p4,sdata.ni.p4,sdata.ns.p4), 
                 add.cell.ids = c("TI-P1", "TM-P1", "NI-P1", "NM-P1","TI-P2", "TM-P2", "NI-P2", "NM-P2",
                                  "TI1-P3",'TM1-P3','NI-P3','NM-P3','TI2-P3','TM2-P3',
                                  "TI-P4",'TS1-P4','TS2-P4','NI-P4','NS-P4'))
rm(lung.ti.p1,lung.tm.p1,lung.ni.p1,lung.nm.p1,lung.ti.p2,lung.tm.p2,lung.ni.p2,lung.nm.p2,
   lung.ni.P3,lung.ni.p4,lung.nm.P3,lung.ns.p4,lung.ti.p4,lung.ti1.P3,lung.ti2.P3,lung.tm1.P3,lung.tm2.P3,
   lung.ts1.p4,lung.ts2.p4,
   sdata.ti.p1,sdata.tm.p1, sdata.ni.p1, sdata.nm.p1,sdata.ti.p2,sdata.tm.p2, sdata.ni.p2, sdata.nm.p2,
   sdata.ti1.P3,sdata.tm1.P3,sdata.ni.P3,sdata.nm.P3,sdata.ti2.P3,sdata.tm2.P3,
   sdata.ti.p4,sdata.ts1.p4,sdata.ts2.p4,sdata.ni.p4,sdata.ns.p4)
gc()

as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
head(alldata@meta.data, 10)
table(alldata$orig.ident)
alldata$orig.ident[alldata$orig.ident == 'P4.1N']='NI_L_P4'
alldata$orig.ident[alldata$orig.ident == 'P4.2N']='NS_L_P4'
alldata$orig.ident[alldata$orig.ident == 'P4.1T']='TI_L_P4'
alldata$orig.ident[alldata$orig.ident == 'P4.2T1']='TS1_L_P4'
alldata$orig.ident[alldata$orig.ident == 'P4.2T2']='TS2_L_P4'
# 
# NI_L_P4  NI_R_P1  NI_R_P2  NI_R_P3  NM_R_P1  NM_R_P2  NM_R_P3  NS_L_P4  TI_L_P4  TI_R_P1  TI_R_P2 TI1_R_P3 TI2_R_P3 
# 5794     8087     9577     9891     9143     9500    12943     6011     4963     8580    10554     9474     6590 
# TM_R_P1  TM_R_P2 TM1_R_P3 TM2_R_P3 TS1_L_P4 TS2_L_P4 
# 9300     9927     9186    11373     3820     5970 
saveRDS(alldata,file='./variables/alldata_merge4.Rds')
#############################################QC: Filter data based on the distributions###########################

#calculate the proportion gene expression that comes from mito proteins.
# Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
#calculate the proportion gene expression that comes from ribosomal proteins.
#\ Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")

#alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
## count log10GenesPerUMI and add it to seurat object metadata
alldata$log10GenesPerUMI <- log10(alldata$nFeature_RNA)/log10(alldata$nCount_RNA)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
# VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 3) + 
#   NoLegend()
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.0, ncol = 3) + 
  NoLegend()
#FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
alldata@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.75)

selected_c <- WhichCells(alldata, expression = nFeature_RNA > 250)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]

data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
# par(mar = c(4, 8, 2, 1))
# C <- data.filt@assays$RNA@counts
# C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
# most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
# boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
#         col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

selected_mito <- WhichCells(data.filt, expression = percent_mito < 25)
#selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 2)

data.filt <- subset(data.filt, cells = selected_mito)
#data.filt <- subset(data.filt, cells = selected_ribo)

dim(data.filt)
data.filt <- subset( data.filt, log10GenesPerUMI > 0.75 )	
table(data.filt$orig.ident)


feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()

# Filter MALAT1
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]

# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
dim(data.filt)
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
data.filt <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]

dim(data.filt)#27350 121354
table(data.filt$orig.ident)

# NI_P1  NI_P2  NI_P3  NM_P1  NM_P2  NM_P3  TI_P1  TI_P2 TI1_P3 TI2_P3  TM_P1 
# 7038   9151   9002   8094   9009  12228   8225   9858   8363   5782   8917 
# TM_P2 TM1_P3 TM2_P3 
# 9465   7873  10061




gc()

table(data.filt$orig.ident)
saveRDS(data.filt,file='./variables/data.filt_P1285.Rds')
#data.filt = readRDS(file='./variables/data.filt_P1285.Rds')
#rm(data.filt)


####################integrate by harmony, 2021-08-06###################
library(harmony)
data.filt = readRDS(file='./variables/data.filt_P1285.Rds')
data.filt = NormalizeData(data.filt)
data.filt = FindVariableFeatures(data.filt, verbose = T)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                      verbose = T)

data.filt = RunPCA(data.filt, verbose = T, npcs = 50)
ElbowPlot(data.filt,ndims = 50)
DimHeatmap(data.filt, dims = 31:35, cells = 500, balanced = TRUE)
DimHeatmap(data.filt, dims = 41:45, cells = 500, balanced = TRUE)

data.filt = RunUMAP(data.filt, dims = 1:50, verbose = T)
data.filt = RunTSNE(data.filt, dims = 1:50, verbose = T)
data.filt$platform = ifelse(grepl('P4',data.filt$orig.ident),'SG','10X')
data.filt$Patient = unlist(lapply(data.filt$orig.ident,function(a){
  strsplit(a,'_')[[1]][3]
}))
DimPlot(data.filt, reduction = "umap",group.by = 'platform')
DimPlot(data.filt, reduction = "umap",group.by = 'orig.ident')
DimPlot(data.filt, reduction = "umap",group.by = 'Patient')

p1.data = data.filt[,data.filt$Patient == 'P4']
DimPlot(p1.data, reduction = "umap",group.by = 'orig.ident')
rm(p1.data)
data.filt = RunHarmony(data.filt,group.by.vars = 'orig.ident')
data.filt = RunUMAP(data.filt,reduction = "harmony", dims = 1:30)
DimPlot(data.filt, reduction = "umap",group.by = 'platform')
DimPlot(data.filt, reduction = "umap",group.by = 'orig.ident')
DimPlot(data.filt, reduction = "umap",group.by = 'Patient')
saveRDS(data.filt,file='./variables/data.filt_harmony.Rds')

data.filt = FindNeighbors(data.filt , reduction = "harmony", dims = 1:30) 
names(data.filt@graphs)
for (res in c( 0.5, 1, 1.5)) {
  data.filt <- FindClusters(data.filt, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
DimPlot(data.filt, label = T,group.by = 'RNA_snn_res.1.5') + NoAxes()

sel.clust = "RNA_snn_res.0.5"

data.filt <- SetIdent(data.filt, value = sel.clust)
table(data.filt@active.ident)
# plot this clustering
DimPlot(data.filt, label = T) + NoAxes()
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(data.filt, features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()

marker.genes = c("EPCAM","KRT19","KRT18","CDH1","DCN","THY1","COL1A1","COL1A2",
                 "PECAM1","CLDN5","FLT1","RAMP2","PTPRC","CD3D","CD3E","CD3G",
                 "NKG7","GNLY","NCAM1","KLRD1","CD79A","IGHM","IGHG3","CD19",'IGHA2',
                 "LYZ","MARCO","CD68","FCGR3A","S100A8","S100A9",
                 "KIT","MS4A2","GATA2")

DotPlot(data.filt, features = marker.genes, group.by = 'RNA_snn_res.0.5',assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() 
saveRDS(data.filt,file='./variables/data.filt_harmony.Rds')

markers_genes <- FindAllMarkers(data.filt,  logfc.threshold = 0.5, 
                                test.use = "wilcox", 
                                min.pct = 0.2,
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_harmony.csv')
markers_genes=read.csv(file = './variables/marker_genes_harmony.csv',row.names = 1)

top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)%>% top_n(25, avg_logFC) 
par(mfrow=c(5,6), mar = c(4, 6, 3, 1))
top25$cluster = factor(top25$cluster,levels = c(0:26))
for (i in 0:26) {
  barplot(sort(setNames(top25$avg_logFC, top25$gene)[top25$cluster == i], F), horiz = T, 
          las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}
DGE_list <- split(markers_genes, markers_genes$cluster)

#Remove null or MT cells
data.filt <- subset(data.filt, idents = c(6,13,16,22), invert = T)
#data.filt$cluster = data.filt$RNA_snn_res.0.5
#data.filt <- SetIdent(data.filt, value = sel.clust)
table(data.filt@active.ident)
dim(data.filt)#33261 142780
###
####Cancer stem cell markers
FeaturePlot(data.filt, reduction = "umap", 
            features = c("ALDH1A1","PROM1","CD44","ALCAM",'NOTCH1','EPCAM'), 
            order = T, slot = "counts", combine = T)

data.filt<-RenameIdents(data.filt, '0'='Epithelial cell','1'='T cell','2'='NK cell','3'='Myeloid cell','4'='Myeloid cell',
                      '5'='T cell','7'='Myeloid cell','8'='Epithelial cell','9'='Mast cell',
                      '10' = 'B cell','11'='Myeloid cell','12'='Myeloid cell','14'='Endothelial cell',
                      '15'='Myeloid cell','17'='Fibroblast cell','18'='Epithelial cell','19'='Epithelial cell',
                      '20'='Epithelial cell','21'='Epithelial cell','23'='B cell','24'='Endothelial cell',
                      '25'='Myeloid cell','26'='Epithelial cell')
data.filt = RenameIdents(data.filt,'T cell'='T&NK cell','NK cell'='T&NK cell')

DimPlot(object = data.filt, label = TRUE) #
data.filt$CellTypeManully = Idents(data.filt)
saveRDS(data.filt,file = './variables/data.filt_harmony_annotated-v2.Rds')
data.filt = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
plot.data = data.frame(Orig.ident = names(table(data.filt$orig.ident)),Number = as.integer(table(data.filt$orig.ident)))
plot.data$Patient =unlist(lapply(as.character(plot.data$Orig.ident),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
ggplot(plot.data,aes(x=Orig.ident,y=Number,fill=Patient))+geom_bar(stat = 'identity',width = 0.65)+
  facet_wrap(vars(Patient),scales = 'free_x')+theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_npg()
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/cell number for orig.idents barplot.pdf
cellType.cols = pal_npg('nrc')(length(unique(data.filt$CellTypeManully)))
names(cellType.cols)=c('Epithelial cell','Endothelial cell','Fibroblast cell',
                       'T&NK cell','B cell','Myeloid cell','Mast cell')
DimPlot(object = data.filt, label = TRUE,group.by = 'CellTypeManully',cols = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Dimplot Celltype annotated-T&NK.pdf

####Re-calculate the cell type markers####

Idents(data.filt) = data.filt$CellTypeManully
markers_genes <- FindAllMarkers(data.filt, logfc.threshold = 0.35, 
                                test.use = "wilcox", 
                                min.pct = 0.1,
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file='Data.filt_harmony_cellTypeMarkers.csv')
markers_genes = read.csv(file='Data.filt_harmony_cellTypeMarkers.csv',row.names = 1)
top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, p_val_adj)%>% top_n(5, avg_logFC) 
DotPlot(data.filt, features = as.character(unique(top5$gene)), group.by = 'CellTypeManully',assay = "RNA",cluster.idents = T) + 
  coord_flip()+ 
  RotatedAxis()+ 
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')+
  scale_size(range = c(1,4))+
  theme(text = element_text(size = 8)) 
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Cell type markers dotplot.pdf
#DoHeatmap(data.filt, features = as.character(unique(top5$gene)), group.by = 'CellTypeManully',assay = "RNA")
#################################Statistics basically###########################################

data.filt$CellTypeManully = unlist(lapply(as.character(data.filt$CellTypeManully),function(a){
  return(strsplit(a,' ')[[1]][1])
}))
annotation = data.filt@meta.data
annotation$Patient = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = rev(c('TS1','TS2','TM','TI','TM1','TM2','TI1','TI2','NS','NM','NI')))
annotation$CellTypeManully = factor(annotation$CellTypeManully,levels = c('T&NK','Epithelial','Myeloid','B',
                                                                          'Mast','Endothelial','Fibroblast'))
names(cellType.cols)=c('Epithelial','Endothelial','Fibroblast',
                       'T&NK','B','Myeloid','Mast')

library(ggplot2)
library(RColorBrewer)
ggplot(data = annotation,aes(x=Patient,fill=type))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = brewer.pal(3,'Paired')[-1])
ggplot(data = annotation,aes(x=Patient,fill=From))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))+theme( strip.background = element_blank())
ggplot(data = annotation,aes(x=Patient))+geom_bar(fill='blue',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+theme( strip.background = element_blank())

ggplot(data = annotation,aes(x=From,fill=CellTypeManully))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free')+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(3,'Set3')))
ggplot(data = annotation,aes(x=From,fill=CellTypeManully))+geom_bar(position = 'stack',width = 0.7)+
  coord_flip()+
  facet_grid( rows= vars(Patient),scales = 'free',space = 'free')+
  theme_cowplot()+
  scale_fill_manual(values = cellType.cols)+
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_blank())
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/cell type sum part4.pdf

#########Correlations between samples################
library(corrplot)
data.filt.i = data.filt[,data.filt$Patient == 'P1']
com2orig = table(data.filt.i$orig.ident,data.filt.i$CellTypeManully)

corrplot(cor(t(com2orig),method='spearman'), method="pie")


data.ss = split(colnames(data.filt.i),data.filt.i$orig.ident)

ss.mexp = matrix(,nrow = length(data.ss),ncol = nrow(data.filt.i@assays$RNA@scale.data),
                     dimnames = list(names(data.ss),rownames(data.filt.i@assays$RNA@scale.data)))

for(i in 1:nrow(ss.mexp)){
  
  ss.mexp[i,] = apply(data.filt.i@assays$RNA@scale.data[,data.ss[[i]]],1,mean)
  
  
}
library(pheatmap)

library(corrplot)
corrplot(cor(t(ss.mexp[,intersect(marker.genes,colnames(ss.mexp))]),method='spearman'), method="pie")
#pheatmap(cor(t(ss.mexp[,intersect(marker.genes,colnames(ss.mexp))]),method='spearman'))

##########################Type0: diff markers for samples, do not considering cell types######################################
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')

patients = as.character(unique(alldata$Patient))

patients.markers = list()
for(i in 1:length(patients)){
  patient = patients[i]
  print(patient)
  sub.data = alldata[,alldata$Patient == patient]
  Idents(sub.data)=sub.data$orig.ident
  orig.idents = as.character(unique(sub.data$orig.ident))
  T.idents = orig.idents[substr(orig.idents,1,1)=='T']
  N.idents = orig.idents[substr(orig.idents,1,1)=='N']
  TOther.idents = T.idents[!grepl('TI',T.idents)]
  NOther.idents = N.idents[!grepl('NI',N.idents)]
  
    res.markers.T= FindMarkers(sub.data,
                               ident.1=orig.idents[grepl('TI',orig.idents)],
                               ident.2=TOther.idents,
                               assay = 'RNA',
                               logfc.threshold=-Inf,
                               min.pct = -Inf)
    res.markers.T$Gene = rownames(res.markers.T)
    
    res.markers.N= FindMarkers(sub.data,
                               ident.1=orig.idents[grepl('NI',orig.idents)],
                               ident.2=NOther.idents,
                               assay = 'RNA',
                               logfc.threshold=-Inf,
                               min.pct = -Inf)
    res.markers.N$Gene = rownames(res.markers.N)
    
    res.markers.T$Compare=paste(c('TI',unique(substr(TOther.idents,1,2))),collapse = 'VS')
    res.markers.N$Compare=paste(c('NI',unique(substr(NOther.idents,1,2))),collapse = 'VS')
    res.markers = rbind(res.markers.T,res.markers.N)
    res.markers$patient = patient
    
  
  patients.markers[[i]]=res.markers
}
saveRDS(patients.markers,file='variables/patients.diff.markers_Type0.Rds')

# top5 = res.markers %>% group_by(Compare) %>% top_n(-5, p_val_adj)  %>% top_n(5, avg_logFC)
# DotPlot(sub.sub.data, features = as.character(top5$Gene),assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() 
patients.markers = readRDS(file='variables/patients.diff.markers_Type0.Rds')
res.markers = Reduce(rbind,patients.markers)
res.markers$tag = paste(res.markers$Gene,res.markers$patient)
res.markers.T = res.markers[grepl('T',res.markers$Compare) ,]
res.markers.N = res.markers[grepl('N',res.markers$Compare) ,]
colnames(res.markers.N)= paste0(colnames(res.markers.N),2)
length(unique(res.markers.T$tag))
length(unique(res.markers.N$tag2))
res.markers.2 = merge(res.markers.T,res.markers.N,by.x='tag',by.y = 'tag2')
nrow(filter(res.markers.2,Gene != Gene2))
res.markers.2.sig = filter(res.markers.2,p_val_adj<0.01|p_val_adj2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & abs(res.markers.2.sig$avg_logFC)>0.25,'SigType']='TOnly'
res.markers.2.sig[res.markers.2.sig$p_val_adj2<0.01 & abs(res.markers.2.sig$avg_logFC2)>0.25,'SigType']='NOnly'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC>0,'SigType']='TNSame'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC<0,'SigType']='TNDiff'
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,Compare) %>% top_n(5, abs(avg_logFC2)+abs(avg_logFC)) 
###Scatter plot to compare Tdiff and NDiff
ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,5,4)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)
###Barplot summary
ggplot(res.markers.2.sig,aes(x=SigType,fill=SigType))+geom_bar(width = 0.6)+
  facet_wrap(~patient,nrow =1)+
  theme_cowplot()+scale_fill_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+scale_y_log10()

res.sig.same = res.markers.2.sig[res.markers.2.sig$SigType == 'TNSame',]
genes.sig.same = as.character(unique(res.sig.same$Gene))
genes.sig.NP = unlist(lapply(genes.sig.same, function(a){
  max(c(length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC>0,'patient'])),
        length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC<0,'patient']))))
  
}))

genes.sig.avgS = unlist(lapply(genes.sig.same, function(a){
  mean(res.sig.same[res.sig.same$Gene==a,'avg_logFC']+res.sig.same[res.sig.same$Gene==a,'avg_logFC2'])
}))
gene.n.sum = data.frame(Gene = genes.sig.same,
                        NP = genes.sig.NP,
                        avgS = genes.sig.avgS)
top.genes = as.character(gene.n.sum[gene.n.sum$NP>=4 ,'Gene']) #None meet the requirement
#res.sig.same.top = res.sig.same[res.sig.same$Gene %in% top.genes,]


res.sig.diff = res.markers.2.sig[res.markers.2.sig$SigType %in% c('TOnly','TNDiff'),]
genes.sig.diff = as.character(unique(res.sig.diff$Gene))
genes.sig.NP = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a,'patient']))
}))

genes.sig.avgS = unlist(lapply(genes.sig.diff, function(a){
  mean(abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC'])+abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC2']))
}))
gene.n.sum.diff = data.frame(Gene = genes.sig.diff,
                             NP = genes.sig.NP,
                             NC = genes.sig.NC,
                             avgS = genes.sig.avgS)
top.genes.diff = as.character(gene.n.sum.diff[gene.n.sum.diff$NP==4,'Gene'])
res.sig.diff.top = res.sig.diff[res.sig.diff$Gene %in% top.genes.diff,]

###
ggplot(res.sig.diff.top,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(res.sig.diff.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

####Type0:Pathway analysis for Type 0 diff Genes, 2021-09-13####
load(file = 'E:/project/metaboliteProteinInteraction/variables/KEGG_path2Gene.RData')
path2cate <- read.csv(file='E:/data/KEGG/KEGG hsa pathway list.csv',row.names = 1)
rownames(path2cate) = path2cate$Name

patients.marker.paths = list()

for(i in 1:length(patients.markers)){
  res.i = patients.markers[[i]]
  patient = unique(res.i$patient)
  cellTypes = as.character(unique(res.i$cellType))
  res.i.TI = res.i[grepl('T',res.i$Compare),]
  res.i.TO = res.i[grepl('N',res.i$Compare),]
 
  print(patient)
  
    
    FC = res.i.TI$avg_logFC
    names(FC)= as.character(res.i.TI$Gene)
    FC = FC[order(FC,decreasing = T)]
   
    path.TI = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
      
    
    
    FC = res.i.TO$avg_logFC
    names(FC)= as.character(res.i.TO$Gene)
    FC = FC[order(FC,decreasing = T)]
    
    path.TO = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
    
    
    if(nrow(path.TI)>0){path.TI$Compare=unique(res.i.TI$Compare)}
    if(nrow(path.TO)>0){path.TO$Compare=unique(res.i.TO$Compare)}
    path.2 = rbind(path.TI,path.TO)
    path.2$patient = patient
    
  
  patients.marker.paths[[i]]=path.2
}
saveRDS(patients.marker.paths,file='variables/patients.diff.marker.pathways_Type0.Rds')

res.markers = Reduce(rbind,patients.marker.paths)
res.markers$tag = paste(res.markers$ID,res.markers$cellType,res.markers$patient)
res.markers.TI = res.markers[grepl('T',res.markers$Compare),]
res.markers.TOthers = res.markers[grepl('N',res.markers$Compare),]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TOthers,res.markers.TI,by.x='tag2',by.y = 'tag')
nrow(filter(res.markers.2,ID != ID))
res.markers.2.sig = filter(res.markers.2,pvalue<0.01|pvalue2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & res.markers.2.sig$pvalue2>=0.01,'SigType']='TOnly'
res.markers.2.sig[res.markers.2.sig$pvalue2<0.01 & res.markers.2.sig$pvalue>=0.01,'SigType']='NOnly'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2>0,'SigType']='TNSame'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2<0,'SigType']='TNDiff'

disease.pathways = path2cate[path2cate$cate=='6. Human Diseases','Name']
res.markers.2.sig = res.markers.2.sig[res.markers.2.sig$ID %in% as.character(disease.pathways)==F,]
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,SigType) %>% top_n(5, abs(NES)+abs(NES2)) 
res.markers.2.sig$cate = path2cate[res.markers.2.sig$ID,'cate']
###Summary bar plot

ggplot(res.markers.2.sig,aes(x=SigType,fill=cate))+geom_bar(width = 0.6)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

####Scattar plot

ggplot(res.markers.2.sig,aes(x=NES,y=NES2,color=SigType))+geom_point(size=3)+
  facet_wrap(~patient,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=NES,y=NES2,color=SigType,label=ID),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)


ggplot(res.markers.2.sig)+geom_point(aes(x=Compare,y=ID,color=NES,size=-log10(pvalue)))+geom_point(aes(x=Compare2,y=ID,color=NES2,size = -log10(pvalue2)))+
  facet_wrap(patient~SigType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_gradient2(low='blue',mid='white',high='red')
#####################################Type1:diff markers for each subtypes#############################
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')

patients = as.character(unique(alldata$Patient))

patients.markers = list()

# top5 = res.markers %>% group_by(Compare) %>% top_n(-5, p_val_adj)  %>% top_n(5, avg_logFC)
# DotPlot(sub.sub.data, features = as.character(top5$Gene),assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() 

for(i in 1:length(patients)){
  patient = patients[i]
  print(patient)
  sub.data = alldata[,alldata$Patient == patient]
  Idents(sub.data)=sub.data$orig.ident
  cellTypes = as.character(unique(sub.data$CellTypeManully))
  orig.idents = as.character(unique(sub.data$orig.ident))
  res.i = data.frame()
  for(j in 1:length(cellTypes)){
    cellType = cellTypes[j]
    print(cellType)
    sub.sub.data = sub.data[,sub.data$CellTypeManully == cellType]
    res.markers.TI= FindMarkers(sub.sub.data,
                                ident.1=orig.idents[grepl('TI',orig.idents)],
                                ident.2=orig.idents[grepl('NI',orig.idents)],
                                assay = 'RNA',
                                logfc.threshold=-Inf,
                                min.pct = -Inf)
    res.markers.TI$Gene = rownames(res.markers.TI)
    other.idents = orig.idents[!grepl('I',orig.idents)]
    res.markers.TOthers= FindMarkers(sub.sub.data,
                                ident.1=other.idents[grepl('T',other.idents)],
                                ident.2=other.idents[grepl('N',other.idents)],
                                assay = 'RNA',
                                logfc.threshold=-Inf,
                                min.pct = -Inf)
    res.markers.TOthers$Gene = rownames(res.markers.TOthers)
    
    res.markers.TI$Compare='TIVSNI'
    res.markers.TOthers$Compare=paste(unique(substr(other.idents,1,2)),collapse = 'VS')
    res.markers = rbind(res.markers.TI,res.markers.TOthers)
    res.markers$cellType = cellType
    res.markers$patient = patient
    res.i=rbind(res.i,res.markers)
  }

  patients.markers[[i]]=res.i
}
saveRDS(patients.markers,file='variables/patients.diff.markers.Rds')

##############Type1: T vs N on different sites##################
patients.markers = readRDS(file='variables/patients.diff.markers.Rds')
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
#####For each patient####
i=1##Change from 1 to 4
res.markers = patients.markers[[i]]
res.markers = filter(res.markers,p_val_adj<0.01, pct.1>0.5 | pct.2>0.5,abs(avg_logFC)>1)
top5 <- res.markers %>%group_by(cellType,Compare) %>% top_n(-5, p_val_adj)  %>% top_n(5, avg_logFC)
top5[top5$p_val_adj <= 1e-200,'p_val_adj']=1e-200

ggplot(top5,aes(x=Compare,y=Gene))+geom_point(aes(size = -log10(p_val_adj),color=avg_logFC))+
  theme_cowplot()+scale_color_gradient2(low='blue',mid = 'white',high='red')+
  facet_wrap(~cellType,scales = 'free',nrow = 1)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='')

res.markers = patients.markers[[i]]
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
res.markers.TI$tag = paste(res.markers.TI$Gene,res.markers.TI$cellType)
res.markers.TOthers$tag2 = paste(res.markers.TOthers$Gene2,res.markers.TOthers$cellType2)

res.markers.2 = merge(res.markers.TI,res.markers.TOthers,by.x='tag',by.y = 'tag2')
nrow(filter(res.markers.2,Gene != Gene2 | cellType != cellType2))
res.markers.2.sig = filter(res.markers.2,p_val_adj<0.01|p_val_adj2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & abs(res.markers.2.sig$avg_logFC)>0.25,'SigType']='TI'
res.markers.2.sig[res.markers.2.sig$p_val_adj2<0.01 & abs(res.markers.2.sig$avg_logFC2)>0.25,'SigType']='TOthers'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC>0,'SigType']='Same'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC<0,'SigType']='Diff'
res.markers.2.sig.top = res.markers.2.sig %>%group_by(cellType,Compare) %>% top_n(5, abs(avg_logFC2)+abs(avg_logFC)) 
  
ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~cellType,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
            mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
            inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

####Four patients together####

res.markers = Reduce(rbind,patients.markers)
res.markers$tag = paste(res.markers$Gene,res.markers$cellType,res.markers$patient)
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TI,res.markers.TOthers,by.x='tag',by.y = 'tag2')
nrow(filter(res.markers.2,Gene != Gene2))
res.markers.2.sig = filter(res.markers.2,p_val_adj<0.01|p_val_adj2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & abs(res.markers.2.sig$avg_logFC)>0.25,'SigType']='TI'
res.markers.2.sig[res.markers.2.sig$p_val_adj2<0.01 & abs(res.markers.2.sig$avg_logFC2)>0.25,'SigType']='TOthers'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC>0,'SigType']='Same'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC<0,'SigType']='Diff'
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,cellType,Compare) %>% top_n(5, abs(avg_logFC2)+abs(avg_logFC)) 
###Scatter plot
sigType.cols = pal_npg('nrc')(5)
names(sigType.cols)=c('Same','NotSig','Diff','TI','TOthers')
res.markers.2.sig$cellType = factor(res.markers.2.sig$cellType,
                                    levels = c('Epithelial cell','Endothelial cell','Fibroblast cell',
                                               'T&NK cell','B cell','Myeloid cell','Mast cell'))
ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_grid(patient~cellType,scales = 'free')+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='TIVSNI',y='TOthersVSNOthers')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_manual(values = sigType.cols)
##file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/ScatterPlot T vs paired N diff genes across cell types_P1234.pdf

####Barplot to show summary about the above scatter plot

ggplot(res.markers.2.sig[res.markers.2.sig$SigType != 'NotSig',],
       aes(x=SigType,fill=SigType))+geom_bar(width = 0.6)+
  facet_grid(patient~cellType)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+scale_y_log10()+
  scale_fill_manual(values = sigType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/Barplot summary about TI VS TO gene number.pdf
res.sig.same = res.markers.2.sig[res.markers.2.sig$SigType == 'Same',]
genes.sig.same = as.character(unique(res.sig.same$Gene))
genes.sig.NP = unlist(lapply(genes.sig.same, function(a){
  max(c(length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC>0,'patient'])),
        length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC<0,'patient']))))
  
}))
genes.sig.NC = unlist(lapply(genes.sig.same, function(a){
  length(unique(res.sig.same[res.sig.same$Gene==a,'cellType']))
}))
genes.sig.avgS = unlist(lapply(genes.sig.same, function(a){
  mean(res.sig.same[res.sig.same$Gene==a,'avg_logFC']+res.sig.same[res.sig.same$Gene==a,'avg_logFC2'])
}))
gene.n.sum = data.frame(Gene = genes.sig.same,
                        NP = genes.sig.NP,
                        NC = genes.sig.NC,
                        avgS = genes.sig.avgS)
top.genes = as.character(gene.n.sum[gene.n.sum$NP==4 & gene.n.sum$NC==1,'Gene'])
res.sig.same.top = res.sig.same[res.sig.same$Gene %in% top.genes,]


res.sig.diff = res.markers.2.sig[res.markers.2.sig$SigType %in% c('Diff','TI','TOthers'),]
genes.sig.diff = as.character(unique(res.sig.diff$Gene))
genes.sig.NP = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a,'patient']))
}))
genes.sig.NC = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a,'cellType']))
}))
genes.sig.avgS = unlist(lapply(genes.sig.diff, function(a){
  mean(abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC'])+abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC2']))
}))
gene.n.sum.diff = data.frame(Gene = genes.sig.diff,
                        NP = genes.sig.NP,
                        NC = genes.sig.NC,
                        avgS = genes.sig.avgS)
top.genes.diff = as.character(gene.n.sum.diff[gene.n.sum.diff$NP==4 & gene.n.sum.diff$NC==1,'Gene'])
res.sig.diff.top = res.sig.diff[res.sig.diff$Gene %in% top.genes.diff,]


ggplot(rbind(res.sig.same.top,res.sig.diff.top),aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(rbind(res.sig.same.top,res.sig.diff.top),
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='TIVSNI',y='TOthersVSNOthers')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_manual(values = sigType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/ScatterPlot T vs paired N diff genes across cell types_P1234_topGenes.pdf


#############Type1 Diff pathway analysis,2021-09-01#################

library(clusterProfiler)
library(org.Hs.eg.db)
gene.mapping = AnnotationDbi::select(org.Hs.eg.db,keys = unique(as.character(res.markers.2$Gene)),keytype = 'SYMBOL',
                                     columns = c('SYMBOL','ENTREZID'))
top.genes.same.id = as.character(unique(gene.mapping[gene.mapping$SYMBOL %in% top.genes,'ENTREZID']))
same.paths = enrichKEGG(top.genes.same.id,pvalueCutoff = 1.2,minGSSize = 2)@result
same.gos = enrichGO(top.genes.same.id,'org.Hs.eg.db',ont = 'BP',pvalueCutoff = 1.2,minGSSize = 2,readable = T)@result

top.genes.diff.id = as.character(unique(gene.mapping[gene.mapping$SYMBOL %in% top.genes.diff,'ENTREZID']))
diff.paths = enrichKEGG(top.genes.diff.id,pvalueCutoff = 1.2,minGSSize = 2)@result
diff.gos = enrichGO(top.genes.diff.id,'org.Hs.eg.db',ont = 'BP',pvalueCutoff = 1.2,minGSSize = 2,readable = T)@result

save.image(file = 'findDiff.RData')

load(file = 'E:/project/metaboliteProteinInteraction/variables/KEGG_path2Gene.RData')
path2cate <- read.csv(file='E:/data/KEGG/KEGG hsa pathway list.csv',row.names = 1)
rownames(path2cate) = path2cate$Name
patients.marker.paths = list()

for(i in 1:length(patients.markers)){
  res.i = patients.markers[[i]]
  patient = unique(res.i$patient)
  cellTypes = as.character(unique(res.i$cellType))
  res.i.TI = res.i[res.i$Compare=='TIVSNI',]
  res.i.TO = res.i[res.i$Compare!='TIVSNI',]
  path.res.i = data.frame()
  print(patient)
  for(j in 1:length(cellTypes)){
    cellType = cellTypes[j]
    print(cellType)
    res.i.TI.cell = res.i.TI[res.i.TI$cellType == cellType,]
    FC = res.i.TI.cell$avg_logFC
    names(FC)= as.character(res.i.TI.cell$Gene)
    FC = FC[order(FC,decreasing = T)]
    if(length(FC[FC==0])>=21261){
      path.TI=data.frame()
    }else{
      path.TI = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
      
    }
    
    res.i.TO.cell = res.i.TO[res.i.TO$cellType == cellType,]
    FC = res.i.TO.cell$avg_logFC
    names(FC)= as.character(res.i.TO.cell$Gene)
    FC = FC[order(FC,decreasing = T)]
    path.TO = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
    
    if(nrow(path.TI!=0)){path.TI$Compare='TIVSNI'}
    path.TO$Compare=unique(res.i.TO$Compare)
    path.2 = rbind(path.TI,path.TO)
    path.2$cellType = cellType
    path.2$patient = patient
    path.res.i=rbind(path.res.i,path.2)
    
  }
  
  patients.marker.paths[[i]]=path.res.i
}
saveRDS(patients.marker.paths,file='variables/patients.diff.marker.pathways.Rds')

patients.marker.paths = readRDS(file='variables/patients.diff.marker.pathways.Rds')
res.markers = Reduce(rbind,patients.marker.paths)
res.markers$tag = paste(res.markers$ID,res.markers$cellType,res.markers$patient)
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TOthers,res.markers.TI,by.x='tag2',by.y = 'tag')
nrow(filter(res.markers.2,ID != ID))
res.markers.2.sig = filter(res.markers.2,pvalue<0.01|pvalue2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & res.markers.2.sig$pvalue2>=0.01,'SigType']='TI'
res.markers.2.sig[res.markers.2.sig$pvalue2<0.01 & res.markers.2.sig$pvalue>=0.01,'SigType']='TOthers'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2>0,'SigType']='Same'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2<0,'SigType']='Diff'

disease.pathways = path2cate[path2cate$cate=='6. Human Diseases','Name']
res.markers.2.sig = res.markers.2.sig[res.markers.2.sig$ID %in% as.character(disease.pathways)==F,]
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,cellType,Compare) %>% top_n(3, abs(NES)+abs(NES2)) 
res.markers.2.sig$cate = path2cate[res.markers.2.sig$ID,'cate']
res.markers.2.sig$cellType = unlist(lapply(res.markers.2.sig$cellType,function(a){
  strsplit(a,' ')[[1]][1]
}))
res.markers.2.sig.top=res.markers.2.sig[res.markers.2.sig$SigType != 'Same',]
res.markers.2.sig.top = res.markers.2.sig.top %>%group_by(patient,cellType,Compare) %>% top_n(3, abs(NES)+abs(NES2)) 

#####Summary plot1
ggplot(res.markers.2.sig,aes(x=SigType,fill=cate))+geom_bar(width = 0.6)+
  facet_grid(patient~cellType,scales = 'free_x',space = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_fill_npg()

####Summary plot 2

res.markers.2.sig.v2 = res.markers.2.sig

res.markers.2.sig.v2$SigType[res.markers.2.sig.v2$SigType != 'Same']='Diff'

res.markers.2.sig.v2 = group_by(res.markers.2.sig.v2,SigType,cellType,patient,cate) %>% summarise(num = n() )

ggbarplot(res.markers.2.sig.v2,x='SigType',y = 'num', add = 'mean_sd',fill='cate',facet.by = c('cate','cellType')) + 
  theme_bw() + 
  scale_fill_manual(values = pal_npg('nrc')(5)) + ylab("Number") +
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/summary barplot pathway cate 2.pdf
ggbarplot(res.markers.2.sig.v2,x='SigType',
          y = 'num', add = 'mean_sd',fill='cate',
          facet.by = c('cate')) + 
  theme_bw() +  stat_compare_means(aes(label = ..p.signif..),label.y = 7,method = 'wilcox.test',label.x = 1.5)+
  scale_fill_manual(values = pal_npg('nrc')(5)) + ylab("Number") +  ylim(c(0,7.5))+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/summary barplot pathway cate summary.pdf

res.markers.2.sig$subcate = path2cate[res.markers.2.sig$ID,'subcate']
res.markers.2.sig$SigType[res.markers.2.sig$SigType!='Same']='Diff'
ggplot(res.markers.2.sig[res.markers.2.sig$cate == '5. Organismal Systems',],
       aes(x=SigType,fill=subcate))+
  geom_bar(width = 0.6)+
  facet_grid(patient~cellType,scales = 'free_x',space = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_fill_npg()
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient&cellType level T vs paired N/Barplot summary about TI VS TO pathways across cell types_only cate5 subcates.pdf
###Scatter plot
ggplot(res.markers.2.sig,aes(x=NES,y=NES2,color=SigType))+geom_point(size=1)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=NES,y=NES2,color=SigType,label=ID),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

####Point plot
res.markers.2.sig.top = res.markers.2.sig.top[order(res.markers.2.sig.top$cate),]
res.markers.2.sig.top$ID = factor(as.character(res.markers.2.sig.top$ID),unique(as.character(res.markers.2.sig.top$ID)))
res.markers.2.sig.top$cate = unlist(lapply(as.character(res.markers.2.sig.top$cate),function(a){
  paste0('C',strsplit(a,'. ',fixed = T)[[1]][1])
}))
ggplot(res.markers.2.sig.top)+
  geom_point(aes(x=Compare,y=ID,color=NES,size=-log10(pvalue)))+
  geom_point(aes(x=Compare2,y=ID,color=NES2,size = -log10(pvalue2)))+
  facet_grid(patient~cate~cellType,scales = 'free',space = 'free')+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8)
        )+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_gradient2(low='blue',mid='white',high='red')+
  scale_size_binned(range = c(0.2,4))
#########Type 2 diff markers: TI VS TOthers, NI VS NOthers 2021-09-02####################
patients = as.character(unique(alldata$Patient))
patients.markers = list()

# top5 = res.markers %>% group_by(Compare) %>% top_n(-5, p_val_adj)  %>% top_n(5, avg_logFC)
# DotPlot(sub.sub.data, features = as.character(top5$Gene),assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() 

for(i in 1:length(patients)){
  patient = patients[i]
  print(patient)
  sub.data = alldata[,alldata$Patient == patient]
  Idents(sub.data)=sub.data$orig.ident
  cellTypes = as.character(unique(sub.data$CellTypeManully))
  orig.idents = as.character(unique(sub.data$orig.ident))
  T.idents = orig.idents[substr(orig.idents,1,1)=='T']
  N.idents = orig.idents[substr(orig.idents,1,1)=='N']
  TOther.idents = T.idents[!grepl('TI',T.idents)]
  NOther.idents = N.idents[!grepl('NI',N.idents)]
  
  res.i = data.frame()
  for(j in 1:length(cellTypes)){
    cellType = cellTypes[j]
    print(cellType)
    sub.sub.data = sub.data[,sub.data$CellTypeManully == cellType]
    res.markers.T= FindMarkers(sub.sub.data,
                                ident.1=orig.idents[grepl('TI',orig.idents)],
                                ident.2=TOther.idents,
                                assay = 'RNA',
                                logfc.threshold=-Inf,
                                min.pct = -Inf)
    res.markers.T$Gene = rownames(res.markers.T)
   
    res.markers.N= FindMarkers(sub.sub.data,
                              ident.1=orig.idents[grepl('NI',orig.idents)],
                              ident.2=NOther.idents,
                              assay = 'RNA',
                              logfc.threshold=-Inf,
                              min.pct = -Inf)
    res.markers.N$Gene = rownames(res.markers.N)
    
    res.markers.T$Compare=paste(c('TI',unique(substr(TOther.idents,1,2))),collapse = 'VS')
    res.markers.N$Compare=paste(c('NI',unique(substr(NOther.idents,1,2))),collapse = 'VS')
    res.markers = rbind(res.markers.T,res.markers.N)
    res.markers$cellType = cellType
    res.markers$patient = patient
    res.i=rbind(res.i,res.markers)
  }
  
  patients.markers[[i]]=res.i
}
saveRDS(patients.markers,file='variables/patients.diff.markers_Type2.Rds')

patients.markers = readRDS(file='variables/patients.diff.markers_Type2.Rds')
####
res.markers = Reduce(rbind,patients.markers)
res.markers$tag = paste(res.markers$Gene,res.markers$cellType,res.markers$patient)
res.markers.T = res.markers[grepl('T',res.markers$Compare) ,]
res.markers.N = res.markers[grepl('N',res.markers$Compare) ,]
colnames(res.markers.N)= paste0(colnames(res.markers.N),2)
length(unique(res.markers.T$tag))
length(unique(res.markers.N$tag2))
res.markers.2 = merge(res.markers.T,res.markers.N,by.x='tag',by.y = 'tag2')
nrow(filter(res.markers.2,Gene != Gene2))
res.markers.2.sig = filter(res.markers.2,p_val_adj<0.01|p_val_adj2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & abs(res.markers.2.sig$avg_logFC)>0.25,'SigType']='TOnly'
res.markers.2.sig[res.markers.2.sig$p_val_adj2<0.01 & abs(res.markers.2.sig$avg_logFC2)>0.25,'SigType']='NOnly'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC>0,'SigType']='TNSame'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC<0,'SigType']='TNDiff'
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,cellType,Compare) %>% top_n(5, abs(avg_logFC2)+abs(avg_logFC)) 
###Scatter plot to compare Tdiff and NDiff
ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,5,4)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)
###Barplot summary
ggplot(res.markers.2.sig,aes(x=SigType,fill=SigType))+geom_bar(width = 0.6)+
  facet_wrap(patient~cellType,nrow =4)+
  theme_cowplot()+scale_fill_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+scale_y_log10()
###Further find
res.sig.same = res.markers.2.sig[res.markers.2.sig$SigType == 'TNSame',]
genes.sig.same = as.character(unique(res.sig.same$Gene))
genes.sig.NP = unlist(lapply(genes.sig.same, function(a){
  max(c(length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC>0,'patient'])),
        length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC<0,'patient']))))
  
}))
genes.sig.NC = unlist(lapply(genes.sig.same, function(a){
  length(unique(res.sig.same[res.sig.same$Gene==a,'cellType']))
}))
genes.sig.avgS = unlist(lapply(genes.sig.same, function(a){
  mean(res.sig.same[res.sig.same$Gene==a,'avg_logFC']+res.sig.same[res.sig.same$Gene==a,'avg_logFC2'])
}))
gene.n.sum = data.frame(Gene = genes.sig.same,
                        NP = genes.sig.NP,
                        NC = genes.sig.NC,
                        avgS = genes.sig.avgS)


res.sig.diff = res.markers.2.sig[res.markers.2.sig$SigType %in% c('TOnly','TNDiff'),]
genes.sig.diff = as.character(unique(res.sig.diff$Gene))
genes.sig.NP = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a,'patient']))
}))
genes.sig.NC = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a & grepl('T',res.sig.diff$Compare),'cellType']))
}))
genes.sig.avgS = unlist(lapply(genes.sig.diff, function(a){
  mean(abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC'])+abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC2']))
}))
gene.n.sum.diff = data.frame(Gene = genes.sig.diff,
                             NP = genes.sig.NP,
                             NC = genes.sig.NC,
                             avgS = genes.sig.avgS)
top.genes.diff = as.character(gene.n.sum.diff[gene.n.sum.diff$NP==4 & gene.n.sum.diff$NC<=2,'Gene'])
res.sig.diff.top = res.sig.diff[res.sig.diff$Gene %in% top.genes.diff,]

###
ggplot(rbind(res.sig.same.top,res.sig.diff.top),aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(res.sig.diff.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)
## Four patients "SOD3" "TKT" 


#############Type2 Diff pathway analysis,2021-09-10#################

patients.marker.paths = list()

for(i in 1:length(patients.markers)){
  res.i = patients.markers[[i]]
  patient = unique(res.i$patient)
  cellTypes = as.character(unique(res.i$cellType))
  res.i.TI = res.i[grepl('T',res.i$Compare),]
  res.i.TO = res.i[grepl('N',res.i$Compare),]
  path.res.i = data.frame()
  print(patient)
  for(j in 1:length(cellTypes)){
    cellType = cellTypes[j]
    print(cellType)
    res.i.TI.cell = res.i.TI[res.i.TI$cellType == cellType,]
    FC = res.i.TI.cell$avg_logFC
    names(FC)= as.character(res.i.TI.cell$Gene)
    FC = FC[order(FC,decreasing = T)]
    if(length(FC[FC==0])>=21261){
      path.TI=data.frame()
    }else{
      path.TI = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
      
    }
    
    res.i.TO.cell = res.i.TO[res.i.TO$cellType == cellType,]
    FC = res.i.TO.cell$avg_logFC
    names(FC)= as.character(res.i.TO.cell$Gene)
    FC = FC[order(FC,decreasing = T)]
    if(length(FC[FC==0])>=20000){
      path.TO=data.frame()
    }else{
      path.TO = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
    }
      
    if(nrow(path.TI)>0){path.TI$Compare=unique(res.i.TI$Compare)}
    if(nrow(path.TO)>0){path.TO$Compare=unique(res.i.TO$Compare)}
    path.2 = rbind(path.TI,path.TO)
    path.2$cellType = cellType
    path.2$patient = patient
    path.res.i=rbind(path.res.i,path.2)
    
  }
  patients.marker.paths[[i]]=path.res.i
}
saveRDS(patients.marker.paths,file='variables/patients.diff.marker.pathways_Type2.Rds')

patients.marker.paths = readRDS(file='variables/patients.diff.marker.pathways_Type2.Rds')
res.markers = Reduce(rbind,patients.marker.paths)
res.markers$tag = paste(res.markers$ID,res.markers$cellType,res.markers$patient)
res.markers.TI = res.markers[grepl('T',res.markers$Compare),]
res.markers.TOthers = res.markers[grepl('N',res.markers$Compare),]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TOthers,res.markers.TI,by.x='tag2',by.y = 'tag')
nrow(filter(res.markers.2,ID != ID))
res.markers.2.sig = filter(res.markers.2,pvalue<0.01|pvalue2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & res.markers.2.sig$pvalue2>=0.01,'SigType']='TOnly'
res.markers.2.sig[res.markers.2.sig$pvalue2<0.01 & res.markers.2.sig$pvalue>=0.01,'SigType']='NOnly'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2>0,'SigType']='TNSame'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2<0,'SigType']='TNDiff'

disease.pathways = path2cate[path2cate$cate=='6. Human Diseases','Name']
res.markers.2.sig = res.markers.2.sig[res.markers.2.sig$ID %in% as.character(disease.pathways)==F,]
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,cellType,Compare) %>% top_n(3, abs(NES)+abs(NES2)) 
res.markers.2.sig$cate = path2cate[res.markers.2.sig$ID,'cate']
###Scatter plot
ggplot(res.markers.2.sig,aes(x=NES,y=NES2,color=SigType))+geom_point(size=1)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=NES,y=NES2,color=SigType,label=ID),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

ggplot(res.markers.2.sig,aes(x=SigType,fill=cate))+geom_bar(width = 0.6)+
  facet_wrap(patient~cellType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

####Type 3, patient specific, TIVSNI, TO vs NO, 2021-09-13####
patients = as.character(unique(alldata$Patient))

patients.markers = list()

for(i in 1:length(patients)){
  patient = patients[i]
  print(patient)
  sub.data = alldata[,alldata$Patient == patient]
  Idents(sub.data)=sub.data$orig.ident
  cellTypes = as.character(unique(sub.data$CellTypeManully))
  orig.idents = as.character(unique(sub.data$orig.ident))
  
   
     res.markers.TI= FindMarkers(sub.data,
                                ident.1=orig.idents[grepl('TI',orig.idents)],
                                ident.2=orig.idents[grepl('NI',orig.idents)],
                                assay = 'RNA',
                                logfc.threshold=-Inf,
                                min.pct = -Inf)
    res.markers.TI$Gene = rownames(res.markers.TI)
    other.idents = orig.idents[!grepl('I',orig.idents)]
    res.markers.TOthers= FindMarkers(sub.data,
                                     ident.1=other.idents[grepl('T',other.idents)],
                                     ident.2=other.idents[grepl('N',other.idents)],
                                     assay = 'RNA',
                                     logfc.threshold=-Inf,
                                     min.pct = -Inf)
    res.markers.TOthers$Gene = rownames(res.markers.TOthers)
    
    res.markers.TI$Compare='TIVSNI'
    res.markers.TOthers$Compare=paste(unique(substr(other.idents,1,2)),collapse = 'VS')
    res.markers = rbind(res.markers.TI,res.markers.TOthers)
  patients.markers[[i]]=res.markers
}
for(i in 1:length(patients)){
  patient = patients[i]
  patients.markers[[i]]$patient = patient
}

saveRDS(patients.markers,file='variables/patients.diff.markers_Type3.Rds')
patients.markers = readRDS(file='variables/patients.diff.markers_Type3.Rds')
res.markers = Reduce(rbind,patients.markers)
res.markers$tag = paste(res.markers$Gene,res.markers$patient)
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TI,res.markers.TOthers,by.x='tag',by.y = 'tag2')
nrow(filter(res.markers.2,Gene != Gene2))
res.markers.2.sig = filter(res.markers.2,p_val_adj<0.01|p_val_adj2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & abs(res.markers.2.sig$avg_logFC)>0.25,'SigType']='TI'
res.markers.2.sig[res.markers.2.sig$p_val_adj2<0.01 & abs(res.markers.2.sig$avg_logFC2)>0.25,'SigType']='TOthers'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC>0,'SigType']='Same'
res.markers.2.sig[res.markers.2.sig$p_val_adj<0.01 & 
                    res.markers.2.sig$p_val_adj2<0.01 & 
                    abs(res.markers.2.sig$avg_logFC)>0.25& 
                    abs(res.markers.2.sig$avg_logFC2)>0.25 &
                    res.markers.2.sig$avg_logFC2*res.markers.2.sig$avg_logFC<0,'SigType']='Diff'
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,Compare) %>% top_n(5, abs(avg_logFC2)+abs(avg_logFC)) 
###Scatter plot
ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='TIVSNI',y='TOthersVSNOthers')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_manual(values = sigType.cols)
####Barplot to show summary about the above scatter plot
ggplot(res.markers.2.sig[res.markers.2.sig$SigType!='NotSig',],aes(x=SigType,fill=SigType))+geom_bar(width = 0.6)+
  facet_wrap(~patient,nrow =1)+
  theme_cowplot()+scale_fill_manual(values=sigType.cols[-2])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+scale_y_log10()
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient level T vs paired N/barplot for T vs paired N diff genes_only patients.pdf
res.sig.same = res.markers.2.sig[res.markers.2.sig$SigType == 'Same',]
genes.sig.same = as.character(unique(res.sig.same$Gene))
genes.sig.NP = unlist(lapply(genes.sig.same, function(a){
  max(c(length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC>0,'patient'])),
        length(unique(res.sig.same[res.sig.same$Gene==a & res.sig.same$avg_logFC<0,'patient']))))
  
}))

genes.sig.avgS = unlist(lapply(genes.sig.same, function(a){
  mean(res.sig.same[res.sig.same$Gene==a,'avg_logFC']+res.sig.same[res.sig.same$Gene==a,'avg_logFC2'])
}))
gene.n.sum = data.frame(Gene = genes.sig.same,
                        NP = genes.sig.NP,
                        avgS = genes.sig.avgS)
top.genes = as.character(gene.n.sum[gene.n.sum$NP>=4 ,'Gene'])
res.sig.same.top = res.sig.same[res.sig.same$Gene %in% top.genes,]



res.sig.diff = res.markers.2.sig[res.markers.2.sig$SigType %in% c('Diff','TI','TOthers'),]
genes.sig.diff = as.character(unique(res.sig.diff$Gene))
genes.sig.NP = unlist(lapply(genes.sig.diff, function(a){
  length(unique(res.sig.diff[res.sig.diff$Gene==a,'patient']))
}))

genes.sig.avgS = unlist(lapply(genes.sig.diff, function(a){
  mean(abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC'])+abs(res.sig.diff[res.sig.diff$Gene==a,'avg_logFC2']))
}))
gene.n.sum.diff = data.frame(Gene = genes.sig.diff,
                             NP = genes.sig.NP,
                             avgS = genes.sig.avgS)
top.genes.diff = as.character(gene.n.sum.diff[gene.n.sum.diff$NP==4 & gene.n.sum.diff$avgS>0.8 ,'Gene'])

res.sig.diff.top = res.sig.diff[res.sig.diff$Gene %in% top.genes.diff,]


ggplot(rbind(res.sig.same.top,res.sig.diff.top),aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(rbind(res.sig.same.top,res.sig.diff.top),
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)


ggplot(res.markers.2.sig,aes(x=avg_logFC,y=avg_logFC2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(rbind(res.sig.same.top,res.sig.diff.top),
                  mapping = aes(x=avg_logFC,y=avg_logFC2,color=SigType,label=Gene),size=2,
                  inherit.aes = F)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='TIVSNI',y='TOthersVSNOthers')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_manual(values = sigType.cols)
##file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/patient level T vs paired N/scatter plot for T vs paired N diff genes_label consis genes.pdf

#CD82 is also known to inhibit cancer invasion and metastasis in non-small-cell lung carcinoma (NSCLC) via 
#multiple mechanisms [ 57 ]. A higher expression of CD82 was reported in tumors which are better differentiated, less likely to metastasize to lymph nodes, and present at an earlier clinical stage in NSCLC [ 58 ].
#RGS1 expression is associated with poor prognosis in multiple myeloma

p1 = FeaturePlot(alldata[,alldata$Patient == 'P2'], features = c("ALOX15B", "CD79B"),split.by = 'orig.ident') + 
  NoLegend()

p2 = DimPlot(alldata[,alldata$Patient == 'P2'], group.by = 'CellTypeManully',split.by = 'orig.ident') 
p1 / p2
####Type3: pathway analysis####

patients.marker.paths = list()

for(i in 1:length(patients.markers)){
  res.i = patients.markers[[i]]
  patient = unique(res.i$patient)
  cellTypes = as.character(unique(res.i$cellType))
  res.i.TI = res.i[res.i$Compare=='TIVSNI',]
  res.i.TO = res.i[res.i$Compare!='TIVSNI',]
  
  
    FC = res.i.TI$avg_logFC
    names(FC)= as.character(res.i.TI$Gene)
    FC = FC[order(FC,decreasing = T)]
    
    path.TI = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
   
    
    
    FC = res.i.TO$avg_logFC
    names(FC)= as.character(res.i.TO$Gene)
    FC = FC[order(FC,decreasing = T)]
    path.TO = GSEA(FC,minGSSize = 3,pvalueCutoff = 1.2,TERM2GENE = path2Gene[,c(1,2)],seed=2223)@result
    
    if(nrow(path.TI!=0)){path.TI$Compare='TIVSNI'}
    path.TO$Compare=unique(res.i.TO$Compare)
    path.2 = rbind(path.TI,path.TO)
    path.2$patient = patient
   
  patients.marker.paths[[i]]=path.2
}
saveRDS(patients.marker.paths,file='variables/patients.diff.marker.pathways_Type3.Rds')

patients.marker.paths = readRDS(file='variables/patients.diff.marker.pathways_Type3.Rds')
res.markers = Reduce(rbind,patients.marker.paths)
res.markers$tag = paste(res.markers$ID,res.markers$patient)
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
colnames(res.markers.TOthers)= paste0(colnames(res.markers.TOthers),2)
length(unique(res.markers.TI$tag))
length(unique(res.markers.TOthers$tag2))
res.markers.2 = merge(res.markers.TOthers,res.markers.TI,by.x='tag2',by.y = 'tag')
nrow(filter(res.markers.2,ID != ID))
res.markers.2.sig = filter(res.markers.2,pvalue<0.01|pvalue2<0.01)
res.markers.2.sig$SigType = 'NotSig'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & res.markers.2.sig$pvalue2>=0.01,'SigType']='TI'
res.markers.2.sig[res.markers.2.sig$pvalue2<0.01 & res.markers.2.sig$pvalue>=0.01,'SigType']='TOthers'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2>0,'SigType']='Same'
res.markers.2.sig[res.markers.2.sig$pvalue<0.01 & 
                    res.markers.2.sig$pvalue2<0.01 & 
                    res.markers.2.sig$NES*res.markers.2.sig$NES2<0,'SigType']='Diff'

disease.pathways = path2cate[path2cate$cate=='6. Human Diseases','Name']
res.markers.2.sig = res.markers.2.sig[res.markers.2.sig$ID %in% as.character(disease.pathways)==F,]
res.markers.2.sig.top = res.markers.2.sig %>%group_by(patient,Compare) %>% top_n(3, abs(NES)+abs(NES2)) 
res.markers.2.sig$cate = path2cate[res.markers.2.sig$ID,'cate']


#####Summary plot
ggplot(res.markers.2.sig,aes(x=SigType,fill=cate))+geom_bar(width = 0.6)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

###Scatter plot
ggplot(res.markers.2.sig,aes(x=NES,y=NES2,color=SigType))+geom_point(size=1)+
  facet_wrap(~patient,scales = 'free',nrow = 1)+
  theme_cowplot()+
  geom_text_repel(res.markers.2.sig.top,
                  mapping = aes(x=NES,y=NES2,color=SigType,label=ID),size=2,
                  inherit.aes = F)+scale_color_manual(values=brewer.pal(5,'Set1')[c(2,3,1,4,5)])+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)

###Point plot
ggplot(res.markers.2.sig)+geom_point(aes(x=Compare,y=ID,color=NES,size=-log10(pvalue)))+geom_point(aes(x=Compare2,y=ID,color=NES2,size = -log10(pvalue2)))+
  facet_wrap(patient~SigType,scales = 'free',nrow = 4)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='')+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)+
  scale_color_gradient2(low='blue',mid='white',high='red')

####PCA analysis of cell type compositions####
com2orig = table(alldata$orig.ident,alldata$CellTypeManully)
com2orig = com2orig/apply(com2orig, 1, sum)
mat = t(com2orig)
mat = mat[apply(mat, 1, sum)>0,]
mat = scale(t(mat))
mat = mat[apply(mat, 1, sd)!=0,]
mat[is.na(mat)]

pcr.int = prcomp(mat)
pcr.int$sdev/sum(pcr.int$sdev)#[1] 2.324118e-01 2.267851e-01 1.662813e-01 1.512029e-01 1.165619e-01 1.067569e-01 4.133088e-17
library(RColorBrewer)
pcr.cluster = data.frame(From = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][3]
})),
name = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][1]
})),
type = unlist(lapply(rownames(mat), function(a){
  substr(a,1,1)
})),
row.names  = rownames(mat),
PC1 = pcr.int$x[,1],
PC2 = pcr.int$x[,2],
PC3 = pcr.int$x[,3])

ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=From,shape=type),size = 3)+
  scale_color_npg()+
  theme_bw()+geom_text_repel(aes(label = name),size=3,nudge_y = 0.1,nudge_x = 0.1)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/PCA plot for cellTypeCom2orig.pdf

# library(gg3D)
# ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2,z=PC3,color=From,shape=type))+
#   scale_color_npg()+
#   theme_void()+
#   axes_3D() +
#   stat_3D()


contributions = abs(pcr.int$rotation[,c(1,2)])
contributions = sweep(contributions, 2, colSums(contributions), "/")

contributions = data.frame(contributions)
contributions$subType = rownames(contributions)
contributions = contributions[order(contributions$PC1,decreasing = T),]
contributions$subType = factor(contributions$subType,levels=as.character(contributions$subType))
p1= ggplot(contributions,aes(x=subType,y=PC1))+geom_bar(stat = 'identity',width = 0.65,fill='blue')+
  theme_cowplot()+
  coord_flip()+ylim(0,0.5)

contributions = contributions[order(contributions$PC2,decreasing = T),]
contributions$subType = factor(contributions$subType,levels=as.character(contributions$subType))

p2= ggplot(contributions,aes(x=subType,y=PC2))+geom_bar(stat = 'identity',width = 0.65,fill='blue')+
  theme_cowplot()+
  coord_flip()+ylim(0,0.5)
p1 | p2
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/PCA contributions PC1 PC2.pdf
####Differetial cell compositions####
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

library(ggpubr)
cellType.cols = pal_npg('nrc')(length(colnames(com2orig)))
names(cellType.cols)=c('Epithelial cell','Endothelial cell','Fibroblast cell',
                       'T&NK cell','B cell','Myeloid cell','Mast cell')
diff.ratios$CellType = factor(diff.ratios$CellType,levels =c('Epithelial cell','Endothelial cell','Fibroblast cell',
                                                             'T&NK cell','B cell','Myeloid cell','Mast cell') )
#ptest <- diff.ratios %>% group_by(CellType) %>% summarize(p.value = kruskal.test(Diff ~Type)$p.value)
ggbarplot(diff.ratios,
          x='CellType',y = 'Diff',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  theme_cowplot() + 
  ylab("Diff Ratio") + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/Differetial ratios of T & N samples barplot.pdf
# ggbarplot(diff.ratios,x='CellType',y = 'Diff', add = 'mean_sd',fill='CellType') + 
#   theme_bw() + 
#   scale_fill_manual(values = cellType.cols) + ylab("Diff Ratio") + 
#   theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8))
# #file:///H:/project/single cell/MPSC&ST/figure results/four samples/Differetial ratios of the samples across patients barplot.pdf

library(reshape2)
library(ggpubr)
ratio.mat = melt(com2orig,measure.vars = colnames(com2orig))
ratio.mat$Type = substr(ratio.mat$Var1,1,1)
ptest2 <- ratio.mat %>% group_by(Var2) %>% summarize(p.value = kruskal.test(value ~Type)$p.value)
ratio.mat$Var2 = factor(ratio.mat$Var2,levels =c('Epithelial cell','Endothelial cell','Fibroblast cell',
                                                          'T&NK cell','B cell','Myeloid cell','Mast cell'))

ggbarplot(ratio.mat,
          x='Var2',y = 'value',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  geom_text(data = ptest2,mapping = aes(x=Var2,y=0.7,label = ifelse(signif(p.value,2)<0.01,'**','')),inherit.aes = F,size=3)+
  theme_cowplot() + 
  ylab("Proportion") + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/com2orig cell composition T vs N barplot.pdf