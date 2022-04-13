##########Test for Spatial transcriptomics£¬R 4.0 2021-07-14######################
####Refer to :https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
########By Di Chen#################

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggsci)
setwd('H:/Project/single cell/scRNAData/')


#############Read input data and images########################
ni.p1.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_1NL/filtered_feature_bc_matrix")
ni.p1.st <- CreateSeuratObject(ni.p1.st, project = "NI_P1",assay = 'Spatial')
ni.p1.st$slice <- 1
ni.p1.st$region <- 'P1_NIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_1NL/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ni.p1.st)]
ni.p1.st[['image']] <- img

ti.p1.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_1TL/filtered_feature_bc_matrix")
ti.p1.st <- CreateSeuratObject(ti.p1.st, project = "TI_P1",assay = 'Spatial')
ti.p1.st$slice <- 1
ti.p1.st$region <- 'P1_TIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_1TL/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ti.p1.st)]
ti.p1.st[['image']] <- img

nm.p1.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_1NM/filtered_feature_bc_matrix")
nm.p1.st <- CreateSeuratObject(nm.p1.st, project = "NM_P1",assay = 'Spatial')
nm.p1.st$slice <- 1
nm.p1.st$region <- 'P1_NMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_1NM/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = nm.p1.st)]
nm.p1.st[['image']] <- img

tm.p1.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_1TM/filtered_feature_bc_matrix")
tm.p1.st <- CreateSeuratObject(tm.p1.st, project = "TM_P1",assay = 'Spatial')
tm.p1.st$slice <- 1
tm.p1.st$region <- 'P1_TMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_1TM/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = tm.p1.st)]
tm.p1.st[['image']] <- img


ni.p2.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_5NL/filtered_feature_bc_matrix")
ni.p2.st <- CreateSeuratObject(ni.p2.st, project = "NI_P2",assay = 'Spatial')
ni.p2.st$slice <- 1
ni.p2.st$region <- 'P2_NIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_5NL/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ni.p2.st)]
ni.p2.st[['image']] <- img

ti.p2.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_5TL/filtered_feature_bc_matrix")
ti.p2.st <- CreateSeuratObject(ti.p2.st, project = "TI_P2",assay = 'Spatial')
ti.p2.st$slice <- 1
ti.p2.st$region <- 'P2_TIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_5TL/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ti.p2.st)]
ti.p2.st[['image']] <- img

nm.p2.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_5NM/filtered_feature_bc_matrix")
nm.p2.st <- CreateSeuratObject(nm.p2.st, project = "NM_P2",assay = 'Spatial')
nm.p2.st$slice <- 1
nm.p2.st$region <- 'P2_NMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_5NM/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = nm.p2.st)]
nm.p2.st[['image']] <- img

tm.p2.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/D_5TM/filtered_feature_bc_matrix")
tm.p2.st <- CreateSeuratObject(tm.p2.st, project = "TM_P2",assay = 'Spatial')
tm.p2.st$slice <- 1
tm.p2.st$region <- 'P2_TMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/D_5TM/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = tm.p2.st)]
tm.p2.st[['image']] <- img



ni.p3.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/P8_I_N/filtered_feature_bc_matrix")
ni.p3.st <- CreateSeuratObject(ni.p3.st, project = "NI_P3",assay = 'Spatial')
ni.p3.st$slice <- 1
ni.p3.st$region <- 'P3_NIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/P8_I_N/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ni.p3.st)]
ni.p3.st[['image']] <- img

ti.p3.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/P8_I_T/filtered_feature_bc_matrix")
ti.p3.st <- CreateSeuratObject(ti.p3.st, project = "TI_P3",assay = 'Spatial')
ti.p3.st$slice <- 1
ti.p3.st$region <- 'P3_TIP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/P8_I_T/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = ti.p3.st)]
ti.p3.st[['image']] <- img


nm.p3.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/P8_M_N/filtered_feature_bc_matrix")
nm.p3.st <- CreateSeuratObject(nm.p3.st, project = "NM_P3",assay = 'Spatial')
nm.p3.st$slice <- 1
nm.p3.st$region <- 'P3_NMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/P8_M_N/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = nm.p3.st)]
nm.p3.st[['image']] <- img

tm.p3.st = Seurat::Read10X(data.dir = "../ourData/stRNA/Count/P8_M_T/filtered_feature_bc_matrix")
tm.p3.st <- CreateSeuratObject(tm.p3.st, project = "TM_P3",assay = 'Spatial')
tm.p3.st$slice <- 1
tm.p3.st$region <- 'P3_TMP'
img <- Seurat::Read10X_Image(image.dir = '../ourData/stRNA/Count/P8_M_T/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = tm.p3.st)]
tm.p3.st[['image']] <- img
#####################             QC         ###########################
st.m = merge(ni.p1.st,c(ti.p1.st,nm.p1.st,tm.p1.st,
                        ni.p2.st,ti.p2.st,nm.p2.st,tm.p2.st,
                        ni.p3.st,ti.p3.st,nm.p3.st,tm.p3.st ),
             add.cell.ids = c('NIP1','TIP1','NMP1','TMP1',
                              'NIP2','TIP2','NMP2','TMP2',
                              "NIP3",'TIP3','NMP3','TMP3'))
names(st.m@images) = c('NIP1','TIP1','NMP1','TMP1',
                       'NIP2','TIP2','NMP2','TMP2',
                       "NIP3",'TIP3','NMP3','TMP3')

st.m <- PercentageFeatureSet(st.m, "^MT-", col.name = "percent_mito")
st.m <- PercentageFeatureSet(st.m, "^HB[^(P)]", col.name = "percent_hb")

P1=VlnPlot(st.m, features = "nCount_Spatial", pt.size = 0) + NoLegend()
P2=VlnPlot(st.m, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
P3=VlnPlot(st.m, features = "percent_mito", pt.size = 0)+ NoLegend()
P4=VlnPlot(st.m, features = "percent_hb", pt.size = 0)+ NoLegend()

P1 / P2 / P3 / P4

ni.p1.st <- PercentageFeatureSet(ni.p1.st, "^MT-", col.name = "percent_mito")
ti.p1.st <- PercentageFeatureSet(ti.p1.st, "^MT-", col.name = "percent_mito")
nm.p1.st <- PercentageFeatureSet(nm.p1.st, "^MT-", col.name = "percent_mito")
tm.p1.st <- PercentageFeatureSet(tm.p1.st, "^MT-", col.name = "percent_mito")
ni.p2.st <- PercentageFeatureSet(ni.p2.st, "^MT-", col.name = "percent_mito")
ti.p2.st <- PercentageFeatureSet(ti.p2.st, "^MT-", col.name = "percent_mito")
nm.p2.st <- PercentageFeatureSet(nm.p2.st, "^MT-", col.name = "percent_mito")
tm.p2.st <- PercentageFeatureSet(tm.p2.st, "^MT-", col.name = "percent_mito")
ni.p3.st <- PercentageFeatureSet(ni.p3.st, "^MT-", col.name = "percent_mito")
ti.p3.st <- PercentageFeatureSet(ti.p3.st, "^MT-", col.name = "percent_mito")
nm.p3.st <- PercentageFeatureSet(nm.p3.st, "^MT-", col.name = "percent_mito")
tm.p3.st <- PercentageFeatureSet(tm.p3.st, "^MT-", col.name = "percent_mito")

ni.p1.st <- PercentageFeatureSet(ni.p1.st, "^HB[^(P)]", col.name = "percent_hb")
ti.p1.st <- PercentageFeatureSet(ti.p1.st, "^HB[^(P)]", col.name = "percent_hb")
nm.p1.st <- PercentageFeatureSet(nm.p1.st, "^HB[^(P)]", col.name = "percent_hb")
tm.p1.st <- PercentageFeatureSet(tm.p1.st, "^HB[^(P)]", col.name = "percent_hb")
ni.p2.st <- PercentageFeatureSet(ni.p2.st, "^HB[^(P)]", col.name = "percent_hb")
ti.p2.st <- PercentageFeatureSet(ti.p2.st, "^HB[^(P)]", col.name = "percent_hb")
nm.p2.st <- PercentageFeatureSet(nm.p2.st, "^HB[^(P)]", col.name = "percent_hb")
tm.p2.st <- PercentageFeatureSet(tm.p2.st, "^HB[^(P)]", col.name = "percent_hb")
ni.p3.st <- PercentageFeatureSet(ni.p3.st, "^HB[^(P)]", col.name = "percent_hb")
ti.p3.st <- PercentageFeatureSet(ti.p3.st, "^HB[^(P)]", col.name = "percent_hb")
nm.p3.st <- PercentageFeatureSet(nm.p3.st, "^HB[^(P)]", col.name = "percent_hb")
tm.p3.st <- PercentageFeatureSet(tm.p3.st, "^HB[^(P)]", col.name = "percent_hb")


ni.p1.st <- subset(ni.p1.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 10 & percent_hb < 1)
ti.p1.st <- subset(ti.p1.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 10 & percent_hb < 1)
nm.p1.st <- subset(nm.p1.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500 & percent_mito < 10 & percent_hb < 1)
tm.p1.st <- subset(tm.p1.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 10 & percent_hb < 1)
ni.p2.st <- subset(ni.p2.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 10 & percent_hb < 1)
ti.p2.st <- subset(ti.p2.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 18 & percent_hb < 1)
nm.p2.st <- subset(nm.p2.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 10 & percent_hb < 1)
tm.p2.st <- subset(tm.p2.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500  & percent_mito < 20 & percent_hb < 1)
ni.p3.st <- subset(ni.p3.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500 & percent_mito < 12 & percent_hb < 1)
ti.p3.st <- subset(ti.p3.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500 & percent_mito < 10 & percent_hb < 1)
nm.p3.st <- subset(nm.p3.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500 & percent_mito < 18 & percent_hb < 1)
tm.p3.st <- subset(tm.p3.st, subset = nFeature_Spatial > 250  & nCount_Spatial > 500 & percent_mito < 10 & percent_hb < 1)



ni.p1.st <- ni.p1.st[!grepl("^MT-", rownames(ni.p1.st)), ]
ni.p1.st <- ni.p1.st[!grepl("^HB[^(P)]", rownames(ni.p1.st)), ]
ni.p1.st <- ni.p1.st[!grepl("^RP[SL]", rownames(ni.p1.st)), ]

ti.p1.st <- ti.p1.st[!grepl("^MT-", rownames(ti.p1.st)), ]
ti.p1.st <- ti.p1.st[!grepl("^HB[^(P)]", rownames(ti.p1.st)), ]
ti.p1.st <- ti.p1.st[!grepl("^RP[SL]", rownames(ti.p1.st)), ]

nm.p1.st <- nm.p1.st[!grepl("^MT-", rownames(nm.p1.st)), ]
nm.p1.st <- nm.p1.st[!grepl("^HB[^(P)]", rownames(nm.p1.st)), ]
nm.p1.st <- nm.p1.st[!grepl("^RP[SL]", rownames(nm.p1.st)), ]

tm.p1.st <- tm.p1.st[!grepl("^MT-", rownames(tm.p1.st)), ]
tm.p1.st <- tm.p1.st[!grepl("^HB[^(P)]", rownames(tm.p1.st)), ]
tm.p1.st <- tm.p1.st[!grepl("^RP[SL]", rownames(tm.p1.st)), ]

ni.p2.st <- ni.p2.st[!grepl("^MT-", rownames(ni.p2.st)), ]
ni.p2.st <- ni.p2.st[!grepl("^HB[^(P)]", rownames(ni.p2.st)), ]
ni.p2.st <- ni.p2.st[!grepl("^RP[SL]", rownames(ni.p2.st)), ]

ti.p2.st <- ti.p2.st[!grepl("^MT-", rownames(ti.p2.st)), ]
ti.p2.st <- ti.p2.st[!grepl("^HB[^(P)]", rownames(ti.p2.st)), ]
ti.p2.st <- ti.p2.st[!grepl("^RP[SL]", rownames(ti.p2.st)), ]

nm.p2.st <- nm.p2.st[!grepl("^MT-", rownames(nm.p2.st)), ]
nm.p2.st <- nm.p2.st[!grepl("^HB[^(P)]", rownames(nm.p2.st)), ]
nm.p2.st <- nm.p2.st[!grepl("^RP[SL]", rownames(nm.p2.st)), ]

tm.p2.st <- tm.p2.st[!grepl("^MT-", rownames(tm.p2.st)), ]
tm.p2.st <- tm.p2.st[!grepl("^HB[^(P)]", rownames(tm.p2.st)), ]
tm.p2.st <- tm.p2.st[!grepl("^RP[SL]", rownames(tm.p2.st)), ]


ni.p3.st <- ni.p3.st[!grepl("^MT-", rownames(ni.p3.st)), ]
ni.p3.st <- ni.p3.st[!grepl("^HB[^(P)]", rownames(ni.p3.st)), ]
ni.p3.st <- ni.p3.st[!grepl("^RP[SL]", rownames(ni.p3.st)), ]

ti.p3.st <- ti.p3.st[!grepl("^MT-", rownames(ti.p3.st)), ]
ti.p3.st <- ti.p3.st[!grepl("^HB[^(P)]", rownames(ti.p3.st)), ]
ti.p3.st <- ti.p3.st[!grepl("^RP[SL]", rownames(ti.p3.st)), ]

nm.p3.st <- nm.p3.st[!grepl("^MT-", rownames(nm.p3.st)), ]
nm.p3.st <- nm.p3.st[!grepl("^HB[^(P)]", rownames(nm.p3.st)), ]
nm.p3.st <- nm.p3.st[!grepl("^RP[SL]", rownames(nm.p3.st)), ]

tm.p3.st <- tm.p3.st[!grepl("^MT-", rownames(tm.p3.st)), ]
tm.p3.st <- tm.p3.st[!grepl("^HB[^(P)]", rownames(tm.p3.st)), ]
tm.p3.st <- tm.p3.st[!grepl("^RP[SL]", rownames(tm.p3.st)), ]


###########################Integrate:remove batch effect#######################
# create a list of the original data that we loaded to start with
st.list = list(ni.p1.st,ti.p1.st,nm.p1.st,tm.p1.st,
               ni.p2.st,ti.p2.st,nm.p2.st,tm.p2.st,
               ni.p3.st,ti.p3.st,nm.p3.st,tm.p3.st)
rm(ni.p1.st,ti.p1.st,nm.p1.st,tm.p1.st,
   ni.p2.st,ti.p2.st,nm.p2.st,tm.p2.st,
   ni.p3.st,ti.p3.st,nm.p3.st,tm.p3.st)
gc()
# sdata.all <- merge(st.list, 
#                  add.cell.ids = c("NI-P1", "TI-P1", "NM-P1", "TM-P1","NI-P2", "TI-P2", "NM-P2", "TM-P2"))

# run SCT on both datasets
st.list = lapply(st.list, SCTransform,
                 assay = "Spatial",
                 variable.features.n=3000,
                 method = "poisson")
#need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 8000 * 1024^2)  # set allowed size to 8K MiB
st.features = SelectIntegrationFeatures(st.list, nfeatures = 2000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
saveRDS(int.anchors,file='../MPSC&ST/variables/int.anchors.stData.Rds')
int.anchors = readRDS(file='../MPSC&ST/variables/int.anchors.stData.biotrainee.Rds')
lung.integrated <- IntegrateData(anchorset = int.anchors,
                                 normalization.method = "SCT",
                                 verbose = FALSE)

names(lung.integrated@images) = c('NIP1','TIP1','NMP1','TMP1',
                                  'NIP2','TIP2','NMP2','TMP2',
                                  'NIP3','TIP3','NMP3','TMP3')

rm(int.anchors, st.list)
gc()


saveRDS(lung.integrated,file='../MPSC&ST/variables/int.stData.Rds')
lung.integrated = readRDS(file='../MPSC&ST/variables/int.stData.Rds')

####Reduction

lung.integrated <- RunPCA(lung.integrated, verbose = FALSE)
DimHeatmap(lung.integrated, dims = 19:27, cells = 500, balanced = TRUE)
DimHeatmap(lung.integrated, dims = 28:36, cells = 500, balanced = TRUE)
DimHeatmap(lung.integrated, dims = 37:45, cells = 500, balanced = TRUE)
DimHeatmap(lung.integrated, dims = 46:50, cells = 500, balanced = TRUE)

lung.integrated <- FindNeighbors(lung.integrated, dims = 1:30)
#lung.integrated <- FindClusters(lung.integrated, verbose = FALSE,resolution = 0.2)
for (res in c( 0.1,0.2,0.3, 0.5,1)) {
  lung.integrated <- FindClusters(lung.integrated, graph.name = "integrated_snn", resolution = res, algorithm = 1)
}
lung.integrated <- RunUMAP(lung.integrated, dims = 1:30)

plot_grid(ncol = 4, 
          DimPlot(lung.integrated, reduction = "umap", group.by = "integrated_snn_res.0.1") + 
            ggtitle("res 0.1"), 
          DimPlot(lung.integrated, reduction = "umap", group.by = "integrated_snn_res.0.2") + 
            ggtitle("res 0.2"), 
          DimPlot(lung.integrated, reduction = "umap", group.by = "integrated_snn_res.0.3") + 
            ggtitle("res 0.3"),
          DimPlot(lung.integrated, reduction = "umap", group.by = "integrated_snn_res.0.5") + 
            ggtitle("res 0.5"))
DimPlot(lung.integrated, reduction = "umap", group.by = "integrated_snn_res.1") + 
  ggtitle("res 1")
DimPlot(lung.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))

SpatialDimPlot(lung.integrated,images = c('NIP1','TIP1'),alpha = c(0.1,0.1))
SpatialDimPlot(lung.integrated,images = c('NIP2','TIP2'),alpha = 1)

sel.clust = "integrated_snn_res.0.1"

lung.integrated <- SetIdent(lung.integrated, value = sel.clust)
table(lung.integrated@active.ident)
SpatialDimPlot(lung.integrated,images = c('NIP1','TIP1'),alpha = c(0.1,0.1))
SpatialPlot(lung.integrated,images = c('NIP2','TIP2'),image.alpha = 0,group.by = 'cellTypeByMarker')
DimPlot(lung.integrated,reduction = 'umap',group.by = 'cellTypeByMarker')
feats <- c("nFeature_Spatial", "nCount_Spatial", "percent_mito")
VlnPlot(lung.integrated, features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()

VlnPlot(lung.integrated, features = feats, pt.size = 0, ncol = 1,group.by = 'orig.ident') + 
  NoLegend()
# differential expression between clusters
de_markers.3 <- FindAllMarkers(lung.integrated, assay = 'integrated',logfc.threshold = 0.1, 
                               test.use = "wilcox", 
                               min.pct = 0.01, 
                               min.diff.pct = 0.2, 
                               only.pos = TRUE)
de_markers_list = split(de_markers.3,de_markers.3$cluster)


##############define spatial clusters based on scRNA data#############################
markers_genes.sc = read.csv(file = '../MPSC&ST/variables/marker_genes_cca_resP12.csv')

top25 <- markers_genes.sc %>% group_by(cluster) %>% top_n(-25, p_val_adj)%>% top_n(25, avg_logFC) 
top25_list = lapply(split(top25,top25$cluster), function(x){
  return(as.character(x$gene))
})
names(top25_list)=paste0('scRNA_C',names(top25_list))
library(fgsea)
# run fgsea for each of the clusters in the list
res <- lapply(de_markers_list, function(x) {
  if(nrow(x)>1){
    gene_rank <- setNames(x$avg_logFC, x$gene)
    fgseaRes <- fgsea(pathways = top25_list, stats = gene_rank,nperm = 1000)
    
  }else{
    fgseaRes = data.frame()
  }
  return(fgseaRes)
})


names(res) <- names(de_markers_list)


#####################################Try to integrate with scRNA data###########################################

scdata = readRDS(file = '../MPSC&ST/variables/data.filt_harmony_annotated-v2.Rds')

# select 500 cells per subclass, first set subclass ass active.ident
Idents(scdata) <- scdata$CellTypeManully
scdata <- subset(scdata, cells = WhichCells(scdata, downsample = 1000))
# check again number of cells per subclass
table(scdata$CellTypeManully)
# First run SCTransform and PCA
scdata <- SCTransform(scdata) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(scdata, label = TRUE)
DimPlot(scdata, group.by = 'orig.ident')

###########Integrate with scRNA
anchors <- FindTransferAnchors(reference = scdata, query = lung.integrated, normalization.method = "SCT",
                               reference.assay = 'RNA',query.assay = 'integrated', reduction = "cca")
predictions.assay <- TransferData(anchorset = anchors, refdata = scdata$CellTypeManully, 
                                  prediction.assay = TRUE, weight.reduction = lung.integrated[["pca"]])
lung.integrated[["predictions"]] <- predictions.assay
dim(GetAssayData(lung.integrated, assay = "predictions"))

lung.integrated <- AddMetaData(lung.integrated, metadata = as.data.frame(t(GetAssayData(lung.integrated, assay = "predictions"))))
head(lung.integrated@meta.data)

prediction.scores <- as.data.frame(t(GetAssayData(lung.integrated, assay = "predictions")))
prediction.scores$max <- NULL
sum(is.na(prediction.scores))
prediction.scores$celltype_prediction <- NA
dim(prediction.scores)
for(i in 1:nrow(prediction.scores)){
  prediction.scores$celltype_prediction[i] <- colnames(prediction.scores)[prediction.scores[i,1:7] == max(as.double(prediction.scores[i,1:7]))]
}

table(prediction.scores$celltype_prediction)
lung.integrated$celltype_prediction <- prediction.scores$celltype_prediction

DimPlot(lung.integrated,group.by ='celltype_prediction')


DefaultAssay(lung.integrated) <- "predictions"
SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('Epithelial cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('Epithelial cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('Fibroblast cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('Fibroblast cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)

SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('Myeloid cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('Myeloid cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)

SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('T&NK cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('T&NK cell'), pt.size.factor = 1.6, ncol = 4, alpha = c(0.1,1),
                   crop = TRUE)

de_markers.pred <- FindAllMarkers(lung.integrated, assay = 'predictions',logfc.threshold = 0.1, 
                                  test.use = "wilcox", 
                                  min.pct = 0.01, 
                                  min.diff.pct = 0.2, 
                                  only.pos = TRUE)
VlnPlot(lung.integrated,features = rownames(lung.integrated@assays$predictions@data),
        assay = 'predictions',pt.size = 0)

Idents(lung.integrated)=lung.integrated$celltype_prediction
DimPlot(object = lung.integrated, label = TRUE) #
DimPlot(object = lung.integrated, group.by = 'orig.ident') #
# library(pals)
# color_pelette <- rev(as.vector(kelly()[4:(length(unique(lung.integrated$celltype_prediction))+3)]))
# names(color_pelette) <- unique(lung.integrated$celltype_prediction)
cellType.cols = pal_npg('nrc')(length(unique(lung.integrated$celltype_prediction)))
names(cellType.cols)=c('Epithelial cell','Endothelial cell','Fibroblast cell',
                       'T&NK cell','B cell','Myeloid cell','Mast cell')

DimPlot(object = lung.integrated,cols = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/Dimplot cellType predictions.pdf
SpatialDimPlot(lung.integrated,images = c('NIP1','TIP1'),cols = cellType.cols)
SpatialDimPlot(lung.integrated,images = c('NMP1','TMP1'),cols = cellType.cols)
SpatialDimPlot(lung.integrated,images = c('NIP2','TIP2'),cols = cellType.cols)
SpatialDimPlot(lung.integrated,images = c('NMP2','TMP2'),cols = cellType.cols)
SpatialDimPlot(lung.integrated,images = c('NIP3','TIP3'),cols = cellType.cols)
SpatialDimPlot(lung.integrated,images = c('NMP3','TMP3'),cols = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/NMP3 TMP3 spatialDimPlot.pdf
lung.integrated$CellTypeManully = Idents(lung.integrated)
com2orig = table(lung.integrated$orig.ident,lung.integrated$CellTypeManully)
#(com2orig)=paste('SNN',colnames(com2orig))
library(reshape2)
com2orig = melt(com2orig,measure.vars = colnames(com2orig))
colnames(com2orig) = c('Orig.ident','Cell.Type','No.')
library(RColorBrewer)
ggplot(com2orig,mapping = aes(x=Cell.Type,y=No.,fill=Orig.ident))+
  geom_bar(stat = 'identity',width = 0.6)+coord_flip()+
  theme_cowplot()+
  scale_fill_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(9,'Set3')))

plot.data = data.frame(Orig.ident = names(table(lung.integrated$orig.ident)),Number = as.integer(table(lung.integrated$orig.ident)))
plot.data$Patient =unlist(lapply(as.character(plot.data$Orig.ident),function(a){
  return(strsplit(a,'_')[[1]][2])
}))
ggplot(plot.data,aes(x=Orig.ident,y=Number,fill=Patient))+geom_bar(stat = 'identity',width = 0.65)+
  facet_wrap(vars(Patient),scales = 'free_x')+theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_npg()

annotation = lung.integrated@meta.data
annotation$Patient = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][2])
}))
annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = c('NI','TI','NM','TM'))
annotation$type = substr(annotation$From,1,1)

library(ggplot2)
ggplot(data = annotation,aes(x=Patient,fill=type))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+theme( strip.background = element_blank())+
  scale_fill_manual(values = brewer.pal(3,'Paired')[-1])
ggplot(data = annotation,aes(x=Patient,fill=From))+geom_bar(position = 'fill',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+
  scale_fill_manual(values = c(brewer.pal(4,'Set2')))+theme( strip.background = element_blank())
ggplot(data = annotation,aes(x=Patient))+geom_bar(fill='blue',width = 0.7)+coord_flip()+
  facet_grid( rows= vars(CellTypeManully))+theme_bw()+theme( strip.background = element_blank())

saveRDS(lung.integrated,file = '../MPSC&ST/variables/spatialDataFilt&Integrate&Annotated_P12_V2.Rds')


####Integrates with scRNA  subtype information####
scdata = readRDS(file = '../MPSC&ST/variables/data.filt_harmony_annotated-With-subtypes.Rds')

# select 500 cells per subclass, first set subclass ass active.ident
Idents(scdata) <- scdata$subType
scdata <- subset(scdata, cells = WhichCells(scdata, downsample = 500))
# check again number of cells per subclass
cells.n = table(scdata$subType)
cells.type = names(cells.n)[cells.n >= 300]
scdata <- scdata[,scdata$subType %in% cells.type]
table(scdata$subType)
# First run SCTransform and PCA
scdata <- SCTransform(scdata) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
save(scdata,file = '../MPSC&ST/variables/scdataWithSubtypes.RData')
load(file = '../MPSC&ST/variables/scdataWithSubtypes.RData')
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(scdata, label = TRUE)
DimPlot(scdata, group.by = 'orig.ident')

###########Integrate
DefaultAssay(lung.integrated)='integrated'
DefaultAssay(scdata)='SCT'
lung.integrated = UpdateSeuratObject(lung.integrated)
lung.integrated = UpdateSCTAssays(lung.integrated)
scdata= UpdateSeuratObject(scdata)
scdata = UpdateSCTAssays(scdata)

#anchors2 <- FindTransferAnchors(reference = scdata, query = lung.integrated, normalization.method = "SCT",
 #                              reference.assay = 'RNA',query.assay = 'integrated', reduction = "cca")
anchors2 <- FindTransferAnchors(reference = scdata, query = lung.integrated, normalization.method = "SCT",
                               reference.assay = 'SCT',query.assay = 'integrated',recompute.residuals = F, reduction = "cca")

predictions.assay <- TransferData(anchorset = anchors2, refdata = scdata$subType, dims = 1:30,
                                  prediction.assay = TRUE, weight.reduction = lung.integrated[["pca"]])

lung.integrated[["predictions2"]] <- predictions.assay
dim(GetAssayData(lung.integrated, assay = "predictions2"))

lung.integrated <- AddMetaData(lung.integrated, metadata = as.data.frame(t(GetAssayData(lung.integrated, assay = "predictions2"))))
head(lung.integrated@meta.data)

prediction.scores <- as.data.frame(t(GetAssayData(lung.integrated, assay = "predictions2")))
prediction.scores$max <- NULL
sum(is.na(prediction.scores))
prediction.scores$celltype_prediction2 <- NA
dim(prediction.scores)
for(i in 1:nrow(prediction.scores)){
  prediction.scores$celltype_prediction2[i] <- colnames(prediction.scores)[prediction.scores[i,1:44] == max(as.double(prediction.scores[i,1:44]))]
}
write.csv(prediction.scores,file='../MPSC&ST/variables/subtypePredictionScores.csv')
table(prediction.scores$celltype_prediction2)
lung.integrated$celltype_prediction2 <- prediction.scores$celltype_prediction2

DimPlot(lung.integrated,group.by ='celltype_prediction2')
DefaultAssay(lung.integrated)='predictions2'
SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('CXCL14+'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('CXCL14+'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP3','TIP3','NMP3','TMP3'), features = c('CXCL14+'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)

SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1','NMP1','TMP1'), features = c('AT2'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2','NMP2','TMP2'), features = c('AT2'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)
SpatialFeaturePlot(lung.integrated,images = c('NIP3','TIP3','NMP3','TMP3'), features = c('AT2'), pt.size.factor = 2,
                   ncol = 4, alpha = c(0.1,1),crop = TRUE)


#####################   This part calclulates the celltype pair neighborhood maps,2021-09-23  ############################
# Example shown for D10 (Run this 12 times for individual sample
lung.integrated = readRDS(file='../MPSC&ST/variables/spatialDataFilt&Integrate&AnnotatedV2.Rds')
prediction.scores <- as.data.frame(t(GetAssayData(lung.integrated, assay = "predictions")))
prediction.scores$max <- NULL
dim(prediction.scores)
prediction.scores.1 <- prediction.scores[colnames(lung.integrated)[lung.integrated$orig.ident == "TM_P3"],]#
prediction.cellTypes.1 = lung.integrated$celltype_prediction[colnames(lung.integrated)[lung.integrated$orig.ident == "TM_P3"]]#
dim(prediction.scores.1)
prediction.scores.1$cellType = prediction.cellTypes.1
spots.info = read.csv(file='H:/project/single cell/ourData/stRNA/Count/P8_M_T/spatial/tissue_positions_list.csv',header = F,row.names = 1)#
colnames(spots.info)=c('in_tissue','array_row','array_col','pxl_col','pxl_row')
spots.info$pairxy=paste(spots.info$array_row,spots.info$array_col)

rownames(prediction.scores.1) = unlist(lapply(rownames(prediction.scores.1), function(a){
  strsplit(a,'_')[[1]][1]
}))


interaction_matrix = matrix(0, ncol = length(unique(lung.integrated$celltype_prediction)), nrow = length(unique(lung.integrated$celltype_prediction)))
rownames(interaction_matrix) <- unique(lung.integrated$celltype_prediction)
colnames(interaction_matrix) <- unique(lung.integrated$celltype_prediction)

for(i in 1:(nrow(prediction.scores.1)-1)){
  set.i = rownames(prediction.scores.1)[i]
  set.i.type = prediction.scores.1[i,'cellType']
  coord.i = spots.info[set.i,]
  row.i = as.integer(coord.i['array_row'])
  col.i = as.integer(coord.i['array_col'])
  candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                           colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
  candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
  spots.info.c.i = spots.info[spots.info$pairxy %in% as.character(candi.coord.i$pairxy),]
  
  retain.j.set = rownames(prediction.scores.1)[(i+1):nrow(prediction.scores.1)]
  retain.j.set.neigh = retain.j.set[retain.j.set %in% rownames(spots.info.c.i)]
  neigh.types = prediction.scores.1[retain.j.set.neigh,'cellType']
  for(neigh.type in neigh.types){
    interaction_matrix[set.i.type,neigh.type]=interaction_matrix[set.i.type,neigh.type]+1
  }
}

interaction_matrix.2 = t(interaction_matrix)
diag(interaction_matrix.2)=0
interaction_matrix <- interaction_matrix + interaction_matrix.2
pheatmap::pheatmap(interaction_matrix)
####Cell heterogeneous scores, 2019-09-26####

hete.scores = numeric(length = nrow(prediction.scores))
for(i in 1:nrow(prediction.scores)){
  set.i = rownames(prediction.scores)[i]
  set.i.pred = prediction.scores[i,c(1:7)]
  hete.scores[i]=length(set.i.pred[set.i.pred>0.1])
}
saveRDS(hete.scores,file='./variables/hete.scores.spatial.Rds')
hete.scores = readRDS(file='../MPSC&ST/variables/hete.scores.spatial.Rds')
lung.integrated$hete.scores = hete.scores
SpatialFeaturePlot(lung.integrated,
                   images = c('NIP1','TIP1','NMP1','TMP1'), 
                   features = c('hete.scores'), 
                   pt.size.factor = 1.6, ncol = 4, 
                   alpha = c(0.5,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,
                   images = c('NIP2','TIP2','NMP2','TMP2'), 
                   features = c('hete.scores'), 
                   pt.size.factor = 1.6, ncol = 4, 
                   alpha = c(0.5,1),
                   crop = TRUE)
SpatialFeaturePlot(lung.integrated,
                   images = c('NIP3','TIP3','NMP3','TMP3'), 
                   features = c('hete.scores'), 
                   pt.size.factor = 1.6, ncol = 4, 
                   alpha = c(0.5,1),
                   crop = TRUE)
lung.integrated$Patient = unlist(lapply(lung.integrated$orig.ident,function(a){
  strsplit(a,'_')[[1]][2]
}))
lung.integrated$Sample = unlist(lapply(lung.integrated$orig.ident,function(a){
  strsplit(a,'_')[[1]][1]
}))
lung.integrated$celltype_prediction = factor(lung.integrated$celltype_prediction,
                                             levels = names(cellType.cols))
ggplot(lung.integrated@meta.data,aes(x=Sample,y=hete.scores,fill = celltype_prediction))+
  facet_grid(Patient~celltype_prediction)+
  geom_boxplot(outlier.size = 1)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Heterogeneity score')+
  scale_fill_manual(values = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/barplot heterogenety scores P123.pdf
####Cell interaction scores, 2021-09-15####
##consis.scores: how many proportions of neighbors share the same cell type with the investigated one
##diff scores: how many proportions of neighbors have different cell types with the investigated one

consis.scores = numeric(length = nrow(prediction.scores.1))
diff.scores = numeric(length = nrow(prediction.scores.1))
for(i in 1:nrow(prediction.scores.1)){
  set.i = rownames(prediction.scores.1)[i]
  set.i.type = prediction.scores.1[i,'cellType']
  coord.i = spots.info[set.i,]
  row.i = as.integer(coord.i['array_row'])
  col.i = as.integer(coord.i['array_col'])
  candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                           colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
  candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
  spots.info.c.i = spots.info[spots.info$pairxy %in% as.character(candi.coord.i$pairxy),]
  retain.j.set.neigh =  rownames(prediction.scores.1)[ rownames(prediction.scores.1) %in% rownames(spots.info.c.i)]
  neigh.types = prediction.scores.1[retain.j.set.neigh,'cellType']
  consis.scores[i]=length(neigh.types[neigh.types == set.i.type])/length(neigh.types)
  diff.scores[i]=length(unique(neigh.types[neigh.types != set.i.type]))/length(neigh.types)
  
}
res.sum = data.frame(row.names = rownames(prediction.scores.1),
                     cellType = prediction.scores.1$cellType,
                     consis.scores,
                     diff.scores)
write.csv(res.sum,file = "../MPSC&ST/variables/spatial/consis&diffScores_TMP3.csv")
p1= ggplot(res.sum,aes(x=cellType,y=consis.scores,fill = cellType))+
  geom_boxplot()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Consistent score')+
  scale_fill_manual(values = brewer.pal(8,'Set2'))
p2= ggplot(res.sum,aes(x=cellType,y=diff.scores,fill = cellType))+
  geom_boxplot()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Diff score')+
  scale_fill_manual(values = brewer.pal(8,'Set2'))
p1/p2



library(pals)
color_pelette <- rev(as.vector(kelly()[2:(length(unique(lung.integrated$celltype_prediction))+1)]))
names(color_pelette) <- unique(lung.integrated$celltype_prediction)

# interaction_matrix <- read.csv("interactions-D10.csv", row.names = 1)
interaction_matrix[lower.tri(interaction_matrix)] <- 0
write.csv(interaction_matrix, file = "../MPSC&ST/variables/spatial/interactions-TMP3.csv")#

library(circlize)
color_used <- color_pelette[colnames(interaction_matrix)]
row_col <- color_used
#row_col[names(row_col) != "TMSB4X high cells"] <- "#cecece"

col <- matrix(rep(color_used, each = ncol(interaction_matrix), T), nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix))
rownames(col) <- rownames(interaction_matrix)
colnames(col) <- colnames(interaction_matrix)
chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = c("name", "grid"),)
circos.clear()
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/chordDiagram  TMP1.pdf

####Consis and diff score compare among samples####
res.sum.1 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TMP1.csv",row.names = 1)
res.sum.2= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TIP1.csv",row.names = 1)
res.sum.3 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NMP1.csv",row.names = 1)
res.sum.4= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NIP1.csv",row.names = 1)

res.sum.1$From = 'TMP'
res.sum.2$From = 'TIP'
res.sum.3$From = 'NMP'
res.sum.4$From = 'NIP'

res.sum1 = rbind(res.sum.1,res.sum.2,res.sum.3,res.sum.4)
res.sum1$cellType = factor(res.sum1$cellType)
res.sum1$Patient = 'P1'


res.sum.1 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TMP2.csv",row.names = 1)
res.sum.2= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TIP2.csv",row.names = 1)
res.sum.3 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NMP2.csv",row.names = 1)
res.sum.4= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NIP2.csv",row.names = 1)

res.sum.1$From = 'TMP'
res.sum.2$From = 'TIP'
res.sum.3$From = 'NMP'
res.sum.4$From = 'NIP'

res.sum2 = rbind(res.sum.1,res.sum.2,res.sum.3,res.sum.4)
res.sum2$cellType = factor(res.sum2$cellType)
res.sum2$Patient = 'P2'


res.sum.1 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TMP3.csv",row.names = 1)
res.sum.2= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_TIP3.csv",row.names = 1)
res.sum.3 = read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NMP3.csv",row.names = 1)
res.sum.4= read.csv(file = "../MPSC&ST/variables/spatial/consis&diffScores_NIP3.csv",row.names = 1)

res.sum.1$From = 'TMP'
res.sum.2$From = 'TIP'
res.sum.3$From = 'NMP'
res.sum.4$From = 'NIP'

res.sum3 = rbind(res.sum.1,res.sum.2,res.sum.3,res.sum.4)
res.sum3$cellType = factor(res.sum3$cellType)
res.sum3$Patient = 'P3'

res.sum = rbind(res.sum1,res.sum2,res.sum3)

test.res = res.sum %>% group_by(Patient,cellType) %>% summarize(p.value = kruskal.test(consis.scores ~From)$p.value)
test.res$p.value = as.character(signif(test.res$p.value,2))
res.sum$cellType = factor(res.sum$cellType,levels = names(cellType.cols))
res.sum$Type = substr(res.sum$From,1,1)

ggplot(res.sum,aes(x=From,y=consis.scores,fill = cellType))+
  facet_grid(Patient~cellType)+
  geom_boxplot()+
  geom_text(data = test.res,mapping = aes(x=2.5,y=1.1,label = p.value),size=3)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Consistent score')+
  scale_fill_manual(values = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/consis&diff scores boxplot all P123.pdf
test.res2 = res.sum %>% group_by(Patient,cellType) %>% summarize(p.value = kruskal.test(consis.scores ~ Type)$p.value)
test.res2$p.value = as.character(signif(test.res2$p.value,2))

ggplot(res.sum,aes(x=Type,y=consis.scores,fill = cellType))+
  facet_grid(Patient~cellType)+
  geom_boxplot(outlier.size = 1)+
  geom_text(data = test.res2,mapping = aes(x=1.5,y=1.1,label = p.value),size=3)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Consistent score')+
  scale_fill_manual(values = cellType.cols)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/consis&diff scores boxplot all P123_N compare T.pdf
res.sum.T = res.sum[res.sum$Type == 'T',]
res.sum.T$location = substr(res.sum.T$From,1,2)
test.res3 = res.sum.T %>% group_by(Patient,cellType) %>% summarize(p.value = kruskal.test(consis.scores ~ location)$p.value)
test.res3$p.value = as.character(signif(test.res3$p.value,2))

ggplot(res.sum.T,aes(x=location,y=consis.scores,fill = cellType))+
  facet_grid(Patient~cellType)+
  geom_boxplot(outlier.size = 1,)+
  geom_text(data = test.res3,mapping = aes(x=1.5,y=1.1,label = p.value),size=3)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Consistent score')+
  scale_fill_manual(values = cellType.cols)
####Check B cell markers
DefaultAssay(lung.integrated)='SCT'
FeaturePlot(lung.integrated, reduction = "umap", 
            features = c("CD79A","IGHM","IGHG3","CD19",'IGHA2'), 
           order = T,  combine = T)

####PCA analysis




####Interactin map compare among samples####

inter.mat.tip1 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TIP1.csv",row.names = 1)
inter.mat.tmp1 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TMP1.csv",row.names = 1)
inter.mat.nip1 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NIP1.csv",row.names = 1)
inter.mat.nmp1 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NMP1.csv",row.names = 1)
inter.mat.tip2 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TIP2.csv",row.names = 1)
inter.mat.tmp2 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TMP2.csv",row.names = 1)
inter.mat.nip2 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NIP2.csv",row.names = 1)
inter.mat.nmp2 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NMP2.csv",row.names = 1)
inter.mat.tip3 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TIP3.csv",row.names = 1)
inter.mat.tmp3 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-TMP3.csv",row.names = 1)
inter.mat.nip3 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NIP3.csv",row.names = 1)
inter.mat.nmp3 = read.csv(file = "../MPSC&ST/variables/spatial/interactions-NMP3.csv",row.names = 1)

library(reshape2)

changeInterMat = function(inter.mat.tip1){
  inter.mat.tip1 = inter.mat.tip1/sum(inter.mat.tip1)
  inter.mat.tip1[lower.tri(inter.mat.tip1)] = NA
  inter.mat.tip1 = data.frame(inter.mat.tip1)
  inter.mat.tip1$source=rownames(inter.mat.tip1)
  inter.mat.tip1 = melt(inter.mat.tip1,
                        measure.vars = c("Epithelial.cell" , "Endothelial.cell" ,
                                         "Myeloid.cell"  ,   "Fibroblast.cell" ,
                                         "B.cell","Mast.cell","T.NK.cell" ))
  inter.mat.tip1 = inter.mat.tip1[!is.na(inter.mat.tip1$value),]
  rownames(inter.mat.tip1) = paste(inter.mat.tip1$source,inter.mat.tip1$variable)
  return(inter.mat.tip1)
}
inter.mat.tip1 = changeInterMat(inter.mat.tip1)
inter.mat.tmp1 = changeInterMat(inter.mat.tmp1)
inter.mat.nip1 = changeInterMat(inter.mat.nip1)
inter.mat.nmp1 = changeInterMat(inter.mat.nmp1)

inter.mat.tip2 = changeInterMat(inter.mat.tip2)
inter.mat.tmp2 = changeInterMat(inter.mat.tmp2)
inter.mat.nip2 = changeInterMat(inter.mat.nip2)
inter.mat.nmp2 = changeInterMat(inter.mat.nmp2)

inter.mat.tip3 = changeInterMat(inter.mat.tip3)
inter.mat.tmp3 = changeInterMat(inter.mat.tmp3)
inter.mat.nip3 = changeInterMat(inter.mat.nip3)
inter.mat.nmp3 = changeInterMat(inter.mat.nmp3)

all.matrix = matrix(nrow=nrow(inter.mat.tip1),ncol=12,
                    dimnames = list(rownames(inter.mat.tip1),
                                    c('TIP1','TMP1','NIP1','NMP1',
                                      'TIP2','TMP2','NIP2','NMP2',
                                      'TIP3','TMP3','NIP3','NMP3')))
all.matrix[,1]=inter.mat.tip1[rownames(all.matrix),'value']
all.matrix[,2]=inter.mat.tmp1[rownames(all.matrix),'value']
all.matrix[,3]=inter.mat.nip1[rownames(all.matrix),'value']
all.matrix[,4]=inter.mat.nmp1[rownames(all.matrix),'value']
all.matrix[,5]=inter.mat.tip2[rownames(all.matrix),'value']
all.matrix[,6]=inter.mat.tmp2[rownames(all.matrix),'value']
all.matrix[,7]=inter.mat.nip2[rownames(all.matrix),'value']
all.matrix[,8]=inter.mat.nmp2[rownames(all.matrix),'value']
all.matrix[,9]=inter.mat.tip3[rownames(all.matrix),'value']
all.matrix[,10]=inter.mat.tmp3[rownames(all.matrix),'value']
all.matrix[,11]=inter.mat.nip3[rownames(all.matrix),'value']
all.matrix[,12]=inter.mat.nmp3[rownames(all.matrix),'value']

pheatmap::pheatmap(all.matrix[])
library(ComplexHeatmap)
ComplexHeatmap::Heatmap(all.matrix,
                        color = colorRamp2(breaks = c(-0.1,0,0.1),colors = c('blue','white','red')),
                        rect_gp = gpar(col = 'black'),
                        cluster_rows = F,
                        clustering_distance_columns = 'spearman',
                        column_names_side  = 'top')

all.df = melt(all.matrix,measure.vars = colnames(all.matrix))
all.df$Type = substr(all.df$Var2,1,1)
all.df$Location= substr(all.df$Var2,1,2)
all.df$Patient= substr(all.df$Var2,3,4)


test.res4 = all.df %>% group_by(Var1) %>% summarize(p.value = kruskal.test(value ~ Type)$p.value)
test.res4$p.value = as.character(signif(test.res4$p.value,2))
test.res4 = as.data.frame(test.res4)
sig.Var1 = as.character(test.res4[test.res4$p.value<0.05,'Var1'])
ggplot(all.df[all.df$Var1%in% sig.Var1,],aes(x=Type,y=value,fill = Type))+
  facet_grid(~Var1)+
  geom_boxplot(outlier.size = 1)+
  geom_text(data = test.res4[test.res4$Var1 %in% sig.Var1,],mapping = aes(x=1.5,y=0.08,label = p.value),size=3,inherit.aes = F)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8,angle = 90),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='Proportion')+
  scale_fill_manual(values= c('blue','red'))


diff.ratio = data.frame(Diff = c(as.double(abs(all.matrix[,'TIP1']-all.matrix[,'TMP1'])),
                               as.double(abs(all.matrix[,'TIP2']-all.matrix[,'TMP2'])),
                               as.double(abs(all.matrix[,'TIP3']-all.matrix[,'TMP3']))
                               ),
                      CellType = rep(rownames(all.matrix),3),
                      Patient = rep(c('P1','P2','P3'),each = nrow(all.matrix)))
diff.ratio.sd = diff.ratio %>% group_by(CellType) %>% summarize(mean = mean(Diff))
diff.ratio.sd = as.data.frame(diff.ratio.sd)
diff.ratio.sd = diff.ratio.sd[order(diff.ratio.sd$mean),]
diff.ratio$CellType = factor(diff.ratio$CellType,levels = rev(as.character(diff.ratio.sd$CellType)))
library(ggpubr)  
  
ggbarplot(diff.ratio,x='CellType',y = 'Diff', add = 'mean_sd',fill='blue') + 
  theme_bw() + 
  ylab("Diff Ratio") + theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/TI TM diff ratio cell type interactions barplot.pdf

#######################Plot features#############

lung.integrated = readRDS(file = '../MPSC&ST/variables/spatialDataFilt&Integrate&AnnotatedV2.Rds')


DefaultAssay(lung.integrated) <- "Spatial"
SpatialFeaturePlot(lung.integrated,images = c('NIP1','TIP1'), 
                   features = c("PTPRC","CD3D","CD3E","CD3G"), pt.size.factor = 1, alpha = c(0.1, 1))
SpatialFeaturePlot(ni.p1.st,
                   features = c("DCN","THY1","COL1A1","COL1A2"), pt.size.factor = 1, alpha = c(0.1, 1))
SpatialFeaturePlot(ti.p1.st,
                   features = c("PTPRC","CD3D","CD3E","CD3G"), pt.size.factor = 1, alpha = c(0.1, 1))

SpatialFeaturePlot(lung.integrated,images = c('TMP1','TIP1'), 
                   features = c("CD82","RGS1"), pt.size.factor = 2)
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2'), 
                   features = c("WIF1",'TNFRSF18'), pt.size.factor = 2, alpha = c(0.5, 1))

#######################Recognize spatially variable features############################################
lung.integrated = readRDS(file='../MPSC&ST/variables/spatialDataFilt&Integrate&AnnotatedV2.Rds')

TIP1 = lung.integrated[,lung.integrated$orig.ident == 'TI_P1']
TIP1@images$NIP1=NULL
TIP1@images$NMP1=NULL
TIP1@images$TMP1=NULL
TIP1@images$NIP2=NULL
TIP1@images$TIP2=NULL
TIP1@images$NMP2=NULL
TIP1@images$TMP2=NULL
TIP1@images$NIP3=NULL
TIP1@images$TIP3=NULL
TIP1@images$NMP3=NULL
TIP1@images$TMP3=NULL
TIP1 <- FindSpatiallyVariableFeatures(TIP1, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(TIP1), r.metric = 5, slot = "data")
TIP1 <- FindSpatiallyVariableFeatures(TIP1, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TIP1@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TIP1,file='spatilTIP1.Rds')
top.clusters <- SpatiallyVariableFeatures(TIP1)[SpatiallyVariableFeatures(TIP1) != 'max']
SpatialPlot(object = TIP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

top.clusters <- head(SpatiallyVariableFeatures(TIP1, assay = "SCT"),8)
DefaultAssay(TIP1)='SCT'
SpatialFeaturePlot(object = TIP1, features = top.clusters,ncol=4, alpha = c(0.1, 1))


NIP1 = lung.integrated[,lung.integrated$orig.ident == 'NI_P1']
NIP1@images$TIP1=NULL
NIP1@images$NMP1=NULL
NIP1@images$TMP1=NULL
NIP1@images$NIP2=NULL
NIP1@images$TIP2=NULL
NIP1@images$NMP2=NULL
NIP1@images$TMP2=NULL
NIP1@images$NIP3=NULL
NIP1@images$TIP3=NULL
NIP1@images$NMP3=NULL
NIP1@images$TMP3=NULL
NIP1 <- FindSpatiallyVariableFeatures(NIP1, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(NIP1), r.metric = 5, slot = "data")
NIP1 <- FindSpatiallyVariableFeatures(NIP1, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NIP1@assays$SCT@scale.data), r.metric = 5, slot = "data")

saveRDS(NIP1,file='spatilNIP1.Rds')


top.clusters <- SpatiallyVariableFeatures(NIP1)[SpatiallyVariableFeatures(NIP1) != 'max']
SpatialPlot(object = NIP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

top.clusters <- head(SpatiallyVariableFeatures(NIP1, assay = "SCT"),8)
DefaultAssay(NIP1)='SCT'
SpatialFeaturePlot(object = NIP1, features = top.clusters,ncol=4, alpha = c(0.1, 1))

TMP1 = lung.integrated[,lung.integrated$orig.ident == 'TM_P1']
TMP1@images$NIP1=NULL
TMP1@images$TIP1=NULL
TMP1@images$NMP1=NULL
TMP1@images$NIP2=NULL
TMP1@images$TIP2=NULL
TMP1@images$NMP2=NULL
TMP1@images$TMP2=NULL
TMP1@images$NIP3=NULL
TMP1@images$TIP3=NULL
TMP1@images$NMP3=NULL
TMP1@images$TMP3=NULL

# TMP1 <- FindSpatiallyVariableFeatures(TMP1, assay = "predictions", selection.method = "markvariogram", 
#                                       features = rownames(TMP1), r.metric = 5, slot = "data")
# top.clusters <- SpatiallyVariableFeatures(TMP1)[SpatiallyVariableFeatures(TMP1) != 'max']
# SpatialPlot(object = TMP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))
TMP1 <- FindSpatiallyVariableFeatures(TMP1, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TMP1@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TMP1,file='spatialTMP1.Rds')

top.clusters <- head(SpatiallyVariableFeatures(TMP1, assay = "SCT"),8)
DefaultAssay(TMP1)='SCT'
SpatialFeaturePlot(object = TMP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


NMP1 = lung.integrated[,lung.integrated$orig.ident == 'NM_P1']
NMP1@images$NIP1=NULL
NMP1@images$TIP1=NULL
NMP1@images$TMP1=NULL
NMP1@images$NIP2=NULL
NMP1@images$TIP2=NULL
NMP1@images$NMP2=NULL
NMP1@images$TMP2=NULL
NMP1@images$NIP3=NULL
NMP1@images$TIP3=NULL
NMP1@images$NMP3=NULL
NMP1@images$TMP3=NULL

NMP1 <- FindSpatiallyVariableFeatures(NMP1, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(NMP1), r.metric = 5, slot = "data")
top.clusters <- SpatiallyVariableFeatures(NMP1)[SpatiallyVariableFeatures(NMP1) != 'max']
SpatialPlot(object = NMP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


NMP1 <- FindSpatiallyVariableFeatures(NMP1, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NMP1@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(NMP1,file='spatialNMP1.Rds')

top.clusters <- head(SpatiallyVariableFeatures(NMP1, assay = "SCT"),8)
DefaultAssay(NMP1)='SCT'
SpatialFeaturePlot(object = NMP1, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


TIP2 = lung.integrated[,lung.integrated$orig.ident == 'TI_P2']
TIP2@images$NIP1=NULL
TIP2@images$TIP1=NULL
TIP2@images$NMP1=NULL
TIP2@images$TMP1=NULL
TIP2@images$NIP2=NULL
TIP2@images$NMP2=NULL
TIP2@images$TMP2=NULL
TIP2@images$NIP3=NULL
TIP2@images$TIP3=NULL
TIP2@images$NMP3=NULL
TIP2@images$TMP3=NULL

TIP2 <- FindSpatiallyVariableFeatures(TIP2, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(TIP2), r.metric = 5, slot = "data")
top.clusters <- SpatiallyVariableFeatures(TIP2)[SpatiallyVariableFeatures(TIP2) != 'max']
SpatialPlot(object = TIP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

TIP2 <- FindSpatiallyVariableFeatures(TIP2, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TIP2@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TIP2,file='spatialTIP2.Rds')
top.clusters <- head(SpatiallyVariableFeatures(TIP2, assay = "SCT"),8)
DefaultAssay(TIP2)='SCT'
SpatialFeaturePlot(object = TIP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


NIP2 = lung.integrated[,lung.integrated$orig.ident == 'NI_P2']
NIP2@images$NIP1=NULL
NIP2@images$TIP1=NULL
NIP2@images$NMP1=NULL
NIP2@images$TMP1=NULL
NIP2@images$TIP2=NULL
NIP2@images$NMP2=NULL
NIP2@images$TMP2=NULL
NIP2@images$NIP3=NULL
NIP2@images$TIP3=NULL
NIP2@images$NMP3=NULL
NIP2@images$TMP3=NULL
NIP2 <- FindSpatiallyVariableFeatures(NIP2, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(NIP2), r.metric = 5, slot = "data")
top.clusters <- SpatiallyVariableFeatures(NIP2)[SpatiallyVariableFeatures(NIP2) != 'max']
SpatialPlot(object = NIP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

NIP2 <- FindSpatiallyVariableFeatures(NIP2, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NIP2@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(NIP2,file='spatialNIP2.Rds')
top.clusters <- head(SpatiallyVariableFeatures(NIP2, assay = "SCT"),8)
DefaultAssay(NIP2)='SCT'
SpatialFeaturePlot(object = NIP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

TMP2 = lung.integrated[,lung.integrated$orig.ident == 'TM_P2']
TMP2@images$NIP1=NULL
TMP2@images$TIP1=NULL
TMP2@images$NMP1=NULL
TMP2@images$TMP1=NULL
TMP2@images$NIP2=NULL
TMP2@images$TIP2=NULL
TMP2@images$NMP2=NULL
TMP2@images$NIP3=NULL
TMP2@images$TIP3=NULL
TMP2@images$NMP3=NULL
TMP2@images$TMP3=NULL

TMP2 <- FindSpatiallyVariableFeatures(TMP2, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(TMP2), r.metric = 5, slot = "data")
top.clusters <- SpatiallyVariableFeatures(TMP2)[SpatiallyVariableFeatures(TMP2) != 'max']
SpatialPlot(object = TMP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

TMP2 <- FindSpatiallyVariableFeatures(TMP2, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TMP2@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TMP2,file = 'spatialTMP2.Rds')
top.clusters <- head(SpatiallyVariableFeatures(TMP2, assay = "SCT"),8)
DefaultAssay(TMP2)='SCT'
SpatialFeaturePlot(object = TMP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


NMP2 = lung.integrated[,lung.integrated$orig.ident == 'NM_P2']
NMP2@images$NIP1=NULL
NMP2@images$TIP1=NULL
NMP2@images$NMP1=NULL
NMP2@images$TMP1=NULL
NMP2@images$NIP2=NULL
NMP2@images$TIP2=NULL
NMP2@images$TMP2=NULL
NMP2@images$NIP3=NULL
NMP2@images$TIP3=NULL
NMP2@images$NMP3=NULL
NMP2@images$TMP3=NULL

NMP2 <- FindSpatiallyVariableFeatures(NMP2, assay = "predictions", selection.method = "markvariogram", 
                                      features = rownames(NMP2), r.metric = 5, slot = "data")
top.clusters <- SpatiallyVariableFeatures(NMP2)[SpatiallyVariableFeatures(NMP2) != 'max']
SpatialPlot(object = NMP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


NMP2 <- FindSpatiallyVariableFeatures(NMP2, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NMP2@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(NMP2,file = 'spatialNMP2.Rds')
top.clusters <- head(SpatiallyVariableFeatures(NMP2, assay = "SCT"),8)
DefaultAssay(NMP2)='SCT'
SpatialFeaturePlot(object = NMP2, features = top.clusters, ncol = 4,alpha = c(0.1, 1))



####P3

TMP3 = lung.integrated[,lung.integrated$orig.ident == 'TM_P3']
TMP3@images$NIP1=NULL
TMP3@images$TIP1=NULL
TMP3@images$NMP1=NULL
TMP3@images$TMP1=NULL
TMP3@images$NIP2=NULL
TMP3@images$TIP2=NULL
TMP3@images$NMP2=NULL
TMP3@images$TMP2=NULL
TMP3@images$NIP3=NULL
TMP3@images$TIP3=NULL
TMP3@images$NMP3=NULL
TMP3 <- FindSpatiallyVariableFeatures(TMP3, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TMP3@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TMP3,file = 'spatialTMP3.Rds')
top.clusters <- head(SpatiallyVariableFeatures(TMP3, assay = "SCT"),8)
DefaultAssay(TMP3)='SCT'
SpatialFeaturePlot(object = TMP3, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

NMP3 = lung.integrated[,lung.integrated$orig.ident == 'NM_P3']
NMP3@images$NIP1=NULL
NMP3@images$TIP1=NULL
NMP3@images$NMP1=NULL
NMP3@images$TMP1=NULL
NMP3@images$NIP2=NULL
NMP3@images$TIP2=NULL
NMP3@images$NMP2=NULL
NMP3@images$TMP2=NULL
NMP3@images$NIP3=NULL
NMP3@images$TIP3=NULL
NMP3@images$TMP3=NULL

NMP3 <- FindSpatiallyVariableFeatures(NMP3, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NMP3@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(NMP3,file = 'spatialNMP3.Rds')
NMP3= readRDS(file = 'spatialNMP3.Rds')
top.clusters <- head(SpatiallyVariableFeatures(NMP3, assay = "SCT"),8)
DefaultAssay(NMP3)='SCT'
SpatialFeaturePlot(object = NMP3, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

TIP3 = lung.integrated[,lung.integrated$orig.ident == 'TI_P3']
TIP3@images$NIP1=NULL
TIP3@images$TIP1=NULL
TIP3@images$NMP1=NULL
TIP3@images$TMP1=NULL
TIP3@images$NIP2=NULL
TIP3@images$TIP2=NULL
TIP3@images$NMP2=NULL
TIP3@images$TMP2=NULL
TIP3@images$NIP3=NULL
TIP3@images$TMP3=NULL
TIP3@images$NMP3=NULL
TIP3 <- FindSpatiallyVariableFeatures(TIP3, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(TIP3@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(TIP3,file = 'spatialTIP3.Rds')
top.clusters <- head(SpatiallyVariableFeatures(TIP3, assay = "SCT"),8)
DefaultAssay(TIP3)='SCT'
SpatialFeaturePlot(object = TIP3, features = top.clusters, ncol = 4,alpha = c(0.1, 1))

NIP3 = lung.integrated[,lung.integrated$orig.ident == 'NI_P3']
NIP3@images$NIP1=NULL
NIP3@images$TIP1=NULL
NIP3@images$NMP1=NULL
NIP3@images$TMP1=NULL
NIP3@images$NIP2=NULL
NIP3@images$TIP2=NULL
NIP3@images$NMP2=NULL
NIP3@images$TMP2=NULL
NIP3@images$NMP3=NULL
NIP3@images$TIP3=NULL
NIP3@images$TMP3=NULL

NIP3 <- FindSpatiallyVariableFeatures(NIP3, assay = "SCT", selection.method = "markvariogram", 
                                      features = rownames(NIP3@assays$SCT@scale.data), r.metric = 5, slot = "data")
saveRDS(NIP3,file = 'spatialNIP3.Rds')
top.clusters <- head(SpatiallyVariableFeatures(NIP3, assay = "SCT"),8)
DefaultAssay(NIP3)='SCT'
SpatialFeaturePlot(object = NIP3, features = top.clusters, ncol = 4,alpha = c(0.1, 1))


cxcl14.markers = c('CXCL14','CLDN2','CEACAM5','CEACAM6','MDK')
SpatialFeaturePlot(object = TIP3, features = cxcl14.markers, ncol = 4,alpha = c(0.1, 1))

DefaultAssay(lung.integrated)='SCT'
SpatialFeaturePlot(lung.integrated,images = c('TIP1','TMP1'), 
                   features = cxcl14.markers, pt.size.factor = 1, alpha = c(0.1, 1))

####selected location analysis, 2021-10-12, R version 3####



Idents(lung.integrated)=lung.integrated$celltype_prediction

SpatialPlot(lung.integrated,images = c('TIP2'),alpha=c(0.1,1),cols = color_pelette)
SpatialDimPlot(lung.integrated,images = c('TMP2'),alpha=c(0.1,1),cols = color_pelette)
SpatialDimPlot(lung.integrated,images = c('TMP3'),alpha=c(0.1,1),cols = color_pelette)
SpatialDimPlot(lung.integrated,images = c('TIP3'),alpha=c(0.1,1),cols = color_pelette)

ISpatialDimPlot(lung.integrated,image = 'TIP2',alpha=c(0.1,1))
LinkedDimPlot(lung.integrated,image = 'TIP2')

TIP2 = readRDS(file = '../MPSC&ST/variables/spatialTIP2.Rds')
Idents(TIP2)=TIP2$celltype_prediction

SpatialDimPlot(TIP2,images = c('TIP2'),alpha=c(0.1,1),cols = color_pelette)
ISpatialDimPlot(TIP2,image = 'TIP2')

####TIP2 locational analysis for selected regions####
TIP2.SRF = read.csv(file='../ourData/stRNA/edit images/TIP2 locations.csv.')
TIP2.SF = read.csv(file='../ourData/stRNA/Count/D_5TL/spatial/tissue_positions_list.csv',header = F,row.names = 1)
TIP2.SF$spot = paste0(rownames(TIP2.SF),'_6')
rownames(TIP2.SF)=paste(TIP2.SF$V3,TIP2.SF$V4)
rownames(TIP2.SRF)=paste(TIP2.SRF$ID_row,TIP2.SRF$ID_col)
TIP2.SRF$spot = TIP2.SF[rownames(TIP2.SRF),'spot']
nrow(TIP2.SRF[is.na(TIP2.SRF$spot),])
rownames(TIP2.SRF)=TIP2.SRF$spot
TIP2.SRF$Region = paste0('SR',TIP2.SRF$Region)
TIP2$SRF = TIP2.SRF[colnames(TIP2),'Region']
TIP2$SRF[is.na(TIP2$SRF)]='NotSelect'

#s.TIP2 = TIP2[,TIP2$SRF != 'NotSelect']
TIP2$SRF[TIP2$SRF == 'SR4']='SR1'
SpatialDimPlot(TIP2,images = c('TIP2'),alpha=c(0.1),
               group.by = 'SRF',
               cols = c('grey',pal_npg('nrc')(3)))

table(TIP2$celltype_prediction,TIP2$SRF)

Idents(TIP2) = TIP2$SRF
de_markers.TIP2 <- FindAllMarkers(TIP2, assay = 'Spatial',logfc.threshold = 0.05, 
                               test.use = "wilcox", 
                               min.pct = 0.01, 
                               min.diff.pct = 0.1, 
                               only.pos = TRUE)
saveRDS(de_markers.TIP2,file = '../MPSC&ST/variables/TIP2 SR markers.Rds')
de_markers.TIP2 = readRDS(file = '../MPSC&ST/variables/TIP2 SR markers.Rds')

SpatialFeaturePlot(lung.integrated,images = c('TIP2','TMP2'),alpha=c(0.1,1),
                   features = c('SLC1A7','MMP7'))
TIP2.svrs = data.frame(Gene = rownames(TIP2@assays[["SCT"]]@meta.features),
                       Order = TIP2@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable.rank"]] ,
                       Score = TIP2@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable"]])
TIP2.svrs = TIP2.svrs[!is.na(TIP2.svrs$Order),]
SR1.markers = filter(de_markers.TIP2,p_val_adj<0.05 & cluster == 'SR1')
rownames(SR1.markers)=SR1.markers$gene
SR2.markers = filter(de_markers.TIP2,p_val_adj<0.05 & cluster == 'SR2')
rownames(SR2.markers)=SR2.markers$gene
SR3.markers = filter(de_markers.TIP2,p_val_adj<0.05 & cluster == 'SR3')
rownames(SR3.markers)=SR3.markers$gene
####MMP7,TB type
TIP2.svrs$SR1 = SR1.markers[as.character(TIP2.svrs$Gene),'avg_logFC']
TIP2.svrs$SR2 = SR2.markers[as.character(TIP2.svrs$Gene),'avg_logFC']
TIP2.svrs$SR3 = SR3.markers[as.character(TIP2.svrs$Gene),'avg_logFC']
####
library(clusterProfiler)
load(file = 'E:/project/metaboliteProteinInteraction/variables/KEGG_path2Gene.RData')
path2cate <- read.csv(file='E:/data/KEGG/KEGG hsa pathway list.csv',row.names = 1)
rownames(path2cate) = path2cate$Name
SR1.paths = enricher(as.character(SR1.markers$gene),pvalueCutoff = 1.2,
                     TERM2GENE = path2Gene[,c(1,2)])@result
SR2.paths = enricher(as.character(SR2.markers$gene),pvalueCutoff = 1.2,
                     TERM2GENE = path2Gene[,c(1,2)])@result

SR3.paths = enricher(as.character(SR3.markers$gene),pvalueCutoff = 1.2,
                     TERM2GENE = path2Gene[,c(1,2)])@result
####TMP2 selected regions####

TMP2 = readRDS(file = '../MPSC&ST/variables/spatialTMP2.Rds')
ISpatialDimPlot(TMP2,image = 'TMP2',alpha=c(0.1,1),group.by = 'celltype_prediction')
TMP2.SRF = read.csv(file='../ourData/stRNA/edit images/TMP2 locations.csv.')
TMP2.SF = read.csv(file='../ourData/stRNA/Count/D_5TM/spatial/tissue_positions_list.csv',header = F,row.names = 1)
TMP2.SF$spot = paste0(rownames(TMP2.SF),'_8')
rownames(TMP2.SF)=paste(TMP2.SF$V3,TMP2.SF$V4)
rownames(TMP2.SRF)=paste(TMP2.SRF$ID_row,TMP2.SRF$ID_col)
TMP2.SRF$spot = TMP2.SF[rownames(TMP2.SRF),'spot']
nrow(TMP2.SRF[is.na(TMP2.SRF$spot),])
rownames(TMP2.SRF)=TMP2.SRF$spot
TMP2.SRF$Region = paste0('SR',TMP2.SRF$Region)
TMP2$SRF = TMP2.SRF[colnames(TMP2),'Region']
TMP2$SRF[is.na(TMP2$SRF)]='NotSelect'

SpatialDimPlot(TMP2,images = c('TMP2'),alpha=c(0.5),
               group.by = 'SRF',
               cols = c('grey',pal_npg('nrc')(6)[c(5,6)]))

table(TMP2$celltype_prediction,TMP2$SRF)

Idents(TMP2) = TMP2$SRF
de_markers.TMP2 <- FindAllMarkers(TMP2, assay = 'Spatial',logfc.threshold = 0.05, 
                                  test.use = "wilcox", 
                                  min.pct = 0.01, 
                                  min.diff.pct = 0.1, 
                                  only.pos = TRUE)
saveRDS(de_markers.TMP2,file = '../MPSC&ST/variables/TMP2 SR markers.Rds')
de_markers.TMP2 = readRDS(file = '../MPSC&ST/variables/TMP2 SR markers.Rds')

SpatialFeaturePlot(lung.integrated,images = c('TMP2','TIP2'),alpha=c(0.1,1),
                   features = c('CHI3L1','SERPINA1'),
                   slot = 'counts')
TMP2.svrs = data.frame(Gene = rownames(TMP2@assays[["SCT"]]@meta.features),
                       Order = TMP2@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable.rank"]] ,
                       Score = TMP2@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable"]])
TMP2.svrs = TMP2.svrs[!is.na(TMP2.svrs$Order),]
SRw.markers.TMP2 = filter(de_markers.TMP2,p_val_adj<0.05 & cluster == 'SRW')
rownames(SRw.markers.TMP2)=SRw.markers.TMP2$gene
SRy.markers.TMP2 = filter(de_markers.TMP2,p_val_adj<0.05 & cluster == 'SRY')
rownames(SRy.markers.TMP2)=SRy.markers.TMP2$gene

####
TMP2.svrs$SRw= SRw.markers.TMP2[as.character(TMP2.svrs$Gene),'avg_logFC']
TMP2.svrs$SRy = SRy.markers.TMP2[as.character(TMP2.svrs$Gene),'avg_logFC']
####
SRw.paths.TMP2 = enricher(as.character(SRw.markers.TMP2$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result

SRy.paths.TMP2 = enricher(as.character(SRy.markers.TMP2$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result
####TIP3 selected regions####

TIP3 = readRDS(file = '../MPSC&ST/variables/spatialTIP3.Rds')
ISpatialDimPlot(TIP3,image = 'TIP3',alpha=c(0.1,1),group.by = 'celltype_prediction')
TIP3.SRF = read.csv(file='../ourData/stRNA/edit images/TIP3 locations.csv.')
TIP3.SF = read.csv(file='../ourData/stRNA/Count/P8_I_T/spatial/tissue_positions_list.csv',header = F,row.names = 1)
TIP3.SF$spot = paste0(rownames(TIP3.SF),'_10')
rownames(TIP3.SF)=paste(TIP3.SF$V3,TIP3.SF$V4)
rownames(TIP3.SRF)=paste(TIP3.SRF$ID_row,TIP3.SRF$ID_col)
TIP3.SRF$spot = TIP3.SF[rownames(TIP3.SRF),'spot']
nrow(TIP3.SRF[is.na(TIP3.SRF$spot),])
rownames(TIP3.SRF)=TIP3.SRF$spot
TIP3.SRF$Region = paste0('SR',TIP3.SRF$Region)
TIP3$SRF = TIP3.SRF[colnames(TIP3),'Region']
TIP3$SRF[is.na(TIP3$SRF)]='NotSelect'

SpatialDimPlot(TIP3,images = c('TIP3'),alpha=c(0.5),
               group.by = 'SRF',
               cols = c('grey',pal_npg('nrc')(3)))

table(TIP3$celltype_prediction,TIP3$SRF)

Idents(TIP3) = TIP3$SRF
de_markers.TIP3 <- FindAllMarkers(TIP3, assay = 'Spatial',logfc.threshold = 0.05, 
                                  test.use = "wilcox", 
                                  min.pct = 0.01, 
                                  min.diff.pct = 0.1, 
                                  only.pos = TRUE)
saveRDS(de_markers.TIP3,file = '../MPSC&ST/variables/TIP3 SR markers.Rds')
de_markers.TIP3 = readRDS(file = '../MPSC&ST/variables/TIP3 SR markers.Rds')

SpatialFeaturePlot(lung.integrated,images = c('TIP3','TIP2'),alpha=c(0.1,1),
                   features = c('MUC1','MAL2'),
                   slot = 'counts')
TIP3.svrs = data.frame(Gene = rownames(TIP3@assays[["SCT"]]@meta.features),
                       Order = TIP3@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable.rank"]] ,
                       Score = TIP3@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable"]])
TIP3.svrs = TIP3.svrs[!is.na(TIP3.svrs$Order),]
SR1.markers.TIP3 = filter(de_markers.TIP3,p_val_adj<0.05 & cluster == 'SR1')
rownames(SR1.markers.TIP3)=SR1.markers.TIP3$gene
SR2.markers.TIP3 = filter(de_markers.TIP3,p_val_adj<0.05 & cluster == 'SR2')
rownames(SR2.markers.TIP3)=SR2.markers.TIP3$gene
SR3.markers.TIP3 = filter(de_markers.TIP3,p_val_adj<0.05 & cluster == 'SR3')
rownames(SR3.markers.TIP3)=SR3.markers.TIP3$gene

length(intersect(rownames(SR1.markers),rownames(SR1.markers.TIP3)))
length(intersect(rownames(SR2.markers),rownames(SR2.markers.TIP3)))
length(intersect(rownames(SR3.markers),rownames(SR3.markers.TIP3)))


####
TIP3.svrs$SR1 = SR1.markers.TIP3[as.character(TIP3.svrs$Gene),'avg_logFC']
TIP3.svrs$SR2 = SR2.markers.TIP3[as.character(TIP3.svrs$Gene),'avg_logFC']
TIP3.svrs$SR3 = SR3.markers.TIP3[as.character(TIP3.svrs$Gene),'avg_logFC']
####
SR1.paths.TIP3 = enricher(as.character(SR1.markers.TIP3$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result

SR2.paths.TIP3 = enricher(as.character(SR2.markers.TIP3$gene),pvalueCutoff = 1.2,
                     TERM2GENE = path2Gene[,c(1,2)])@result
SR3.paths.TIP3 = enricher(as.character(SR3.markers.TIP3$gene),pvalueCutoff = 1.2,
                     TERM2GENE = path2Gene[,c(1,2)])@result


####TMP3, selected regions####

TMP3 = readRDS(file = '../MPSC&ST/variables/spatialTMP3.Rds')
ISpatialDimPlot(TMP3,image = 'TMP3',group.by = 'celltype_prediction')
TMP3.SRF = read.csv(file='../ourData/stRNA/edit images/TMP3 locations.csv.')
TMP3.SF = read.csv(file='../ourData/stRNA/Count/P8_M_T/spatial/tissue_positions_list.csv',header = F,row.names = 1)
TMP3.SF$spot = paste0(rownames(TMP3.SF),'_12')
rownames(TMP3.SF)=paste(TMP3.SF$V3,TMP3.SF$V4)
rownames(TMP3.SRF)=paste(TMP3.SRF$ID_row,TMP3.SRF$ID_col)
TMP3.SRF$spot = TMP3.SF[rownames(TMP3.SRF),'spot']
nrow(TMP3.SRF[is.na(TMP3.SRF$spot),])
rownames(TMP3.SRF)=TMP3.SRF$spot
TMP3.SRF$Region = paste0('SR',TMP3.SRF$Region)
TMP3$SRF = TMP3.SRF[colnames(TMP3),'Region']
TMP3$SRF[is.na(TMP3$SRF)]='NotSelect'

SpatialDimPlot(TMP3,images = c('TMP3'),alpha=c(0.5),
               group.by = 'SRF',
               cols = c('grey',pal_npg('nrc')(3)),
               stroke = 0)

table(TMP3$celltype_prediction,TMP3$SRF)

Idents(TMP3) = TMP3$SRF
#TMP3[['integrated']]=NULL
de_markers.TMP3 <- FindAllMarkers(TMP3, assay = 'Spatial',
                                  logfc.threshold = 0.05, 
                                  test.use = "wilcox", 
                                  min.pct = 0.01, 
                                  min.diff.pct = 0.1, 
                                  only.pos = TRUE)
saveRDS(de_markers.TMP3,file = '../MPSC&ST/variables/TMP3 SR markers.Rds')
de_markers.TMP3 = readRDS(file = '../MPSC&ST/variables/TMP3 SR markers.Rds')

SpatialFeaturePlot(lung.integrated,images = c('TMP3','TIP3'),alpha=c(0.1,1),
                   features = c('PGC','SLPI'),
                   slot = 'counts')
TMP3.svrs = data.frame(Gene = rownames(TMP3@assays[["SCT"]]@meta.features),
                       Order = TMP3@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable.rank"]] ,
                       Score = TMP3@assays[["SCT"]]@meta.features[["markvariogram.spatially.variable"]])
TMP3.svrs = TMP3.svrs[!is.na(TMP3.svrs$Order),]
SR1.markers.TMP3 = filter(de_markers.TMP3,p_val_adj<0.05 & cluster == 'SR1')
rownames(SR1.markers.TMP3)=SR1.markers.TMP3$gene
SR2.markers.TMP3 = filter(de_markers.TMP3,p_val_adj<0.05 & cluster == 'SR2')
rownames(SR2.markers.TMP3)=SR2.markers.TMP3$gene
SR3.markers.TMP3 = filter(de_markers.TMP3,p_val_adj<0.05 & cluster == 'SR3')
rownames(SR3.markers.TMP3)=SR3.markers.TMP3$gene

####
TMP3.svrs$SR1 = SR1.markers.TMP3[as.character(TMP3.svrs$Gene),'avg_logFC']
TMP3.svrs$SR2 = SR2.markers.TMP3[as.character(TMP3.svrs$Gene),'avg_logFC']
TMP3.svrs$SR3 = SR3.markers.TMP3[as.character(TMP3.svrs$Gene),'avg_logFC']
####
SR1.paths.TMP3 = enricher(as.character(SR1.markers.TMP3$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result

SR2.paths.TMP3 = enricher(as.character(SR2.markers.TMP3$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result
SR3.paths.TMP3 = enricher(as.character(SR3.markers.TMP3$gene),pvalueCutoff = 1.2,
                          TERM2GENE = path2Gene[,c(1,2)])@result

length(intersect(rownames(SR1.markers),rownames(SR1.markers.TMP3)))
length(intersect(rownames(SR2.markers),rownames(SR2.markers.TMP3)))
length(intersect(rownames(SR3.markers),rownames(SR3.markers.TMP3)))

####Merge TIP2, TIP3, TMP3 regional features####

####Cell types among differernt regions####
TIP2.com2Type = table(TIP2$celltype_prediction,TIP2$SRF)[,-1]
TIP3.com2Type = table(TIP3$celltype_prediction,TIP3$SRF)[,-1]
TMP3.com2Type = table(TMP3$celltype_prediction,TMP3$SRF)[,-1]
TMP2.com2Type = table(TMP2$celltype_prediction,TMP2$SRF)[,-1]

# TIP2.com2Type = table(TIP2$celltype_prediction,TIP2$SRF)[,-1]
# TIP3.com2Type = table(TIP3$celltype_prediction,TIP3$SRF)[,-1]
# TMP3.com2Type = table(TMP3$celltype_prediction,TMP3$SRF)[,-1]

library(reshape2)
TIP2.com2Type = t(TIP2.com2Type)/apply(TIP2.com2Type,2,sum)
TIP2.com2Type = data.frame(TIP2.com2Type)
TIP2.com2Type$From = 'TIP2'

TIP3.com2Type = t(TIP3.com2Type)/apply(TIP3.com2Type,2,sum)
TIP3.com2Type = data.frame(TIP3.com2Type)
TIP3.com2Type$From = 'TIP3'

TMP3.com2Type = t(TMP3.com2Type)/apply(TMP3.com2Type,2,sum)
TMP3.com2Type = data.frame(TMP3.com2Type)
TMP3.com2Type$From = 'TMP3'

TMP2.com2Type = t(TMP2.com2Type)/apply(TMP2.com2Type,2,sum)
TMP2.com2Type = data.frame(TMP2.com2Type)
TMP2.com2Type$From = 'TMP2'

merge.com2Type = rbind(TIP2.com2Type,TIP3.com2Type,TMP3.com2Type,TMP2.com2Type)

merge.com2Type$Var2 = factor(merge.com2Type$Var2,
                             levels=c("Epithelial cell","Endothelial cell",
                                      "Fibroblast cell","Myeloid cell",
                                      "B cell","Mast cell" , "T&NK cell"))
ggplot(data = merge.com2Type,mapping = aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat = 'identity',width = 0.65)+
  facet_grid(~From,scales = 'free',space = 'free')+scale_fill_manual(values = cellType.cols)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='Ratio')

ptest <- merge.com2Type %>% filter(Var1 %in% c('SR1','SR2','SR3')) %>% group_by(Var2) %>% summarize(p.value = kruskal.test(Freq ~ Var1)$p.value)
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/P2 P3 regional cell type compositions.pdf
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/P2 P3 regional epi subtype compositions.pdf

####Marker genes
nrow(SR1.markers.TIP3)#46
nrow(SR2.markers.TIP3)#4102
nrow(SR3.markers.TIP3)#2319

nrow(SR1.markers.TMP3)#50
nrow(SR2.markers.TMP3)#12
nrow(SR3.markers.TMP3)#132


nrow(SR1.markers)#793
nrow(SR2.markers)#25
nrow(SR3.markers)#3172


nrow(SRw.markers.TMP2)#3447
nrow(SRy.markers.TMP2)#83

num.df = data.frame(Number = c(46,4102,2319,50,12,132,793,25,3172,3447,84),
                    PR = c(rep( c('PR1','PR2','PR3'),3),'PRm','PRp'),
                    From =c( rep(c('TIP3','TMP3','TIP2'),each = 3),'TMP2','TMP2'))
num.df$From = factor(num.df$From,levels = c('TIP2','TMP2','TIP3','TMP3'))
ggplot(num.df,aes(x=PR,y=Number,fill=From))+geom_bar(stat='identity',width = 0.6)+
  facet_grid(~From,scales = 'free',space = 'free')+scale_fill_npg()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='Number')+scale_y_log10()

length(intersect(rownames(SR1.markers),rownames(SR1.markers.TMP3)))#8

length(intersect(rownames(SR1.markers.TIP3),rownames(SR1.markers.TMP3)))#4
intersect(rownames(SR2.markers.TIP3),rownames(SR2.markers.TMP3))
length(intersect(rownames(SR3.markers.TIP3),rownames(SR3.markers.TMP3)))#53
length(intersect(rownames(SR3.markers),rownames(SR3.markers.TMP3)))#63
length(intersect(rownames(SR3.markers),intersect(rownames(SR3.markers.TIP3),rownames(SR3.markers.TMP3))))#28
intersect(rownames(SR1.markers),intersect(rownames(SR1.markers.TIP3),rownames(SR1.markers.TMP3)))
intersect(rownames(SR2.markers),intersect(rownames(SR2.markers.TIP3),rownames(SR2.markers.TMP3)))


p1=VlnPlot(TIP2, group.by = "SRF", features = 'PDZK1IP1', pt.size = 0, ncol = 1) + 
  NoLegend()
p2=VlnPlot(TIP3, group.by = "SRF", features = 'PDZK1IP1', pt.size = 0, ncol = 1) + 
  NoLegend()
p3=VlnPlot(TMP3, group.by = "SRF", features = 'PDZK1IP1', pt.size = 0, ncol = 1) + 
  NoLegend()
p1|p2|p3
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/vlnplot PDZK1IP1 in TIP2 TIP3 TMP3.pdf

all.SR3 = intersect(rownames(SR3.markers),intersect(rownames(SR3.markers.TIP3),rownames(SR3.markers.TMP3)))
all.SR3 = data.frame(Gene = all.SR3,row.names = all.SR3)
all.SR3$TIP2.S = SR3.markers[rownames(all.SR3),'avg_logFC']
all.SR3$TIP3.S = SR3.markers.TIP3[rownames(all.SR3),'avg_logFC']
all.SR3$TMP3.S = SR3.markers.TMP3[rownames(all.SR3),'avg_logFC']
all.SR3$meanS = apply(as.matrix(all.SR3[,-1]),1,mean)
all.SR3$sd = apply(as.matrix(all.SR3[,-1]),1,sd)
all.SR3 = arrange(all.SR3,-meanS,-sd)
all.SR3$Gene = factor(all.SR3$Gene,levels = as.character(all.SR3$Gene))
ggplot(all.SR3,mapping = aes(x=Gene,y=meanS))+geom_bar(stat = 'identity',width = 0.6,fill='blue')+
  geom_errorbar(aes(ymin=meanS-sd,ymax=meanS+sd),width = 0.1)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='mean(FC)')
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/consisSR3Markers meanFC barplot.pdf
write.csv(all.SR3,file = '../MPSC&ST/variables/all three samples regional features for SR3.csv')

p1= SpatialDimPlot(TIP2,images = c('TIP2'),alpha=c(0.5),
                   group.by = 'SRF',
                   cols = c('grey','yellow','green','red'),pt.size.factor = 1.2,stroke = 0)
p2= SpatialDimPlot(TIP3,images = c('TIP3'),alpha=c(0.5),
                   group.by = 'SRF',
                   cols = c('grey','yellow','green','red'),pt.size.factor = 1.2,stroke = 0)
p3= SpatialDimPlot(TMP3,images = c('TMP3'),alpha=c(0.5),
               group.by = 'SRF',
               cols = c('grey','yellow','green','red'),pt.size.factor = 1.2,stroke = 0)
p1 | p2 | p3
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TIP3','TMP3'),alpha = c(0.2,1),
                   pt.size.factor = 1.2,stroke = 0,crop=F,
                   features = c('PRDX4','ILF2','ELOC'),
                   slot = 'counts')
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/TIP2 TIP3 TMP3 three consis markers for SR3.pdf
#SR1 shared marker
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TIP3','TMP3'),
                   alpha=c(0.4,1),
                   pt.size.factor = 1.5,stroke = 0,crop=F,
                   features = c('PDZK1IP1'),
                   slot = 'counts')
SpatialFeaturePlot(lung.integrated,images = c("TMP2"),
                   alpha=c(0.4,1),
                   pt.size.factor = 1.5,stroke = 0,crop=F,
                   features = c('PDZK1IP1'),
                   slot = 'counts')
library(VennDiagram)
library(grid)
venn.p <- venn.diagram(x=list(TIP3.SR1=rownames(SR1.markers.TIP3),
                              TMP3.SR1=rownames(SR1.markers.TMP3)),
                       filename=NULL,
                       fill = c('red','blue'))
grid.draw(venn.p)

venn.p <- venn.diagram(x=list(TIP3.SR2=rownames(SR2.markers.TIP3),
                              TMP3.SR2=rownames(SR2.markers.TMP3)),
                       filename=NULL,
                       fill = c('red','blue'))
grid.draw(venn.p)
venn.p <- venn.diagram(x=list(TIP3.SR3=rownames(SR3.markers.TIP3),
                              TMP3.SR3=rownames(SR3.markers.TMP3)),
                       filename=NULL,
                       fill = c('red','blue'))
grid.draw(venn.p)


venn.p <- venn.diagram(x=list(TIP3.SR1=rownames(SR1.markers.TIP3),
                              TMP3.SR1=rownames(SR1.markers.TMP3),
                              TIP2.SR1 = rownames(SR1.markers)),
                       filename=NULL,
                       fill = c('red','blue','pink'))
grid.draw(venn.p)

venn.p <- venn.diagram(x=list(TIP3.SR2=rownames(SR2.markers.TIP3),
                              TMP3.SR2=rownames(SR2.markers.TMP3),
                              TIP2.SR2 = rownames(SR2.markers)),
                       filename=NULL,
                       fill = c('red','blue','pink'))
grid.draw(venn.p)
venn.p <- venn.diagram(x=list(TIP3.SR3=rownames(SR3.markers.TIP3),
                              TMP3.SR3=rownames(SR3.markers.TMP3),
                              TIP2.SR3 = rownames(SR3.markers)),
                       filename=NULL,
                       fill = c('red','blue','pink'))
grid.draw(venn.p)

venn.p <- venn.diagram(x=list(TIP2.SR1=rownames(SR1.markers),
                              TIP2.SR2=rownames(SR2.markers),
                              TIP2.SR3 = rownames(SR3.markers),
                              TMP2.SRW = rownames(SRw.markers.TMP2),
                              TMP2.SRY = rownames(SRy.markers.TMP2)),
                       filename=NULL,
                       fill = c('red','blue','pink','yellow','blue'))
grid.draw(venn.p)
venn.p <- venn.diagram(x=list(TIP2.SR1=rownames(SR1.markers),
                              TIP2.SR2=rownames(SR2.markers),
                              TIP2.SR3 = rownames(SR3.markers)),
                       filename=NULL,
                       fill = c('red','blue','pink'))
grid.draw(venn.p)
####Specific markers
##SR1 TIP2
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('MMP7','TM4SF4','CD55'),
                   slot = 'counts')
##SR3 TIP2
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('TFF3','SLC1A7','HPGD'),
                   slot = 'counts')
##SR2 TIP3
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('C16orf89','GDF15','MALL'),
                   slot = 'counts')
##SR3 TIP3
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('IGLC7','MIA','PGC'),
                   slot = 'counts')
##SR2 TMP3
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('CDC42EP1','CES1','RNF145'),
                   slot = 'counts')#poor effects in image
##SR3 TMP3
SpatialFeaturePlot(lung.integrated,images = c("TIP2",'TMP2','TIP3','TMP3'),
                   pt.size.factor = 1.5,stroke = 0,crop=F,alpha = c(0.1,1),
                   features = c('DDIT3','HOPX','TSTD1'),
                   slot = 'counts')

de_markers.TIP3$From = 'TIP3'
de_markers.TMP3$From = 'TMP3'
de_markers.TIP2$From = 'TIP2'

de_markers.m = rbind(de_markers.TIP3,de_markers.TMP3,de_markers.TIP2)
de_markers.top = de_markers.m %>% filter(cluster != 'NotSelect') %>% 
  group_by(From,cluster) %>% top_n(-5, p_val_adj)%>% top_n(5, avg_logFC) 
top.genes = unique(as.character(de_markers.top$gene))

de_markers.m = de_markers.m[de_markers.m$gene %in% top.genes & de_markers.m$cluster != 'NotSelect',]
de_markers.top$cluster = factor(de_markers.top$cluster,levels = c('SR1','SR2','SR3'))
ggplot(de_markers.top,aes(x=gene,y=avg_logFC,color=From))+geom_point(aes(size = -log10(p_val_adj)))+
  facet_grid(cluster ~ From,scales = 'free')+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_color_npg()+
  scale_size_continuous(range = c(1,3))
#file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/specific markers point plot.pdf
####Markers pathways
SR1.paths$From='TIP2'
SR1.paths$SRF = 'SR1'
SR1.paths.TMP3$From = 'TMP3'
SR1.paths.TMP3$SRF = 'SR1'
SR1.paths.TIP3$From = 'TIP3'
SR1.paths.TIP3$SRF = 'SR1'

SR3.paths$From = 'TIP2'
SR3.paths$SRF = 'SR3'
SR2.paths.TIP3$From = 'TIP3'
SR2.paths.TIP3$SRF = 'SR2'
SR2.paths.TMP3$From = 'TMP3'
SR2.paths.TMP3$SRF = 'SR2'
SR2.paths$From = 'TIP2'
SR2.paths$SRF = 'SR2'

SR3.paths.TIP3$From = 'TIP3'
SR3.paths.TIP3$SRF = 'SR3'
SR3.paths.TMP3$From = 'TMP3'
SR3.paths.TMP3$SRF='SR3'

SRy.paths.TMP2$From = 'TMP2'
SRy.paths.TMP2$SRF='SRy'
SRw.paths.TMP2$From = 'TMP2'
SRw.paths.TMP2$SRF='SRw'

SR.paths = rbind(SR1.paths,SR1.paths.TMP3,SR1.paths.TIP3,
                 SR2.paths,SR2.paths.TIP3,SR2.paths.TMP3,
                 SR3.paths,SR3.paths.TIP3,SR3.paths.TMP3)
SR.paths = rbind(SRy.paths.TMP2,
                 SRw.paths.TMP2)
SR.paths.sig = SR.paths[SR.paths$pvalue<0.05 & SR.paths$Count>=2,]
SR.paths.sig$cate = path2cate[as.character(SR.paths.sig$Description),'cate']
SR.paths.sig = SR.paths.sig[SR.paths.sig$cate != '6. Human Diseases',]
SR.paths.sig = arrange(SR.paths.sig,cate,desc(p.adjust))
SR.paths.sig = group_by(SR.paths.sig,SRF,From) %>% top_n(-10, pvalue) 
SR.paths.sig$Description = factor(SR.paths.sig$Description,
                                  levels = unique(SR.paths.sig$Description))
ggplot(data = SR.paths.sig,mapping = aes(x=Description,y=-log10(pvalue),fill=cate))+
  geom_bar(stat = 'identity',width = 0.6)+
  facet_grid(From~SRF,scales = 'free_y',space = 'free')+
  coord_flip()+
  theme_bw()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='-log10(p.adjust)')+
  scale_fill_manual(values=brewer.pal(5,'Dark2'))
##file:///H:/project/single cell/MPSC&ST/figure results/spatial/three samples/locational/barplot selected regions enriched pathways.pdf

SR.GN = data.frame(N = c(nrow(SR1.markers.TIP3),#46
                         nrow(SR2.markers.TIP3),#4102
                         nrow(SR3.markers.TIP3),#2319
                         
                         nrow(SR1.markers.TMP3),#50
                         nrow(SR2.markers.TMP3),#12
                         nrow(SR3.markers.TMP3),#132
                         
                         
                         nrow(SR1.markers),#793
                         nrow(SR2.markers),#25
                         nrow(SR3.markers)#3172
                         ),
                   Regions = c(rep(c('SR1','SR2','SR3'),3)),
                   From = c(rep(c('TIP3','TMP3','TIP2'),each = 3))
)
ggplot(data = SR.GN,mapping = aes(x=From,y=N,fill=Regions))+geom_bar(stat='identity',width = 0.65)+
  facet_grid(~Regions)+
  theme_bw()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='Number')+scale_y_log10()+
  scale_fill_manual(values=c('yellow','green','red'))
p1=VlnPlot(TIP2, group.by = "SRF", features = c("LSM5","MDH2","MRFAP1","RANBP1"), pt.size = 0, ncol = 1) + 
  NoLegend()
p2=VlnPlot(TIP3, group.by = "SRF", features = c("LSM5","MDH2","MRFAP1","RANBP1"), pt.size = 0, ncol = 1) + 
  NoLegend()
p3=VlnPlot(TMP3, group.by = "SRF", features = c("LSM5","MDH2","MRFAP1","RANBP1"), pt.size = 0, ncol = 1) + 
  NoLegend()
p1|p2|p3
####Notes
####Compare to TCGA-LUAD
# intersect(int.genes,fav.genes)
# [1] "AK2"    "PMM1"   "PRR15L"
# > intersect(int.genes,unfav.genes)
# [1] "LSM5"   "MDH2"   "MRFAP1" "RANBP1"

####Survival, see file:///E:/codes/R codes/survivalPlot.R

VlnPlot(lung.integrated, group.by = "orig.ident", features = c("TNFRSF18"), pt.size = 0, ncol = 1) + 
  NoLegend()
VlnPlot(lung.integrated, group.by = "orig.ident", features = c("CD79B"), pt.size = 0, ncol = 1) + 
  NoLegend()
VlnPlot(lung.integrated, group.by = "orig.ident", features = c("ALOX15B", "LCN2", "PRSS12" , "WIF1"), pt.size = 0, ncol = 1) + 
  NoLegend()

VlnPlot(lung.integrated, group.by = "orig.ident", features = c('RGS1','CD82','IGHM'), pt.size = 0, ncol = 1) + 
  NoLegend()
SpatialFeaturePlot(lung.integrated,images = c('NIP2','TIP2'),features = 'WIF1',alpha = c(0.1,1))
SpatialFeaturePlot(lung.integrated,images = c('NMP1','TMP1'),features = 'TNFRSF18',alpha = c(0.1,1))
SpatialFeaturePlot(lung.integrated,images = c('NMP1','TMP1','NIP1','TIP1'),features = 'TNFRSF18',pt.size.factor = 2,alpha = c(0.5,1))
SpatialFeaturePlot(lung.integrated,images = c('NMP1','TMP1','NIP1','TIP1'),features = 'TNFRSF18',pt.size.factor = 2,alpha = c(0.5,1))

####Table S2
de_markers.m = rbind(de_markers.TIP3,de_markers.TMP3,de_markers.TIP2,de_markers.TMP2)
de_markers.top = de_markers.m %>% filter(cluster != 'NotSelect') %>% 
  group_by(From,cluster) %>% filter(p_val_adj<0.01 & abs(avg_logFC)>1)
write.table(de_markers.top,file = '../MPSC&ST/variables/selectedRegion markers.txt',sep = '\t')
