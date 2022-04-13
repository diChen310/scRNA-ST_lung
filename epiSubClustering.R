####Sub-clustering of the Epithelial cell, 2021-08-17####
library(Seurat)
library(SeuratObject)
library(cowplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(harmony)
library(ggsci)
setwd('H:/project/single cell/MPSC&ST/')
#
####Sub-clustering of the Epithelial cell, 2021-08-17####

#
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
epi.data = alldata[,alldata$CellTypeManully == 'Epithelial cell']
dim(epi.data)# 33261 34161
# selected_c <- WhichCells(epi.data, expression = nFeature_RNA > 200)
# selected_f <- rownames(epi.data)[Matrix::rowSums(epi.data) > 10]
# epi.data <- subset(epi.data, features = selected_f, cells = selected_c)
# epi.data@active.assay
epi.data = NormalizeData(epi.data)
epi.data = FindVariableFeatures(epi.data)
epi.data = ScaleData(epi.data, vars.to.regress = c("nFeature_RNA", "percent_mito"))
epi.data = RunPCA(epi.data, verbose = T, npcs = 50)
epi.data = RunUMAP(epi.data, dims = 1:50, verbose = T)
epi.data = RunHarmony(epi.data,c( "orig.ident" ))

epi.data <- RunUMAP(epi.data,reduction = "harmony", dims = 1:30)
epi.data <- FindNeighbors(epi.data, reduction = "harmony",dims = 1:30, verbose = FALSE)
names(epi.data@graphs)

for (res in c(0.1,0.2,0.3,0.5)) {
  epi.data <- FindClusters(epi.data, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
plot_grid(ncol = 3, DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.3",pt.size=0.01) + 
            ggtitle("louvain_0.3"), DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
            ggtitle("louvain_0.5"), DimPlot(epi.data, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
            ggtitle("louvain_0.1"))
sel.clust = "RNA_snn_res.0.1"
epi.data <- SetIdent(epi.data, value = sel.clust)
table(epi.data@active.ident)
DimPlot(epi.data, reduction = "umap", label = T)
VlnPlot(epi.data, features = c("G2M.Score", "percent_hb",'percent_mito','nFeature_RNA'), 
        group.by = "RNA_snn_res.0.1", 
        ncol = 2, pt.size = 0)
epi.data$orig.ident = factor(epi.data$orig.ident,
                             levels = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
                                        'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
                                        'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
                                        'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4'))

DimPlot(epi.data,group.by = 'RNA_snn_res.0.1',split.by = 'orig.ident',ncol = 4)

markers_genes <- FindAllMarkers(epi.data, logfc.threshold = 0.25, 
                                test.use = "wilcox", 
                                min.pct = 0.1, 
                                only.pos = TRUE, 
                                assay = "RNA")
write.csv(markers_genes,file = './variables/marker_genes_epiSubClustering-v2.csv')
markers_genes = read.csv(file = './variables/marker_genes_epiSubClustering-v2.csv',row.names = 1)
markers_genes = split(markers_genes,markers_genes$cluster)
library(ggplot2)
library(cowplot)
top.mark = markers_genes[['1']][1:15,]
top.mark$gene = factor(top.mark$gene,levels = as.character(top.mark$gene[15:1]))
ggplot(top.mark,mapping = aes(x=avg_log2FC,y=gene))+geom_point(mapping = aes(size = pct.1,color=pct.1),fill='grey')+
  theme_cowplot()+scale_color_gradient2(low='white',high = 'red')
markers <- read.delim("Human_cell_markers.txt")
markers <- markers[markers$speciesType == "Human", ]
markers <- markers[markers$tissueType == "Lung", ]
View(markers)
View(markers_genes)
#View(markers_genes[["3"]])
# remove strange characters etc.
celltype_list <- lapply(unique(markers$cellName), function(x) {
  x <- paste(markers$geneSymbol[markers$cellName == x], sep = ",")
  #x <- gsub("[[]|[]]| |-", ",", x)
  x <- unlist(strsplit(x, split = ","))
  x <- gsub(' ','',x)
  x <- unique(x[!x %in% c("", "NA", "family")])
  x <- casefold(x, upper = T)
})
names(celltype_list) <- unique(markers$cellName)
library(fgsea)
DGE_list = markers_genes
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
  if(nrow(x)>1){
    gene_rank <- setNames(x$avg_log2FC, x$gene)
    fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank,nperm = 500)
  }else{
    fgseaRes = data.frame()
  }
  return(fgseaRes)
})

FeaturePlot(epi.data, reduction = "umap", 
            features = c("SFTPD","ETV5","SFTPA1","MUC1"), 
            order = T, slot = "counts", combine = T)


epi.data<-RenameIdents(epi.data, '0'='AT2','1'='CXCL14+','2'='Ciliated',
                       '3'='T&NK cell','4'='Club','5'='Myeloid cell',
                       '6'='Basal','7'='AT1','8'='AT2')
epi.data <- subset(epi.data, idents = c('T&NK cell','Myeloid cell'), invert = T)

epi.data$subType = Idents(epi.data)
#AT2:SFTPD, ETV5,SFTPA1,MUC1
#Ciliated Cell:CAPS, TPPP3,FOXJ1,TP73,CCDC78,The ciliary membrane-associated proteome revealsactin-binding proteins as key components of cilia
#AT1:AGER,CLIC5,PDPN
#Basal cell: KRT5, KRT14
#Tuft Cell: ASCL2
#Club cell : SCGB1A1,TSPAN8

# alldata$SubClustering = 'Others'
# alldata$SubClustering[colnames(epi.data)] = paste('Epi',epi.data$CCA_snn_res.0.2)
# Idents(alldata)=alldata$SubClustering
# markers_test = FindMarkers(alldata,ident.1 = 'Epi 9',
#                            logfc.threshold = 0.25,
#                            test.use = "wilcox", 
#                            min.pct = 0.1,
#                            only.pos = TRUE,
#                            assay = "RNA")
# 
marker.genes = c("SFTPD","ETV5","SFTPA1","MUC1","KRT5",'KRT17',
                 "AGER","CLIC5","C8orf4",'CAPS','TPPP3','SCGB1A1',"CXCL14","CLDN2",'CEACAM5','CEACAM6','MDK')
#marker.genes = intersect(marker.genes,as.character(markers_genes$gene))
# alldata <- ScaleData(alldata, features = marker.genes, assay = "RNA")
# DoHeatmap(alldata, features = marker.genes, group.by = 'CellTypeManully',assay = "CCA")
# #Dot plot for the features of each clusters
DotPlot(epi.data, features = marker.genes,assay = "RNA",cluster.idents = T) + coord_flip()+ RotatedAxis() 

DotPlot(epi.data, features = marker.genes, group.by = 'subType',assay = "RNA",cluster.idents = T) + 
  coord_flip()+ 
  RotatedAxis() +
  scale_color_gradient2(low='blue',mid = 'white',high = 'red')+
  scale_size(range = c(1,4))+
  theme(text = element_text(size = 8)) 
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/epi/V2/dotplot celltype markers.pdf

DimPlot(epi.data,group.by = 'subType',reduction = 'umap',label = T)
plot_grid(ncol = 3, 
          DimPlot(epi.data, group.by = 'Patient') + NoAxes(), 
          DimPlot(epi.data, group.by = "orig.ident") + NoAxes(), 
          DimPlot(epi.data, group.by = "type") + NoAxes())
FeaturePlot(epi.data, reduction = "umap", 
            features = c("ALDH1A1","PROM1","CD44","ALCAM",'NOTCH1','EPCAM'), 
            order = T, slot = "counts", combine = T)

subtype.cols = pal_npg('nrc')(length(unique(epi.data$subType)))
names(subtype.cols)=c('AT2','CXCL14+','Ciliated',
                      'Club','Basal','AT1')

epi.data$subType = factor(epi.data$subType,levels =c('AT2','CXCL14+','Ciliated',
                                                     'Club','Basal','AT1') )
DimPlot(epi.data,group.by = 'subType',cols = subtype.cols,label = T)
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/epi/V2/Dimplot subtype anno labeled.pdf
epi.data.tumor = epi.data[,epi.data$type == 'TUMOR']
DimPlot(epi.data.tumor,group.by = 'Patient')
library(pals)
epi.data.tumor$orig.ident2 = unlist(lapply(as.character(epi.data.tumor$orig.ident),function(a){
  paste0(substr(a,1,2),strsplit(a,'_',fixed = T)[[1]][3])
}))
color_pelette <- rev(as.vector(kelly()[2:(length(unique(epi.data.tumor$orig.ident2))+1)]))
names(color_pelette) <- unique(epi.data.tumor$orig.ident2)

DimPlot(epi.data.tumor[,epi.data.tumor$subType == 'CXCL14+'],group.by = 'orig.ident2',
        cols = brewer.pal(9,'Set1')[-1],split.by = 'Patient')

DimPlot(epi.data.tumor,group.by = 'subType',split.by = 'Patient')
DimPlot(epi.data.tumor,group.by = 'orig.ident',split.by = 'subType',
        cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))

DimPlot(epi.data[,epi.data$Patient == 'P4'],group.by = 'subType',split.by = 'orig.ident',ncol = 2)

p1 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P1'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p2 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P2'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p3 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P3'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p4 = DimPlot(epi.data.tumor[,epi.data.tumor$Patient == 'P4'],group.by = 'orig.ident',split.by = 'subType',
             cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set3')))
p1/p2/p3/p4




epi.data.normal = epi.data[,epi.data$type != 'TUMOR']

epi.data.normal$orig.ident2 = unlist(lapply(as.character(epi.data.normal$orig.ident),function(a){
  paste0(substr(a,1,2),strsplit(a,'_',fixed = T)[[1]][3])
}))

DimPlot(epi.data.normal,group.by = 'orig.ident2',
        cols = brewer.pal(9,'Set1')[-1],split.by = 'Patient')

####

annotation = epi.data@meta.data

annotation$From = unlist(lapply(as.character(annotation$orig.ident),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
annotation$From = factor(annotation$From ,levels = c('TI','TM','TI1','TI2','TM1','TM2','TS1','TS2','NI','NM','NS'))
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


saveRDS(epi.data,file = './variables/epi.data.V2.Rds')
epi.data = readRDS(file = './variables/epi.data.V2.Rds')
# markers_genes = read.csv(file = './variables/marker_genes_epiSubClustering.csv',row.names = 1)
##############################################################################################################################################
###############################################monocle trajectory analysis,2021-08-17#############################################################

library(monocle)
epi.data.i = epi.data[,epi.data$Patient == 'P4']#Change P1, P2,P3,P4
#all together
#epi.data.i = epi.data
data <- as(as.matrix(epi.data.i@assays$RNA@counts), 'sparseMatrix') 
gc()
pd <-  epi.data.i@meta.data
pd$cellType = pd$subType
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#Construct monocle cds
#Monocle works well with both relative expression data and count-based measures (e.g. UMIs). 
#In general, it works best with transcript count data, especially UMI data. 
#Whatever your data type, it is critical that specify the appropriate distribution for it. 
#FPKM/TPM values are generally log-normally distributed, while UMIs or read counts are better modeled with the negative binomial. 
#To work with count data, specify the negative binomial distribution as the expressionFamily argument to newCellDataSet
cds <- newCellDataSet(data,
                      phenoData = new('AnnotatedDataFrame',data=pd),
                      featureData = new('AnnotatedDataFrame',data=fData),
                      expressionFamily=negbinomial.size())
rm(epi.data.i,data,pd)
gc()
#Pre-process the data
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
print(Sys.time())
#fData = new('AnnotatedDataFrame',data=fData)

###Determine the genes
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

###Cluster cells by monocle################
markers_genes = read.csv(file = 'marker_genes_epiSubClustering.csv',row.names = 1)
marker_diff <- markers_genes[markers_genes$p_val_adj<0.01,]
top100 <- marker_diff %>% group_by(cluster) %>% top_n(-100, p_val_adj)%>% top_n(100, avg_logFC)
semisup_clustering_genes <- unique(top100$gene)
cds <- setOrderingFilter(cds, semisup_clustering_genes)
plot_ordering_genes(cds)
# plot_pc_variance_explained(cds, return_all = F)
# 
# cds <- reduceDimension(cds, max_components = 2, num_dim = 3,
#                         norm_method = 'log',
#                         reduction_method = 'tSNE',
#                         residualModelFormulaStr = "~orig.ident + num_genes_expressed",
#                         verbose = T)
# cds <- clusterCells(cds, num_clusters = 6)
# plot_cell_clusters(cds, 1, 2, color = "cellType")+facet_wrap(~orig.ident,nrow = 4)
# plot_cell_clusters(cds, 1, 2, color = "Patient")
####################end of clustering###############################
# disp_table <- dispersionTable(cds)
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) 
# cds_myo <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
# plot_ordering_genes(cds)

# diff_test_res <- differentialGeneTest(cds[expressed_genes,],
#                                       fullModelFormulaStr = "~orig.ident")
# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# cds <- setOrderingFilter(cds, ordering_genes)
# plot_ordering_genes(cds)


###reduce dimensions
cds <- reduceDimension(cds,max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by="orig.ident", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="State", size=1,show_backbone=TRUE)
plot_cell_trajectory(cds,color_by="cellType", size=1,show_backbone=TRUE)+facet_wrap(~orig.ident,nrow = 5)

# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning", 
#so we often have to call orderCells again using the root_state argument to specify the beginning. 
#
cds <- orderCells(cds, root_state =2)
plot_cell_trajectory(cds, color_by = "Pseudotime")

# BEAM_res=BEAM(cds,branch_point = 1,cores = 6)
# BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
# saveRDS(BEAM_res,file = 'variables/P2_beamRes.Rds')
# tmp1=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval==0))[1:25],],
#                                  branch_point = 1,
#                                  num_clusters = 4, #??§»??????????????group
#                                  cores = 4,
#                                  branch_labels = c("Cell fate 1", "Cell fate 2"),
#                                  #hmcols = NULL, #??????
#                                  hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
#                                  branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2??????????????
#                                  use_gene_short_name = T,
#                                  show_rownames = T,
#                                  return_heatmap = T #??????????§»????????
# )
# tmp1$ph_res
cds@phenoData@data$orig.ident = factor(cds@phenoData@data$orig.ident,
                                       levels = c("TI_R_P1",  "TM_R_P1",  "NI_R_P1",  "NM_R_P1",  "TI_R_P2",  "TM_R_P2",  "NI_R_P2",  "NM_R_P2", 
                                                  "TI1_R_P3", "TM1_R_P3", "NI_R_P3",  "NM_R_P3",  "TI2_R_P3", "TM2_R_P3", "TI_L_P4",  "TS1_L_P4",
                                                  "TS2_L_P4", "NI_L_P4",  "NS_L_P4" ))
p1= ggplot(cds@phenoData@data,aes(x=orig.ident,fill=State))+geom_bar(width = 0.7)+theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(1,'Set1')))

p2= ggplot(cds@phenoData@data,aes(x=subType,fill=State))+
  geom_bar(position = 'fill',width = 0.7)+
  theme_cowplot()+
  facet_wrap(~orig.ident,nrow = 1)+
  scale_fill_manual(values = c(brewer.pal(8,'Set2'),brewer.pal(1,'Set1')))+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
p1
p2
p1/p2
saveRDS(cds,file = 'variables/P1234.cds.Rds')
cds = readRDS(file = 'variables/P4.cds.Rds')

#####################################sample specific markers for each subtypes#############################
subtypes = as.character(unique(epi.data$subType))

subtype.sample.markers = list()

for(i in 1:length(subtypes)){
  subtype = subtypes[i]
  epi.sub = epi.data[,epi.data$subType == subtype]
  Idents(epi.sub)=epi.sub$orig.ident
  epi.sub.P1 = epi.sub[,epi.sub$Patient == 'P1']
  res.markers.P1 = FindAllMarkers(epi.sub.P1,assay = 'RNA',logfc.threshold=-Inf,min.pct = -Inf)
  epi.sub.P2 = epi.sub[,epi.sub$Patient == 'P2']
  res.markers.P2 = FindAllMarkers(epi.sub.P2,assay = 'RNA',logfc.threshold=-Inf,min.pct = -Inf)
  epi.sub.P3 = epi.sub[,epi.sub$Patient == 'P3']
  res.markers.P3 = FindAllMarkers(epi.sub.P3,assay = 'RNA',logfc.threshold=-Inf,min.pct = -Inf)
  epi.sub.P4 = epi.sub[,epi.sub$Patient == 'P4']
  res.markers.P4 = FindAllMarkers(epi.sub.P4,assay = 'RNA',logfc.threshold=-Inf,min.pct = -Inf)
  if(nrow(res.markers.P1)>0){res.markers.P1$Patient = 'P1'}
  if(nrow(res.markers.P2)>0){res.markers.P2$Patient = 'P2'}
  if(nrow(res.markers.P3)>0){res.markers.P3$Patient = 'P3'}
  if(nrow(res.markers.P4)>0){res.markers.P4$Patient = 'P4'}
  
  res.markers = rbind(res.markers.P1,
                      res.markers.P2,
                      res.markers.P3,
                      res.markers.P4)
  res.markers$Subtype = subtype
  subtype.sample.markers[[i]]=res.markers
}

i=7##Change from 1 to 7
res.markers = subtype.sample.markers[[i]]
top5 <- res.markers %>% group_by(cluster) %>% top_n(-5, p_val_adj)  %>% top_n(5, avg_log2FC)
top5[top5$p_val_adj <= 1e-200,'p_val_adj']=1e-200

ggplot(top5,aes(x=cluster,y=gene))+geom_point(aes(size = -log10(p_val_adj),color=avg_log2FC))+
  theme_cowplot()+scale_color_gradient2(low='blue',mid = 'white',high='red')+facet_wrap(Patient~Subtype,scales = 'free',nrow = 1)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='')

##pathway analysis##
#library(org.Hs.eg.db)
library(clusterProfiler)
load("~/scRNA/path.info.RData")
#path.info = path.info[c(-)]
#paths = clusterProfiler::download_KEGG('hsa')
# gene.map=AnnotationDbi::select(org.Hs.eg.db,keys = unique(res.markers$gene),columns = c('SYMBOL','ENTREZID'),keytype = 'SYMBOL')
# gene.map = gene.map[!duplicated(gene.map$SYMBOL),]
# rownames(gene.map)=gene.map$SYMBOL
subtype.sample.paths = data.frame()
for(k in 1:length(subtypes)){
  res.markers = subtype.sample.markers[[k]]
  res.markers = split(res.markers,res.markers$cluster)
  
  for(i in 1:length(res.markers)){
    sample.res = res.markers[[i]]
    FC=sample.res$avg_log2FC
    names(FC)=as.character(sample.res$gene)
    FC = FC[!is.na(names(FC))]
    FC = FC[order(FC,decreasing = T)]
    sample.path = GSEA(FC,minGSSize = 5,pvalueCutoff = 1.2,seed = 111,TERM2GENE = path.info)
    path.res = sample.path@result
    if(nrow(path.res)>0){
      path.res$Patient = as.character(unique(sample.res$Patient))
      path.res$Subtype = as.character(unique(sample.res$Subtype))
      path.res$Sample = as.character(unique(sample.res$cluster))
      subtype.sample.paths = rbind(subtype.sample.paths,path.res)
    }
    
  }
  
}
library(reshape2)
library(pheatmap)

i=7
res.paths = subtype.sample.paths[subtype.sample.paths$Subtype == subtypes[i],]
top5 <- res.paths %>% group_by(Sample) %>% top_n(-5, p.adjust)  %>% top_n(5,abs(NES) )
top5[top5$p.adjust <= 1e-5,'p.adjust']=1e-5
top5 = res.paths[res.paths$ID %in% top5$ID,]
top5 = top5[top5$p.adjust<0.1,]

ggplot(top5,aes(x=Sample,y=ID))+geom_point(aes(size = -log10(p.adjust),color=NES))+
  theme_cowplot()+scale_color_gradient2(low='blue',mid = 'white',high='red')+facet_wrap(Patient~Subtype,scales = 'free',nrow = 1)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  labs(x='',y='')
mat = acast(res.paths,ID ~ Sample,value.var = 'NES')
mat[is.na(mat)]=0

pheatmap(mat[as.character(unique(top5$ID)),],show_rownames = T,
         clustering_distance_cols='correlation')
sel.cols = c('TM_R_P1','TI_R_P1','NM_R_P1','NI_R_P1',
             'TM_R_P2','TI_R_P2','NM_R_P2','NI_R_P2',
             'TM1_R_P3','TM2_R_P3','TI1_R_P3','TI2_R_P3','NM_R_P3','NI_R_P3',
             'TS1_L_P4','TS2_L_P4','TI_L_P4','NS_L_P4','NI_L_P4')
pheatmap(mat[as.character(unique(top5$ID)),intersect(sel.cols,colnames(mat))],show_rownames = T,cluster_cols = F,
         clustering_distance_cols='correlation',gaps_col = c(4))

######PCA######

top5 <- subtype.sample.paths %>% group_by(Sample,Subtype) %>% top_n(-5, p.adjust)  %>% top_n(5,abs(NES) )
top5[top5$p.adjust <= 1e-5,'p.adjust']=1e-5
top5 = subtype.sample.paths[subtype.sample.paths$ID %in% top5$ID,]
top5 = top5[top5$p.adjust<0.1,]

mat = acast(subtype.sample.paths,ID ~ Sample+Subtype,value.var = 'NES')
mat[is.na(mat)]=0
mat = scale(t(mat))
mat = mat[apply(mat, 1, sd)!=0,]
mat[is.na(mat)]
pcr.int = prcomp(mat)
library(RColorBrewer)
library(ggrepel)
pcr.cluster = data.frame(Patient = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][3]
})),
Subtype = unlist(lapply(rownames(mat), function(a){
  strsplit(a,'_')[[1]][4]
})),
Orig.ident = unlist(lapply(rownames(mat), function(a){
  paste0(strsplit(a,'_')[[1]][c(1,2,3)],collapse = '_')
})),
Type = unlist(lapply(rownames(mat), function(a){
  substr(a,1,1)
})),
row.names  = rownames(mat),
PC1 = pcr.int$x[,1],
PC2 = pcr.int$x[,2])
ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=Type),size = 3)+
  scale_color_manual(values = c(brewer.pal(8,'Set1'),brewer.pal(9,'Set3')[4:9]))+
  theme_cowplot()+facet_wrap(vars(Subtype),scales = 'free')+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggplot(pcr.cluster,mapping = aes(x=PC1,y=PC2))+
  geom_point(mapping = aes(color=Patient),size = 3)+
  scale_color_manual(values = c(brewer.pal(8,'Set1')))+
  theme_cowplot()+facet_wrap(vars(Subtype),scales = 'free',nrow = 2)+
  geom_text_repel(aes(label=Orig.ident,color=Patient),size=2,nudge_x = .3,nudge_y = .3,max.overlaps = 100)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))


####PCA analysis####
epi.data = readRDS(file = './variables/epi.data.V2.Rds')
com2orig = table(epi.data$orig.ident,epi.data$subType)
com2orig = com2orig/apply(com2orig, 1, sum)
mat = t(com2orig)
mat = mat[apply(mat, 1, sum)>0,]
mat = scale(t(mat))
mat = mat[apply(mat, 1, sd)!=0,]
mat[is.na(mat)]

pcr.int = prcomp(mat)
pcr.int$sdev/sum(pcr.int$sdev)#[1] 2.324118e-01 2.267851e-01 1.662813e-01 1.512029e-01 1.165619e-01 1.067569e-01 4.133088e-17
library(RColorBrewer)
library(ggrepel)
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



####epi.data cxcl14+ subtype, 2021-09-26####
epi.data = readRDS(file = './variables/epi.data.V2.Rds')
cxcl14.markers = c('CXCL14','CLDN2','CEACAM5','CEACAM6','MDK','TMSB4X')
VlnPlot(epi.data, features = cxcl14.markers, 
        group.by = "RNA_snn_res.0.1", 
        ncol = 2, pt.size = 0)
FeaturePlot(epi.data, reduction = "umap", 
            features = cxcl14.markers, 
            order = T, slot = "counts", combine = T)

####Tumor epi data analysis, 2021-10-08####

epi.data = readRDS(file = './variables/epi.data.V2.Rds')
#epi.data.tumor =  epi.data[,epi.data$subType == 'CXCL14+']
#epi.data.tumor =  epi.data.tumor[,epi.data.tumor$type == 'TUMOR']
epi.data.tumor =  epi.data[,epi.data$type == 'TUMOR']
epi.data.tumor = NormalizeData(epi.data.tumor)
epi.data.tumor = FindVariableFeatures(epi.data.tumor,nfeatures = 500)
epi.data.tumor = ScaleData(epi.data.tumor, vars.to.regress = c("nFeature_RNA", "percent_mito"))
epi.data.tumor = RunPCA(epi.data.tumor, verbose = T, npcs = 20)
#epi.data.tumor = RunUMAP(epi.data.tumor, dims = 1:50, verbose = T)
epi.data.tumor = RunHarmony(epi.data.tumor,c( "orig.ident" ))

epi.data.tumor <- RunUMAP(epi.data.tumor,reduction = "harmony", dims = 1:20)
epi.data.tumor <- RunTSNE(epi.data.tumor,reduction = "harmony", dims = 1:20)

DimPlot(epi.data.tumor, reduction = "umap", group.by = 'orig.ident')
epi.data.tumor$location = unlist(lapply(epi.data.tumor$orig.ident,function(a){
  substr(a,1,2)
}))
epi.data.tumor$facet = paste(epi.data.tumor$Patient,epi.data.tumor$subType)
DimPlot(epi.data.tumor, reduction = "tsne", group.by = 'location',split.by = 'Patient',
        cols = c(brewer.pal(8,'Set1')))
DimPlot(epi.data.tumor, reduction = "umap", 
        group.by = 'location',split.by = 'facet',ncol = 4,
        cols = c(brewer.pal(8,'Set1'),brewer.pal(8,'Set2')))

DimPlot(epi.data.tumor, reduction = "umap", group.by = 'Patient')
DimPlot(epi.data.tumor, reduction = "umap", group.by = 'subType')
DimPlot(epi.data.tumor, reduction = "tsne", group.by = 'subType')
DimPlot(epi.data.tumor, reduction = "tsne", group.by = 'Patient')


epi.data.tumor.i = epi.data.tumor[,epi.data.tumor$Patient =="P3"]
epi.data.tumor.i = NormalizeData(epi.data.tumor.i)
epi.data.tumor.i = FindVariableFeatures(epi.data.tumor.i)
epi.data.tumor.i = ScaleData(epi.data.tumor.i, vars.to.regress = c("nFeature_RNA", "percent_mito"))
epi.data.tumor.i = RunPCA(epi.data.tumor.i, verbose = T, npcs = 20)#npcs = 50 for P1,2,4, npcs = 20 for P3
#epi.data.tumor = RunUMAP(epi.data.tumor, dims = 1:50, verbose = T)
epi.data.tumor.i = RunHarmony(epi.data.tumor.i,c( "orig.ident" ))

epi.data.tumor.i <- RunUMAP(epi.data.tumor.i,reduction = "harmony", dims = 1:10)

DimPlot(epi.data.tumor.i, reduction = "umap", group.by = 'orig.ident')


####check int genes####
VlnPlot(epi.data,group.by = 'subType',
        cols = pal_npg('nrc')(10),
        features = c("ALOX15B", "LCN2", "PRSS12" , "WIF1"),sort = 'counts',
        ncol = 1,
        pt.size = 0)


com2orig = table(epi.data$orig.ident,epi.data$subType)
com2orig = com2orig/apply(com2orig, 1, sum)
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
diff.ratios$CellType = factor(diff.ratios$CellType,levels = rev(c('CXCL14+','AT2','Ciliated',
                                                                  'Club','AT1','Basal')))
library(ggpubr)
                                                           
ptest <- diff.ratios %>% group_by(CellType) %>% summarize(p.value = kruskal.test(Diff ~Type)$p.value)
ggbarplot(diff.ratios,
          x='CellType',y = 'Diff',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  theme_cowplot() + 
  ylab("Diff Ratio") + 
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))+
  coord_flip()
###file:///H:/project/single cell/MPSC&ST/figure results/four samples/epi/V2/diff cell composition/Differetial ratios of T & N samples barplot.pdf

library(reshape2)
ratio.mat = melt(com2orig,measure.vars = colnames(com2orig))
ratio.mat$Type = substr(ratio.mat$Var1,1,1)
ptest2 <- ratio.mat %>% group_by(Var2) %>% summarize(p.value = kruskal.test(value ~Type)$p.value)
# CXCL14+
#   0.01663868
# 2
# Club
# 0.09864761
# 3
# Ciliated
# 0.11667745
# 4
# AT2
# 0.24767627
# 5
# Basal
# 0.36372233
# 6
# AT1
# 0.40896134
#write.csv(ptest2,file = './variables/ptest for immu subtype diff T&N.csv')
ratio.mat$Var2 = factor(ratio.mat$Var2,levels = rev(c('CXCL14+','AT2','Ciliated',
                                                                  'Club','AT1','Basal')))

ggbarplot(ratio.mat,
          x='Var2',y = 'value',
          add = 'mean_sd',fill='Type',
          position=position_dodge(0.7)) + 
  geom_text(data = ptest2,mapping = aes(x=Var2,y=0.7,label = ifelse(signif(p.value,2)<0.05,'*','')),inherit.aes = F,size=3)+
  coord_flip()+
  theme_cowplot() + 
  ylab("Proportion") + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_fill_manual(values = c('blue','red'))
##file:///H:/project/single cell/MPSC&ST/figure results/four samples/epi/V2/com2orig TVSN compare barplot.pdf
