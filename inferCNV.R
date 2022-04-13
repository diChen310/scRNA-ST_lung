####Infer cnv based on the scRNA data####
#Prepare gene order file
library(infercnv)
library(Seurat)
setwd('H:/Project/single cell/MPSC&ST/')

alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
epi.data = readRDS(file = './variables/epi.data.Rds')

# genes.loc = read.delim(file='H:/project/single cell/documents/infercnv/gencode_v21_gen_pos.complete.txt',
#                        header = F)
# genes = unlist(lapply(genes.loc$V1, function(a){
#   strsplit(a,'\\|')[[1]][1]
# }))
# genes.loc$V1 = genes
# genes.loc = genes.loc[!duplicated(genes.loc$V1),]
# write.table(genes.loc,file = './documents/inferCNV/genes location.txt',row.names = F,sep='\t',quote = F)

lung.p1 = alldata[, alldata$orig.ident %in% c('TI_R_P2')]
epi.p1 = epi.data[, epi.data$orig.ident %in% c('TI_R_P2')]
#lung.p1 = alldata[, alldata$orig.ident %in% c('TI_P2','NI_P2')]
lung.p1 = lung.p1[,lung.p1$CellTypeManully %in% c('Epithelial cell','T&NK cell')]
lung.p1@active.assay = 'RNA'
selected_c <- WhichCells(lung.p1, expression = nFeature_RNA > 200)
selected_f <- rownames(lung.p1)[Matrix::rowSums(lung.p1) > 3]
lung.p1 <- subset(lung.p1, features = selected_f, cells = selected_c)
epi.p1.samples = intersect(colnames(epi.p1),colnames(lung.p1))
length(epi.p1.samples)
ncol(epi.p1)
#
annotation = lung.p1@meta.data
annotation$CellTypeInput = as.character(annotation$CellTypeManully)
annotation$CellName = rownames(annotation)
annotation[epi.p1.samples,'CellTypeInput']=as.character(epi.p1[,epi.p1.samples]$subType)
annotation = annotation[annotation$CellTypeInput != 'Epithelial cell',]

counts_matrix = as.matrix(lung.p1@assays$RNA@counts[,rownames(annotation)])

write.table(annotation[,c('CellName','CellTypeInput')],
            file = './documents/inferCNV/res2/P2_TI_anno.txt',row.names = F,sep='\t',quote = F,col.names = F)
#nm.type = unique(as.character(annotation$CellTypeInput))[grepl('NM_P',unique(as.character(annotation$CellTypeInput)))]
#ni.type = unique(as.character(annotation$CellTypeInput))[grepl('NI_P',unique(as.character(annotation$CellTypeInput)))]
ref.type = unique(as.character(annotation$CellTypeInput))[grepl('T&NK cell',unique(as.character(annotation$CellTypeInput)))]

infercnv_ti = CreateInfercnvObject(raw_counts_matrix = counts_matrix[,annotation$CellName],
                                   gene_order_file = './documents/inferCNV/genes location.txt',
                                   annotations_file = './documents/inferCNV/res2/P2_TI_anno.txt',
                                   ref_group_names = ref.type)
infercnv_obj = infercnv::run(infercnv_ti,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="./documents/inferCNV/res3/p2_TI_refT/",  
                             cluster_by_groups=T,   # 聚类
                             denoise=T, #去噪
                             HMM=T) # 是否基于HMM预测CNV


