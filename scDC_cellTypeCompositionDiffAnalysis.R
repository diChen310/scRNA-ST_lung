####Differential cell type composition analysis,2021-09-09####
####Based on R version 4, 32-bit, for all cells####
# devtools::install_github("SydneyBioX/scDC",force=T)

setwd('H:/project/single cell/MPSC&ST/')
library(scDC)
library(Seurat)
alldata = readRDS(file = './variables/data.filt_harmony_annotated-v2.Rds')
#data.i = alldata[,alldata$Patient == 'P1']

res_BCa = scDC_noClustering(alldata$CellTypeManully, 
alldata$orig.ident, 
calCI = TRUE, 
calCI_method = "BCa",
nboot = 1000)
saveRDS(res_BCa,file='./variables/scDCRes.Rds')
res_BCa = readRDS(file='./variables/scDCRes.Rds')
library(ggplot2)
library(ggpubr)
library(ggsci)
df_toPlot <- res_BCa$results
df_toPlot$median <- apply(res_BCa$thetastar, 1, median)
df_toPlot$cond <- unlist(lapply(as.character(df_toPlot$subject),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
n_celltype = length(unique(df_toPlot$cellTypes))
ggplot(df_toPlot, aes(x = subject, y = median, fill = cond)) + 
  geom_bar(stat = "identity",position = "dodge", alpha = 0.8) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  geom_errorbar(aes(ymin = conf_low,ymax = conf_high), color = 'black', width = 0.1, lwd = 1, 
                                                       position = position_dodge(width = 0.5)) + 
  facet_wrap(cond~cellTypes, 
             ncol = n_celltype, nrow=4,scales='free',
             labeller = labeller(cellTypes = label_wrap_gen(width = 10, multi_line = TRUE))) + 
  coord_flip() + 
  ylim(c(0,1)) + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
                                                                                                                                                     
# barplotCI(res_BCa, condition = unlist(lapply(as.character(res_BCa$results$subject),function(a){
#   return(strsplit(a,'_')[[1]][3])
# })))+facet_wrap(cond~cellTypes,nrow=4,scales='free')

res_BCa$thetastar[1:5,1:10]
colnames(res_BCa$thetastar)=paste('B',1:ncol(res_BCa$thetastar))
library(reshape2)
df_box = cbind(res_BCa$info,res_BCa$thetastar)
df_box$cond <- unlist(lapply(as.character(df_box$subject),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
df_box = melt(df_box,measure.vars =colnames(res_BCa$thetastar) )
df_box = df_box[,c(-3)]
library(dplyr)
df_box.i = filter(df_box,cond == 'P1',cellTypes == 'B cell')
ggplot(df_box.i, aes(x=subject,y = value, fill = cond)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  coord_flip() + 
  ylim(c(0,1)) + stat_compare_means(aes(label = ..p.signif..),ref.group = 'TM_R_P1',label.y = 0.8)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
df_box$facet = paste(df_box$cond,df_box$cellTypes)
df_box$cellTypes = unlist(lapply(as.character(df_box$cellTypes),function(a){
  return(strsplit(a,' ')[[1]][1])
}))
df_box$cellTypes = factor(df_box$cellTypes,levels =c('Epithelial','Endothelial','Fibroblast',
                                                     'T&NK','B','Myeloid','Mast'))
cellType.cols = pal_npg('nrc')(length(unique(df_box$cellTypes)))
names(cellType.cols)=c('Epithelial','Endothelial','Fibroblast',
                       'T&NK','B','Myeloid','Mast')
ggbarplot(df_box,x='subject',y = 'value', fill = 'cellTypes',add = 'mean_sd',
          facet.by = c('cond','cellTypes'),scales = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = cellType.cols) + ylab("Proportion") + 
  coord_flip()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
#file:///H:/project/single cell/MPSC&ST/figure results/four samples/diff/diff cell compostion/cell-type composition diff barplot.pdf           

library(dplyr)
df_box$Type = substr(df_box$subject,1,1)

ptest <- df_box %>% filter(Type == 'T') %>% group_by(cond,cellTypes) %>% summarize(p.value = kruskal.test(value ~subject)$p.value)
ptest.n <- df_box %>% filter(Type == 'N') %>% group_by(cond,cellTypes) %>% summarize(p.value = kruskal.test(value ~subject)$p.value)

ggbarplot(df_box[df_box$Type == 'T',],x='subject',y = 'value', fill = 'cond',add = 'mean_sd',facet.by = c('cond','cellTypes'),scales='free_y') + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  coord_flip() + 
  ylim(c(0,1)) + stat_compare_means(aes(label = ..p.signif..),label.y = 0.8)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggbarplot(df_box[df_box$Type == 'N',],x='subject',y = 'value', fill = 'cond',add = 'mean_sd',facet.by = c('cond','cellTypes'),scales='free_y') + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  coord_flip() + 
  ylim(c(0,1)) + stat_compare_means(aes(label = ..p.signif..),label.x = 0.8)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

####Based on R version 4, for epi cells, 2021-09-13,for B cells,2021-09-27 ####
#### t&NK,2021-09-29
setwd('H:/project/single cell/MPSC&ST/')
#install.packages('ModelMetrics')
library(scDC)
library(Seurat)
#alldata = readRDS(file = './variables/epi.data.V2.Rds')#for epi cells

alldata = readRDS(file = './variables/tcell.data.Rds')

res_BCa = scDC_noClustering(alldata$subType, 
alldata$orig.ident, 
calCI = TRUE, 
calCI_method = "BCa",
nboot = 1000)
saveRDS(res_BCa,file='./variables/scDCRes_tCells.Rds')
res_BCa = readRDS(file='./variables/scDCRes_endoCells.Rds')
library(ggplot2)
library(ggpubr)

df_toPlot <- res_BCa$results
df_toPlot$median <- apply(res_BCa$thetastar, 1, median)
df_toPlot$cond <- unlist(lapply(as.character(df_toPlot$subject),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
n_celltype = length(unique(df_toPlot$cellTypes))
ggplot(df_toPlot, aes(x = subject, y = median, fill = cond)) + 
  geom_bar(stat = "identity",position = "dodge", alpha = 0.8) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  geom_errorbar(aes(ymin = conf_low,ymax = conf_high), color = 'black', width = 0.1, lwd = 1, 
                                                       position = position_dodge(width = 0.5)) + 
  facet_wrap(cond~cellTypes, 
             ncol = n_celltype, nrow=4,scales='free',
             labeller = labeller(cellTypes = label_wrap_gen(width = 10, multi_line = TRUE))) + 
  coord_flip() + 
  ylim(c(0,1)) + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
           
                                                                                                                                                       
# barplotCI(res_BCa, condition = unlist(lapply(as.character(res_BCa$results$subject),function(a){
#   return(strsplit(a,'_')[[1]][3])
# })))+facet_wrap(cond~cellTypes,nrow=4,scales='free')

res_BCa$thetastar[1:5,1:10]
colnames(res_BCa$thetastar)=paste('B',1:ncol(res_BCa$thetastar))
library(reshape2)
df_box = cbind(res_BCa$info,res_BCa$thetastar)
df_box$cond <- unlist(lapply(as.character(df_box$subject),function(a){
  return(strsplit(a,'_')[[1]][3])
}))
df_box = melt(df_box,measure.vars =colnames(res_BCa$thetastar) )
df_box = df_box[,c(-3)]

df_box$facet = paste(df_box$cond,df_box$cellTypes)
df_box$subject = unlist(lapply(as.character(df_box$subject),function(a){
  return(strsplit(a,'_')[[1]][1])
}))
subtype.cols=pal_npg('nrc')(8)
names(subtype.cols)=rev(names(sort(table(fib.data$subType))))#manully
df_box$cellTypes = factor(df_box$cellTypes,levels = rev(names(sort(table(endo.data$subType)))))
ggbarplot(df_box,x='subject',y = 'value', fill = 'cellTypes',add = 'mean_sd',facet.by = c('cond','cellTypes'),scales='free') + 
  theme_bw() + 
  scale_fill_manual(values=subtype.cols) + ylab("Proportion") + 
  coord_flip() + 
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
#epi: file:///H:/project/single cell/MPSC&ST/figure results/four samples/epi/V2/diff cell composition/cell composition ggbar.pdf
#endo:file:///H:/project/single cell/MPSC&ST/figure results/four samples/endo/cell composition ggbar.pdf
library(dplyr)
df_box$Type = substr(df_box$subject,1,1)

ptest <- df_box %>% filter(Type == 'T') %>% group_by(cond,cellTypes) %>% summarize(p.value = kruskal.test(value ~subject)$p.value)
ptest.n <- df_box %>% filter(Type == 'N') %>% group_by(cond,cellTypes) %>% summarize(p.value = kruskal.test(value ~subject)$p.value)

ggbarplot(df_box[df_box$Type == 'T',],x='subject',y = 'value', fill = 'cond',add = 'mean_sd',facet.by = c('cond','cellTypes'),scales='free_y') + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  coord_flip() + 
  ylim(c(0,1)) + stat_compare_means(aes(label = ..p.signif..),label.x = 0.8)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

ggbarplot(df_box[df_box$Type == 'N',],x='subject',y = 'value', fill = 'cond',add = 'mean_sd',facet.by = c('cond','cellTypes'),scales='free_y') + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + ylab("Proportion") + 
  coord_flip() + 
  ylim(c(0,1)) + stat_compare_means(aes(label = ..p.signif..),label.x = 0.8)+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))
