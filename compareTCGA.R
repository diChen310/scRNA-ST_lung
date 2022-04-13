##########Interesting genes comparing to TCGA, 2019-09-16####################
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ComplexHeatmap)
library(dplyr)
setwd('H:/project/single cell/MPSC&ST/')
patients.markers = readRDS(file='variables/patients.diff.markers.Rds')
res.markers = Reduce(rbind,patients.markers)
res.markers$tag = paste(res.markers$Gene,res.markers$patient)
res.markers.TI = res.markers[res.markers$Compare == 'TIVSNI',]
res.markers.TOthers = res.markers[res.markers$Compare != 'TIVSNI',]
View(res.markers.TI)
View(res.markers.TOthers)
tumor = 'LUAD' # LUSC or LUAD
exp.tcga<-read.delim(file=paste0("E:/data/TCGA/FRData/",tumor,'ExpressionNormalized.txt'),sep='\t',check.names = FALSE)
samples.tp = colnames(exp.tcga)[substr(colnames(exp.tcga),14,15) == '01']
samples.nt = colnames(exp.tcga)[substr(colnames(exp.tcga),14,15) == '11']

samples.nt.short = substr(samples.nt,1,12)
samples.tp.short = substr(samples.tp,1,12)
paired.tp = samples.tp[samples.tp.short %in% samples.nt.short]
paired.nt = samples.nt[samples.nt.short %in% samples.tp.short]

paired.tp = paired.tp[order(substr(paired.tp,1,12))]
paired.nt = paired.nt[order(substr(paired.nt,1,12))]
substr(paired.nt,1,12) == substr(paired.tp,1,12)
exp.m = log2(exp.tcga+1)
exp.m = exp.m[apply(exp.m, 1, sd)!=0,]
exp.m = exp.m[apply(exp.m, 1, function(a){length(a[a==0])/length(a)})<0.6,]
diff.genes = apply(exp.m[,c(paired.tp,paired.nt)],1,function(a){wilcox.test(as.double(a[1:(length(a)/2)]),as.double(a[((length(a)/2)+1):length(a)]))$p.value})
fold.genes = apply(exp.m[,c(paired.tp,paired.nt)],1,function(a){(mean(as.double(a[1:(length(a)/2)])))-(mean(as.double(a[((length(a)/2)+1):length(a)])))})
tumor.diff = data.frame(PValue = diff.genes,LogFC=fold.genes,row.names = rownames(exp.m))
tumor.diff$FDR=p.adjust(tumor.diff$PValue,'fdr')
tumor.diff.genes = rownames(tumor.diff[tumor.diff$FDR<0.01 & abs(tumor.diff$LogFC)>1,])

p2.ti.sig.genes = filter(res.markers.TI, p_val_adj<0.05 ) %>% group_by(patient,cellType) %>% summarise(N=n())
p2.ti.sig.genes.in = filter(res.markers.TI, p_val_adj<0.05 & Gene %in% tumor.diff.genes) %>% group_by(patient,cellType) %>% summarise(N=n())
ratio= p2.ti.sig.genes.in[,3]/p2.ti.sig.genes[,3]
names(ratio)='ratio'
p2.ti.sig.genes[,4]=ratio

p2.tO.sig.genes = filter(res.markers.TOthers, p_val_adj<0.05 ) %>% group_by(patient,cellType) %>% summarise(N=n())
p2.tO.sig.genes.in = filter(res.markers.TOthers, p_val_adj<0.05 & Gene %in% tumor.diff.genes) %>% group_by(patient,cellType) %>% summarise(N=n())
ratio = p2.tO.sig.genes.in[,3]/p2.tO.sig.genes[,3]
names(ratio)='ratio'
p2.tO.sig.genes[,4]=ratio

p2.ti.sig.genes = data.frame(p2.ti.sig.genes)
p2.ti.sig.genes$From = 'TIVSNI'
p2.tO.sig.genes = data.frame(p2.tO.sig.genes)
p2.tO.sig.genes$From = 'TOVSNO'
colnames(p2.tO.sig.genes)=colnames(p2.ti.sig.genes)
sig.sum = rbind(p2.ti.sig.genes,p2.tO.sig.genes)
ggplot(sig.sum,
       aes(x=From,y=ratio,fill=cellType))+geom_bar(width = 0.6,stat = 'identity')+
  facet_grid(patient~cellType)+
  theme_cowplot()+scale_fill_manual(values=brewer.pal(7,'Set2'))+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='ratio')
ggplot(sig.sum,
       aes(x=From,y=N,fill=cellType))+geom_bar(width = 0.6,stat = 'identity')+
  facet_grid(patient~cellType)+
  theme_cowplot()+scale_fill_manual(values=brewer.pal(7,'Set2'))+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='ratio')+scale_y_log10()
##################This part help analyze the differential expressions between tumor and matched normal samples########
#int.genes = c('RGS1','CD82','IGHM','HSPA1A','DSTN','ANXA1','MT1X','NNMT')
#int.genes = c("ALOX15B", "CD79B" ,  "ITSN2"  , "LCN2"   , "PRSS12" , "TNFRSF18" ,"WIF1"    )
int.genes = c('CXCL14','CLDN2','CEACAM5', 'CEACAM6', 'MDK')
#int.genes = as.character(all.SR3$Gene)
int.genes = intersect(rownames(exp.m),int.genes)
exp.int = t(exp.m[int.genes,c(paired.tp,paired.nt)])

exp.int.m = melt(exp.int,measure.vars = int.genes)
colnames(exp.int.m)=c('Sample','Gene','log2mRNA')
exp.int.m$Type = ifelse(exp.int.m$Sample %in% paired.tp,'T','N')
test.res = exp.int.m %>% group_by(Gene) %>% summarize(p.value = kruskal.test(log2mRNA ~Type)$p.value)
test.res$p.value = as.character(signif(test.res$p.value,2))
ggplot(exp.int.m,mapping=aes(x=Type,y=log2mRNA,fill = Type))+
  facet_wrap(vars(Gene),nrow =1,scales = 'free')+
  geom_boxplot()+
  geom_text(data = test.res,mapping = aes(x=1.5,y=16.3,label = p.value),size=3,inherit.aes = F)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='log2mRNA')+
  scale_fill_manual(values = c('blue','red'))

#########################This part help find clustering of the input samples and TCGA samples#################
# 
# logFC.TCGA = matrix(0,nrow = nrow(exp.m),ncol=length(paired.tp),
#                     dimnames = list(rownames(exp.m),
#                                     substr(paired.tp,1,12)))
# for(i in 1:ncol(logFC.TCGA)){
#   logFC.TCGA[,i]= as.double(exp.m[,paired.tp[i]])-as.double(exp.m[,paired.nt[i]])
# }
# 
# input.sc = matrix(0,nrow=length(unique(res.markers.TI$Gene)),ncol = 8,
#                   dimnames = list(unique(res.markers.TI$Gene),
#                                   c('TIP1','TIP2','TIP3','TIP4',
#                                     'TMP1','TMP2','TMP3','TSP4')))
# input.sc[,1]= res.markers.TI[res.markers.TI$patient == 'P1',][rownames(input.sc),'avg_logFC']
# input.sc[,2]= res.markers.TI[res.markers.TI$patient == 'P2',][paste0(rownames(input.sc),2),'avg_logFC']
# input.sc[,3]= res.markers.TI[res.markers.TI$patient == 'P3',][paste0(rownames(input.sc),3),'avg_logFC']
# input.sc[,4]= res.markers.TI[res.markers.TI$patient == 'P4',][paste0(rownames(input.sc),4),'avg_logFC']
# 
# input.sc[,5]= res.markers.TOthers[res.markers.TOthers$patient == 'P1',][paste0(rownames(input.sc),1),'avg_logFC']
# input.sc[,6]= res.markers.TOthers[res.markers.TOthers$patient == 'P2',][paste0(rownames(input.sc),11),'avg_logFC']
# input.sc[,7]= res.markers.TOthers[res.markers.TOthers$patient == 'P3',][paste0(rownames(input.sc),12),'avg_logFC']
# input.sc[,8]= res.markers.TOthers[res.markers.TOthers$patient == 'P4',][paste0(rownames(input.sc),13),'avg_logFC']
# 
# 
# share.genes = intersect(rownames(logFC.TCGA),rownames(input.sc))
# 
# both.mat = cbind(logFC.TCGA[share.genes,],
#                  input.sc[share.genes,])
# ha = HeatmapAnnotation(df = data.frame(row.names  = c(colnames(both.mat)),
#                                        Source=c(rep('TCGA',ncol(logFC.TCGA)),
#                                                 rep('Input',ncol(input.sc)))),
#                        col=list(Source = c('TCGA'='blue','Input'='red')))
# Heatmap(both.mat,top_annotation = ha,show_row_names = F)
# 
# res.sig.TI = res.markers.TI[res.markers.TI$p_val_adj<0.01,]
# res.sig.TO = res.markers.TOthers[res.markers.TOthers$p_val_adj<0.01,]
# sig.genes = union(unique(res.sig.TI$Gene),unique(res.sig.TO$Gene))
# input.sc.short = input.sc[sig.genes,]
# 
# Heatmap(input.sc.short,show_row_names = F)

##############################This part try to find prognosis relevances###################################
library(survival)
library(survminer)
surv_coxph <- function(gene,exp.tumor,clin){
  if(length(exp.tumor[gene,][exp.tumor[gene,]==0])/ncol(exp.tumor)<0.6){
    clin <- clin[rownames(clin) %in% colnames(exp.tumor),]
    notDead <- is.na(clin$days_to_death)
    if (length(notDead) > 0) {
      clin[notDead,]$days_to_death <-
        clin[notDead,]$days_to_last_follow_up
    }
    one.data <- data.frame(overall_time = as.numeric(as.character(clin$days_to_death)),status = ifelse(clin$vital_status == 'dead',1,0),mRNA = as.double(exp.tumor[gene,rownames(clin)]))
    
    cox.res <- coxph(Surv(overall_time, status) ~ mRNA,data=one.data)
    cox.sum <- summary(cox.res)
    p <- cox.sum$logtest['pvalue']
    HR <- cox.sum$coefficients[1,'exp(coef)']
    
  }else{
    p <- 1
    HR <-'NA'
  }
  return(list(p,HR))
}
clin<-read.delim(file = paste0("E:/data/TCGA/GDCData/",tumor,"_clin.txt"),sep='\t')
print(paste(nrow(clin[clin$vital_status=='dead',]),tumor))
rownames(clin)=clin$bcr_patient_barcode

clin.sub = clin
tumor.patients <- substr(samples.tp,1,12)
clin.sub <- clin.sub[rownames(clin.sub) %in% tumor.patients,]
exp.tumor <- exp.m[,samples.tp]
colnames(exp.tumor)<-substr(colnames(exp.tumor),1,12)
p.values <- c()
HRs <- c()
for(gene in rownames(exp.tumor)){
  
  surv.res <- surv_coxph(gene,exp.tumor,clin.sub)
  p.values <- append(p.values,surv.res[[1]])
  HRs <- append(HRs, surv.res[[2]])
  
  
}

write.table(data.frame(p.values,HRs,row.names = rownames(exp.tumor)),file=paste0('variables/',tumor,'_allGene_coxph.txt'),sep='\t')

#surv.res = read.delim(file=paste0('../MPSC&ST/variables/',tumor,'_allGene_coxph.txt'))

surv.res = read.delim(file=paste0('variables/',tumor,'_allGene_coxph.txt'))
fav.genes = rownames(surv.res[surv.res$p.values<0.05 & surv.res$HRs<1,])
unfav.genes =  rownames(surv.res[surv.res$p.values<0.05 & surv.res$HRs>1,])
intersect(int.genes,fav.genes)
intersect(int.genes,unfav.genes)
gN = length(unique(res.markers.TI$Gene))
res.markers.TI = res.markers.TI[res.markers.TI$avg_logFC<0,]
res.markers.TOthers = res.markers.TOthers[res.markers.TOthers$avg_logFC<0,]

p2.ti.sig.genes = as.data.frame(filter(res.markers.TI, p_val_adj<0.05 ) %>% group_by(patient,cellType) %>% summarise(N=n()))
fav.ratios = unlist( lapply(1:nrow(p2.ti.sig.genes), function(i){
  pi = p2.ti.sig.genes[i,'patient']
  celli = p2.ti.sig.genes[i,'cellType']
  Ni = p2.ti.sig.genes[i,'N']
  genes.i = res.markers.TI[res.markers.TI$patient==pi & res.markers.TI$cellType == celli & res.markers.TI$p_val_adj<0.05,'Gene']
  
  p= fisher.test(matrix(c(length(intersect(genes.i,fav.genes)),length(genes.i)-length(intersect(genes.i,fav.genes)),
                      length(fav.genes)-length(intersect(genes.i,fav.genes)),
                      gN-length(fav.genes)-length(genes.i)+length(intersect(genes.i,fav.genes))),nrow = 2))$p.value
  return(p)
}))
unfav.ratios = unlist( lapply(1:nrow(p2.ti.sig.genes), function(i){
  pi = p2.ti.sig.genes[i,'patient']
  celli = p2.ti.sig.genes[i,'cellType']
  Ni = p2.ti.sig.genes[i,'N']
  genes.i = res.markers.TI[res.markers.TI$patient==pi & res.markers.TI$cellType == celli & res.markers.TI$p_val_adj<0.05,'Gene']
  p= fisher.test(matrix(c(length(intersect(genes.i,unfav.genes)),length(genes.i)-length(intersect(genes.i,unfav.genes)),
                          length(unfav.genes)-length(intersect(genes.i,unfav.genes)),
                          gN-length(unfav.genes)-length(genes.i)+length(intersect(genes.i,unfav.genes))),nrow = 2))$p.value
  return(p)
}))
p2.ti.sig.genes$favR=fav.ratios
p2.ti.sig.genes$unfavR=unfav.ratios

p2.to.sig.genes =  as.data.frame(filter(res.markers.TOthers, p_val_adj<0.05 ) %>% group_by(patient,cellType) %>% summarise(N=n()))
fav.ratios = unlist(lapply(1:nrow(p2.to.sig.genes), function(i){
  pi = p2.to.sig.genes[i,'patient']
  celli = p2.to.sig.genes[i,'cellType']
  Ni = p2.to.sig.genes[i,'N']
  genes.i = res.markers.TOthers[res.markers.TOthers$patient==pi & res.markers.TOthers$cellType == celli & res.markers.TOthers$p_val_adj<0.05,'Gene']
  p= fisher.test(matrix(c(length(intersect(genes.i,fav.genes)),length(genes.i)-length(intersect(genes.i,fav.genes)),
                          length(fav.genes)-length(intersect(genes.i,fav.genes)),
                          gN-length(fav.genes)-length(genes.i)+length(intersect(genes.i,fav.genes))),nrow = 2))$p.value
  return(p)
}))
unfav.ratios = unlist( lapply(1:nrow(p2.to.sig.genes), function(i){
  pi = p2.to.sig.genes[i,'patient']
  celli = p2.to.sig.genes[i,'cellType']
  Ni = p2.to.sig.genes[i,'N']
  genes.i = res.markers.TI[res.markers.TOthers$patient==pi & res.markers.TOthers$cellType == celli & res.markers.TOthers$p_val_adj<0.05,'Gene']
  p= fisher.test(matrix(c(length(intersect(genes.i,unfav.genes)),length(genes.i)-length(intersect(genes.i,unfav.genes)),
                          length(unfav.genes)-length(intersect(genes.i,unfav.genes)),
                          gN-length(unfav.genes)-length(genes.i)+length(intersect(genes.i,unfav.genes))),nrow = 2))$p.value
  return(p)
}))
p2.to.sig.genes$favR=fav.ratios
p2.to.sig.genes$unfavR=unfav.ratios

p2.ti.sig.genes$From = 'TIVSNI'
p2.to.sig.genes$From = 'TOVSNO'
colnames(p2.to.sig.genes)=colnames(p2.ti.sig.genes)
sig.sum = rbind(p2.ti.sig.genes,p2.to.sig.genes)
library(reshape2)
sig.sum= melt(sig.sum,measure.vars = c('favR','unfavR'))
sig.sum$cellType = unlist(lapply(sig.sum$cellType,function(a){
  strsplit(a,' ')[[1]][1]
}))
ggplot(sig.sum,
       aes(x=From,y=-log10(value),fill=variable))+geom_bar(width = 0.6,stat = 'identity',position = 'stack')+
  facet_grid(patient~cellType)+
  theme_cowplot()+scale_fill_manual(values=c('purple','cyan'))+
  theme(axis.text = element_text(size=8),axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        text = element_text(size = 8))+
  labs(x='',y='ratio')
