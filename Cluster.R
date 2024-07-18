
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(infercnv);library(CaSpER);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(Startrac);library(scRepertoire)
library(SCEVAN)
library(GseaVis)
library(irGSEA)
library(scMetabolism)
library(ChIPseeker)
library(org.Hs.eg.db)
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(VennDiagram)
library(SCPA)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

my36colors <-c('#3638a9','#e8e967','#f39562','#c7eedb','#c496d4','#7e7913','#00223d','#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")

mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1)),
        axis.text.x = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5),
        text=element_text(family="sans",face = 'bold'))

set.seed(1234)


#-----++++MEL167-CTC--------
counts1=readRDS(paste0(save.data,'NCG_MEL167_Kidney_2.human.rds'))
colnames(counts1) = paste0('MEL167_Kidney_M2_',colnames(counts1))
counts2=readRDS(paste0(save.data,'NCG_MEL167_Kidney.human.rds'))
colnames(counts2) = paste0('MEL167_Kidney_M1_',colnames(counts2))
counts3=readRDS(paste0(save.data,'NCG_MEL167_Lung.human.rds'))
colnames(counts3) = paste0('MEL167_Lung_M1_',colnames(counts3))
counts4=readRDS(paste0(save.data,'NCG_MEL167_Lung_2.human.rds'))
colnames(counts4) = paste0('MEL167_Lung_M2_',colnames(counts4))
counts5=readRDS(paste0(save.data,'NCG_MEL167_Ovary.human.rds'))
colnames(counts5) = paste0('MEL167_Ovary_M1_',colnames(counts5))
counts6=Read10X(paste0(save.data,'/MEL167'))
colnames(counts6) = paste0('MEL167_Cellline_M1_',colnames(counts6))
counts7=readRDS(paste0(save.data,'NCG_MEL167_LungNormal_and_CTC.human.rds'))
colnames(counts7) = paste0('MEL167_2ndCellline_M1_',colnames(counts7))



counts=Reduce(cbind,list(counts1,counts2,counts3,counts4,counts5,counts6,counts7))
counts <- CreateSeuratObject(counts =counts)
counts[["percent.mt"]] = PercentageFeatureSet(counts, pattern = "^[Mm]T-")
counts[["percent.RPL"]] = PercentageFeatureSet(counts, pattern = "^RP[SL].+")
counts =counts %>% subset(percent.mt<10 & nFeature_RNA>200)
counts <- NormalizeData(object = counts, normalization.method = "LogNormalize", scale.factor = 10000)
counts <- FindVariableFeatures(object = counts, selection.method = "vst", nfeatures = 2000)
counts <- ScaleData(counts, vars.to.regress = c("nCount_RNA"))
counts <- RunPCA(object =counts, features = VariableFeatures(object = counts))


message=rownames(counts@meta.data) %>% str_split_fixed('_',n=6) 
counts@meta.data[,'Name']=message[,1];counts@meta.data[,'Tissue']=message[,2]
counts@meta.data[,'Patient']=message[,3]

counts <- RunHarmony(counts, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T)
counts <- FindNeighbors(object = counts, dims = 1:35, reduction = "harmony")#
counts <- FindClusters(object = counts, resolution = 0.5)
counts<- RunUMAP(object = counts, dims = 1:30, reduction = "harmony")#

#----clean-------
counts= counts %>% subset(seurat_clusters %in% c(0:7,9:11))
counts=Delete.recluster(counts,FindNei.dim=35)

counts= counts %>% subset(seurat_clusters %in% c(0:9))
counts=Delete.recluster(counts,FindNei.dim=35)

name=c('0'='C1_LGALS3','1'='M1_NRG3','2'='M2_HES1','3'='C2_STMN1','4'='M3_APOE','5'='M4_JUN','6'='M5_CD36','7'='M6_PTTG1','8'='M7_VEGFA','9'='C1_LGALS3')
counts$recluster=name[counts$seurat_clusters]

saveRDS(counts,paste0(save.data,'MEL167.PM.rds'))


counts=Reduce(cbind,list(counts1,counts2,counts4,counts6))
counts <- CreateSeuratObject(counts =counts)
counts[["percent.mt"]] = PercentageFeatureSet(counts, pattern = "^[Mm]T-")
counts[["percent.RPL"]] = PercentageFeatureSet(counts, pattern = "^RP[SL].+")
counts =counts %>% subset(percent.mt<10 & nFeature_RNA>200)
counts <- NormalizeData(object = counts, normalization.method = "LogNormalize", scale.factor = 10000)
counts <- FindVariableFeatures(object = counts, selection.method = "vst", nfeatures = 2000)
counts <- ScaleData(counts, vars.to.regress = c("nCount_RNA"))
counts <- RunPCA(object =counts, features = VariableFeatures(object = counts))


message=rownames(counts@meta.data) %>% str_split_fixed('_',n=6) 
counts@meta.data[,'Name']=message[,1];counts@meta.data[,'Tissue']=message[,2]
counts@meta.data[,'Patient']=message[,3]

counts <- FindNeighbors(object = counts, dims = 1:35)#, reduction = "harmony"
counts <- FindClusters(object = counts, resolution = 0.5)
counts<- RunUMAP(object = counts, dims = 1:30)#, reduction = "harmony"


counts$recluster=counts$seurat_clusters
Print.report(counts,c('MLANA','PMEL','MITF','TYR','SOX10','PLP1'), 'All.cell.Report' )

markers=wilcoxauc(counts,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
top2heatmap=top2heatmap[!duplicated(top2heatmap$feature),]
marker_heatmap2(top2heatmap,counts,save.pic,'recluster','heatmap.all', levels(as.factor(top2heatmap$group)) ,height=25)

