
dyn.load('/data1/anaconda/lib/libhdf5_hl.so.100')
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
library(SeuratDisk)
library(CellChat)
library(schex)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

my36colors <-c('#3638a9','#e8e967','#f39562','#c7eedb','#c496d4','#7e7913','#00223d','#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")
my36colors1<-c('#594C57',"#8F6D5F", '#EAE5E3','#DB9C5E',"#420B2F","#C25160","#CC5D20",'#2E59A7','#80A492','#32788A','#E1D279','#A64036')
my36colors2 <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"))  
my36colors3 <- c('#279e6895','#d6272895','#aa40fc95','#8c564b95','#b5bd6195','#17becf95','#e377c295','#aec7e895','#ffbb7895','#98df8a95','#ff989695','#c5b0d595','#1f77b495','#ff7f0e95')  



MEL167=readRDS(paste0(save.data,'MEL167.PM.rds'))
MEL167= MEL167 %>% subset(Tissue !='2ndCellline')
MEL167 <- RunUMAP(object = MEL167, dims = 1:30, reduction = "harmony")

MEL167$Tissue1=MEL167$Tissue
MEL167@meta.data[MEL167$Tissue %in% c('Kidney','Lung','Ovary'),'Tissue1']='Meta'

# MEL167@meta.data[MEL167$recluster %in% c('1st_IL24','1st_POSTN','1st_prolif_CDK1','1st_prolif_PTTG1','1st_TF','1st_SPP1') ,'recluster1']='1st'
# #MEL167@meta.data[MEL167$recluster %in% c('2nd_SOX5') ,'recluster1']='2nd'
# MEL167@meta.data[MEL167$recluster %in% c('KidneyMeta_SERPINE2','LungMeta_ANXA1','OvaryMeta_PMEL','OvaryMeta_TTC14') ,'recluster1']='Meta'
# 
# PEM22@meta.data[PEM22$recluster %in% c('1st_ACAT2','1st_CCNB1','1st_LGALS3') ,'recluster1']='1st'
# PEM22@meta.data[PEM22$recluster %in% c('KidneyMeta_FAM172A','KidneyMeta_SOX5','LungMeta_KLF6','OvaryMeta_ISG15','OvaryMeta_JUN','OvaryMeta_STMN1') ,'recluster1']='Meta'
# 

#-----++++Fig1---------


p=DimPlot(MEL167,group.by = 'recluster')+mytheme+scale_color_manual(values = my36colors3)
ggsave(paste0(save.pic,'MEL167-cluster.pdf'),width = 6.2,height = 5,p)
p=DimPlot(MEL167,group.by = 'Tissue1')+mytheme+scale_color_manual(values = my36colors2)
ggsave(paste0(save.pic,'MEL167-Tissue1.pdf'),width = 6.2,height = 5,p)
p=DimPlot(MEL167,group.by = 'Tissue')+mytheme+scale_color_manual(values = my36colors1)
ggsave(paste0(save.pic,'MEL167-Tissue.pdf'),width = 6,height = 5,p)



markers=wilcoxauc(MEL167,group_by = 'recluster')  %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}


features=c('LGALS3','SCD','CD36','STMN1','CCNB1','TOP2A','NRG3','SOX5','PDE4D','HES1','SOD3','EGR1','APOE','IGFBP7','BCAN','JUN','TAGLN2','FOS','SERPINE2','EEF2','MT2A','PTTG1','NUSAP1','UBE2C','VEGFA','IGFBP5','SLC2A1')

p=DotPlot(MEL167,group.by = 'recluster' ,features = features[length(features):1] ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#17becf95')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=1, vjust=1,size = rel(1),color = 'black'))
ggsave(paste0(save.pic,'MEL167-Dotplot.pdf'),width =8,height = 4,p)



# ----TF Expression-----

kegmt=read.table(paste0(save.data,'tf-target-infomation.txt'),header=T,sep='\t')
markers=wilcoxauc(MEL167,group_by = 'Tissue1')  #%>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
markers=markers[markers$group=='Meta',]
markers = markers[markers$feature %in% unique(kegmt$TF) , ] 
markers=markers[order(markers$logFC,decreasing = T),]

genes.choose=markers[1:10,'feature']
p=DotPlot(MEL167,group.by = 'Tissue1' ,features = c(genes.choose) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#17becf95')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=1, vjust=1,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=1, vjust=1,size = rel(1),color = 'black'))
ggsave(paste0(save.pic,'MEL167-Dotplot-TF.pdf'),width = 6,height = 2.5,p)


#------JUN and HES1-------
p=FeaturePlot(MEL167,c('HES1','JUN'),blend = T)+mytheme
ggsave(paste0(save.pic,'MEL167-HES1-JUN-exp.pdf'),width = 20,height = 5,p)

#-----#-----Monocle2--------

set.seed(1300)
Mono= monoclepic(MEL167,num=300,onlyTumor = F,use.col='recluster')
p=plot_cell_trajectory(Mono, color_by = "recluster",cell_size = 2)+mytheme+scale_color_manual(values =my36colors3)
ggsave(paste0(save.pic,'MEL167-Monocle2.pdf'),width = 6,height = 5,p)

test=orderCells(Mono,root_state = 3) 
p=plot_cell_trajectory(test, color_by = 'Pseudotime',cell_size = 2)+mytheme+scale_color_gradientn(colors = c('#cacaca30','#fef4d2','#17becf95'))
ggsave(paste0(save.pic,'MEL167-persudotime.pdf'),width = 5.7,height = 5,p)


gene.choose=genes.choose
Pseudotime=data.frame(row.names = colnames(test),Pseudotime=test$Pseudotime,barcode=colnames(test)) %>% {.[order(.$Pseudotime),] }
plots <- list()
x_axis=1:dim(Pseudotime)[1]
for(gene_id in gene.choose){
  y_axis <- as.numeric( test@assayData$exprs[gene_id,rownames(Pseudotime)] )
  lo <- loess(y_axis~x_axis)
  xl <- seq(min(x_axis),max(x_axis),(max(x_axis) - min(x_axis))/1000)
  y_axis <- predict(lo,xl)
  #y_axis <- rescale(y_axis, newrange = c(0,1))
  df <- data.frame(cells=1:length(y_axis),Expression=as.numeric(y_axis))
  df$gene <- gene_id
  df$cells <- factor(df$cells, levels=df$cells)
  num_subpops <- length(unique(df$population))
  plots[[gene_id]] <- df
}
tables <- do.call(rbind, plots)
# tables$trajectory=1
# tables[tables$gene %in% c('KRT18','EPCAM'),'trajectory']='1'
# tables[tables$gene %in% c('GMNN'),'trajectory']='2'
#
alldfs <- tables
p=ggplot(alldfs, aes(x=cells, y=Expression, group=gene, colour=gene)) +  #,linetype=trajectory
  theme_bw() + geom_line(size=1) + xlab('Pseudotime') +mytheme +
  guides(col=guide_legend(direction="vertical")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5))+
  scale_color_manual("", values=my36colors3)

ggsave(paste0(save.pic,'MEL167-TF-persudotime.pdf'),width = 5.7,height = 5,p)



#----HES1/JUN expression -----
scRNA <- make_hexbin(MEL167, nbins = 20, dimension_reduction = "UMAP")
p <- plot_hexbin_gene(scRNA,type = "data",gene = "HES1", action = "mean") + scale_fill_viridis_c(option = "A")+mytheme
ggsave(paste0(save.pic,'MEL167-HES1.umap-exp.pdf'),width =6,height = 5,p)

p <- plot_hexbin_gene(scRNA,type = "data",gene = "JUN", action = "mean") + scale_fill_viridis_c(option = "A")+mytheme
ggsave(paste0(save.pic,'MEL167-JUN.umap-exp.pdf'),width =6,height = 5,p)


markers=wilcoxauc(MEL167,group_by = 'Tissue1')  %>% { .[.$logFC>0.58 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
markers=markers[markers$group=='Meta',]


Go.analysis=enrichGO(gene=markers$feature,OrgDb = org.Hs.eg.db,ont='ALL',pvalueCutoff = 1,keyType = 'SYMBOL')
Pathway=Go.analysis@result
Pathway=Pathway[Pathway$ID%in% c('GO:0006979','GO:0000302','GO:0034599','GO:0005925','GO:0051881','GO:0009612'),]
Pathway[,'logP']= -log2(Pathway$pvalue)
Pathway$Description=as.factor(Pathway$Description) %>% factor(levels = Pathway[order(Pathway$logP),'Description'])
p=ggplot(Pathway,aes(x=logP,y=Description))+mytheme+geom_bar(stat='identity')+scale_fill_manual(values = my36colors)
ggsave(paste0(save.pic,'Pathway-SC-Meta.pdf'),width =6,height = 6.3,p)


p <- SCpubr::do_NebulosaPlot(sample = MEL167,features = 'PIEZO1',viridis.palette = "D")  
ggsave(paste0(save.pic,'MEL167-PIEZO2-exp.pdf'),width =5,height = 5,p)





#------++++Fig2  vivo vitro RNA-------
CTCvitro=readRDS(paste0(save.data,'MEL167_KO_in_vitro_TPM_matrix.RDS'))
rownames(CTCvitro)=CTCvitro$gene_id

CTCvitro.hes1=CTCvitro[,c(4:9 , 16:18)]
colnames(CTCvitro.hes1)=c('KO_1_1','KO_1_2','KO_1_3','KO_2_1','KO_2_2','KO_2_3','VEC_1','VEC_2','VEC_3')


gmt <- readLines(paste0(save.data,'h.all.v7.0.symbols.gmt'))
gmt <- strsplit(gmt, "\t")
names(gmt) <- vapply(gmt, function(y) y[1], character(1))
geneset = lapply(gmt, "[", -c(1:2))
 

ssGSEA_matrix <- GSVA::gsva(as.matrix(CTCvitro.hes1),geneset,method='plage',ssgsea.norm=T) 

ac=data.frame(row.names = colnames(CTCvitro.hes1),V1=c(rep('KO_1',3) ,rep('KO_2',3),rep('VEC',3) ))
bk = unique(c(seq(-2,2, length=100)))
pdf(paste0(save.pic,'Vitro.heatmap.pathway.pdf'),width = 6,height = 3)
pheatmap(ssGSEA_matrix[c(29,9,8,27,17,43,23),],show_colnames =F,show_rownames = T,
         cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = bk,color=colorRampPalette(c("blue", "white", "orange"))(100),
         annotation_col=ac)

dev.off()


Human.RNA <- read.table(paste0(save.data,'/RNA-seq/MEL167_KO_HES1_CD298_expected_count_output.matrix'))
colnames(Human.RNA)=c('KO_2_1','KO_2_2','KO_2_3','KO_2_4','VEC_1','VEC_2','KO_1_1','KO_1_2','KO_1_3','KO_1_4','KO_1_5','KO_1_6','VEC_3','VEC_4','VEC_5','VEC_6')
Human.TPM <- readRDS(paste0(save.data,'/RNA-seq/MEL167_KO_HES1_CD298_TPM_matrix.RDS'))
rownames(Human.TPM)=Human.TPM$gene_id
Human.TPM=Human.TPM[,-1:-3]
colnames(Human.TPM)=c('KO_2_1','KO_2_2','KO_2_3','KO_2_4','VEC_1','VEC_2','KO_1_1','KO_1_2','KO_1_3','KO_1_4','KO_1_5','KO_1_6','VEC_3','VEC_4','VEC_5','VEC_6')
colnames(Human.RNA)=c('KO_2_1','KO_2_2','KO_2_3','KO_2_4','VEC_1','VEC_2','KO_1_1','KO_1_2','KO_1_3','KO_1_4','KO_1_5','KO_1_6','VEC_3','VEC_4','VEC_5','VEC_6')

CTCvivo.hes1=Human.TPM

res.new=data.frame(row.names = colnames(CTCvitro.hes1) , group=c( rep('KO',6) ,rep('Vector' , 3)    ) )
design <- model.matrix(~0+factor(res.new$group))
colnames(design)=levels(factor(res.new$group))

contrast.matrix<-makeContrasts("KO-Vector",levels = design)
fit <- lmFit( log2(CTCvitro.hes1[,rownames(res.new)]+1),design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
et<-topTable(fit2, coef=1, n=Inf) %>% na.omit()
res <- et %>% dplyr::mutate( name=rownames(.) ) %>% dplyr::mutate(logpvalue= -log10(.$P.Value) ) %>%  { .[! .$name %>% str_count('^MT-|RPL|IGHV|LOC|TRAV|AL|AF|WI2|SNORA|CTC-|XXyac|TRBV|RPS|IGK|LINC|^AC|^AP|^CTD|XXbac|RP.+') %>% as.logical(),]} %>%
  { .[is.infinite(.$logpvalue),'logpvalue'] = max(.$logpvalue[.$logpvalue!=max(.$logpvalue)]) ;. } %>%  {.[.$logpvalue>50,'logpvalue']=50 ;.}



log2FC=0.58
pvalue=0.05
f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(pvalue))
  dff <- rbind(data.frame(x = inputx + log2FC, y = y),
               data.frame(x = -(inputx + log2FC), y = y))
  return(dff)
}
dff_curve <- f(2)
res$curve_y <- case_when(
  res$logFC > 0 ~ 1/(res$logFC-log2FC) + (-log10(pvalue)),
  res$logFC <= 0 ~ 1/(-res$logFC-log2FC) + (-log10(pvalue))
  
)
  
res$group2 <- case_when(
  res$logpvalue > res$curve_y & res$logFC >= log2FC ~ 'up',
  res$logpvalue > res$curve_y & res$logFC <= -log2FC ~ 'down',
  TRUE ~ 'none'
)

res$group2 <- factor(res$group2, levels = c("up","down","none")) #指定顺序
mycol2 <- c("#F8B60670","#4A198570","#d8d8d870")
p4 <- ggplot(data = res,
             aes(x = logFC, y = logpvalue, color = group2)) + #使用新的分组
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2)) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0)),
                     limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
  scale_colour_manual(name = "", values = alpha(mycol2, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y), #曲线坐标
            color = "black",lty = "dashed", size = 0.7) +
  mytheme

tophit= res#[res$group2 %in% c('up','down'),] %>% dplyr::group_by(group2) %>% top_n(100,wt=abs(logFC)) %>% data.frame()
tophit=tophit[tophit$name %in% intersect.gene,]


p5 <- p4 +
  geom_text_repel(data = tophit,
                  aes(x = logFC, y = logpvalue, label = name),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)
res.vitro=res

ggsave(paste0(save.pic,'Vitro.huoshan.pdf'),width =6,height = 6.3,p5)


res.new=data.frame(row.names = colnames(CTCvivo.hes1) , group=c( rep('KO',4) ,rep('Vector' , 2), rep('KO',6) ,rep('Vector' , 4)   ) )
design <- model.matrix(~0+factor(res.new$group))
colnames(design)=levels(factor(res.new$group))

contrast.matrix<-makeContrasts("KO-Vector",levels = design)
fit <- lmFit( log2(CTCvivo.hes1[,rownames(res.new)]+1),design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
et<-topTable(fit2, coef=1, n=Inf) %>% na.omit()
res <- et %>% dplyr::mutate( name=rownames(.) ) %>% dplyr::mutate(logpvalue= -log10(.$P.Value) ) %>%  { .[! .$name %>% str_count('^MT-|RPL|LOC|IGHV|TRAV|AL|AF|WI2|SNORA|CTC-|XXyac|TRBV|RPS|IGK|LINC|^AC|^AP|^CTD|XXbac|RP.+') %>% as.logical(),]} %>%
  { .[is.infinite(.$logpvalue),'logpvalue'] = max(.$logpvalue[.$logpvalue!=max(.$logpvalue)]) ;. } %>%  {.[.$logpvalue>50,'logpvalue']=50 ;.}



log2FC=0.58
pvalue=0.05
f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(pvalue))
  dff <- rbind(data.frame(x = inputx + log2FC, y = y),
               data.frame(x = -(inputx + log2FC), y = y))
  return(dff)
}
res$curve_y <- case_when(
  res$logFC > 0 ~ 1/(res$logFC-log2FC) + (-log10(pvalue)),
  res$logFC <= 0 ~ 1/(-res$logFC-log2FC) + (-log10(pvalue))
  
)

res$group2 <- case_when(
  res$logpvalue > res$curve_y & res$logFC >= log2FC ~ 'up',
  res$logpvalue > res$curve_y & res$logFC <= -log2FC ~ 'down',
  TRUE ~ 'none'
)

res$group2 <- factor(res$group2, levels = c("up","down","none")) #指定顺序
mycol2 <- c("#F8B60670","#4A198570","#d8d8d870")
p4 <- ggplot(data = res,
             aes(x = logFC, y = logpvalue, color = group2)) + #使用新的分组
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 2)) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0)),
                     limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
  scale_colour_manual(name = "", values = alpha(mycol2, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y), #曲线坐标
            color = "black",lty = "dashed", size = 0.7) +
  mytheme

tophit= res#[res$group2 %in% c('up','down'),] %>% dplyr::group_by(group2) %>% top_n(100,wt=abs(logFC)) %>% data.frame()
tophit=tophit[tophit$name %in% intersect.gene,]


p5 <- p4 +
  geom_text_repel(data = tophit,
                  aes(x = logFC, y = logpvalue, label = name),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)

res.vivo=res

intersect.gene=intersect(res.vivo[res.vivo$group2=='up','name'] , res.vitro[res.vitro$group2=='up','name'] )

ggsave(paste0(save.pic,'Vivo.huoshan.pdf'),width =6,height = 6.3,p5)


Go.analysis=enrichGO(gene=intersect.gene,OrgDb = org.Hs.eg.db,ont='ALL',pvalueCutoff = 1,keyType = 'SYMBOL')
Pathway=Go.analysis@result
Pathway=Pathway[Pathway$ID%in% c('GO:0030318','GO:0050931','GO:0048066','GO:0018212','GO:0051591','GO:0022407','GO:0002693'),]
Pathway[,'logP']= -log2(Pathway$pvalue)
Pathway$Description=as.factor(Pathway$Description) %>% factor(levels = Pathway[order(Pathway$logP),'Description'])
p=ggplot(Pathway,aes(x=logP,y=Description))+mytheme+geom_bar(stat='identity')+scale_fill_manual(values = my36colors)
ggsave(paste0(save.pic,'Pathway-vivo-vitro.pdf'),width =6,height = 6.3,p)



#韦恩图
pdf(paste0(save.pic,'Venn.pdf'),width = 4,height=4)
venn_ploy <- venn.diagram(
  x = list(SC.cTEC = res.vivo[res.vivo$group2=='up','name'],T2 = res.vitro[res.vitro$group2=='up','name']),
  filename = NULL,fill = c('#31B1FF80','#FF559680') )
grid.draw(venn_ploy)
dev.off()




#-----对一些基因基因绘制------
CTC.bulk= cbind(CTCvitro.hes1 %>% {colnames(.)=paste0('Vitro_',colnames(.)) ;.} , CTCvivo.hes1 %>% {colnames(.)=paste0('Vivo_',colnames(.)) ;.}  )
gene.data=melt(data.frame(gene=c('MITF','CITED1'),CTC.bulk[c('MITF','CITED1'),]))
gene.data$variable=  gene.data$variable %>% str_extract('.+_KO|.+_VEC')
gene.data$variable = gene.data$variable %>% as.factor() %>% factor(levels = c('Vitro_VEC','Vitro_KO','Vivo_VEC','Vivo_KO'))

pic=list()
for(i in c('MITF','CITED1')){
  pic[[i]]=ggplot(gene.data[gene.data$gene==i,],aes(variable,value))+
    stat_summary(mapping=aes(fill=variable),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.8,linewidth=0.2,color='black')+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
    labs(x = "",y = "Mean (log2(TPM+1))")+mytheme+scale_fill_manual(values = c('#31B1FF80','#FF559680','#31B1FF80','#FF559680') )+theme(
      axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(0.8)))+theme(legend.position = 'none')+ggtitle(i)+
    theme(plot.title = element_text(size=rel(0.8)))+
    theme(axis.title.y = element_text(size=10))+theme(
      axis.text.y = element_text(size = rel(0.8)))+geom_point(size=2.5)+stat_compare_means(comparisons = list(c('Vitro_KO','Vitro_VEC') , c('Vivo_KO','Vivo_VEC')), label = "p.signif")
}

p=CombinePlots(pic,ncol=2)
ggsave(paste0(save.pic,'Bulk.barplot-MITF_CITED1.pdf'),width =5,height = 4,p)


CTC.bulk= cbind(CTCvitro.hes1 %>% {colnames(.)=paste0('Vitro_',colnames(.)) ;.} , CTCvivo.hes1 %>% {colnames(.)=paste0('Vivo_',colnames(.)) ;.}  )
gene.data=melt(data.frame(gene=c('VEGFA','MMP2'),CTC.bulk[c('VEGFA','MMP2'),]))
gene.data$variable=  gene.data$variable %>% str_extract('.+_KO|.+_VEC')
gene.data$variable = gene.data$variable %>% as.factor() %>% factor(levels = c('Vitro_VEC','Vitro_KO','Vivo_VEC','Vivo_KO'))

pic=list()
for(i in c('VEGFA','MMP2')){
  pic[[i]]=ggplot(gene.data[gene.data$gene==i,],aes(variable,value))+
    stat_summary(mapping=aes(fill=variable),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.8,linewidth=0.2,color='black')+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
    labs(x = "",y = "Mean (log2(TPM+1))")+mytheme+scale_fill_manual(values = c('#31B1FF80','#FF559680','#31B1FF80','#FF559680') )+theme(
      axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(0.8)))+theme(legend.position = 'none')+ggtitle(i)+
    theme(plot.title = element_text(size=rel(0.8)))+
    theme(axis.title.y = element_text(size=10))+theme(
      axis.text.y = element_text(size = rel(0.8)))+geom_point(size=2.5)+stat_compare_means(comparisons = list(c('Vitro_KO','Vitro_VEC') , c('Vivo_KO','Vivo_VEC')), label = "p.signif")
}

p=CombinePlots(pic,ncol=2)
ggsave(paste0(save.pic,'Bulk.barplot-VEGFA_MMP2.pdf'),width =5,height = 4,p)



CTC.bulk= cbind(CTCvitro.hes1 %>% {colnames(.)=paste0('Vitro_',colnames(.)) ;.} , CTCvivo.hes1 %>% {colnames(.)=paste0('Vivo_',colnames(.)) ;.}  )
gene.data=melt(data.frame(gene=c('PODXL','EDNRB'),CTC.bulk[c('PODXL','EDNRB'),]))
gene.data$variable=  gene.data$variable %>% str_extract('.+_KO|.+_VEC')
gene.data$variable = gene.data$variable %>% as.factor() %>% factor(levels = c('Vitro_VEC','Vitro_KO','Vivo_VEC','Vivo_KO'))

pic=list()
for(i in c('PODXL','EDNRB')){
  pic[[i]]=ggplot(gene.data[gene.data$gene==i,],aes(variable,value))+
    stat_summary(mapping=aes(fill=variable),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.8,linewidth=0.2,color='black')+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
    labs(x = "",y = "Mean (log2(TPM+1))")+mytheme+scale_fill_manual(values = c('#31B1FF80','#FF559680','#31B1FF80','#FF559680') )+theme(
      axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(0.8)))+theme(legend.position = 'none')+ggtitle(i)+
    theme(plot.title = element_text(size=rel(0.8)))+
    theme(axis.title.y = element_text(size=10))+theme(
      axis.text.y = element_text(size = rel(0.8)))+geom_point(size=2.5)+stat_compare_means(comparisons = list(c('Vitro_KO','Vitro_VEC') , c('Vivo_KO','Vivo_VEC')), label = "p.signif")
}

p=CombinePlots(pic,ncol=2)
ggsave(paste0(save.pic,'Bulk.barplot-PODXL_EDNRB.pdf'),width =5,height = 4,p)


CTC.bulk= cbind(CTCvitro.hes1 %>% {colnames(.)=paste0('Vitro_',colnames(.)) ;.} , CTCvivo.hes1 %>% {colnames(.)=paste0('Vivo_',colnames(.)) ;.}  )
gene.data=melt(data.frame(gene=c('ITGB1','ITGAV'),CTC.bulk[c('ITGB1','ITGAV'),]))
gene.data$variable=  gene.data$variable %>% str_extract('.+_KO|.+_VEC')
gene.data$variable = gene.data$variable %>% as.factor() %>% factor(levels = c('Vitro_VEC','Vitro_KO','Vivo_VEC','Vivo_KO'))

pic=list()
for(i in c('ITGB1','ITGAV')){
  pic[[i]]=ggplot(gene.data[gene.data$gene==i,],aes(variable,value))+
    stat_summary(mapping=aes(fill=variable),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.8,linewidth=0.2,color='black')+
    stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
    labs(x = "",y = "Mean (log2(TPM+1))")+mytheme+scale_fill_manual(values = c('#31B1FF80','#FF559680','#31B1FF80','#FF559680') )+theme(
      axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(0.8)))+theme(legend.position = 'none')+ggtitle(i)+
    theme(plot.title = element_text(size=rel(0.8)))+
    theme(axis.title.y = element_text(size=10))+theme(
      axis.text.y = element_text(size = rel(0.8)))+geom_point(size=2.5)+stat_compare_means(comparisons = list(c('Vitro_KO','Vitro_VEC') , c('Vivo_KO','Vivo_VEC')), label = "p.signif")
}

p=CombinePlots(pic,ncol=2)
ggsave(paste0(save.pic,'Bulk.barplot-ITGB1_ITGAV.pdf'),width =5,height = 4,p)



CTCpatients=read.table(paste0(save.data,'RSEM.tpm.PRJNA662599.csv'),header=T,sep=',')
CTCpatients[,1]=(CTCpatients[,1] %>% str_split_fixed('_',2))[,2] 
CTCpatients=CTCpatients[!duplicated(CTCpatients$X),]
rownames(CTCpatients) = CTCpatients$X
CTCpatients=CTCpatients[,-1]

CTCpatients = CTCpatients[,colnames(CTCpatients) %>% str_count('BloodCTC|Primary|WhiteCell') %>% as.logical()]
CTCpatients=log2(CTCpatients+1)

tmp.pic=data.frame(t(CTCpatients[c('HES1','JUN'),]) , (colnames(CTCpatients) %>% str_split_fixed('_',5)))
tmp.pic[85:90,'X3']='WhiteCell'
tmp.pic[tmp.pic$X2=='PKPEM22','X3']='Meta'

p=ggplot(tmp.pic,aes(x=HES1,y=JUN,color=X3))+geom_point(size=4)+mytheme+scale_color_manual(values = my36colors3)
ggsave(paste0(save.pic,'Pat-CTC-HES-FOS.pdf'),width =5,height = 4,p)


tmp.pic=data.frame(t(CTCpatients[c('HES1','MMP2'),]) , (colnames(CTCpatients) %>% str_split_fixed('_',5)))
tmp.pic[85:90,'X3']='WhiteCell'
tmp.pic[tmp.pic$X2=='PKPEM22','X3']='Meta'

p=ggplot(tmp.pic,aes(x=HES1,y=MMP2,color=X3))+geom_point(size=4)+mytheme+scale_color_manual(values = my36colors3)
ggsave(paste0(save.pic,'Pat-CTC-HES-MMP2.pdf'),width =5,height = 4,p)


tmp.pic=data.frame(t(CTCpatients[c('HES1','PODXL'),]) , (colnames(CTCpatients) %>% str_split_fixed('_',5)))
tmp.pic=tmp.pic[tmp.pic$X3=='BloodCTC',]
pdf(paste0(save.pic,'Pat-CTC-HES-PODXL.pdf'),height=4.5,width = 4)
smoothScatter(tmp.pic$HES1, tmp.pic$PODXL)
dev.off()

tmp.pic=data.frame(t(CTCpatients[c('HES1','MMP2'),]) , (colnames(CTCpatients) %>% str_split_fixed('_',5)))
tmp.pic=tmp.pic[tmp.pic$X3=='BloodCTC',]
pdf(paste0(save.pic,'Pat-CTC-HES-MMP2.pdf'),height=4.5,width = 4)
smoothScatter(tmp.pic$HES1, tmp.pic$MMP2)
dev.off()



#----------+++++endothelial talk 2 Celline --------

Endo=readRDS(paste0(save.data,'Endo.all.rds'))
Myeloid=readRDS(paste0(save.data,'Myeloid.all.rds'))
Caf=readRDS(paste0(save.data,'CAF.all.rds'))
p=DimPlot(Endo,group.by = 'recluster1')+scale_color_manual(values = my36colors3)+mytheme
ggsave(paste0(save.pic,'Endo.umap.pdf'),width = 6,height = 5,p)


pic=densitypic(Endo,'Tissue1')
p=  cowplot::plot_grid(plotlist=pic, ncol=2, nrow=1) 
ggsave(paste0(save.pic,'Endo-Tissue-density.pdf'),width = 10,height = 5,p) 


p=DimPlot(Endo,group.by = 'Tissue')+scale_color_manual(values = my36colors2)+mytheme
ggsave(paste0(save.pic,'Endo.-Tissue.umap.pdf'),width = 5.8,height = 5,p)


Meta.bar=Endo@meta.data  %>% {.$Tissue %>% as.character(); . } 
Meta.bar=Meta.bar %>% data.table() %$% .[,.N,.(Tissue,recluster1)] %>% data.frame() %>% {ddply(.,'recluster1','transform',per=N/sum(N))}


Ptumor.pie <- Meta.bar[which(Meta.bar$recluster1=="Endo_Fabp4"),]
rownames(Ptumor.pie) <- Ptumor.pie$Tissue

labs <- paste0(Ptumor.pie$Tissue," \n(", round(Ptumor.pie$N/sum(Ptumor.pie$N)*100,2), "%)")

pdf(paste0(save.pic,'Endo.-Tissue-pie-FABP4.pdf'),height=4,width = 4)
pie(Ptumor.pie$N, labels = labs, init.angle=90, border="white",col=my36colors2)
dev.off()

markers=wilcoxauc(Endo,group_by = 'Tissue1')
markers=markers[markers$group=='Tumor',]
rownames(markers)=markers$feature
kegmt<-read.gmt(paste0(main.path,'/data/c5.go.v7.4.symbols.gmt'))  
kegmt[,1]=kegmt[,1] %>% str_remove('GO_')  %>% tolower()%>% capitalize()


humangene=homologene::homologene(rownames(Endo), inTax = 10090, outTax = 9606)
humangene=humangene[humangene$`10090`%in% rownames(Endo),]
humangene=humangene[!duplicated(humangene$`10090`),]
rownames(humangene)=humangene$`10090`
markers[,'humangene']=humangene[markers$feature,2]
markers=markers[!is.na(markers$humangene),]

hark<-GSEA(markers$logFC %>% { names(.)=markers$humangene;.} %>% sort(decreasing = T),TERM2GENE = kegmt) #GSEA分析
GSEA=hark@result
#GSEA=GSEA[GSEA$NES>0,]


mygene=GSEA[1,'core_enrichment'] %>% str_split('/') %>% unlist()
p=gseaNb(object = hark,geneSetID = 'Gobp_response_to_oxidative_stress',addPval = T,pvalX = 0.9,pvalY = 0.8,subPlot = 1)
ggsave(paste0(save.pic,'Endo.response_to_oxidative_stress.pdf'),p,width = 5,height = 2.5)

p=gseaNb(object = hark,geneSetID = 'Gocc_adherens_junction',addPval = T,pvalX = 0.9,pvalY = 0.8,subPlot = 1)
ggsave(paste0(save.pic,'Endo.Gocc_adherens_junction.pdf'),p,width = 5,height = 2.5)


source('/data2/liuxf/Mela_mouse_human_New/program/Vec.vectore.R')

VEC = Endo@reductions$umap@cell.embeddings
rownames(VEC) = colnames(Endo)
PCA =Endo@reductions$pca@cell.embeddings
PCA=vector.rankPCA(PCA)

OUT=vector.buildGrid(VEC, N=30,SHOW=T)
# 构建网络
OUT=vector.buildNet(OUT, CUT=1, SHOW=T)
# 计算量子极化 (QP) 分数
OUT=vector.getValue(OUT, PCA, SHOW=T)
# 获取像素的 QP 分数
OUT=vector.gridValue(OUT,SHOW=T)
# 找起始点
OUT=vector.autoCenter(OUT,SHOW=T)
OUT=vector.selectCenter(OUT)
# 推断向量
pdf(paste0(save.pic,'Endo.-Vector.pdf'),height=6,width = 6)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=T, COL=OUT$COL, SHOW.SUMMIT=T)
dev.off()


#-----Moncle3-Endo-------
set.seed(1300)

cds=Endo

cds$recluster1=as.factor(cds$recluster1)
#subset(CD4.ex,cells=  (doBy::sampleBy(formula = ~ recluster,frac = 0.5,data =CD4.ex@meta.data) %>% {.$barcode} )   )

data <- GetAssayData(cds, assay = 'RNA', slot = 'counts')
gene_annotation <- data.frame(gene_short_name = rownames(cds)) %>% {rownames(.)=rownames(data);.}
cds <- new_cell_data_set(data,
                         cell_metadata = cds@meta.data,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- reduce_dimension(cds, preprocess_method = "PCA")#umap降维
cds@int_colData$reducedDims$UMAP <- Embeddings(Endo, reduction = "umap")[rownames(cds@int_colData$reducedDims$UMAP),]##从seurat导入整合过的umap坐标
cds <- cluster_cells(cds,k = 20)
cds <- learn_graph(cds)

cds@colData$recluster1=as.factor(cds@colData$recluster1)

p=plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, trajectory_graph_color = "grey20",label_roots = F,#label_cell_groups = F,
             label_branch_points = FALSE, color_cells_by="recluster1")+mytheme+scale_color_manual(values=my36colors3)+mytheme
ggsave(paste0(save.pic,'Epi.Monocle3.recluster.pdf'),p,height = 5,width=6.5)
ggsave(paste0(save.pic,'Epi.Monocle3.recluster.jpg'),p,height = 5,width=6.5)


cds <- order_cells(cds, root_cells = rownames(cds@colData[cds@colData$recluster1=='Endo_Col4a1',]) )
p=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE,label_roots = F)+mytheme+theme(text=element_text(family="sans",face = 'bold',color="black"))+scale_color_gradientn(colors = c('#fef4d2','green','blue'))
ggsave(paste0(save.pic,'Epi.Monocle3.pseudotime.pdf'),p,height = 5,width=6)


#-----Monocle2-------
set.seed(1234)
nouse=Endo@meta.data
nouse=nouse[nouse$recluster1 %in% c('Endo_Col4a1','Endo_Fabp4') ,]
nouse=nouse[!duplicated(nouse$nFeature_RNA),]

Mono= monoclepic(Endo %>% subset(cells=rownames(nouse)),num=400,onlyTumor = F,use.col='recluster1')


bks = unique(c(seq(-2,2, length=100)))
pse.gene=c('Tjp1','Tek','Cdh5','Cldn5','Adgrl2','Pcdh17','Ramp2','Jund','Junb','Jun','Fos','Fosb','Zfp36','Egr1','Hspa1a','Nr4a1','Atf3','Selp','Vcam1','Jam2')
m <- genSmoothCurves(Mono[pse.gene,], new_data = data.frame(Pseudotime = seq(min(Mono@phenoData@data[,'Pseudotime']),max(Mono@phenoData@data[,'Pseudotime']), length.out = 100)))
m=m[pse.gene,] %>% { .[!apply(., 1, sd) == 0, ]} %>% {Matrix::t(scale(Matrix::t(.), center = TRUE))} 

pdf(file=paste0(save.pic,'monocle.heatmap.markers.pdf'),height = 4,width=6)
p=pheatmap::pheatmap(m,show_colnames = F,show_rownames = T,cluster_rows = F,breaks = bks,cluster_cols = F,color=blue2green2red(length(bks) - 1), border_color = NA)
dev.off()


markers=wilcoxauc(Endo %>% subset(recluster1 %in%c('Endo_Col4a1','Endo_Fabp4') ),group_by = 'recluster1') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
intersect(markers[markers$group=='Endo_Fabp4','feature'] , kegmt[kegmt$term=='Gobp_cellular_extravasation' ,2]%>% tolower()%>% capitalize())
intersect(markers[,'feature'] , kegmt[kegmt$term=='Gobp_cellular_extravasation' ,2]%>% tolower()%>% capitalize())





#-----step1: CellChat cellline_Endo-------
MEL167@meta.data[,'HES1.exp']=MEL167@assays$RNA@data['HES1',] %>% as.numeric()
now.meta=MEL167@meta.data
now.meta=now.meta[now.meta$Tissue=='Cellline',]
now.meta=now.meta[now.meta$HES1.exp>0 ,]
now.meta[,'group']='Low'
now.meta[now.meta$HES1.exp >  quantile(now.meta$HES1.exp,0.5),'group']='High'
# Endo1=Endo
# del.cell=Endo1 %>% subset(recluster1=='Endo_Fabp4')
# del.cell=del.cell@meta.data[del.cell@assays$RNA@counts['Sele',]==0 , ]
# Endo1=Endo %>% subset(cells= setdiff(rownames(Endo@meta.data) , rownames(del.cell) ) )
# 
now.meta1=Endo@meta.data
now.meta1[,'group']=now.meta1$recluster1

now.count=MEL167@assays$RNA@data[,rownames(now.meta)]



now.count1=Endo@assays$RNA@data[,rownames(now.meta1)]
humangene=homologene::homologene(rownames(Endo), inTax = 10090, outTax = 9606)
humangene=humangene[humangene$`9606`%in% rownames(MEL167),]
now.count=now.count[humangene$`9606`,]
now.count1=now.count1[humangene$`10090`,]


new_count2=cbind(now.count, now.count1 )


now.meta=now.meta[,c('Patient','group')]
now.meta1=now.meta1[,c('Patient','group')]

cellchat <- createCellChat(object = new_count2, meta = rbind(now.meta,now.meta1), group.by = "group")
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) %>% identifyOverExpressedGenes( only.pos = T) %>% identifyOverExpressedInteractions() %>% projectData(PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) %>% filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()


pdf(paste0(save.pic,'/Endo2HES1-heatmap.pdf'),width = 5,height = 5)
p=netVisual_heatmap(cellchat,color.heatmap = c("white", "orange", "red") )
p
dev.off()



pdf(paste0(save.pic,'/Endo2HES1.pdf'),width = 5,height = 3)
p=netVisual_bubble(cellchat,sources.use =c(1:5),targets.use = c(6,7),color.heatmap = c("Spectral"),n.colors = 11 )
#p$data=p$data[(p$data$source %in% c('C1','C2','C3') & p$data$target %in% c('Mac_CCL3' ,'Mac_CXCL10' ,'Mac_MMP9' ,'Mac_SLC40A1')) ,]
p
dev.off()

pdf(paste0(save.pic,'/HES12Endo.pdf'),width = 5,height = 4)
p=netVisual_bubble(cellchat,sources.use =c(6:7),targets.use = c(1:5),color.heatmap = c("Spectral"),n.colors = 11 )
#p$data=p$data[(p$data$source %in% c('C1','C2','C3') & p$data$target %in% c('Mac_CCL3' ,'Mac_CXCL10' ,'Mac_MMP9' ,'Mac_SLC40A1')) ,]
p
dev.off()



#--------体外HES1 高低分组-------

tophit= res.CTC.HES1ko[res.CTC.HES1ko$name %in% c('ANGPT1','ANGPT2','APP','CADM1','LAMA4','LAMB1','LAMC1','MIF','MPZL1','VCAM1','KIT','EDNRB','EPHA3','ITGAV','ITGB1','CD44'),]
tophit=tophit[tophit$color=='KO',]

xmax=max(abs(min(res.CTC.HES1ko$logFC)),max(res.CTC.HES1ko$logFC))
p <- ggplot(res.CTC.HES1ko, aes(x = logFC, y = logpvalue)) +xlim(-xmax,xmax)+geom_point(aes(color = color),size=3)+
  scale_color_manual(values = c('#F8323A40',"grey",'#0183CB40')) +mytheme+
  geom_text_repel(data = tophit,aes(label = name),size = 3,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"),segment.color='black',colour = "black")+ 
  geom_vline(aes(xintercept=-0.58),colour="darkgrey", linetype="dashed")+geom_vline(aes(xintercept=0.58),colour="darkgrey", linetype="dashed") +geom_hline(aes(yintercept=1.3),colour="darkgrey", linetype="dashed")

Gene.choose=c('ANGPT1','ANGPT2','APP','CADM1','LAMA4','LAMB1','LAMC1','MIF','MPZL1','VCAM1','KIT','EDNRB','EPHA3','ITGAV','ITGB1','CD44')
p=StackedVlnPlot(MEL167, Gene.choose, 'recluster',color=my36colors3,pt.size=0)

