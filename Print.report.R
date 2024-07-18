Print.report=function(Proseurat,signatures,file.name,show.legend=F,show.TCR=F){
  all.cell.result=list()
  
  if(show.legend==T){
    markers=wilcoxauc(Proseurat) %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% group_by(group) %>% top_n(5,wt=logFC) %>% data.frame()
    now.repel=markers %>% data.table() %>%   { .[,.(This.features=paste(feature,collapse = ' ~ ')   ),by='group' ]  } %>% data.frame() %>% { colnames(.)[1]='seurat_clusters'; rownames(.)= .[,1]; .}
    nouse=cbind(Proseurat@meta.data,Proseurat@reductions$umap@cell.embeddings) %>%  {dplyr::mutate(., now.repel=now.repel[as.character(.$seurat_clusters),'This.features'])} 
    all.cell.result[[1]]=ggplot(nouse, aes( x=UMAP_1 ,y= UMAP_2, color= seurat_clusters ))+geom_point(size=0.4)+mytheme+scale_color_manual(values = my36colors)+
      geom_text_repel(data = nouse[!duplicated(nouse$seurat_clusters),],aes(label = now.repel ),size = 3,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))
  }else{
    all.cell.result[[1]]=DimPlot(Proseurat,label=T)+mytheme
  }
  
  all.cell.result[[2]]=DimPlot(Proseurat,group.by = 'recluster')+mytheme+scale_color_manual(values = c(my36colors,my36colors) )
  all.cell.result[[3]]=DimPlot(Proseurat,group.by = 'Tissue')+mytheme+scale_color_manual(values = my36colors)
  all.cell.result[[4]]=DimPlot(Proseurat,label=T,group.by = 'Patient')+mytheme
  if(show.TCR==T){
    
    Proseurat@meta.data[is.na(Proseurat@meta.data)]=0
    Proseurat@meta.data[,'TCRnum.log']=log10(as.integer(Proseurat@meta.data[,'TCRnum'])+1)
    all.cell.result[[6]]=FeaturePlot(subset(Proseurat,subset=TCRnum.log>0),features='TCRnum.log',cols = c('#f0ff8945','#69ffba','#3600ff'),pt.size =0.5)+mytheme+labs(title = element_blank())
    Proseurat@meta.data[Proseurat$TCRnum>0,'TCRnum']='Have'
    Proseurat@meta.data[Proseurat$TCRnum==0,'TCRnum']='No-Have'
    all.cell.result[[7]]=DimPlot(Proseurat,group.by = 'TCRnum')+mytheme+scale_color_manual(values = c('red','#d5d5d501'))
    
  }
  
  for(i in signatures) all.cell.result[[i]]= FeaturePlot(Proseurat,i,max.cutoff = 3,label=T)+mytheme
  pic=CombinePlots(all.cell.result,ncol=4)
  ggsave(paste0(save.pic,file.name,'.jpg'),pic,width = 20,height = ceiling(length(all.cell.result)/4)*4)
}