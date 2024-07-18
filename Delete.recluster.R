Delete.recluster=function(Proseurat,FindNei.dim=35){
  Proseurat <- NormalizeData(object = Proseurat, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(vars.to.regress = c("nCount_RNA"))  %>%  { RunPCA(. , features = VariableFeatures(object = .))  ; . }
  Proseurat  <- RunHarmony(Proseurat, group.by.vars = c("Patient"), dims.use = 1:33, verbose = T) %>%  FindNeighbors(dims = 1:FindNei.dim, reduction = "harmony") %>%  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30, reduction = "harmony")
  return(Proseurat)
}