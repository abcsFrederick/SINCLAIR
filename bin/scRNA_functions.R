##################################################################
# Seurat Pre-processing
##################################################################
SEURAT_CLUSTERING = function(so){
  
  so$Sample = tail(strsplit(h5,"/")[[1]],1)
  
  so <- RunPCA(object = so, 
               features = VariableFeatures(object = so), 
               do.print = TRUE, 
               pcs.print = 1:5,
               genes.print = 0,verbose=F,npcs = 30)
  npcs = 30
  so <- FindNeighbors(so,dims = 1:npcs)
  so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.1,algorithm = 3)
  so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.2,algorithm = 3)
  so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.5,algorithm = 3)
  so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.6,algorithm = 3)
  so <- FindClusters(so,dims = 1:npcs, print.output = 0, resolution = 0.8,algorithm = 3)
  so <- RunUMAP(so,dims = 1:npcs,n.components = 3L)
  return(so)
}

CONVERT_TO_HUMAN_GENELIST  <- function(gns){
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  return(as.character(unlist(mapped$MUS )))
}

MAIN_PROCESS_SO<-function(so_in, ref){
  # so_in=so_filt
  
  # assign genes depending on ref input
  if(ref=="hg38" || ref == "hg19"){
    print("--proccesing human data")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  } else if(ref=="mm10"){
    print("--proccesing mouse data")
    s.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$s.genes)
    g2m.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$g2m.genes)
  }
  
  # process
  so_1 = SCTransform(so_in)
  so_2 = NormalizeData(so_1,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000,
                       assay = "RNA")
  so_3 = ScaleData(so_2, assay = "RNA")
  so_4 = CellCycleScoring(so_3,
                          s.features = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident = TRUE)
  so_out = SEURAT_CLUSTERING(so_4)
  return(so_out)
}

##################################################################
#
##################################################################
RUN_SINGLEr = function(obj,refFile,fineORmain){
  obj = DietSeurat(obj, graphs = "umap")
  sce = as.SingleCellExperiment(obj,assay = "SCT")
  ref = refFile
  s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]])
  return(s$pruned.labels)
}

MAIN_LABEL_SO<-function(so_in, ref){
  if(ref == "hg38" || ref == "hg19"){
    so_in$HPCA_main <- RUN_SINGLEr(so_in,celldex::HumanPrimaryCellAtlasData(),"label.main")
    so_in$HPCA <-  RUN_SINGLEr(so_in,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    so_in$BP_encode_main <-  RUN_SINGLEr(so_in,celldex::BlueprintEncodeData(),"label.main")
    so_in$BP_encode <-  RUN_SINGLEr(so_in,celldex::BlueprintEncodeData(),"label.fine")
    so_in$monaco_main <-  RUN_SINGLEr(so_in,celldex::MonacoImmuneData(),"label.main")
    so_in$monaco <-     RUN_SINGLEr(so_in,celldex::MonacoImmuneData(),"label.fine")
    so_in$immu_cell_exp_main <-  RUN_SINGLEr(so_in,celldex::DatabaseImmuneCellExpressionData(),
                                            "label.main")
    so_in$immu_cell_exp <- RUN_SINGLEr(so_in,celldex::DatabaseImmuneCellExpressionData(),
                                      "label.fine")
    so_in$annot = so_in$HPCA_main
  } else if(ref == "mm10"){
    so_in$immgen_main <-  RUN_SINGLEr(so_in,celldex::ImmGenData(),"label.main")
    so_in$immgen <- RUN_SINGLEr(so_in,celldex::ImmGenData(),"label.fine")
    so_in$mouseRNAseq_main <-  RUN_SINGLEr(so_in,celldex::MouseRNAseqData(),"label.main")
    so_in$mouseRNAseq <- RUN_SINGLEr(so_in,celldex::MouseRNAseqData(),"label.fine")
    so_in$annot = so_in$immgen_main
  }
  return(so_in)
}

##################################################################
#
##################################################################
PROCESS_DOUBLETS <-function(dfso){
  dfso <- SCTransform(dfso)
  dfso <- RunPCA(dfso, pc.genes = dfso@var.genes, pcs.print = 0,verbose = F,npcs =10)
  npcs = 10
  dfso <- RunUMAP(dfso, verbose=TRUE,dims = 1:npcs)
  
  sweep.res.list_kidney <- paramSweep_v3(dfso,PCs = 1:10, sct = T)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## Homotypic Doublet Proportion Estimate 
  homotypic.prop <- modelHomotypic(dfso$annot)
  perc = 0.005 * (length(colnames(dfso))/1000)
  nExp_poi <- round(perc*length(colnames(dfso)))#dfso@cell.names
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies 
  dfso <- doubletFinder_v3(dfso, 
                           pN = 0.25, pK = 0.09, 
                           nExp = nExp_poi, 
                           reuse.pANN = FALSE,PCs = 1:10,sct = T)
  
  pAAN=tail(names(dfso@meta.data),2)[1]
  dfso <- doubletFinder_v3(dfso, 
                           pN = 0.25, pK = 0.09, 
                           nExp = nExp_poi.adj, 
                           reuse.pANN = pAAN,PCs = 1:10,sct = T)
  return(dfso)
}

##################################################################
#
##################################################################
RUN_SINGLEr_seurat = function(obj,refFile,fineORmain){
  # find av expressions
  avg = AverageExpression(obj,assays = "SCT")
  avg = as.data.frame(avg)
  
  s = SingleR(test = as.matrix(avg),ref = refFile,labels = refFile[[fineORmain]])
  
  clustAnnot = s$labels
  names(clustAnnot) = colnames(avg)
  names(clustAnnot) = gsub("SCT.","",names(clustAnnot))
  
  obj$clustAnnot = clustAnnot[match(obj$seurat_clusters,names(clustAnnot))]
  return(obj$clustAnnot)
}

RUN_CORRECTION_SEURAT = function(obj,npcs,ref,resolution){
  # what is this catching? why another scale?
  if(ncol(obj@assays$integrated@scale.data)!=ncol(obj@assays$integrated@data)){
    obj = ScaleData(obj)
  }
  
  # run PCA, FindNeighbors
  obj <- RunPCA(object = obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj,dims = 1:npcs)
  
  # for each input resolution, find clusters
  for (res in resolution){
    obj <- FindClusters(obj,dims = 1:npcs, 
                        print.output = 0, resolution = res,algorithm = 3)
  }
  
  # rename cols
  colnames(obj@meta.data) = gsub("integrated_snn_res","SLM_int_snn_res",colnames(obj@meta.data))
  colnames(obj@meta.data) = gsub("SCT_snn_res","SLM_SCT_snn_res",colnames(obj@meta.data))
  
  # run UMAP
  obj <- RunUMAP(object = obj, reduction = "pca",
                 dims = 1:npcs,n.components = 3)
  
  # this is essentially the same labeling as was done in preprocessing, except now you're adding
  # "clustAnnot" in the front - is it ncessary? can we keep the obj$HPCA_main, or is the prefix needed
  # for all of these?
  if(ref == "hg38" || ref == "hg19"){
    obj$clustAnnot_HPCA_main <- RUN_CORRECTION_SEURAT(obj,celldex::HumanPrimaryCellAtlasData(),"label.main")
    obj$clustAnnot_HPCA <-  RUN_CORRECTION_SEURAT(obj,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    obj$clustAnnot_BP_encode_main <-  RUN_CORRECTION_SEURAT(obj,celldex::BlueprintEncodeData(),"label.main")
    obj$clustAnnot_BP_encode <-  RUN_CORRECTION_SEURAT(obj,celldex::BlueprintEncodeData(),"label.fine")
    obj$clustAnnot_monaco_main <-  RUN_CORRECTION_SEURAT(obj,celldex::MonacoImmuneData(),"label.main")
    obj$clustAnnot_monaco <- RUN_CORRECTION_SEURAT(obj,celldex::MonacoImmuneData(),"label.fine")
    obj$clustAnnot_immu_cell_exp_main <-  RUN_CORRECTION_SEURAT(obj,celldex::DatabaseImmuneCellExpressionData(),"label.main")
    obj$clustAnnot_immu_cell_exp <- RUN_CORRECTION_SEURAT(obj,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
  } else if(ref == "mm10"){
    obj$clustAnnot_immgen_main <-  RUN_CORRECTION_SEURAT(obj,celldex::ImmGenData(),"label.main")
    obj$clustAnnot_immgen <- RUN_CORRECTION_SEURAT(obj,celldex::ImmGenData(),"label.fine")
    obj$clustAnnot_mouseRNAseq_main <-  RUN_CORRECTION_SEURAT(obj,celldex::MouseRNAseqData(),"label.main")
    obj$clustAnnot_mouseRNAseq <- RUN_CORRECTION_SEURAT(obj,celldex::MouseRNAseqData(),"label.fine")
  }
  return(obj)
}

##################################################################
#
##################################################################
RUN_SINGLEr_AVERAGE = function(obj,refFile,fineORmain){
  avg = AverageExpression(obj,assays = "SCT")
  avg = as.data.frame(avg)
  ref = refFile
  s = SingleR(test = as.matrix(avg),ref = ref,labels = ref[[fineORmain]])
  
  clustAnnot = s$labels
  names(clustAnnot) = colnames(avg)
  names(clustAnnot) = gsub("SCT.","",names(clustAnnot))
  
  obj$clustAnnot = clustAnnot[match(obj$seurat_clusters,names(clustAnnot))]
  return(obj$clustAnnot)
}

FIND_CLUSTERS_BY_RES<-function(resolution,obj,algorithm_val){
  for (res in resolution){
    obj <- FindClusters(obj,dims = 1:npcs, 
                        print.output = 0, resolution = res,algorithm = algorithm_val)
  }
  return(obj)
}

RUN_CORRECTION_HARMONY = function(obj,npcs,harm, resolution){
  obj <- RunPCA(object = obj, npcs = npcs, verbose = FALSE)
  
  if(harm=="Yes") {
    obj <- FindNeighbors(obj,dims = 1:npcs,reduction = "harmony")
  } else if(harm=="No") {
    obj <- FindNeighbors(obj,dims = 1:npcs)
  }
  
  # for each input resolution, find clusters
  FIND_CLUSTERS_BY_RES(obj,res,3)
  colnames(obj@meta.data) = gsub("integrated_snn_res","SLM_int_snn_res",colnames(obj@meta.data))
  colnames(obj@meta.data) = gsub("SCT_snn_res","SLM_SCT_snn_res",colnames(obj@meta.data))
  
  FIND_CLUSTERS_BY_RES(obj,res,4)
  colnames(obj@meta.data) = gsub("integrated_snn_res","Leiden_int_snn_res",colnames(obj@meta.data))
  colnames(obj@meta.data) = gsub("SCT_snn_res","Leiden_SCT_snn_res",colnames(obj@meta.data))
  
  if(harm=="Yes") {
    obj <- RunUMAP(object = obj,dims = 1:npcs, n.components = 3,reduction = "harmony")
  } else if(harm=="No") {
    obj <- RunUMAP(object = obj, reduction = "pca",dims = 1:npcs,n.components = 3)
  }

  # this is essentially the same labeling as was done in preprocessing, except now you're adding
  # "clustAnnot" in the front - is it ncessary? can we keep the obj$HPCA_main, or is the prefix needed
  # for all of these?
  if(ref == "hg38" || ref == "hg19"){
    obj$clustAnnot_HPCA_main <- RUN_SINGLEr_AVERAGE(obj,celldex::HumanPrimaryCellAtlasData(),"label.main")
    obj$clustAnnot_HPCA <-  RUN_SINGLEr_AVERAGE(obj,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    obj$clustAnnot_BP_encode_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::BlueprintEncodeData(),"label.main")
    obj$clustAnnot_BP_encode <-  RUN_SINGLEr_AVERAGE(obj,celldex::BlueprintEncodeData(),"label.fine")
    obj$clustAnnot_monaco_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::MonacoImmuneData(),"label.main")
    obj$clustAnnot_monaco <- RUN_SINGLEr_AVERAGE(obj,celldex::MonacoImmuneData(),"label.fine")
    obj$clustAnnot_immu_cell_exp_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::DatabaseImmuneCellExpressionData(),"label.main")
    obj$clustAnnot_immu_cell_exp <- RUN_SINGLEr_AVERAGE(obj,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
  }
  
  if(ref == "mm10"){
    obj$clustAnnot_immgen_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::ImmGenData(),"label.main")
    obj$clustAnnot_immgen <- RUN_SINGLEr_AVERAGE(obj,celldex::ImmGenData(),"label.fine")
    obj$clustAnnot_mouseRNAseq_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::MouseRNAseqData(),"label.main")
    obj$clustAnnot_mouseRNAseq <- RUN_SINGLEr_AVERAGE(obj,celldex::MouseRNAseqData(),"label.fine")
  }
  return(obj)
}

RUN_CORRECTION_RPCA = function(obj,npcs){
  obj <- RunPCA(object = obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj,dims = 1:npcs)
  
  FIND_CLUSTERS_BY_RES(obj,res,3)
  colnames(obj@meta.data) = gsub("integrated_snn_res","SLM_int_snn_res",colnames(obj@meta.data))
  colnames(obj@meta.data) = gsub("SCT_snn_res","SLM_SCT_snn_res",colnames(obj@meta.data))
  
  FIND_CLUSTERS_BY_RES(obj,res,4)
  colnames(obj@meta.data) = gsub("integrated_snn_res","Leiden_int_snn_res",colnames(obj@meta.data))
  colnames(obj@meta.data) = gsub("SCT_snn_res","Leiden_SCT_snn_res",colnames(obj@meta.data))
  
  obj <- RunUMAP(object = obj, reduction = "pca", 
                 dims = 1:npcs,n.components = 3)

  if(ref == "hg38" || ref == "hg19"){
    obj$clustAnnot_HPCA_main <- RUN_SINGLEr_AVERAGE(obj,celldex::HumanPrimaryCellAtlasData(),"label.main")
    obj$clustAnnot_HPCA <-  RUN_SINGLEr_AVERAGE(obj,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    obj$clustAnnot_BP_encode_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::BlueprintEncodeData(),"label.main")
    obj$clustAnnot_BP_encode <-  RUN_SINGLEr_AVERAGE(obj,celldex::BlueprintEncodeData(),"label.fine")
    obj$clustAnnot_monaco_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::MonacoImmuneData(),"label.main")
    obj$clustAnnot_monaco <- RUN_SINGLEr_AVERAGE(obj,celldex::MonacoImmuneData(),"label.fine")
    obj$clustAnnot_immu_cell_exp_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::DatabaseImmuneCellExpressionData(),"label.main")
    obj$clustAnnot_immu_cell_exp <- RUN_SINGLEr_AVERAGE(obj,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
  } else if(ref == "mm10"){
    obj$clustAnnot_immgen_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::ImmGenData(),"label.main")
    obj$clustAnnot_immgen <- RUN_SINGLEr_AVERAGE(obj,celldex::ImmGenData(),"label.fine")
    obj$clustAnnot_mouseRNAseq_main <-  RUN_SINGLEr_AVERAGE(obj,celldex::MouseRNAseqData(),"label.main")
    obj$clustAnnot_mouseRNAseq <- RUN_SINGLEr_AVERAGE(obj,celldex::MouseRNAseqData(),"label.fine")
  }
  return(obj)
}