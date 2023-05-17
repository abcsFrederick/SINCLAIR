###Run Seurat Clustering 
seuratClustering = function(so){
  
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

doublets <-function(dfso){
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

convertHumanGeneList  <- function(gns){
  egs <- mapIds(org.Hs.eg.db, gns, "ENTREZID","SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Mus.musculus","Homo.sapiens")
  mapped$MUS <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  return(as.character(unlist(mapped$MUS )))
}

runSingleR = function(obj,refFile,fineORmain){
  obj = DietSeurat(obj, graphs = "umap")
  sce = as.SingleCellExperiment(obj,assay = "SCT")
  ref = refFile
  s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]])
  return(s$pruned.labels)
}

PROCESS_SO<-function(so_in){
  # assign genes depending on ref input
  if(ref=="hg38" || ref == "hg19"){
    print("--proccesing human data")
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  } else if(ref=="mm10"){
    print("--proccesing mouse data")
    s.genes <- convertHumanGeneList(cc.genes$s.genes)
    g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
  }
  
  # process
  so_1 = SCTransform(so_in)
  so_2 = NormalizeData(so_1,
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000,
                       assay = "RNA")
  so_3 = ScaleData(so_2,assay = "RNA")
  so_4 = CellCycleScoring(so_3,
                          s.features = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident = TRUE)
  so_out = seuratClustering(so_4)
  return(so_out)
}

LABEL_SO<-function(so_in){
  if(ref == "hg38" || ref == "hg19"){
    so_in$HPCA_main <- runSingleR(so_in,celldex::HumanPrimaryCellAtlasData(),"label.main")
    so_in$HPCA <-  runSingleR(so_in,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    so_in$BP_encode_main <-  runSingleR(so_in,celldex::BlueprintEncodeData(),"label.main")
    so_in$BP_encode <-  runSingleR(so_in,celldex::BlueprintEncodeData(),"label.fine")
    so_in$monaco_main <-  runSingleR(so_in,celldex::MonacoImmuneData(),"label.main")
    so_in$monaco <-     runSingleR(so_in,celldex::MonacoImmuneData(),"label.fine")
    so_in$immu_cell_exp_main <-  runSingleR(so_in,celldex::DatabaseImmuneCellExpressionData(),
                                            "label.main")
    so_in$immu_cell_exp <- runSingleR(so_in,celldex::DatabaseImmuneCellExpressionData(),
                                      "label.fine")
    so_in$annot = so_in$HPCA_main
  } else if(ref == "mm10"){
    so_in$immgen_main <-  runSingleR(so_in,celldex::ImmGenData(),"label.main")
    so_in$immgen <- runSingleR(so_in,celldex::ImmGenData(),"label.fine")
    so_in$mouseRNAseq_main <-  runSingleR(so_in,celldex::MouseRNAseqData(),"label.main")
    so_in$mouseRNAseq <- runSingleR(so_in,celldex::MouseRNAseqData(),"label.fine")
    so_in$annot = so_in$immgen_main
  }
  return(so_in)
}
