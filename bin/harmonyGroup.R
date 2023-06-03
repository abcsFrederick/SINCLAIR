.libPaths( c("/data/CCBR_Pipeliner/db/PipeDB/scrna5",.libPaths()))


library(Seurat)
library(stringr)
library(reticulate)
library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library(celldex)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(SeuratWrappers)

source("source.R")
args <- commandArgs(trailingOnly = T)


so <- as.character(args[1])
outDirSeurat = as.character(args[2])
ref = as.character(args[3])
pcs = as.character(args[4])
resolutions = as.character(strsplit(gsub(",+",",",as.character(args[5])),split=",")[[1]])
resolutions = as.numeric(resolutions)
so = readRDS(so)



DefaultAssay(so) = "RNA"


so[["RNA"]] <- as(so[["RNA"]], Class = "Assay5")
so[["RNA"]] <- split(so[["RNA"]], f = so$Sample)

so <- SCTransform(so)
so <- RunPCA(so)


so <- IntegrateLayers(object = so, method = HarmonyIntegration, normalization.method = "SCT", verbose = F,new.reduction = "harmony")

so = runClustIntegrated(so,pcs,resolutions,"harmony")

if(ref == "hg38"){
    so$clustAnnot_HPCA_main <- singleRClusters(so,celldex::HumanPrimaryCellAtlasData(),"label.main")
    so$clustAnnot_HPCA <-  singleRClusters(so,celldex::HumanPrimaryCellAtlasData(),"label.fine")
    so$clustAnnot_BP_encode_main <-  singleRClusters(so,celldex::BlueprintEncodeData(),"label.main")
    so$clustAnnot_BP_encode <-  singleRClusters(so,celldex::BlueprintEncodeData(),"label.fine")
    so$clustAnnot_monaco_main <-  singleRClusters(so,celldex::MonacoImmuneData(),"label.main")
    so$clustAnnot_monaco <- singleRClusters(so,celldex::MonacoImmuneData(),"label.fine")
    so$clustAnnot_immu_cell_exp_main <-  singleRClusters(so,celldex::DatabaseImmuneCellExpressionData(),"label.main")
    so$clustAnnot_immu_cell_exp <- singleRClusters(so,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
  }

if(ref == "mm10"){

    so$clustAnnot_immgen_main <-  singleRClusters(so,celldex::ImmGenData(),"label.main")
    so$clustAnnot_immgen <- singleRClusters(so,celldex::ImmGenData(),"label.fine")
    so$clustAnnot_mouseRNAseq_main <-  singleRClusters(so,celldex::MouseRNAseqData(),"label.main")
    so$clustAnnot_mouseRNAseq <- singleRClusters(so,celldex::MouseRNAseqData(),"label.fine")

}
