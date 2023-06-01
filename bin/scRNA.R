.libPaths( c("/data/CCBR_Pipeliner/db/PipeDB/scrna5",.libPaths()))


library(Seurat)
library(stringr)
library("DoubletFinder")
library(SingleR)
library(scRNAseq)
library(SingleCellExperiment)
library(celldex)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(djvdj)
library(miQC)
library(Routliers)
library(flexmix)
library(SeuratWrappers)

source("source.R")


options(future.globals.maxSize = 1e12)
#options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = T)


CRfolder =  as.character(args[1]) 
h5 = as.character(paste0("cellrangerOut/",CRfolder,"/outs/filtered_feature_bc_matrix.h5") )
ref =  as.character(args[2])  
outFile = as.character(args[3])
qcDefault = as.character(args[4]) #manual/miqc/mads
manualParms = as.vector(args[5])
doubletFinder = as.character(args[6])

manualParms = strsplit(manualParms,",")
print(manualParms)

rnaCounts = Read10X_h5(h5)

if (class(rnaCounts) == "list") {
rnaCounts = rnaCounts$'Gene Expression'
}

so <- CreateSeuratObject(rnaCounts)



groupFile = read.delim("groups.tab",header=F,stringsAsFactors = F)
sample = groupFile$V3[groupFile$V1 == CRfolder & groupFile$V4 == "gex"]

so$Sample = sample
#so$groups =  groupFile$V2[groupFile$V3 == sample & groupFile$V4 == "gex"]
print(head(so$Sample))

if (nrow(groupFile[groupFile$V3 == sample & groupFile$V4 == "vdj",]) > 0 ) {

tcrSamples = groupFile$V1[groupFile$V3 == sample & groupFile$V4 == "vdj"]
tcrSamples = paste0("cellrangerOut/",tcrSamples,"/outs")

so = import_vdj(input = so, vdj_dir = tcrSamples,  filter_paired = FALSE  )

}




if(ref=="hg38"){so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")}
if(ref=="mm10"){so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")}

if(ref=="hg38"){so[["percent.RP"]] <- PercentageFeatureSet(so, pattern = "RPL|RPS")}
if(ref=="mm10"){so[["percent.RP"]] <- PercentageFeatureSet(so, pattern = "Rpl|Rps")}



so <- subset(so, subset = nFeature_RNA > 200)

so <- RunMiQC(so, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.7,model.slot = "flexmix_model") 
nCount_out = outliers_mad(so$nCount_RNA,threshold = 3)$LL_CI_MAD#
nFeature_out = outliers_mad(so$nFeature_RNA,threshold = 3)$LL_CI_MAD
mt_out = outliers_mad(so$percent.mt,threshold = 3)$UL_CI_MAD


if(qcDefault == "manual") {so =  subset(so,nCount_RNA > manualParms[[1]][1] & nCount_RNA < manualParms[[1]][2] & nFeature_RNA > manualParms[[1]][3] & nFeature_RNA < manualParms[[1]][4] & percent.mt > manualParms[[1]][5] & percent.mt < manualParms[[1]][6])}
if(qcDefault == "miqc") {so = subset(so, miQC.keep == "keep")}
if(qcDefault == "mads") {so = subset(so, subset = nFeature_RNA > nFeature_out & nCount_RNA > nFeature_out & percent.mt < mt_out)}


if(ref=="hg38"){
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
}

if(ref=="mm10"){
s.genes <- convertHumanGeneList(cc.genes$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)
}

so = seuratSample(so,30)

if(ref == "hg38"){
so$HPCA_main <- runSingleRCell(so,celldex::HumanPrimaryCellAtlasData(),"label.main")
so$HPCA <-  runSingleRCell(so,celldex::HumanPrimaryCellAtlasData(),"label.fine")
so$BP_encode_main <-  runSingleRCell(so,celldex::BlueprintEncodeData(),"label.main")
so$BP_encode <-  runSingleRCell(so,celldex::BlueprintEncodeData(),"label.fine")
so$monaco_main <-  runSingleRCell(so,celldex::MonacoImmuneData(),"label.main")
so$monaco <-     runSingleRCell(so,celldex::MonacoImmuneData(),"label.fine")
so$immu_cell_exp_main <-  runSingleRCell(so,celldex::DatabaseImmuneCellExpressionData(),"label.main")
so$immu_cell_exp <- runSingleRCell(so,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
so$annot = so$HPCA_main
}

if(ref == "mm10"){

so$immgen_main <-  runSingleRCell(so,celldex::ImmGenData(),"label.main")
so$immgen <- runSingleRCell(so,celldex::ImmGenData(),"label.fine")
so$mouseRNAseq_main <-  runSingleRCell(so,celldex::MouseRNAseqData(),"label.main")
so$mouseRNAseq <- runSingleRCell(so,celldex::MouseRNAseqData(),"label.fine")

so$annot = so$immgen_main

}


if(doubletFinder == "Y"){
dfso = doublets(so)
so$DF_hi.lo = dfso[[tail(names(dfso@meta.data),1)]]
so=subset(so,cells=names(so$DF_hi.lo)[so$DF_hi.lo =="Singlet"])
}



so = UpdateSeuratObject(so)
saveRDS(so,outFile)
