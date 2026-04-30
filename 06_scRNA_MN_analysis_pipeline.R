# ==============================================================================
# Single-Cell RNA-seq Analysis Pipeline for Motor Neuron Characterization
# ==============================================================================
# 
# Description:
#   This script performs a comprehensive single-cell RNA-seq (scRNA-seq) 
#   analysis workflow to characterize motor neuron subpopulations. 
#   The pipeline includes:
#     - Data preprocessing and Quality Control (QC)
#     - Normalization and highly variable feature selection
#     - Dimensionality reduction (PCA, UMAP, t-SNE)
#     - Clustering and cell-type identification
#     - Dataset integration (Control vs. Training conditions)
#     - Motor neuron (MN) subpopulation extraction
#     - Cell-type annotation using SingleR and marker genes
#
# Study:
#   "Agility Training Enhances Motor Temporal Precision by Reweighting
#    Spinal Phase-Locked Commissural Inhibition"
#
# Dependencies:
#   Seurat, dplyr, ggplot2, ggpubr, SingleR, scater,
#   SingleCellExperiment, patchwork, tidyverse, textshape, etc.
#
# Author:
#   Wei Liu, 2026
#
# Notes:
#   - Raw data: 10X Genomics format
#   - Reference dataset used for annotation: GSE243993_MN.rds
#   - Integration performed using Seurat anchors
# ==============================================================================

rm(list = ls()); gc()

# Install and load necessary packages
if(!require(Seurat)) install.packages("Seurat")
if(!require(dplyr)) install.packages("dplyr")
if(!require(patchwork)) install.packages("patchwork")
if(!require(R.utils)) install.packages("R.utils")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(devtools)) install.packages("devtools")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(Hmisc)) install.packages("Hmisc")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(cowplot)) install.packages("cowplot")
if(!require(matrixStats)) install.packages("matrixStats")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(textshape)) install.packages('textshape')
library(textshape)
if(!require(scater)) BiocManager::install('scater')
library(scater)
if(!require(SingleCellExperiment)) BiocManager::install('SingleCellExperiment')
library(dplyr)
library(SingleR)

# Load reference dataset
myref <- readRDS("GSE243993_MN.rds")
# Store current labels in celltype
myref$celltype <- Idents(myref)
table(Idents(myref))

# Read in reference dataset -------
Refassay <- log1p(AverageExpression(myref, verbose = FALSE)$RNA) # Calculate log1p average expression
# Ref <- textshape::column_to_rownames(Ref, loc = 1) # Alternative method to get reference matrix
head(Refassay) # Check the expression matrix structure

# Reference dataset needs to be constructed into a SingleCellExperiment object
ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=Refassay))
ref_sce = scater::logNormCounts(ref_sce)
logcounts(ref_sce)[1:4,1:4]
colData(ref_sce)$Type = colnames(Refassay)
ref_sce # Construction complete

# ==============================================================================
# HUC Control Data Processing
# ==============================================================================
HUCControl.data <- Read10X(data.dir = "matrix") 
HUCControl <- CreateSeuratObject(counts = HUCControl.data, project = "Control", assay = "RNA", min.cells = 3, min.features = 200)
HUCControl[["percent.mt"]] <- PercentageFeatureSet(HUCControl, pattern = "mt-")

plot1 <- FeatureScatter(HUCControl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HUCControl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(HUCControl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident') 
HUCControl <- subset(HUCControl, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 ) 

# Normalize expression data
HUCControl <- NormalizeData(HUCControl, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes
HUCControl <- FindVariableFeatures(HUCControl, selection.method = "vst", nfeatures = 2000)

# View the top 10 most highly variable genes
top10 <- head(VariableFeatures(HUCControl), 10) 

plot1 <- VariableFeaturePlot(HUCControl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2 # Plot variable features with and without labels

# Centering and scaling data matrix
HUCControl <- ScaleData(HUCControl, features = rownames(HUCControl))

HUCControl <- RunPCA(HUCControl, features = VariableFeatures(object = HUCControl))
print(HUCControl[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HUCControl, dims = 1:2, reduction = "pca")
DimPlot(HUCControl, reduction = "pca")
DimHeatmap(HUCControl, dims = 1:15, cells = 500, balanced = TRUE)

HUCControl <- JackStraw(HUCControl, num.replicate = 100)
HUCControl <- ScoreJackStraw(HUCControl, dims = 1:20)
JackStrawPlot(HUCControl, dims = 1:20)
ElbowPlot(HUCControl) 

HUCControl <- FindNeighbors(HUCControl, dims = 1:16)
HUCControl <- FindClusters(HUCControl, resolution = 0.6)
HUCControl <- RunUMAP(HUCControl, dims = 1:16)
HUCControl <- RunTSNE(HUCControl, dims = 1:16)
DimPlot(object = HUCControl, label = T)

# Subset clusters with high expression of snap25a, elavl4, elavl3, rbfox1
HUCControlNeuron <- subset(x = HUCControl, idents=c('Neuron')) 
HUCControlNeuron <- NormalizeData(HUCControlNeuron, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes for Neuron subset
HUCControlNeuron <- FindVariableFeatures(HUCControlNeuron, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HUCControlNeuron), 10) # View top 10 variable genes

plot1 <- VariableFeaturePlot(HUCControlNeuron)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 # Plot variable features with and without labels

HUCControlNeuron <- ScaleData(HUCControlNeuron, features = rownames(HUCControlNeuron))
HUCControlNeuron <- RunPCA(HUCControlNeuron, features = VariableFeatures(object = HUCControlNeuron))
print(HUCControlNeuron[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HUCControlNeuron, dims = 1:2, reduction = "pca")
DimPlot(HUCControlNeuron, reduction = "pca")
DimHeatmap(HUCControlNeuron, dims = 1:15, cells = 500, balanced = TRUE)

HUCControlNeuron <- JackStraw(HUCControlNeuron, num.replicate = 100)
HUCControlNeuron <- ScoreJackStraw(HUCControlNeuron, dims = 1:20)
JackStrawPlot(HUCControlNeuron, dims = 1:20)
ElbowPlot(HUCControlNeuron) 

HUCControlNeuron <- FindNeighbors(HUCControlNeuron, dims = 1:20)
HUCControlNeuron <- FindClusters(HUCControlNeuron, resolution = 0.6)
HUCControlNeuron <- RunUMAP(HUCControlNeuron, dims = 1:20)
DimPlot(object = HUCControlNeuron, label = T)

# ==============================================================================
# HUC Training Data Processing
# ==============================================================================
HUCTraining.data <- Read10X(data.dir = "matrix") 
HUCTraining <- CreateSeuratObject(counts = HUCTraining.data, project = "Training", assay = "RNA", min.cells = 3, min.features = 200)
HUCTraining[["percent.mt"]] <- PercentageFeatureSet(HUCTraining, pattern = "mt-")

plot1 <- FeatureScatter(HUCTraining, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HUCTraining, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(HUCTraining, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident') 
HUCTraining <- subset(HUCTraining, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 ) 

# Normalize expression data
HUCTraining <- NormalizeData(HUCTraining, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes
HUCTraining <- FindVariableFeatures(HUCTraining, selection.method = "vst", nfeatures = 2000)

# View the top 10 most highly variable genes
top10 <- head(VariableFeatures(HUCTraining), 10) 

plot1 <- VariableFeaturePlot(HUCTraining)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
## When using repel, set xnudge and ynudge to 0 for optimal results
plot1 + plot2 # Plot variable features with and without labels

# Centering and scaling data matrix
HUCTraining <- ScaleData(HUCTraining, features = rownames(HUCTraining))

HUCTraining <- RunPCA(HUCTraining, features = VariableFeatures(object = HUCTraining))
print(HUCTraining[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HUCTraining, dims = 1:2, reduction = "pca")
DimPlot(HUCTraining, reduction = "pca")
DimHeatmap(HUCTraining, dims = 1:15, cells = 500, balanced = TRUE)

HUCTraining <- JackStraw(HUCTraining, num.replicate = 100)
HUCTraining <- ScoreJackStraw(HUCTraining, dims = 1:20)
JackStrawPlot(HUCTraining, dims = 1:20)
ElbowPlot(HUCTraining) 

HUCTraining <- FindNeighbors(HUCTraining, dims = 1:16)
HUCTraining <- FindClusters(HUCTraining, resolution = 0.6)
HUCTraining <- RunUMAP(HUCTraining, dims = 1:16)
HUCTraining <- RunTSNE(HUCTraining, dims = 1:16)
DimPlot(object = HUCTraining, label = T)

# Subset clusters with high expression of snap25a, elavl4, elavl3, rbfox1
HUCTrainingNeuron <- subset(x = HUCTrainingNeuron, idents=c('Neuron')) 
HUCTrainingNeuron <- NormalizeData(HUCTrainingNeuron, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes for Training Neuron subset
HUCTrainingNeuron <- FindVariableFeatures(HUCTrainingNeuron, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HUCTrainingNeuron), 10) # View top 10 variable genes

plot1 <- VariableFeaturePlot(HUCTrainingNeuron)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 # Plot variable features with and without labels

HUCTrainingNeuron <- ScaleData(HUCTrainingNeuron, features = rownames(HUCTrainingNeuron))
HUCTrainingNeuron <- RunPCA(HUCTrainingNeuron, features = VariableFeatures(object = HUCTrainingNeuron))
print(HUCTrainingNeuron[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HUCTrainingNeuron, dims = 1:2, reduction = "pca")
DimPlot(HUCTrainingNeuron, reduction = "pca")
DimHeatmap(HUCTrainingNeuron, dims = 1:15, cells = 500, balanced = TRUE)

HUCTrainingNeuron <- JackStraw(HUCTrainingNeuron, num.replicate = 100)
HUCTrainingNeuron <- ScoreJackStraw(HUCTrainingNeuron, dims = 1:20)
JackStrawPlot(HUCTrainingNeuron, dims = 1:20)
ElbowPlot(HUCTrainingNeuron) 

HUCTrainingNeuron <- FindNeighbors(HUCTrainingNeuron, dims = 1:19)
HUCTrainingNeuron <- FindClusters(HUCTrainingNeuron, resolution = 0.6)
HUCTrainingNeuron <- RunUMAP(HUCTrainingNeuron, dims = 1:19)
DimPlot(object = HUCTrainingNeuron, label = T)

# ==============================================================================
# Merge and Integrate Datasets
# ==============================================================================
cord.grouped <-
  merge(HUCControlNeuron,
        y = c(HUCTrainingNeuron),
        add.cell.ids = c(
          'C Neuron',
          'T Neuron'
        ),
        project = "HUC Neuron"
  )
head(colnames(cord.grouped))
tail(colnames(cord.grouped))
unique(sapply(X = strsplit(colnames(cord.grouped), split = "_"), FUN = "[", 1))
table(cord.grouped$orig.ident)

combinedgroup.list <- SplitObject(cord.grouped, split.by = "orig.ident")
combinedgroup.list <- lapply(X = combinedgroup.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

combine.anchors <- FindIntegrationAnchors(object.list = combinedgroup.list, dims = 1:20)
combine.combined <- IntegrateData(anchorset = combine.anchors, dims = 1:20)
DefaultAssay(combine.combined) <- "integrated"

combine.combined <- ScaleData(combine.combined, verbose = FALSE)
combine.combined <- RunPCA(combine.combined, npcs = 30, verbose = FALSE)
combine.combined <- RunUMAP(combine.combined, reduction = "pca", dims = 1:20)
combine.combined <- FindNeighbors(combine.combined, reduction = "pca", dims = 1:20)
combine.combined <- FindClusters(combine.combined, resolution = 1)
DefaultAssay(combine.combined) <- "RNA"
DimPlot(object = combine.combined, label = T)

# ==============================================================================
# Motor Neuron (MN) Analysis
# ==============================================================================
# Subset clusters with high expression of chata, slc18a3a
MN <- subset(x = combine.combined, idents=c('MNs')) 
MN <- NormalizeData(MN, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes for MN subset
MN <- FindVariableFeatures(MN, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(MN), 10) # View top 10 variable genes

plot1 <- VariableFeaturePlot(MN)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 # Plot variable features with and without labels

MN <- ScaleData(MN, features = rownames(MN))
MN <- RunPCA(MN, features = VariableFeatures(object = MN))
print(MN[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MN, dims = 1:2, reduction = "pca")
DimPlot(MN, reduction = "pca")
DimHeatmap(MN, dims = 1:15, cells = 500, balanced = TRUE)

MN <- JackStraw(MN, num.replicate = 100)
MN <- ScoreJackStraw(MN, dims = 1:20)
JackStrawPlot(MN, dims = 1:20)
ElbowPlot(MN) 

MN <- FindNeighbors(MN, dims = 1:20)
MN <- FindClusters(MN, resolution = 0.4)
MN <- RunUMAP(MN, dims = 1:20, min.dist = 1)
DimPlot(object = MN, label = T, pt.size = 2)

# ==============================================================================
# Preliminary Cell Annotation using SingleR
# ==============================================================================
bfreaname.pbmc <- MN
testdata <- GetAssayData(bfreaname.pbmc)
pred <- SingleR(test=testdata, ref=ref_sce, 
                labels=ref_sce$Type)
head(pred) 
table(pred$labels)
as.data.frame(table(pred$labels))

cellType=data.frame(seurat=bfreaname.pbmc@meta.data$seurat_clusters,
                    predict=pred$labels)
sort(table(cellType[,1]))
table(cellType[,1:2])

lalala <- as.data.frame(table(cellType[,1:2]))
finalmap <- lalala %>% group_by(seurat) %>% top_n(n = 1, wt = Freq)
finalmap <- finalmap[order(finalmap$seurat),]$predict
print(finalmap)

MN <- UpdateSeuratObject(object = bfreaname.pbmc)
new.cluster.ids <- as.character(finalmap)
names(new.cluster.ids) <- levels(MN)
MN <- RenameIdents(MN, new.cluster.ids)
DimPlot(object = MN, split.by ="orig.ident", label = T) + NoLegend()

# ==============================================================================
# Further Annotation Referencing Pallucchi et al., Nat Neurosci. 2024
# ==============================================================================
VlnPlot(object = MN, features = c('itga3a','bmp16','ret','pcdh11','oxria','grik2','
                                 grintb','esrrga'), slot = 'scale.data')
VlnPlot(object = MN, features = c('tmprss4b','scn4bb','pvalbo','atpla3b','apooa',
                                  'glrba','glrba','stmn3'), slot = 'scale.data')
VlnPlot(object = MN, features = c('chrna2b','calb1','fndc4a','pcdh10a','sox2',
                                  'necab2','ramp1','calb2b','pcp4a'), slot = 'scale.data')
VlnPlot(object = MN, features = c('nes','pcna','mki67','tubb5','tuba1a',
                                  'hnrnpa0l','rps6','nhlh2','olig4','ebf2'
                                  ,'tuba8l5','elavl3','rassf1','nkx6.1'), slot = 'scale.data')

new.cluster.ids <- c('Slow MNs','Immature MNs','Fast MNs','Slow MNs')
names(new.cluster.ids) <- levels(MN)
MN1 <- RenameIdents(MN, new.cluster.ids)

MN1$celltype <- Idents(MN1)
MN1$group <- MN1$orig.ident
MN1 <- ScaleData(MN1, features = rownames(MN1))

Idents(MN1) <- factor(Idents(MN1), levels = c('Slow MNs','Fast MNs','Immature MNs'))
DimPlot(object = MN1, label = T, pt.size = 3, cols = c('#CE0800','#2222A8',"grey"), split.by = 'orig.ident')
saveRDS(object = MN1, file = 'scobj.rds')