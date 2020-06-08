library(Seurat)
library(dplyr)
library(patchwork)

# Loading subsetted embryo cell counts and metadata

load("data/experimentation/main/all_embryo_counts.rda")

# Initialize the Seurat object with the raw (non-normalized data). 

all_embryo_cell <- CreateSeuratObject(counts = all_embryo_cell_counts, project = "all_embryo_cell", min.cells = 5, min.features = 5)

# Add meta data to the Seurat object

all_embryo_cell <- AddMetaData(object = all_embryo_cell, all_embryo_meta)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

all_embryo_cell[["percent.mt"]] <- PercentageFeatureSet(all_embryo_cell, pattern = "^MT-")

# Visualize QC metrics as a violin plot

VlnPlot(all_embryo_cell, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all_embryo_cell, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_embryo_cell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Histogram distribution: frequency (number of cells) and nFeature, 
# set up threshold (500, 7200), if it looks Gaussian, take log2

hist(all_embryo_cell@meta.data$nFeature_RNA,100)

# Histogram distribution: frequency (number of cells) and nCount, set up threashold 
# in log2, cut the left tail at 17, no need to set up max since no right tail
# calculate 2^17 = 131072, it is the counts we need to delete

hist(all_embryo_cell@meta.data$nCount_RNA,100)
hist(log2(all_embryo_cell@meta.data$nCount_RNA),100)
2^17

# QC

all_embryo_cell_QC <- subset(all_embryo_cell, subset = 500 < nFeature_RNA & nFeature_RNA < 7200 & nCount_RNA > 131072)
dim(all_embryo_cell_QC@meta.data)

# Quickly visualize data distribution as a violin plot

VlnPlot(all_embryo_cell_QC, features = c("total_counts", "num_detect_genes", "percent_ribo"), ncol = 3)
FeatureScatter(all_embryo_cell_QC, feature1 = "total_counts", feature2 = "num_detect_genes")

# normalizing the data. Default is log, 10000X)

all_embryo_cell_QC <- NormalizeData(all_embryo_cell_QC, normalization.method = "LogNormalize", scale.factor = 10000)

# identification of highly variable features (feature selection)

all_embryo_cell_QC <- FindVariableFeatures(all_embryo_cell_QC, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(all_embryo_cell_QC), 10)
plot1 <- VariableFeaturePlot(all_embryo_cell_QC)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

# Scaling the data. Shifts the expression of each gene, so that the mean expression across cells is 0. 
# Scales the expression of each gene, so that the variance across cells is 1. 
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# May check SCTransform as a new method

all.genes <- rownames(all_embryo_cell_QC)
all_embryo_cell_QC <- ScaleData(all_embryo_cell_QC, features = all.genes)
all_embryo_cell_QC <- ScaleData(all_embryo_cell_QC, vars.to.regress = "total_counts")


# Perform linear dimension reduction (PCA)

all_embryo_cell_QC <- RunPCA(all_embryo_cell_QC, features = VariableFeatures(object = all_embryo_cell_QC))
DimHeatmap(all_embryo_cell_QC, dims = 1:20, cells = 500, balanced = TRUE)

# Determine the "dimensionality" of the dataset (choosing components)
# checking PCA genes, run statistical test, elbow plots.
# advise users to err on the higher side when choosing this parameter.

all_embryo_cell_QC <- JackStraw(all_embryo_cell_QC, num.replicate = 100)
all_embryo_cell_QC <- ScoreJackStraw(all_embryo_cell_QC, dims = 1:20)

save(all_embryo_cell_QC, file="data/experimentation/main/all_embryo_cell_QC.rda")