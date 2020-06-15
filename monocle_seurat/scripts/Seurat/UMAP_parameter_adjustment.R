library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(seqtools)

# Loading quality controlled embryo cell data

load("data/experimentation/main/all_embryo_cell_QC.rda")

# Plotting JackStraw and Elbow Plots

JackStrawPlot(all_embryo_cell_QC, dims = 1:100)
ElbowPlot(all_embryo_cell_QC)

# Container for umaps of different parameters

trials <- c()

# ============================= DIMENSION TESTING ===============================


dimensions <- seq(10, 35, by = 5)
num_trials = length(dimensions)

# Including 1st 15 dimensions in PCA using KNN

for (dim in dimensions){
  seurat_dim <- FindNeighbors(all_embryo_cell_QC, dims = 1:dim)
  trials <- c(trials, seurat_dim)
}

for (i in c(1:num_trials)){
  trials[[i]] = FindClusters(trials[[i]], resolution = 1.4)
}
# ============================= RESOLUTION TESTING ===============================

# Resolutions to be tested

resolutions = seq(.5, 2.0, by=.3)

# Clustering same umap with different resolutions and storing in umap_res

for (i in c(1:length(umaps))){
  umaps <- c(umaps, FindNeighbors(all_embryo_cell_QC, dims = 1:15))
}


for (i in c(1:length(umaps))){
  umaps[[i]] <- FindClusters(umaps[[i]], resolution = resolutions[[i]])
}

# ================================================================================

# Processing umaps

for (i in c(1:num_trials)){
  trials[[i]] = RunUMAP(trials[[i]], dims = 1:dimensions[[i]], seed.use=42L, n.neighbors = 30L)
}

# Plotting processed umaps

p <- DimPlot(trials[[1]], reduction = "umap", pt.size = .5) + ggtitle(dimensions[1])
for (i in c(2:num_trials)){
  p <- p + (DimPlot(trials[[i]], reduction = "umap", pt.size = .5) + ggtitle(dimensions[i]))
}

p

for (i in c(1:num_trials)){
  trials[[i]][["ClusterIdents"]] <- Idents(object = trials[[i]])
  trials[[i]]@assays$RNA = NULL
}

# Saving processed umaps with resolutions .6 to 1.8

save(trials, data_RNA, file="data/experimentation/umap_parameter_adjustment/trial_dims_10-35.rda")
