library(Seurat)
library(dplyr)
library(patchwork)

# Loading quality controlled embryo cell data

load("data/experimentation/main/all_embryo_cell_QC.rda")

# Plotting JackStraw and Elbow Plots

JackStrawPlot(all_embryo_cell_QC, dims = 1:20)
ElbowPlot(all_embryo_cell_QC)

# Including 1st 15 dimensions in PCA using KNN

all_embryo_cell_QC_15N <- FindNeighbors(all_embryo_cell_QC, dims = 1:15)

# Container for umaps of different resolutions

umaps_res <- c();

# Resolutions to be tested

resolutions = c(.6, 1.0, 1.4, 1.8)

# Clustering same umap with different resolutions and storing in umap_res

for (i in c(1:4)){
  umaps_res <- c(umaps_res, FindClusters(all_embryo_cell_QC_15N, resolution = resolutions[i]))
}

# Processing umaps in umap_res

umaps_processed = c();
for (umap in umaps_res){
  umaps_processed <- c(umaps_processed, RunUMAP(umap, dims = 1:15, seed.use=42L, n.neighbors = 30L))
}

# Plotting processed umaps

p <- DimPlot(umaps_processed[[1]], reduction = "umap", pt.size = .5) + ggtitle(resolutions[1])
for (i in c(2:length(umaps_processed))){
  p <- p + (DimPlot(umaps_processed[[i]], reduction = "umap", pt.size = .5) + ggtitle(resolutions[i]))
}

p

# Saving processed umaps with resolutions .6 to 1.8

save(umaps_processed, file="data/experimentation/umap_parameter_adjustment/umaps_res_.6-1.8.rda")
