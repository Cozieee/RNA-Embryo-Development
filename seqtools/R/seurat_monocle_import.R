
seurat_to_cds <- function (seurat) {
  dataUMAP <- data_RNA@data

  cell_metadata_UMAP <- new('AnnotatedDataFrame', data = seurat@meta.data)

  fData_UMAP <- data.frame(gene_short_name = row.names(dataUMAP), row.names = row.names(dataUMAP))

  gene_metadata_UMAP <- new('AnnotatedDataFrame', data = fData_UMAP)
  gene_metadata_UMAP <- as(gene_metadata_UMAP, "data.frame")
  cell_metadata_UMAP <- as(cell_metadata_UMAP, "data.frame")

  #Construct monocle cds

  cds_from_seurat <- new_cell_data_set(dataUMAP,
                                       cell_metadata = cell_metadata_UMAP,
                                       gene_metadata = gene_metadata_UMAP)

  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)

  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

  ### Assign the cluster info

  list_cluster <- seurat@meta.data[[sprintf("ClusterIdents", 1.4, dim)]]

  names(list_cluster) <- data_RNA@data@Dimnames[[2]]

  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

  ### Could be a space-holder, but essentially fills out louvain parameters

  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

  ### Assign UMAP coordinate

  cds_from_seurat@int_colData@listData[["reducedDims"]][["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings

  ### Assign feature loading for downstream module analysis

  cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

  return (cds_from_seurat)
}

fill_rna <- function(object){

  ret <- object
  obj_class = class(object)

  if (obj_class == "Seurat"){
    ret@assays$RNA <- data_RNA
  } else if (obj_class == "cell_data_set") {
    ret@assays@data@listData["counts"] <- data_RNA
  }

  return (ret)
}


del_rna <- function(object){

  ret <- object
  obj_class = class(object)

  if (obj_class == "Seurat"){
    ret@assays$RNA <- NULL
  } else if (obj_class == "cell_data_set") {
    ret@assays@data@listData <- list()
  }

  return (ret)
}
