gen_cds <- function(dataRNA, metaRNA) {
  dataUMAP <- dataRNA@data

  cell_metadata_UMAP <- new('AnnotatedDataFrame', data = metaRNA)

  fData_UMAP <- data.frame(gene_short_name = row.names(dataUMAP), row.names = row.names(dataUMAP))

  gene_metadata_UMAP <- new('AnnotatedDataFrame', data = fData_UMAP)
  gene_metadata_UMAP <- as(gene_metadata_UMAP, "data.frame")
  cell_metadata_UMAP <- as(cell_metadata_UMAP, "data.frame")

  #Construct monocle cds

  return (new_cell_data_set(dataUMAP,
                                   cell_metadata = cell_metadata_UMAP,
                                   gene_metadata = gene_metadata_UMAP))
}

seurat_to_cds <- function (seurat, dataRNA) {
  cds_from_seurat <- gen_cds(dataRNA)

  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)

  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

  ### Assign the cluster info

  list_cluster <- seurat@meta.data[[sprintf("ClusterIdents", 1.4, dim)]]

  names(list_cluster) <- dataRNA@data@Dimnames[[2]]

  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

  ### Could be a space-holder, but essentially fills out louvain parameters

  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

  ### Assign UMAP coordinate

  cds_from_seurat@int_colData@listData[["reducedDims"]][["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings

  ### Assign feature loading for downstream module analysis

  cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

  return (cds_from_seurat)
}

add_rna <- function(object, dataRNA){

  ret <- object
  obj_class = class(object)

  if (obj_class == "Seurat"){
    ret@assays$RNA <- dataRNA
  } else if (obj_class == "cell_data_set") {
    ret@assays@data@listData["counts"] <- dataRNA
  }else if (obj_class == "gg") {
    ret$plot_env$cds <- dataRNA
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
  } else if (obj_class == "gg") {
    ret$plot_env$cds = NULL
    ret$layers$mapping = NULL
  }

  return (ret)
}

get_rna <- function(object){

  ret <- object
  obj_class = class(object)

  if (obj_class == "Seurat"){
    return (ret@assays$RNA)
  } else if (obj_class == "cell_data_set") {
    return (ret@assays@data@listData)
  }else if (obj_class == "gg") {
    return (ret$plot_env$cds)
  }

}
