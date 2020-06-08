library(dplyr)
library(patchwork)
library(ggplot2)
library(monocle3)

############################### Monocle3 Trajectory#################################

# Loading processed umaps
load("data/experimentation/umap_parameter_adjustment/umaps_res_.6-1.8.rda")

# converts Seurat umap object to monocle object

umap_to_monocle <- function (umap) {
  dataUMAP <- umap@assays$RNA@data
  cell_metadata_UMAP <- new('AnnotatedDataFrame', data = umap@meta.data)
  
  fData_UMAP <- data.frame(gene_short_name = row.names(dataUMAP), row.names = row.names(dataUMAP))
  
  gene_metadata_UMAP <- new('AnnotatedDataFrame', data = fData_UMAP)
  gene_metadata_UMAP <- as(gene_metadata_UMAP, "data.frame")
  cell_metadata_UMAP <- as(cell_metadata_UMAP, "data.frame")
  
  #Construct monocle cds
  return (new_cell_data_set(dataUMAP, 
                            cell_metadata = cell_metadata_UMAP,
                            gene_metadata = gene_metadata_UMAP))
}

# Processes monocle object; can change parameter

process_monocle <- function (cell_set) {

  #PreProcess
  cell_set <- preprocess_cds(cell_set, num_dim = 50)
  
  # Reduce dimensionality and visualize the results
  cell_set <- reduce_dimension(cell_set)
  
  # Cluster your cells
  
  cell_set <- cluster_cells(cell_set)
  
  # learn the trajectory graph
  
  cell_set <- learn_graph(cell_set)
  
  # tell Monocle where the "beginning" of the biological process is. 
  # We do so by choosing regions of the graph that we mark as "roots" of the trajectory.
  
  cell_set <- order_cells(cell_set)
  
  return (cell_set)
}

# plots a single tissue colored by stage; all other tissues are colored grey

plot_single_tissue <- function (cell_set, tissue_name) {
  
  # Retrieving tissue and stages metadata columns
  
  tissues <- cell_set@colData$tissue
  stages <- cell_set@colData$stage
  
  # color palette to distinguish cell stages; grey for cell types that are not the tissue provided
  
  c_palette = c("grey",  "#FAE755", "#BADB52", "#6BBD76", "#48958B", "#3D688B", "#43377C")
  
  # reference labels for all stages 
  
  stage_labels = c(0, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5)
  
  # create new column call t_stage

  col_new = c()
  for (i in c(1:nrow(cell_set@colData)) ){
    
      # cells that are the tissue type provided are recolored based on stage 
    
      if (tissues[i] == tissue_name){
        s1 = toString(stages[i])
        col_new <- c(col_new, as.numeric(substr(s1,2, nchar(s1))))
        
      # other cells are colored with grey
        
      } else {
        col_new <- c(col_new, 0.0)
      }
  }
  cell_set@colData$t_stage = as.factor(col_new)
  
  # creating ggplot colored by t_stage column
  
  ggplot <- plot_cells(cell_set,
                       color_cells_by = "t_stage",
                       label_cell_groups=FALSE,
                       label_leaves=TRUE,
                       label_branch_points=TRUE,
                       graph_label_size=1.5,
                       cell_size=.75
  )
  
  # assigning correct colors to t_stage labels
  
  color_vals = c()
  stage_levels = levels(cell_set@colData$t_stage)
  j = 1
  
  for (i in c(1:length(stage_labels))) {
     if (j <= length(stage_levels) && stage_levels[j] == stage_labels[i]){
        color_vals <- c(color_vals, verdis[i])
        j = j + 1
     }
  }
  
  # returns manually colored ggplot 
  
  return (ggplot + scale_color_manual(values = color_vals) + ggtitle(tissue_name))
}

plot_single_stage <- function (cell_set, stage_name) {
  
  # Retrieving tissue and stages metadata columns
  
  tissues <- cell_set@colData$tissue
  stages <- cell_set@colData$stage
  
  # color palette to distinguish cell stages; grey for cell types that are not the tissue provided
  
  c_palette = c("grey",  "#FAE755", "#BADB52", "#6BBD76", "#48958B", "#3D688B", "#43377C")
  
  # reference labels for all tissues 
  
  tissue_labels = c("Other", "Body", "Brain", "Head", "Limb", "Liver", "YS")
  
  # create new column call t_stage
  
  col_new = c()
  for (i in c(1:nrow(cell_set@colData)) ){
    
    # cells that are the tissue type provided are recolored based on tissue 
    
    if (stages[i] == stage_name){
      col_new <- c(col_new, toString(tissues[i]))
      
      # other cells are colored with grey
      
    } else {
      col_new <- c(col_new, "Other")
    }
  }
  cell_set@colData$t_stage = as.factor(col_new)
  
  tissue_levels = levels(cell_set@colData$t_stage)
  
  other_index = which(tissue_levels == "Other")
  
  tissue_levels <- c("Other", tissue_levels[- other_index])

  cell_set@colData$t_stage <- factor(cell_set@colData$t_stage, levels = tissue_levels)
  
  # creating ggplot colored by t_stage column
  
  ggplot <- plot_cells(cell_set,
                       color_cells_by = "t_stage",
                       label_cell_groups=FALSE,
                       label_leaves=TRUE,
                       label_branch_points=TRUE,
                       graph_label_size=1.5,
                       cell_size=.75
  )
  
  # assigning correct colors to t_stage tissues
  
  color_vals = c()
  tissue_levels = levels(cell_set@colData$t_stage)
  j = 1
  
  for (i in c(1:length(tissue_labels))) {
    if (j <= length(tissue_levels) && tissue_levels[j] == tissue_labels[i]){
      color_vals <- c(color_vals, c_palette[i])
      j = j + 1
    }
  }
  
  # returns manually colored ggplot 
  
  return (ggplot + scale_color_manual(values = color_vals) + ggtitle(stage_name))
}

# converting all umaps into monocle cell sets

cell_sets <- c()

for (umap in umaps_processed) {
  cell_sets <- c(tdcds_umaps, umap_to_monocle(umap))
}

# Retrieving monocle objects w/ res = 1.4C

cell_set_1.4C <- process_monocle(cell_sets[[3]])

# Reordering stages from alphabetical to chronological

cell_set_1.4C@colData$stage <- factor(cell_set_1.4C@colData$stage, levels = c("E7.5", "E8.5", "E9.5", "E10.5", "E12.5", "E14.5"))

# Plotting single tissue stages for all tissues present

p <- plot_single_stage(cell_set_1.4C, "YS")

tissue_levels = c("Liver", "Limb", "Head", "Brain", "Body")
for (tissue_name in tissue_levels){
   p <- p + plot_single_tissue(cell_set_1.4C, tissue_name)
}

p

# Coloring single stage tissues for all stages

p <- plot_single_stage(cell_set_1.4C, "E7.5")

stage_levels = c("E8.5", "E9.5", "E10.5", "E12.5", "E14.5")
for (stage_name in stage_levels){
  p <- p + plot_single_stage(cell_set_1.4C, stage_name)
}

p

# ggplot <- plot_cells(tdcds_TSNE,
#              color_cells_by = "stage",
#              label_cell_groups=FALSE,
#              label_leaves=TRUE,
#              label_branch_points=TRUE,
#              graph_label_size=1.5,
#              cell_size=.5
#              )

# plot_cells(umap_1.4C, color_cells_by = "partition")

####################Order the cells in pseudotime

tdcds_TSNE <- order_cells(tdcds_TSNE)
plot_cells(tdcds_TSNE,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#####################Top 10 genes#########################################

# plot_cells() to visualize how individual genes vary along the trajectory
top10_genes <- c("Col4a1" , "Ngp" ,    "Col4a2" , "Lamb1" ,  "Mpo" ,    "Ltf"  ,   "Sparc" ,  "Lama1" ,  "Cd74" ,   "Hbb-bh1")
plot_cells(tdcds_TSNE,    genes=top10_genes,       label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


################### can save your calculated R object#######################

# saveRDS(all_1816_cell, file = "all_1816_cell.rds")
# save(list=c("all_1816_cell", "xxx","xxx"), file="single_cell_objects_04092020.rda")

#try2 (reference)
saveRDS(all_embryo_cell_QC_15N_1.4C_tSNE, file = "all_embryo_cell_QC_15N_1.4C_tSNE.rds")
save(list=c("all_embryo_cell_QC_15N_1.4C_tSNE", "all_embryo_cell_QC"), file="single_cell_objects_05272020.rda")

#Trajecotry
saveRDS(tdcds_TSNE, file = "15N1.4C_Trajectory_tdcds_TSNE_05272020.rds")
save(list=c("tdcds_TSNE", "?"), file="Embryo_15N1.4C_Trajectory_objects_05272020.rda")

# save the workspace in a session (all objects are saved automatically)
save.image(file="filename.RData")
savehistory()

#try2 (reference)
save.image(file="Embryo_scRNAseq_15N1.4C_05272020.RData")
savehistory()

#Trajectory
save.image(file="Embryo_scRNAseq_15N1.4C_Trajectory_05272020.RData")
savehistory()

