library(dplyr)
library(patchwork)
library(ggplot2)
library(monocle3)

detach("package:seqtools", unload=TRUE)
library(seqtools)

############################### Monocle3 Trajectory#################################

# Loading processed umaps
load("data/experimentation/umap_parameter_adjustment/trial_dims_10-35.rda")

for(i in c(1:6)){
  
  # Convert Seurat objects to monocle object
  
  trials[[i]] <- seurat_to_cds(trials[[i]]) 
  
  # Learn UMAP graph for trajectory
  
  trials[[i]] <- learn_graph(trials[[i]])
  
  # Reorder stage factors
  
  trials[[i]]@colData$stage <- factor(trials[[i]]@colData$stage, levels = c("E7.5","E8.5","E9.5","E10.5","E12.5","E14.5"))
}

# Convenience function to add dimension to multi-plot title

dim_title <- function(levels, dim){
  title <- levels
  title[1] = paste(title[1], " (Dim ", dim, ")", sep ="")
  return (title)
}

# levels to be plotted

tissue_levels = c("YS", "Body", "Head", "Liver", "Limb", "Brain")
stage_levels = c("E7.5", "E8.5", "E9.5", "E10.5", "E12.5", "E14.5")

for (i in c(1:length(trials))){
  
  # plotting trial for each tissue
  
  p <- multi_plot(plot_single_tissue, trials[[i]], tissue_levels, dim_title(tissue_levels, dimensions[i]))
  print(p)
  readline()
  
  # plotting trial for each stage
  
  p <- multi_plot(plot_single_stage, trials[[i]], stage_levels, dim_title(stage_levels, dimensions[i]))
  print(p)
  readline()
}

####################Order the cells in pseudotime

for (i in c(2:6)){
  p <- plot_cells(cell_set_24N,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  trajectory_graph_segment_size = .5,
                  cell_size = .75,
                  graph_label_size=1.5) + ggtitle("Dim 24")
}

p

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

