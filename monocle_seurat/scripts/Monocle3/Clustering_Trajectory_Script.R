library(dplyr)
library(patchwork)
library(ggplot2)
library(monocle3)

detach("package:seqtools", unload=TRUE)
library(seqtools)

############################### Monocle3 Trajectory#################################

# Loading processed umaps
load("data/experimentation/umap_param_adj/mono_autoN_32-37C.rds")
load("data/experimentation/umap_param_adj/mono_autoN_38-43C.rds")
load("data/experimentation/main/cds_data_RNA.rds")

trials <- c()

low = 32
high = 37
dimensions = seq(low, high, by=1)

IMGPATH = paste("data/images/plots/mono_autoN_", low,"-", high,"C/",sep="")
ggsave(paste(IMGPATH, "label/", "pseudotime", ".jpg", sep = ""), plot = p, width = 18, height = 10)

for (i in c(1:6)){
  trials32_37C[[i]] <- del_rna(trials32_37C[[i]])
}

for (i in c(1:6)){
  trials38_43C[[i]] <- order_cells(trials38_43C[[i]])
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
  
  # 6-panel plot for each tissue
  
  plots[["all_tissues"]][[i]] <- multi_plot(plot_single_tissue, trials[[i]], tissue_levels, dim_title(tissue_levels, dimensions[i]))
  
  # 6-panel plot for each stage
  
  plots[["all_stages"]][[i]] <- multi_plot(plot_single_stage, trials[[i]], stage_levels, dim_title(stage_levels, dimensions[i]))
}

####################Order the cells in pseudotime
p <- plot_time(t1[[1]], 1)
for (i in c(2:6)){
  p <- p + plot_time(t1[[i]], i)
}

p

plot_time <- function(cds, i){
  return (plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  trajectory_graph_segment_size = .5,
                  cell_size = .5,
                  graph_label_size=1.5) + ggtitle(dimensions[[i]]))
}

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

############################################################################
# for (i in c(1:6)){
#   trials <- c(trials, control_cds)
#   trials[[i]] <- preprocess_cds(trials[[i]], num_dim = dimensions[i])
#   trials[[i]] <- reduce_dimension(trials[[i]])
#   trials[[i]] <- cluster_cells(trials[[i]])
#   trials[[i]] <- learn_graph(trials[[i]])
#   trials[[i]]@colData$stage <- factor(trials[[i]]@colData$stage, levels = c("E7.5","E8.5","E9.5","E10.5","E12.5","E14.5"))
# }

# for(i in c(1:6)){
#   
#   # Convert Seurat objects to monocle object
#   
#   trials[[i]] <- seurat_to_cds(trials[[i]]) 
#   
#   # Learn UMAP graph for trajectory
#   
#   trials[[i]] <- learn_graph(trials[[i]])
#   
#   # Reorder stage factors
#   
#   trials[[i]]@colData$stage <- factor(trials[[i]]@colData$stage, levels = c("E7.5","E8.5","E9.5","E10.5","E12.5","E14.5"))
# }

# plots <- list()
# 
# plots[["stage"]] <- multi_plot(plot_by_label, trials, "stage", dimensions)
# plots[["tissue"]] <- multi_plot(plot_by_label, trials, "tissue", dimensions)
# plots[["all_tissues"]] <- list()
# plots[["all_stages"]] <- list()

# p <- multi_plot(plot_by_label, trials38_43C, "tissue", dimensions)
# 
# plot_clusters <- function(cds, i){
#   return (plot_cells(cds,
#                      label_cell_groups=FALSE,
#                      label_leaves=TRUE,
#                      label_branch_points=TRUE,
#                      graph_label_size=1.5,
#                      cell_size=.5,
#                      trajectory_graph_segment_size = .5)
#           + ggtitle(dimensions[i]))
# }
# 
# p <- plot_clusters(trials32_37C[[1]], 1)
# for (i in c(2:6)){
#   p <- p + plot_clusters(trials32_37C[[i]], i)
# }
# 
# p
# 
# ggsave(paste(IMGPATH, "label/", "cluster", ".jpg", sep = ""), plot = p, width = 18, height = 10)
