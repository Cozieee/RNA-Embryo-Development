library(Seurat)
library(dplyr)
library(patchwork)

##################################GENE EXPRESSION###############################

# Finding differentially expressed features (cluster biomarkers)
# find markers of cluster 1

#try2
cluster1.markers <- FindMarkers(all_embryo_cell_QC_15N_1.4C_UMAP, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
write.csv(cluster1.markers, "cluster1.markers.csv")



# find markers for every cluster compared to all remaining cells, report only the positive ones

#try2 
all_embryo_cell_QC_15N_1.4C_tSNE.markers <- FindAllMarkers(all_embryo_cell_QC_15N_1.4C_tSNE, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#try2
all_embryo_cell_QC_15N_1.4C_tSNE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
View(all_embryo_cell_QC_15N_1.4C_tSNE.markers)
write.csv(all_embryo_cell_QC_15N_1.4C_tSNE.markers, "all_embryo_cell_QC_15N_1.4C_tSNE.markers.csv")



# Visualize gene expression
#try2
VlnPlot(all_embryo_cell_QC_15N_1.4C_tSNE, features = c("Tmem119", "Clec7a"))


# plot raw counts as well
#try2
VlnPlot(all_embryo_cell_QC_15N_1.4C_tSNE, features = c("Tmem119", "Clec7a"), slot = "counts", log = TRUE)
FeaturePlot(all_embryo_cell_QC_15N_1.4C_tSNE, features = c("Tmem119", "Clec7a"))



# use heatmap to visualize gene expression
#try2
top10 <- all_embryo_cell_QC_15N_1.4C_tSNE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(all_embryo_cell_QC_15N_1.4C_tSNE, features = top10$gene) + NoLegend()
write.csv(top10, "top10.csv")


# Assigning cell type identity to clusters
#try2
new.cluster.ids <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x")
names(new.cluster.ids) <- levels(all_embryo_cell_QC_15N_1.4C_tSNE)
all_embryo_cell_QC_15N_1.4C_tSNE <- RenameIdents(all_embryo_cell_QC_15N_1.4C_tSNE, new.cluster.ids)
DimPlot(all_embryo_cell_QC_15N_1.4C_tSNE, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23")
names(new.cluster.ids) <- levels(all_embryo_cell_QC_15N_1.4C_tSNE)
all_embryo_cell_QC_15N_1.4C_tSNE <- RenameIdents(all_embryo_cell_QC_15N_1.4C_tSNE, new.cluster.ids)
DimPlot(all_embryo_cell_QC_15N_1.4C_tSNE, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# can save your calculated R object

# saveRDS(all_1816_cell, file = "all_1816_cell.rds")
# save(list=c("all_1816_cell", "xxx","xxx"), file="single_cell_objects_04092020.rda")
#try2
saveRDS(all_embryo_cell_QC_15N_1.4C_tSNE, file = "all_embryo_cell_QC_15N_1.4C_tSNE.rds")
save(list=c("all_embryo_cell_QC_15N_1.4C_tSNE", "all_embryo_cell_QC"), file="single_cell_objects_05272020.rda")



# save the workspace in a session (all objects are saved automatically)
save.image(file="filename.RData")
savehistory()

#try2
save.image(file="Embryo_scRNAseq_15N1.4C_05272020.RData")
savehistory()