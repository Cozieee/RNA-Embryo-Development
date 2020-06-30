tissue_labels = c("Other", "YS", "Body", "Head", "Liver", "Limb", "Brain")
stage_labels = c("Other", "E7.5", "E8.5", "E9.5", "E10.5", "E12.5", "E14.5")

c_palette = c("grey", "#F33B5A", "#FF8614", "#F6E702", "#75D21B", "#22AAEB", "#884AF8")

tissue_palette = c("Other" = "grey", "YS" = "#F33B5A", "Body"="#FF8614",
                   "Head"="#F6E702", "Liver"="#75D21B", "Limb"="#22AAEB", "Brain"="#884AF8")

stage_palette = c("Other" = "grey", "E7.5" = "#F33B5A", "E8.5"="#FF8614",
                  "E9.5"="#F6E702", "E10.5"="#75D21B", "E12.5"="#22AAEB", "E14.5"="#884AF8")

# Only labels c2 of cells whose c1 label matches name. All else labeled 'Other'

gen_tissue_stage <- function (c1, c2, name) {
  col_new = c()

  for (i in c(1:length(c1)) ){

    # cells with the name provided, in are labeled based off c2

    if (c1[i] == name){
      col_new <- c(col_new, toString(c2[i]))

      # other cells are colored with grey

    } else {
      col_new <- c(col_new, "Other")
    }
  }

  return (col_new)
}

plot_single_tissue <- function (cds, tissue_name, title = NULL) {

  # Retrieving tissue and stages metadata columns

  tissues <- cds@colData$tissue
  stages <- cds@colData$stage

  # create new column call t_stage

  col_new <- gen_tissue_stage(tissues, stages, tissue_name)

  cds@colData$t_stage = as.factor(col_new)

  cds@colData$t_stage <- factor(cds@colData$t_stage, levels = stage_labels)

  # creating ggplot colored by t_stage column

  ggplot <- del_rna(plot_cells(cds,
                               color_cells_by = "t_stage",
                               label_cell_groups=FALSE,
                               label_leaves=TRUE,
                               label_branch_points=TRUE,
                               graph_label_size=1.5,
                               cell_size=.75,
                               trajectory_graph_segment_size = .5
  ))

  # assigning correct colors to t_stage labels

  return (ggplot + scale_color_manual(values = stage_palette) + ggtitle(title))
}

plot_single_stage <- function (cds, stage_name, title = NULL) {

  # Retrieving tissue and stages metadata columns

  tissues <- cds@colData$tissue
  stages <- cds@colData$stage

  # create new column call t_stage

  col_new <- gen_tissue_stage(stages, tissues, stage_name)

  cds@colData$t_stage = as.factor(col_new)

  cds@colData$t_stage <- factor(cds@colData$t_stage, levels = tissue_labels)

  # creating ggplot colored by t_stage column

  ggplot <- del_rna(plot_cells(cds,
                               color_cells_by = "t_stage",
                               label_cell_groups=FALSE,
                               label_leaves=TRUE,
                               label_branch_points=TRUE,
                               graph_label_size=1.5,
                               cell_size=.75,
                               trajectory_graph_segment_size = .5
  ))

  # returns manually colored ggplot

  return (ggplot + scale_color_manual(values = tissue_palette) + ggtitle(title))
}

plot_by_label <- function(cds, label = NULL, title = NULL) {
  if (is.null(label)){
    return (del_rna(plot_cells(cds,
                       label_cell_groups=FALSE,
                       label_leaves=TRUE,
                       label_branch_points=TRUE,
                       graph_label_size=1.5,
                       cell_size=.5,
                       trajectory_graph_segment_size = .5)
            + scale_color_manual(values = tail(c_palette, -1))
            + ggtitle(title)))
    } else {
      return (del_rna(plot_cells(cds,
                         color_cells_by = label,
                         label_cell_groups=FALSE,
                         label_leaves=TRUE,
                         label_branch_points=TRUE,
                         graph_label_size=1.5,
                         cell_size=.5,
                         trajectory_graph_segment_size = .5)
              + scale_color_manual(values = tail(c_palette, -1))
              + ggtitle(title)))
    }
}

multi_plot <- function(plot_func, cds_set, labels, titles = NULL) {

  if (length(labels) == 1 || is.null(labels) ){
    labels = rep(labels, length(cds_set))
  } else if(!is.vector(cds_set)){
    temp <- c()
    for (i in c(1:length(labels))){
      temp <- c(temp, cds_set)
    }
    cds_set <- temp
  }

  p <- plot_func(cds_set[[1]], labels[1], titles[1])

  for (i in c(2:length(cds_set))){
    p <- p + plot_func(cds_set[[i]], labels[i], titles[i])
  }

  return (p)
}
