
plot_single_tissue <- function (cds, tissue_name, title = NULL) {

  # Retrieving tissue and stages metadata columns

  tissues <- cds@colData$tissue
  stages <- cds@colData$stage

  # color palette to distinguish cell stages; grey for cell types that are not the tissue provided

  c_palette = c("grey",  "#B83679", "#EB7184", "#B3DE45", "#53B3AE", "#3A8AB9", "#192965")

  # reference labels for all stages

  stage_labels = c(0, 7.5, 8.5, 9.5, 10.5, 12.5, 14.5)

  # create new column call t_stage

  col_new = c()
  for (i in c(1:nrow(cds@colData)) ){

    # cells that are the tissue type provided are recolored based on stage

    if (tissues[i] == tissue_name){
      s1 = toString(stages[i])
      col_new <- c(col_new, as.numeric(substr(s1,2, nchar(s1))))

      # other cells are colored with grey

    } else {
      col_new <- c(col_new, 0.0)
    }
  }
  cds@colData$t_stage = as.factor(col_new)

  # creating ggplot colored by t_stage column

  ggplot <- plot_cells(cds,
                       color_cells_by = "t_stage",
                       label_cell_groups=FALSE,
                       label_leaves=TRUE,
                       label_branch_points=TRUE,
                       graph_label_size=1.5,
                       cell_size=.75,
                       trajectory_graph_segment_size = .5
  )

  # assigning correct colors to t_stage labels

  color_vals = c()
  stage_levels = levels(cds@colData$t_stage)
  j = 1

  for (i in c(1:length(stage_labels))) {
    if (j <= length(stage_levels) && stage_levels[j] == stage_labels[i]){
      color_vals <- c(color_vals, c_palette[i])
      j = j + 1
    }
  }

  # returns manually colored ggplot

  return (ggplot + scale_color_manual(values = color_vals) + ggtitle(title))
}

plot_single_stage <- function (cds, stage_name, title = NULL) {

  # Retrieving tissue and stages metadata columns

  tissues <- cds@colData$tissue
  stages <- cds@colData$stage

  # color palette to distinguish cell stages; grey for cell types that are not the tissue provided

  c_palette = c("grey",  "#B83679", "#EB7184", "#B3DE45", "#53B3AE", "#3A8AB9", "#192965")

  # reference labels for all tissues

  tissue_labels = c("Other", "Body", "Brain", "Head", "Limb", "Liver", "YS")

  # create new column call t_stage

  col_new = c()
  for (i in c(1:nrow(cds@colData)) ){

    # cells that are the tissue type provided are recolored based on tissue

    if (stages[i] == stage_name){
      col_new <- c(col_new, toString(tissues[i]))

      # other cells are colored with grey

    } else {
      col_new <- c(col_new, "Other")
    }
  }
  cds@colData$t_stage = as.factor(col_new)

  tissue_levels = levels(cds@colData$t_stage)

  other_index = which(tissue_levels == "Other")

  tissue_levels <- c("Other", tissue_levels[- other_index])

  cds@colData$t_stage <- factor(cds@colData$t_stage, levels = tissue_levels)

  # creating ggplot colored by t_stage column

  ggplot <- plot_cells(cds,
                       color_cells_by = "t_stage",
                       label_cell_groups=FALSE,
                       label_leaves=TRUE,
                       label_branch_points=TRUE,
                       graph_label_size=1.5,
                       cell_size=.75,
                       trajectory_graph_segment_size = .5
  )

  # assigning correct colors to t_stage tissues

  color_vals = c()
  tissue_levels = levels(cds@colData$t_stage)
  j = 1

  for (i in c(1:length(tissue_labels))) {
    if (j <= length(tissue_levels) && tissue_levels[j] == tissue_labels[i]){
      color_vals <- c(color_vals, c_palette[i])
      j = j + 1
    }
  }

  # returns manually colored ggplot

  return (ggplot + scale_color_manual(values = color_vals) + ggtitle(title))
}

plot_by_label <- function(cds, label, title = NULL) {

  return (plot_cells(cds,
                     color_cells_by = label,
                     label_cell_groups=FALSE,
                     label_leaves=TRUE,
                     label_branch_points=TRUE,
                     graph_label_size=1.5,
                     cell_size=.5,
                     trajectory_graph_segment_size = .5)
          + scale_color_manual(values = c("#B83679", "#EB7184", "#B3DE45", "#53B3AE", "#3A8AB9", "#192965"))
          + ggtitle(title))
}

multi_plot <- function(plot_func, cds_set, labels, titles = NULL) {

  if (length(labels) == 1){
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
