# Make sure the working directory is set to monocle3 folder.

load("data/unprocessed/single_cell_required_objects_06032020.rda")

# get adult cell names (excluding WM)

E7.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E7.5"] 

E8.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E8.5"] 

E9.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E9.5"] 

E10.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E10.5"] 

E12.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E12.5"] 

E14.5_all_names <- rownames(my_data_design.new2) [my_data_design.new2$stage == "E14.5"] 

# check total number of cells

length(E7.5_all_names) + length(E8.5_all_names) + length(E9.5_all_names) + length(E10.5_all_names) + length(E12.5_all_names) + length(E14.5_all_names)

# take a quick look at the cell names

head(E7.5_all_names)

# subset raw read counts (no ERCC)

all_embryo_cell_counts <- my_data_cleaned_no_ERCC_raw [ , c(E7.5_all_names,E8.5_all_names,
                                                            E9.5_all_names,E10.5_all_names,E12.5_all_names,E14.5_all_names) ]

# subset meta data

all_embryo_meta <- my_data_design.new2 [c(E7.5_all_names,E8.5_all_names,E9.5_all_names,E10.5_all_names,E12.5_all_names,E14.5_all_names), ]

# check dimension

dim(all_embryo_cell_counts)
dim(all_embryo_meta)

# make sure cell names are identical

identical (rownames(all_embryo_meta), names(all_embryo_cell_counts))

save(all_embryo_cell_counts, all_embryo_meta, file = "experimentation/main/all_embryo_counts.rda")