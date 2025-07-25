library(Seurat)
library(dplyr)
library(tibble)
library(PANDA)

random_cluster_select = function(sc_count, sc_meta, prop)
{
    # For the train dataset
    sc_meta_train = sc_meta %>%
      mutate(rowname = rownames(.)) %>%  # Add row names as a new column
      group_by(major_cell_type) %>%
      slice_sample(prop = prop) %>%       # Sample 20% of rows per group
      ungroup() %>%
      column_to_rownames(var = "rowname")  # Restore row names after sampling
    
    sc_count_train = sc_count[, rownames(sc_meta_train)]
    
    # For the test dataset
    sc_meta_test_barcode = setdiff(rownames(sc_meta), rownames(sc_meta_train))
    sc_meta_test = sc_meta[sc_meta_test_barcode,]
    sc_count_test = sc_count[, rownames(sc_meta_test)]
    
    result = list(sc_count_train, sc_meta_train, sc_count_test, sc_meta_test)
    
    return(result)
}

decon_PANDA = function(sc_count, sc_meta, st_count, dir)
{
    sc_labels = sc_meta$major_cell_type
    sc_count = t(as.matrix(sc_count))
    st_count = t(as.matrix(st_count))
    sc_results = sc_train(sc_count, sc_labels, n_archetypes_vec = 10, save_dir = paste0(dir, "/sc_results"))
    st_results = st_train(st_count, sc_results = sc_results, save_dir = paste0(dir, "/st_result"))
    saveRDS(st_results, paste0(dir, "/result.RDS"))
}


# Read the spatial data
merged = readRDS("../../../../Result/Region_annotation/processed_merged_spatial.rds")

# Read sc data
sc_count = readRDS("../../Reference/sc_count_with_VAS.rds")
sc_meta = readRDS("../../Reference/sc_meta_with_VAS.rds")


for (sample in c("0309", "0310", "0618", "0623", "0802", "0915"))
{
    res = random_cluster_select(sc_count, sc_meta, prop = 0.5)
    tmp = subset(merged, subset = ID == sample)
    st_count = tmp@assays$SCT@counts
    decon_PANDA(res[[1]], res[[2]], st_count, dir = paste0(sample, "/50/train"))
    decon_PANDA(res[[3]], res[[4]], st_count, dir = paste0(sample, "/50/test"))
}
