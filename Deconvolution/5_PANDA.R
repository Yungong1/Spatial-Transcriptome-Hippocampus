# PANDA for all spatial samples
library(PANDA)
library(dplyr)
library(Seurat)

decon_PANDA = function(sc_count, sc_meta, st_count, dir)
{
    sc_labels = sc_meta$major_cell_type
    sc_count = t(as.matrix(sc_count))
    st_count = t(as.matrix(st_count))
    sc_results = sc_train(sc_count, sc_labels, n_archetypes_vec = 10, n_cores = 48, save_dir = paste0(dir, "/sc_results"))
    st_results = st_train(st_count, sc_results = sc_results, save_dir = paste0(dir, "/st_result"))
    saveRDS(st_results, paste0(dir, "/result.RDS"))
}

# Read the spatial data
merged = readRDS("../../Result/Region_annotation/processed_merged_spatial.rds")

# Read sc data
sc_count = readRDS("Reference/sc_count_with_VAS.rds")
sc_meta = readRDS("Reference/sc_meta_with_VAS.rds")

# For all samples together
st_count = merged@assays$SCT@counts

overlap_gene = intersect(rownames(st_count), rownames(sc_count))

st_count = st_count[overlap_gene,]
sc_count = sc_count[overlap_gene,]

decon_PANDA(sc_count, sc_meta, st_count, "Decon_result/PANDA_results/All_cell_ref/")