# The Result for this method is not good, just abandon this
library(Seurat)
library(dplyr)
library(CARD)
library(Matrix)
library(tidyr)
library(ggplot2)

# Perform the Card spatial
CARD_perform = function(sc_count, spatial_data, sample)
{
    tmp = subset(spatial_data, subset = ID == sample)
    spatial_count = tmp@assays$SCT@counts
    spatial_location = GetTissueCoordinates(tmp)
    spatial_location = spatial_location[,1:2]

    # Create CARD object
    CARD_obj = createCARDObject(
        sc_count = sc_count,
        sc_meta = sc_meta,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        ct.varname = "major_cell_type",
        ct.select = unique(sc_meta$major_cell_type),
        sample.varname = "projid",
        minCountGene = 100,
        minCountSpot = 5)
    # Perform CARD deconvolution
    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

    result_prop = data.frame(CARD_obj@Proportion_CARD)
    result_prop$x = spatial_location$x
    result_prop$y = spatial_location$y
    write.csv(result_prop, paste0("Decon_result//CARD/", sample,"/proportion.csv"))
}

# Read the reference data
sc_meta = readRDS("Reference/origin_sc_meta.rds")
sc_count = readRDS("Reference/origin_sc_count.rds")

# Read the spatial data
merged = readRDS("../../Result/Region_annotation/processed_merged_spatial.rds")

# read the meta
meta = read.csv("../../Result/Raw_data_for_R/meta_data.csv")
merged[["annotation"]] = meta$annotation

# Create the file
for (sample in c("0309", "0310", "0618", "0623", "0802", "0915"))
{
    CARD_perform(sc_count, merged, sample)
}
