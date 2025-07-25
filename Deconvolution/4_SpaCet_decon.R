library(SpaCET)
library(Seurat)
library(Matrix)

Decon = function(sc_count, sc_meta, spatial_data, sample)
{
    # Data clearning for spatial transcriptome
    tmp = subset(spatial_data, subset = ID == sample)
    spatial_count = tmp@assays$SCT@counts
    #spatial_count = merged@assays$SCT@counts
    spatial_location = GetTissueCoordinates(tmp)
    spatial_location = spatial_location[,1:2]
    
    # Data cleaning for snRNA reference
    meta_sub = data.frame(cbind(rownames(sc_meta), sc_meta$major_cell_type)) 
    colnames(meta_sub) = c("cellID", "cellType")
    rownames(meta_sub) = meta_sub$cellID
    
    # Perform the Decon
    SpaCET_obj <- create.SpaCET.object(
      counts= spatial_count,
      spotCoordinates=spatial_location,
      imagePath=NA,
      platform = "Visium"
    )
        
    sc_lineageTree = list()
    for (i in unique(meta_sub$cellType))
    {
        sc_lineageTree[i] = i
    }
    
    SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
      SpaCET_obj, 
      sc_counts=sc_count, 
      sc_annotation=meta_sub, 
      sc_lineageTree=sc_lineageTree, 
      coreNo=48
    )
    
    # Extract the proportion mat
    proportion_mat = SpaCET_obj@results$deconvolution$propMat
    proportion_mat = data.frame(t(proportion_mat))
    proportion_mat$x = spatial_location$x
    proportion_mat$y = spatial_location$y
    
    # Write the csv file for the proportion mat
    write.csv(proportion_mat, paste0("Decon_result/SpaCet/", sample, "/proportion.csv"))
    
    # Write the RDS for the object
    saveRDS(SpaCET_obj, paste0("Decon_result/SpaCet/", sample, "/SpaCET_obj.RDS"))
}

# Read the spatial data
merged = readRDS("../../Result/Region_annotation/processed_merged_spatial.rds")

# read the count and meta data for sc ref
sc_meta = readRDS("Reference/origin_sc_meta.rds")
sc_count = readRDS("Reference/origin_sc_count.rds")

for (sample in unique(merged@meta.data$ID))
{
    Decon(sc_count, sc_meta, merged, sample)
}