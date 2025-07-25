library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)


read_Visium = function(dir, ID, diagnosis, sc_transform = FALSE)
# read the 10X Viisum data

{
    spatial_dir = paste0(dir, "spatial/")
    
    # if tissue_positions_list.csv exist, change the the name into tissue1_positions_list.csv
    for (file in list.files(spatial_dir))
    {
        if (file == "tissue_positions_list.csv")
            {
                current_name = paste0(spatial_dir, "tissue_positions_list.csv")
                new_name = paste0(spatial_dir, "tissue1_positions_list.csv")
                
                file.rename(current_name, new_name)

            }
    }

    # load the data
    data = Load10X_Spatial(
      dir,
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial",
      slice = "slice1",
      bin.size = NULL,
      filter.matrix = TRUE,
      to.upper = FALSE,
      image = NULL
    )
    
    # change the name back
    file.rename(new_name, current_name)
    
    data$ID = ID
    data$diagnosis = diagnosis
    
    if (sc_transform == TRUE)
    {
        data =  SCTransform(data, assay = "Spatial", verbose = FALSE)
    }
    
    
    return(data)

}