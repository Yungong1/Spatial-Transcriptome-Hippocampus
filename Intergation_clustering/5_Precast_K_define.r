library(PRECAST)
library(Seurat)
source("utils.R")
library(SingleCellExperiment)

# load the data
data0309 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//Processed_data/M0309.Result/",
    ID = "0309",
    diagnosis = "Control"
)

data0310 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//Processed_data/M0310.Result/",
    ID = "0310",
    diagnosis = "Control"
)

data0618 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//BMK_DATA_20240111093424_1/BD288-0618-matrix/outs/",
    ID = "0618",
    diagnosis = "Part"
)

data0623 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//BMK_DATA_20240111093424_1/BD288-0623-matrix/outs/",
    ID = "0623",
    diagnosis = "Part"
)

data0802 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//Processed_data/M0802/",
    ID = "0802",
    diagnosis = "AD"
)

data0915 = read_Visium(
    "../Data/10X_Hippocampus_Spatial_transcriptome//spatial_data//Processed_data/M0915//",
    ID = "0915",
    diagnosis = "AD"
)

seuList = list(data0309, data0310, data0618, data0623, data0802, data0915)

# add the position info into the meta table for each sample
for (i in 1:length(seuList))
{
    seuList[[i]]$col = GetTissueCoordinates(seuList[[i]])$x
    seuList[[i]]$row = GetTissueCoordinates(seuList[[i]])$y
}

metadataList <- lapply(seuList, function(x) x@meta.data)

for (r in seq_along(metadataList)) {
    meta_data <- metadataList[[r]]
    cat(all(c("row", "col") %in% colnames(meta_data)), "\n")  ## the names are correct!
}
                       
# Prepare PRECAST object
set.seed(2023)
preobj <- CreatePRECASTObject(seuList = seuList, selectGenesMethod = "SPARK-X", gene.number = 2000)
                       
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 1, maxIter = 30, verbose = TRUE)
                       
# simulate the data from 2 to 20 to find the best K
K_number = c()
sequence = seq(2, 30, by = 1)
for (num in sequence)
{
    PRECASTObj <- PRECAST(PRECASTObj, K = num)
    ## backup the fitting results in resList
    resList <- PRECASTObj@resList
    PRECASTObj <- SelectModel(PRECASTObj)

    print(PRECASTObj@seuList)
    #> NULL
    seuInt <- IntegrateSpaData(PRECASTObj)
    
    BIC = as.numeric(PRECASTObj@resList$icMat[1,2])
    K_number = c(K_number, BIC)
}
print(paste0("K_number:", K_number))
#print(paste0("BCI:", BCI))
                      
saveRDS(K_number, "../Result/Clustering_result/BCI.rds")
