library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
source("utils.R")
library(CellChat)

# Load the data
merge = readRDS("../Result/Region_annotation/processed_merged_spatial.rds")

# Read the meta data from python result (annotation)
meta = read.csv("../Result/Raw_data_for_R/meta_data.csv")

merge$annotation = meta$annotation

# Check the result
tmp = subset(merge, subset = ID == "0618")

Idents(tmp) = "annotation"
SpatialDimPlot(tmp, label = TRUE, label.size = 3)

# Prepare input data for cellchat
data.input = Seurat::GetAssayData(tmp, slot = "data", assay = "SCT")

meta = data.frame(
    labels = Seurat::Idents(tmp), 
    samples = "sample1", 
    row.names = names(Seurat::Idents(tmp))
)

meta$samples <- factor(meta$samples)

# Check the label
unique(meta$labels)

tmp_loc = Seurat::GetTissueCoordinates(tmp, scale = NULL, cols = c("imagerow", "imagecol"))
spatial.locs = data.frame(cbind(tmp_loc[["x"]], tmp_loc[["y"]]))
rownames(spatial.locs) = rownames(spatial.locs)

scalefactors = jsonlite::fromJSON(
    txt = file.path(
        "../Data/10X_Hippocampus_Spatial_transcriptome/spatial_data/BMK_DATA_20240111093424_1/BD288-0618-matrix/outs/spatial/", 
        'scalefactors_json.json')
)

spot.size = 55 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(
    coordinates = spatial.locs, 
    ratio = spatial.factors$ratio, 
    tol = spatial.factors$tol
)
min(d.spatial[d.spatial!=0])

# Create the object
cellchat <- createCellChat(
    object = data.input, 
    meta = meta, 
    group.by = "labels",
    datatype = "spatial", 
    coordinates = spatial.locs, 
    spatial.factors = spatial.factors
)

# Set the ligand receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(
    CellChatDB, 
    search = c("Secreted Signaling"), 
    key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 48) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
options(future.globals.maxSize = 1024 * 1024 * 1024) # Set to 1 GiB, or adjust as needed
#> The number of highly variable ligand-receptor pairs used for signaling inference is 648

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.5,
                              contact.dependent = TRUE, contact.range = 100)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

# Save the data
saveRDS(cellchat, "cell_chat_test_result.RDS")
