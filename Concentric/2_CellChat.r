library(CellChat)
library(dplyr)

cell_cell_interaction = function(cellchat)
{ 
    
    CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # set the used database in the object
    cellchat@DB <- CellChatDB.use
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 48) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    options(future.globals.maxSize = 2 * 1024^3) # 2 GB

    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250,  
                              contact.dependent = TRUE, contact.range = 100,
                              scale.distance = FALSE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    return(cellchat)
}


pseudo_mat = readRDS("../Deconvolution2/BayesPrism_pseudo_matrix/pseudo_spot_concentric.rds")

# combined the layer 0 and layer 1
meta_data = pseudo_mat@meta.data

# Add the region_celltype column
#pseudo_mat$region_celltype = paste0(pseudo_mat$major_annotation, "_", pseudo_mat$celltype) 

# Replace 'your_dataframe' with the actual name of your dataframe
meta_data = meta_data %>%
  mutate(distance = ifelse(distance == "layer0", "layer1", distance))
pseudo_mat$distance1 = meta_data$distance

# Separate them based on the diagnosis
Control = subset(pseudo_mat, subset = diagnosis == "Control")
Part = subset(pseudo_mat, subset = diagnosis == "Part")
AD = subset(pseudo_mat, subset = diagnosis == "AD")

# Separate based on the distance
Control_layer1 = subset(Control, subset = distance == "layer1")
Control_layer3 = subset(Control, subset = distance == "layer3")

Part_layer1 = subset(Part, subset = distance == "layer1")
Part_layer3 = subset(Part, subset = distance == "layer3")

AD_layer1 = subset(AD, subset = distance == "layer1")
AD_layer3 = subset(AD, subset = distance == "layer3")

########## For Control_layer1
data_input = Control_layer1@assays$RNA$data
meta = Control_layer1@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/Control_layer1_noregion.RDS")

########## For Control_layer3
data_input = Control_layer3@assays$RNA$data
meta = Control_layer3@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/Control_layer3_noregion.RDS")

########## For Part_layer1
data_input = Part_layer1@assays$RNA$data
meta = Part_layer1@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/Part_layer1_noregion.RDS")

########## For Part_layer3
data_input = Part_layer3@assays$RNA$data
meta = Part_layer3@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/Part_layer3_noregion.RDS")

########## For AD_layer1
data_input = AD_layer1@assays$RNA$data
meta = AD_layer1@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/AD_layer1_noregion.RDS")

########## For AD_layer3
data_input = AD_layer3@assays$RNA$data
meta = AD_layer3@meta.data

# Create the Cellchat object
cellchat = createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype",                           
    datatype = "RNA"
)

result = cell_cell_interaction(cellchat)
saveRDS(result, "Results/CellChat/AD_layer3_noregion.RDS")

