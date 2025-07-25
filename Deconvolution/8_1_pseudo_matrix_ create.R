library(Seurat)
suppressWarnings(library(BayesPrism))
library(Matrix)
library(dplyr)

# Read the Spatial meta data
# Check the location of the End cells
meta_info = read.csv("../../Result/Region_annotation/meta_info.csv")

meta_info$individual = as.character(meta_info$individual)
meta_info$individual = paste0("0",meta_info$individual)
unique(meta_info$individual)

# Read the result of the bayePrism
decon_all_data = readRDS("Decon_result/BayesPrism/decon_result_allgenes.rds")

# Convert the data into pseudo single cell
meta_combined = data.frame()
exp_mat = data.frame()
for (celltype in rownames(decon_all_data@prism@phi_cellType@phi))
{ 
    tmp = get.exp(
            bp=decon_all_data,
            state.or.type="type",
            cell.name=celltype
            )
    
    meta_info_ordered = meta_info[match(rownames(tmp), meta_info$X), ]
    meta_info_ordered$celltype = celltype
    rownames(meta_info_ordered) = paste0(rownames(meta_info_ordered), "_", celltype)
    rownames(tmp) = paste0(rownames(tmp), "_", celltype)
    exp_mat = rbind(exp_mat, tmp)
    meta_combined = rbind(meta_combined, meta_info_ordered)
}

exp_mat = round(exp_mat)

# Save the pseudo matrix
# Save the data
sparse_matrix <- as(as.matrix(exp_mat), "sparseMatrix")

writeMM(obj = sparse_matrix, file = "BayesPrism_pseudo_matrix/pseudo_sn_matrix.mtx")

# Write the gene names
write.csv(data.frame(colnames(exp_mat)), "BayesPrism_pseudo_matrix/varname.csv")

# Write the barcode names
write.csv(data.frame(rownames(exp_mat)), "BayesPrism_pseudo_matrix/barcode.csv")

# Write the meta info
write.csv(meta_combined, "BayesPrism_pseudo_matrix/meta_data.csv")

# Convert the "End" into VAS
meta_combined = read.csv("BayesPrism_pseudo_matrix/meta_data.csv")

meta_combined = meta_combined %>%
  mutate(celltype = ifelse(celltype == "End", "VAS", celltype))
write.csv(meta_combined, "BayesPrism_pseudo_matrix/meta_data.csv")

# Check the results
matrix = readMM("BayesPrism_pseudo_matrix/pseudo_sn_matrix.mtx")
barcode = read.csv("BayesPrism_pseudo_matrix/barcode.csv")
gene = read.csv("BayesPrism_pseudo_matrix/varname.csv")

meta_data = read.csv("BayesPrism_pseudo_matrix/meta_data.csv")
meta_data$X = paste0(meta_data$X, "_", meta_data$celltype)

rownames(matrix) = meta_data$X
colnames(matrix) = gene$colnames.exp_mat.

# Creat the Seurat object
pbmc <- CreateSeuratObject(counts = t(matrix), project = "pbmc3k", min.cells = 3, min.features = 100)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:20)

meta_data = subset(meta_data, subset = X %in% colnames(pbmc))

for (col in colnames(meta_data))
{
    pbmc[[col]] = meta_data[[col]]
}

DimPlot(pbmc, reduction = "umap", group.by = "celltype")

# SaveRDS
saveRDS(pbmc, "BayesPrism_pseudo_matrix/pseudo_spot.rds")