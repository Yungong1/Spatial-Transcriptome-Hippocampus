library(Seurat)
suppressWarnings(library(BayesPrism))
library(dplyr)
library(tidyverse)


decon = function(sc_count_filter, sc_meta_filter, merged)
{
    sc.dat = t(sc_count_filter)

    cell.type.labels = sc_meta_filter$major_cell_type

    #cell.state.labels = c(rep("state1", 1000), rep("state2", nrow(sc.dat)-1000))
    cell.state.labels = sc_meta_filter$major_cell_type

    sc.stat <- plot.scRNA.outlier(
        input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
        cell.type.labels=cell.type.labels,
        species="hs", #currently only human(hs) and mouse(mm) annotations are supported
        return.raw=TRUE #return the data used for plotting.
    #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
    )

    bk.stat <- plot.bulk.outlier(
        bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
        sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
        cell.type.labels=cell.type.labels,
        species="hs", #currently only human(hs) and mouse(mm) annotations are supported
        return.raw=TRUE
        #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
    )
    
    sc.dat.filtered <- cleanup.genes(
        input=sc.dat,
        input.type="count.matrix",
        species="hs",
        gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
        exp.cells=5
    )
    
    #plot.bulk.vs.sc(
    #    sc.input = sc.dat.filtered,
    #    bulk.input = bk.dat
        #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
    #)
    
    #sc.dat.filtered.pc <- select.gene.type(
    #    sc.dat.filtered,
    #    gene.type = "protein_coding"
    #)

    print(length(cell.type.labels) == length(cell.state.labels))
    bk.dat = as.matrix(bk.dat)
    sc.dat = as.matrix(sc.dat)
    sc.dat.filtered = as.matrix(sc.dat.filtered)
    
    myPrism <- new.prism(
        reference=sc.dat.filtered,
        mixture=bk.dat,
        input.type="count.matrix",
        cell.type.labels = cell.type.labels,
        cell.state.labels = cell.state.labels,
        key= NULL,
        outlier.cut=1,
        outlier.fraction=1,
    )

    bp.res <- run.prism(prism = myPrism, n.cores=48)
    return(bp.res)

}

random_select = function(sc_count, sc_meta, prop)
{
    # Random select the cells in each cluster
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

# Read the spatial data

merged = readRDS("../../../Result/Region_annotation/processed_merged_spatial.rds")

# read the meta
meta = read.csv("../../../Result/Raw_data_for_R/meta_data.csv")

merged[["annotation"]] = meta$annotation

bk.dat = t(merged@assays$SCT@counts)

# read the count and meta data for sc ref
sc_meta = readRDS("../Reference/origin_sc_meta.rds")
sc_count = readRDS("../Reference/origin_sc_count.rds")

# Random select the data
result = random_select(sc_count = sc_count, sc_meta = sc_meta, prop = 0.5)

# Perform the bayesPrism for two data
decon_result_test = decon(result[[1]], result[[2]], merged)
decon_result_val = decon(result[[3]], result[[4]], merged)

# Save the data
saveRDS(decon_result_test, "Test1/50_decon_result_test.RDS")
saveRDS(decon_result_val, "Test1/50_decon_result_val.RDS")
