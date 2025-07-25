library(Seurat)
suppressWarnings(library(BayesPrism))

decon = function(sc_count_filter, sc_meta_filter, bk.dat)
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

    bp.res <- run.prism(prism = myPrism, n.cores=64)
    return(bp.res)

}

# Read the spatial data
merged = readRDS("../../Result/Region_annotation/processed_merged_spatial.rds")

# read the meta
meta = read.csv("../../Result/Raw_data_for_R/meta_data.csv")

merged[["annotation"]] = meta$annotation

bk.dat = t(merged@assays$SCT@counts)

# read the count and meta data for sc ref
sc_meta = readRDS("Reference/origin_sc_meta.rds")
sc_count = readRDS("Reference/origin_sc_count.rds")

overlap = intersect(colnames(bk.dat), rownames(sc_count))
sc_count = sc_count[overlap, ]
bk.dat = bk.dat[,overlap]

#result = random_select(sc_count = sc_count, sc_meta = sc_meta, per = 0.3)
decon_result = decon(sc_count, sc_meta, bk.dat)
path = "Decon_result/BayesPrism/decon_result_allgenes.rds"
saveRDS(decon_result, path)