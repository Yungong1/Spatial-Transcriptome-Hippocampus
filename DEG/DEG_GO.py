import scanpy as sc
#import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import utils
import importlib
importlib.reload(utils)
import numpy as np
import enrichrpy
#import enrichrpy.enrichr as een
#import enrichrpy.plotting as epl
import os

# Read the data
merge = sc.read_h5ad("../Result/Region_annotation/annotated_hippo.h5ad")

# Change the annotation
annotation_map = {
    'SUB': 'SUB',
    'CA1.SP': 'CA1',
    'SO': 'SO',
    'CA1.SR': 'SR',
    'ML': 'ML',
    'CA4': 'CA4',
    'VAS': 'VAS',
    'DG.GCL': 'DG',
    'CA3': 'CA3',
    'CA2.SP': 'CA2'
}
# Apply the mapping
merge.obs['major_annotation'] = merge.obs['major_annotation'].replace(annotation_map)

# Check the DE genes
sc.tl.rank_genes_groups(merge, groupby="major_annotation", method = "wilcoxon")

DE_table = utils.DE_table(merge, gene_numbers = 3000)
DE_table = DE_table[(DE_table["Value"] > 0.5) & (DE_table["P_adjusted"] < 0.05)]

# Save the GO result
sub_table = DE_table[DE_table["Group_key"] == "VAS"]
result = utils.GO_enrichment(sub_table)

utils.GO_plot(result[:10])
