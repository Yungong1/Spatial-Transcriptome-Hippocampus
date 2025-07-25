def gene_mean_exp(table, data):
    list1 = table["genes"].tolist()
    list2 = data.gene_names.tolist()
    indices = [list2.index(item) for item in list1 if item in list2]
    table["mean_exp"] = data.exp_matrix[:, indices].mean(axis=0).tolist()[0]                  

    # print the gene expression
def DE_table(data, gene_numbers):
    import pandas as pd
    import numpy as np
    result = data.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    table = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']}).head(gene_numbers)
    # Assuming 'table' is your DataFrame
    # Separating gene names, values, and P-values
    gene_names_df = table.iloc[:, ::3]  # Select columns with gene names
    values_df = table.iloc[:, 1::3]     # Select columns with values
    p_values_df = table.iloc[:, 2::3]   # Select columns with P-values

    # Reset index to keep the index column for melting
    gene_names_df = gene_names_df.reset_index()
    values_df = values_df.reset_index()
    p_values_df = p_values_df.reset_index()

    # Melting the DataFrames
    melted_gene_names = pd.melt(gene_names_df, id_vars=['index'], var_name='Group_key', value_name='Gene_name')
    melted_values = pd.melt(values_df, id_vars=['index'], var_name='Group_key', value_name='Value')
    melted_p_values = pd.melt(p_values_df, id_vars=['index'], var_name='Group_key', value_name='P_adjusted')

    # Adjusting the 'Group_key' to have consistent group names in all DataFrames
    melted_gene_names['Group_key'] = melted_gene_names['Group_key'].str.replace('_n', '')
    melted_values['Group_key'] = melted_values['Group_key'].str.replace('_l', '')
    melted_p_values['Group_key'] = melted_p_values['Group_key'].str.replace('_p', '')
    #print(melted_p_values['Group_key'])
    # Combining the melted DataFrames
    long_table = pd.DataFrame({
        'index': melted_gene_names['index'],
        'Group_key': melted_gene_names['Group_key'],
        'Gene_name': melted_gene_names['Gene_name'],
        'Value': melted_values['Value'],
        'P_adjusted': melted_p_values['P_adjusted']
    })
    
    matrix = data.to_df()
    long_table["mean_exp"] = np.mean(matrix[long_table["Gene_name"].tolist()]).tolist()

    #mean_exp = []
    #for genes in long_table["Gene_name"]:
    #    tmp = np.mean(data[:, genes].X)
    #    mean_exp.append(tmp)
    #long_table["Mean_exp"] = mean_exp
    
    return long_table

def feature_plot(data, gene, dot_size = 35):
    """
    Plot the spatial gene expression 
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    df = data.obs.copy()
        
    
    matrix = data.to_df().copy()
    df["gene_exp"] = matrix[gene].tolist()

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        df["col"], 
        df["row"], 
        c=df["gene_exp"], 
        cmap="viridis", 
        s=dot_size, 
        edgecolor="k", 
        alpha=0.8
    )
        
        
        
    # Add a colorbar
    cbar = plt.colorbar(sc)
    cbar.set_label("Gene Expression Level", fontsize=12)

    # Add labels and title
    plt.title("Scatter Plot of Rows vs Columns Colored by Gene Expression", fontsize=14)
    plt.xlabel("col (X-axis)", fontsize=12)
    plt.ylabel("row (Y-axis)", fontsize=12)

    # Show the plot
    plt.grid(True, linestyle='--', alpha=0.7)
    # Remove axes (coordinates and labels)
    plt.axis("off")


def cluster_plot(data, group_by="cluster", dot_size=35):
    import pandas as pd
    import matplotlib.pyplot as plt

    # Define a distinct color palette
    color_palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2",
        "#c7c7c7", "#dbdb8d", "#9edae5", "#1f77b4", "#aec7e8"
    ]
    unique_groups = data.obs[group_by].unique()

    # Copy data
    df = data.obs.copy()
    df[group_by] = df[group_by].astype("category")  # Ensure the selected group_by column is categorical

    # Create the scatter plot
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        df["col"],
        df["row"],
        c=df[group_by].cat.codes.map(lambda x: color_palette[x % len(color_palette)]),
        s=dot_size,
        edgecolor="k",
        alpha=0.8,
    )

    # Add a legend
    handles = [
        plt.Line2D([0], [0], marker='o', color=color_palette[i % len(color_palette)], 
                   linestyle='', markersize=8, label=str(cat))
        for i, cat in enumerate(df[group_by].cat.categories)
    ]
    plt.legend(handles=handles, title=group_by.capitalize(), bbox_to_anchor=(1.05, 1), loc="upper left")

    # Customize layout
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.axis("off")


    
def DE_subregion(data, annotation, region, group, ref):
    import scanpy as sc
    import squidpy as sq
    import pandas as pd
    import matplotlib.pyplot as plt
    import anndata as ad

    
    # region specific DE in part and control
    sub_region = data[data.obs[annotation] == region]
    
    sub_region = sub_region[sub_region.obs["diagnosis"].isin([group, ref])]
    # upregulate
    sc.tl.rank_genes_groups(sub_region, groupby="diagnosis", method="wilcoxon")
    result_table = DE_table(sub_region, gene_numbers = 30000)
    result_table = result_table[(result_table["Value"]>0.25) & (result_table["P_adjusted"]<0.05) & (result_table["mean_exp"]>0.25)]

    return result_table

def gene_detected_plot(data):
    
    """
    Plot the genes detected in each spot
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    
    df = pd.DataFrame(data.obs[["n_genes_by_counts", "col", "row"]])

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(df['col'], df['row'], c=df['n_genes_by_counts'], cmap='viridis', s=35, edgecolor='k')
    plt.colorbar(scatter, label='n_genes_by_counts')
    plt.xlabel('Column (col)')
    plt.ylabel('Row (row)')
    plt.title('Scatter Plot Colored by n_genes_by_counts')
    plt.gca().invert_yaxis()  # Invert y-axis to match common image orientation

def GO_enrichment(gene_table):
    import enrichrpy
    import enrichrpy.enrichr as een
    import enrichrpy.plotting as epl

    GO_result = een.get_pathway_enrichment(gene_table["Gene_name"], gene_set_library='GO_Biological_Process_2023')
    GO_result = GO_result[GO_result["Adjusted p-value"] < 0.05]
    return GO_result

def GO_plot(GO_result):
    import pandas as pd
    import matplotlib.pyplot as plt

    # Calculate gene count for each GO term
    GO_result["Gene count"] = [len(genes) for genes in GO_result["Overlapping genes"]]

    # Create a DataFrame
    df = pd.DataFrame(GO_result)

    # Remove the (GO:XXXX) portion from each term name
    df["Term name"] = df["Term name"].str.replace(r" \(GO:\d+\)", "", regex=True)

    # Sort the DataFrame by Gene count in ascending order
    df = df.sort_values(by="Gene count", ascending=True)

    # Set up the figure
    plt.figure(figsize=(7, 4))  # Adjust the figure size

    # Create a scatter plot (dot plot)
    scatter = plt.scatter(
        df["Gene count"], 
        df["Term name"], 
        c=df["Adjusted p-value"], 
        cmap="viridis", 
        s=50  # Adjust the dot size
    )

    # Add a color bar to indicate Adjusted p-value
    cbar = plt.colorbar(scatter)
    cbar.set_label("Adjusted P-value")

    # Add axis labels and title
    plt.xlabel("Number of Overlapping Genes")
    plt.ylabel("GO Term")
    plt.title("GO Term Enrichment - Dot Plot")
    plt.tight_layout()

def overlap_genes(file_path, region_list, comparison, group):
    
    import pandas as pd
    # Calculate the overlap genes

    up_gene_dict = {}
    for region in region_list:
        tmp = pd.read_csv("{}/{}/{}.csv".format(file_path, region, comparison))
        tmp1 = tmp[tmp["Group_key"] == group]["Gene_name"]
        up_gene_dict[region] = tmp1

    common_genes = set.intersection(*(set(genes) for genes in up_gene_dict.values()))
    return common_genes

def unique_up_regulate(data, file_path, group_list, group_target, comparison_group, check_region):
    import pandas as pd
    merge = data
    # Calculate the unique upregulated genes in each group
    up_gene_dict = {}
    for group in group_list:
        tmp = pd.read_csv("{}/{}/{}.csv".format(file_path, group, comparison_group))
        tmp = tmp[tmp["Group_key"] == group_target]["Gene_name"]
        up_gene_dict[group] = tmp
    pyn_layer = group_list

    overlap_genes = None  # Start with None to initialize with the first set
    all_genes = set()  # To collect all unique genes

    # Dictionary to store unique genes for each group
    unique_genes_by_group = {}

    # Iterate through each group in pyn_layer
    for group in pyn_layer:
        if group in up_gene_dict:
            # Extract the gene list for this table (assuming 'Gene' column contains gene names)
            current_gene_set = set(up_gene_dict[group])  # Replace "Gene" with the actual column name

            # Collect genes from other groups
            other_genes = set()
            for other_group in pyn_layer:
                if other_group != group and other_group in up_gene_dict:
                    other_genes |= set(up_gene_dict[other_group])  # Union of other group genes

            # Find unique genes for the current group
            unique_genes_by_group[group] = current_gene_set - other_genes


    up_gene_table = {}
    for group in merge.obs["major_annotation"].unique():
        tmp = pd.read_csv("{}/{}/{}.csv".format(file_path, group, comparison_group))
        tmp = tmp[tmp["Group_key"] == group_target]
        up_gene_table[group] = tmp

    # Dictionary to store ordered unique genes for each group
    ordered_unique_genes_by_group = {}

    # Iterate through each group in unique_genes_by_group
    for group, unique_genes in unique_genes_by_group.items():
        if group in up_gene_table:  # Check if the group exists in the table list
            # Filter the up_gene_table for the current group
            table = up_gene_table[group]

            # Filter rows where 'Gene_name' is in the unique genes set
            filtered_table = table[table["Gene_name"].isin(unique_genes)]

            # Sort by 'Value' column (descending order)
            sorted_table = filtered_table.sort_values(by="Value", ascending=False)

            # Store the ordered list of genes
            ordered_unique_genes_by_group[group] = sorted_table

    unique_gene_tmp = ordered_unique_genes_by_group.get(check_region)
    
    return unique_gene_tmp

def GO_plot(GO_result, term_number, size = (7,5)):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    # Add number of overlapping genes for plotting
    df = pd.DataFrame(GO_result[:term_number])
    df['Number of Overlapping Genes'] = df['Overlapping genes'].apply(len)
    # Sort data for better visualization
    df = df.sort_values(by="Adjusted p-value", ascending=True)
    df = df.sort_values(by="Number of Overlapping Genes", ascending=False)
    
    # Create the plot
    plt.figure(figsize=size)
    barplot = sns.barplot(
        data=df,
        x="Number of Overlapping Genes",
        y="Term name",
        hue="Adjusted p-value",
        dodge=False,
        palette="viridis",
    )

    # Add color bar for continuous values
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=df["Adjusted p-value"].min(), vmax=df["Adjusted p-value"].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=barplot)
    cbar.set_label("Adjusted p-value", fontsize=12)

    # Customize the plot
    plt.title("Gene Ontology Term Enrichment", fontsize=14)
    plt.xlabel("Number of Overlapping Genes", fontsize=12)
    plt.ylabel("Term name", fontsize=12)
    plt.legend([], [], frameon=False)  # Remove the default legend
    plt.tight_layout()