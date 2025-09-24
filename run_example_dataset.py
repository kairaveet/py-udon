# load udon scripts
import pyudon

# load packages for satay-udon
import scanpy as sc
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt


sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")

os.chdir("/data/salomonis2/LabFiles/Kairavee/udon-python/") # this is the folder that contains the udon scripts including this script 

folder_path = "/data/salomonis2/LabFiles/Kairavee/udon_data/GEO_Submission_Processed_Files_06-29-2022"  # Path to the folder containing the .h5 files
file_extension = ".h5"

# Get a list of all files in the folder
files_in_folder = os.listdir(folder_path)

# Filter the list to include only files ending with .h5
file_list = [file for file in files_in_folder if file.endswith(file_extension)]
file_list.sort()


# Initialize an empty list to store the loaded data
adata_list = []

# Iterate over the file list with tqdm for the progress bar
for filename in tqdm(file_list, desc="Loading Data"):
    # Load data from each file and append it to the adata_list
    adata_list.append(sc.read_10x_h5(os.path.join(folder_path, filename)))
              
donor_list = ['PID201056', 'PID201087', 'PID201295', 'PID200673', 
              'PID200704', 'PID200794', 'PID200805', 'PID200609', 
              'PID200637', 'PID200657', 'PID200967', 'PID200974',
              'PID200353', 'PID200593', 'PID200594', 'PID200641',
              'PID200812', 'PID201208', 'PID200396', 'PID200622',
              'PID2007581','PID2007582', 'PID200965', 'PID200970',
              'PID201020','PID200785']

# Assuming adata_list is your list of anndata.AnnData objects
for i, i_adata in enumerate(adata_list):
    # Update the barcodes to include donor id
    i_adata.X = csr_matrix(i_adata.X)
    i_adata.obs_names = i_adata.obs_names + f".{donor_list[i]}"

    # standard scanpy/seurat like QC steps
    i_adata.var_names_make_unique()
    sc.pp.filter_cells(i_adata, min_genes=200)
    sc.pp.filter_genes(i_adata, min_cells=3)
    i_adata.var["mt"] = i_adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(i_adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    i_adata = i_adata[i_adata.obs.n_genes_by_counts < 2500, :]
    i_adata = i_adata[i_adata.obs.pct_counts_mt < 5, :]

# Concatenate the anndata objects in adata_list along the observations axis
adata = ad.concat(adata_list, axis=0, join='outer', merge='first', uns_merge='same', label='donor_ids', keys=donor_list)
print(adatas)

'''
AnnData object with n_obs × n_vars = 229707 × 22914
    obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'donor_ids'
    var: 'gene_ids', 'feature_types', 'genome', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'
'''

metadata = pd.read_csv("/data/salomonis2/LabFiles/Kairavee/udon_data/sjia_metadata.txt", sep="\t")  # 209955 cells
metadata.index = metadata.cell_id

metadata['basic_disease_type'] = metadata['disease_subtype']
metadata.loc[metadata['basic_disease_type'] == "Active", 'basic_disease_type'] = "Disease"
metadata.loc[metadata['basic_disease_type'] == "Inactive", 'basic_disease_type'] = "Disease"
metadata.loc[metadata['basic_disease_type'] == "LungDisease", 'basic_disease_type'] = "Disease"
metadata.loc[metadata['basic_disease_type'] == "MAS", 'basic_disease_type'] = "Disease"

# find common cell ids between metadata and adatas
common_cell_ids = list(set(metadata.index).intersection(adata.obs_names))  # 209425 cells
metadata_f = metadata.loc[common_cell_ids, :]

adata = adata[common_cell_ids, :]

# add metadata to adatas.obs
adata.obs = adata.obs.join(metadata_f)

print(adata)
'''
AnnData object with n_obs × n_vars = 209425 × 22914
    obs: 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'donor_ids', 'cell_id', 'cell_type', 'disease_subtype', 'patient_id', 'UMAP_1', 'UMAP_2', 'basic_disease_type'
    var: 'gene_ids', 'feature_types', 'genome', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'

'''

sc.pp.normalize_total(adata, target_sum=1e4)  # this is a necessary step
sc.pp.log1p(adata)


######################################### UDON #########################################

# add udon object -- read below before running the UDON pipeline from here on.
# calculate_pseudobulk_folds: Calculate pseudobulk fold values for each gene in the adata object

# :param adata: input ann data which contains cell by gene gene expression in adata.X
# :param groups: the column name in adata.obs that contains the cell type cluster
# :param donors: the column name in adata.obs that contains the donor id of the cell
# :param disease_type_col: the column name in adata.obs that contains which cells are “Disease” and which cells are “Healthy” for fold calculation
# :param min_cells: minimum number of cells present for calculation of any pseudobulk vector
# :param baseline:in the disease_type_col, the character string that refers to the baseline (denominator) for pseudobulk fold calculation. For udon, the baseline are cells labeled as “Control” in the disease_type_col
# :param query: : in the disease_type_col, the character string that refers to the query (numerator) for pseudobulk fold calculation. For udon, the baseline are cells labeled as “Disease” in the disease_type_col
# :param control_donor_specificity: how the baseline for pseudobulk-fold should be calculated? If “collective” (default), then the baseline (typically the control samples for UDON) will be calculated agnostic of donors.
# :return: adata with slot adata.varm['pseudobulk_folds'] containing the pseudobulk fold values

adata = calculate_pseudobulk_folds(uadata=adata, groups='cell_type', donors='patient_id', disease_type_col='basic_disease_type', min_cells=5, baseline='Control', query='Disease')

# perform feature selection on the adata.uns['pseudobulk_folds'] 
adata = feature_selection_wrapper(adata, species='human', fold_threshold=0.5, samples_differing=3, intercorr_threshold=0.3, corr_n_events=10, pca_corr_threshold=0.3, n_components=30)


print("Writing the feature metadata to a file...")
feature_metadata = adata.var
feature_metadata.to_csv("/data/salomonis2/LabFiles/Kairavee/udon_example_run/udon_feature_metadata.txt", sep="\t")

writing_dir = "/data/salomonis2/LabFiles/Kairavee/udon_example_run/"

# Define the parameter sets you want to run
parameter_sets = [
	{'rank': 15},
    {'rank': 20},
    {'rank': 25},
    {'rank': 30}
]

# Loop over each parameter set
for params in parameter_sets:
    # Create the output filename based on parameter values
    output_filename = os.path.join(writing_dir, f"udon_k{params['rank']}")

    # Run the udon clustering_wrapper function with the current parameters
    adata = clustering_wrapper(
        adata=adata, 
        output_filename=output_filename,
        rank=params['rank'], 
        min_group_size=3, 
        top_n_genes=100, 
        rho_threshold=0.2, 
        marker_finder_rho=0.2, 
        write_files=True)

    # perform pathway enrichment
    adata = enrichment_wrapper(adata=adata, species="Hs", p_val=0.1)
    
    # write the enrichment files 
    all_pathways = adata.uns['significant_pathway_enrichment'] 
    top_5_pathways = adata.uns['significant_pathway_enrichment_top_n']
    
    all_pathways.to_csv(os.path.join(output_filename, "/all_pathways.txt"), sep="\t")
    top_5_pathways.to_csv(os.path.join(output_filename, "/top_5_pathways.txt"), sep="\t")
    
    # visualize the results
    clusters_clean = adata.uns['udon_clusters']
    markers_heatmap = adata.uns['marker_heatmap']
    markers_top100 = adata.uns['udon_marker_genes_top_n']
    udon_metadata = clusters_clean
    udon_metadata['cell_type'] = udon_metadata.index.str.split('__').str[0]
    udon_metadata['donor'] = udon_metadata.index.str.split('__').str[1]
    
    plot_metadata_heatmap(udon_clusters=clusters_clean, udon_metadata=udon_metadata, metadata_col='cell_type', path_to_save_figure=os.path.join(output_filename, "/udon_cell_type.pdf"))
    plot_metadata_heatmap(udon_clusters=clusters_clean, udon_metadata=udon_metadata, metadata_col='donor', path_to_save_figure=os.path.join(output_filename, "/udon_donor.pdf"))
    plot_markers_df(marker_heatmap=markers_heatmap, markers_df=markers_top100, clusters=clusters_clean, path_to_save_figure=os.path.join(output_filename, "/udon_marker_heatmap.pdf"))
    
    
    # save the final anndata object
    adata.write_h5ad(os.path.join(output_filename, "/adata_udon.h5ad"), compression='lzf')


######################################### SATAY-UDON #########################################

runs = os.listdir(writing_dir)
runs.remove("udon_feature_metadata.txt")

# requires binary clinical data so we can convert categorical data into binary using pd.get_dummies function or use make_mean_based_binary from satay_udon.py
clinical_metadata = pd.read_excel("//data/salomonis2/LabFiles/Kairavee/udon_data/donor_metadata.xlsx",index_col=0)
df_binary = pd.get_dummies(clinical_metadata, columns=['Disease status','Sex','Systemic Features','Active arthritis','Fever'])
df_binary = df_binary.drop(columns=['Sample','Age_yrs'])
df_binary = df_binary.astype(int)
df_numeric_binary = make_mean_based_binary(df=clinical_metadata, cols_to_convert=['Age_yrs'])

clinical_metadata_binary = pd.concat([df_numeric_binary, df_binary], axis=1)

for r in runs:
    os.chdir(os.path.join(writing_dir, r))
    udon_clusters = pd.read_csv("udon_clusters.txt", sep="\t", index_col=0)
    print(r)

    p_vals_dict = dict()

    for c in clinical_metadata_binary.columns:
        p_val_matrix_c = fishers_clinical_feats(udon_clusters=udon_clusters, clinical_metadata_df=clinical_metadata_binary, clinical_measure_key=c, p_val=0.1, n_samples=3)
        p_val_matrix_c = p_val_matrix_c['p_val']

        p_val_matrix_c[p_val_matrix_c > 0.1] = 1

        mask_red = p_val_matrix_c < 0.1

        if mask_red.sum().sum() == 0:
            print(f"No significant p-values for {c}")
            continue

        p_vals_dict[c] = p_val_matrix_c

        # Create the heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(p_val_matrix_c,  annot=True, fmt=".2f",
                    cbar=False, linewidths=.5, linecolor='black',
                    cmap='rocket_r')

        plt.title(f'Heatmap for {c}')
        plt.tight_layout()

        # Save the heatmap as a PDF
        plt.savefig(f'{c}_pval_heatmap_robust_colored.pdf')
        plt.close()
        
        