import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
#from cmh import CMH


def make_mean_based_binary(df, cols_to_convert, direction='greater'):
    """
    Convert a dataframe to binary based on the mean value of the variable.
    :param df: pandas DataFrame
    :param cols_to_convert: list of columns to convert to binary
    :param direction: if mean above threshold, set to 1, else 0; setting to "lesser" does the opposite --
    use 'lesser' option when clinical var is more interesting when value is lower
    :return: binary pandas DataFrame
    """

    df = df[cols_to_convert]

    # Clean up the numeric clinical data to make it binary based on average values of clinical parameters
    numeric_data = df.apply(pd.to_numeric, errors='coerce')

    # Calculate column means
    colmeans_clinical_data = numeric_data.mean(axis=0, skipna=True)

    if direction == 'lesser':
        df_binary = (numeric_data <= colmeans_clinical_data).astype(int)
    else:
        df_binary = (numeric_data >= colmeans_clinical_data).astype(int)

    return df_binary


### HOPE's version
def fishers_clinical_feats(udon_clusters, clinical_metadata_df, clinical_measure_key,  p_val=0.1, n_samples=3, patient_key="HTB_ID"):
    """
    Runs the fisher's test on binarized clinical data focusing on the clinical_measure_key

    udon_clusters: Udon clusters that were output to udon_clusters.txt
    clincial_metadata_df: already binarized metadata
    clinical_measure_key: ONE variable to from the metadata consider for testing
    p_val: p-value cutoff to consider something significant (default 0.1)
    n_samples: minimum number of samples needed to consider something singificant (default 3)
    patient_key: column for patient (as linking to udon_clusters)

    Returns: p-value matrix for each of the clusters-celltypes where a p-value indicates if that cluster is
        enriched in the given clinical_measure_key=1. P-values are changed to NaN if # samples < n_samples
    """

    ## PROJECT the donor metadata to the udon clusters
    donor_to_clinical_measure = clinical_metadata_df[clinical_measure_key].to_dict()
    udon_clusters['clinical_measure'] = udon_clusters['donors'].map(donor_to_clinical_measure)

    # Only consider where we have patient information
    udon_clusters_updated = udon_clusters.dropna(subset="clinical_measure", axis=0)
    print(clinical_measure_key, ":", len(udon_clusters["donors"].unique()) - len(udon_clusters_updated["donors"].unique()), "Patients could not be considered due to lack of clinical data, leaving", 
        len(udon_clusters_updated["donors"].unique()))

    cell_types = udon_clusters_updated['cell_types'].unique()
    cell_types.sort()#
    clusters = udon_clusters_updated['cluster'].unique()


    # Initialize matrices for p-values and odds ratios
    p_val_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)
    odds_ratio_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)

    # Loop through each cluster and cell type
    for cluster in clusters:
        for cell_type in cell_types:
            # Subset data for the current cluster and other clusters
            pseudobulks_in_cluster = udon_clusters_updated[udon_clusters_updated['cluster'] == cluster]
            pseudobulks_in_other_clusters = udon_clusters_updated[udon_clusters_updated['cluster'] != cluster]

            # Compute the counts for the Fisher's test
            n_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 1) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
            n_not_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 0) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
            n_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 1) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])
            n_not_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 0) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])

            # Create a contingency table for Fisher's test
            fisher_table = np.array([[n_cm_ct_in_cluster, n_not_cm_ct_in_cluster],
                                        [n_cm_ct_in_other_clusters, n_not_cm_ct_in_other_clusters]])

            # Perform Fisher's exact test
            odds_ratio, p_value = fisher_exact(fisher_table, alternative='greater')

            # Store the results
            p_val_matrix.at[cell_type, cluster] = p_value
            odds_ratio_matrix.at[cell_type, cluster] = odds_ratio

            # Apply additional conditions based on p-value and sample size
            if p_value < p_val and fisher_table[0, 0] < n_samples:
                p_val_matrix.at[cell_type, cluster] = np.nan # Set p-value to NaN if sample size is less than threshold
        
    return p_val_matrix



### HOPE's version
def prep_metadata_fortest(clinical_metadata_df, clinical_measure_key, categorical_variables, numerical_variables, patient_key="HTB_ID"):
    """
    Gets the metadata properly prepped for analysis for the given item. 
    Specifically it removes cases where the patient has an NA value or other notes 
    for unknown. It then binarizes the variable.

    Final dataframe should only include the columns of interest

    WARNING: dummy variables mean that a single variable --> more than one
    """
    # Remove rows with NA in the specified column
    metadata_df = clinical_metadata_df.dropna(subset=clinical_measure_key, axis=0)
    metadata_df.index = metadata_df[patient_key]
    # Remove cases where not sure of info
    metadata_df = metadata_df[~metadata_df[clinical_measure_key].isin(["NA/Unk", "Not Assigned", "UNKNOWN", "unknown", "INDETERMINANT"])]
    # only consider the variable of interest so can consider all columns in case dummy
    metadata_df = metadata_df[clinical_measure_key].to_frame()
    # Get the binary variables
    metadata_df = get_binary_variable(metadata_df, clinical_measure_key, categorical_variables, numerical_variables)
    return(metadata_df)


### HOPE's version
def get_binary_variable(metadata_df, clinical_measure_key, categorical_variables, numerical_variables):
    num_options = len(metadata_df[clinical_measure_key].unique())
    # if only one or 0 options then skip and return an empty dataframe
    if num_options < 2:
        print(clinical_measure_key, "has only zero or one option after cleaning:", metadata_df[clinical_measure_key].unique())
        return pd.DataFrame()
    elif set(metadata_df[clinical_measure_key].unique()) == {0, 1}:
        next
    elif (num_options == 2) and (clinical_measure_key in categorical_variables):
        # if binary categorical variable have the Yes/Y be a 1 and No/N be a 0
        metadata_df[clinical_measure_key] = metadata_df[clinical_measure_key].replace(["Y", "Yes", "YES", True, "TRUE", "True"], 1)
        metadata_df[clinical_measure_key] = metadata_df[clinical_measure_key].replace(["N", "No", "NO", False, "FALSE", "False"], 0)
    elif (num_options > 2) and (clinical_measure_key in categorical_variables):
        # if nonbinary categorical variable, the use pandas to get dummies
        metadata_df = pd.get_dummies(metadata_df, columns=[clinical_measure_key])
        for col in metadata_df.columns:
            if clinical_measure_key in col:
                metadata_df[col] = metadata_df[col].astype(int)
    elif clinical_measure_key in numerical_variables:
        # if numerical variable, use the UDON function to split into binary based on means
        metadata_df = make_mean_based_binary(df=metadata_df, cols_to_convert=[clinical_measure_key])
    else:
        raise Exception("The variable does not reach any requirements:", clinical_measure_key, "and has num options being", num_options)
    # check that final variable is binary
    #if set(metadata_df[clinical_measure_key].unique()) != {0, 1}:
        #raise Exception("Unsuccessful conversion to a binary variable:", metadata_df[clinical_measure_key].unique())
    #else:
    return(metadata_df)


    

# def cmh_clinical_feats(clinical_metadata_df, clinical_measure_key, batch_key, udon_clusters=None, adata=None, udon_clusters_key='udon_clusters', n_samples=3):

#     # Check input validity
#     if udon_clusters is None and adata is None:
#         raise ValueError("Either udon_clusters or adata must be provided.")
#     if udon_clusters is not None and adata is not None:
#         raise ValueError("Provide only one of udon_clusters or adata, not both.")
#     if adata is not None and udon_clusters_key is None:
#         raise ValueError("udon_cluster_key must be provided when adata is used.")

#     # Get udon clusters depending on the input type
#     if adata is not None:
#         udon_clusters = adata.uns[udon_clusters_key]

#     udon_clusters['cell_types'] = udon_clusters.index.str.split('__').str[0]
#     cell_types = udon_clusters['cell_types'].unique()
#     cell_types.sort()
#     clusters = udon_clusters['cluster'].unique()
#     udon_clusters['donors'] = udon_clusters.index.str.split('__').str[1]

#     # Remove rows with NA in the specified column
#     metadata_df = clinical_metadata_df.dropna(subset=clinical_measure_key, axis=0)

#     # project the donor metadata to the udon clusters
#     # create donor to clinical measure dictionary
#     donor_to_clinical_measure = metadata_df[clinical_measure_key].to_dict()
#     donor_to_batch = metadata_df[batch_key].to_dict()
#     udon_clusters['clinical_measure'] = udon_clusters['donors'].map(donor_to_clinical_measure)
#     udon_clusters['batches'] = udon_clusters['donors'].map(donor_to_batch)
#     donor_to_batch = metadata_df[batch_key].to_dict()
#     udon_clusters['batches'] = udon_clusters['donors'].map(donor_to_batch)

#     batches = udon_clusters['batches'].unique()

#     # Initialize matrices for p-values and odds ratios
#     p_val_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)

#     udon_clusters_og = udon_clusters.copy()

#     for cluster in clusters:
#         for cell_type in cell_types:
#             print(cluster, cell_type)
#             contigencies = np.ndarray(shape=(2,2,2), dtype=int)

#             for batch_component in batches:
#                 udon_clusters = udon_clusters_og.copy()
#                 print(batch_component)

#                 # limit analysis to the individual batch
#                 udon_clusters = udon_clusters[udon_clusters['batches'] == batch_component]

#                 # Subset data
#                 pseudobulks_in_cluster = udon_clusters[udon_clusters['cluster'] == cluster]
#                 pseudobulks_in_other_clusters = udon_clusters[udon_clusters['cluster'] != cluster]

#                 # Create contingency table
#                 # Compute the counts for the Fisher's test
#                 n_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 1) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
#                 n_not_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 0) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
#                 n_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 1) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])
#                 n_not_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 0) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])

#                 # Create a contingency table for Fisher's test
#                 contigency_mat = np.array([[n_cm_ct_in_cluster, n_not_cm_ct_in_cluster],
#                                          [n_cm_ct_in_other_clusters, n_not_cm_ct_in_other_clusters]])

#                 contigencies[batch_component] = contigency_mat

#             # Determine if CMH or Fisher's test should be used
#             if np.shape(contigencies)[2] < 2 or np.any(contigencies[0, 0, :] < n_samples):
#                 p_value = np.nan
#             else:
#                 udon_clusters['cm_ct_binary'] = 'no'
#                 udon_clusters.loc[(udon_clusters['clinical_measure'] == 1 & udon_clusters['cell_types'] == cell_type), 'cm_ct_binary'] = 'yes'
#                 udon_clusters['cluster_binary'] = 'no'
#                 udon_clusters.loc[udon_clusters['cluster'] == cluster, 'cluster_binary'] = 'yes'

#                 _, _, p_value = CMH(udon_clusters, 'cm_ct_binary', 'cluster_binary', stratifier='batches', raw=True)

#             # Store results in matrices
#             p_val_matrix.at[cell_type, cluster] = p_value

#     return p_val_matrix



def fdr_correction(p_val_matrix, alpha=0.05, method='fdr_bh'):
    """
    Perform multiple hypothesis correction using the Benjamini-Hochberg method.
    :param p_val_matrix: pandas DataFrame of p-values
    :param alpha: significance level
    :param method: method for multiple hypothesis correction
    :return: pandas DataFrame of corrected p-values
    """
    # Flatten the p-value matrix
    p_vals_flat = p_val_matrix.values.flatten()

    # Perform FDR correction
    _, p_vals_corrected, _, _ = multipletests(p_vals_flat, alpha=alpha, method=method)

    # Reshape the corrected p-values
    p_vals_corrected_reshaped = p_vals_corrected.reshape(p_val_matrix.shape)

    # Convert the corrected p-values to a DataFrame
    p_vals_corrected_df = pd.DataFrame(p_vals_corrected_reshaped, index=p_val_matrix.index, columns=p_val_matrix.columns)

    return p_vals_corrected_df


def satay_udon(clinical_metadata_df, clinical_measure_keys, batch_key, udon_clusters=None, adata=None, udon_clusters_key='udon_clusters', p_val=0.1, n_samples=3):

    if batch_key is None:
        print("No batch key provided. Running Fisher's exact test.")
        for clinical_measure_key in clinical_measure_keys:

            stats = fishers_clinical_feats(clinical_metadata_df, clinical_measure_key, udon_clusters=udon_clusters, adata=adata,udon_clusters_key=udon_clusters_key, p_val=p_val, n_samples=n_samples)
            p_val_matrix = stats['p_val']

    # cmh_clinical_feats

    # fdr_correction

    return adata


## Testing out CODE
# metadata_practice_df = pd.DataFrame({"HTB_ID":[1, 2, 3, 4, 5, 6, ], 
#                             "Binary":["Yes", "No", "No", "Yes", "Yes", "Yes"], 
#                             "Binary_Num":[1, 0, 0, 1, 1, 1], 
#                             "Cat_Three": ["Adverse", "Favorable", "Intermediate", "Intermediate", 
#                             "Favorable", "Adverse"], 
#                             "Cat_Four": ["Option1", "Option2", "Option3", "Option4", "Option1", "Option2"], 
#                             "Numerical": [0.4, 0.5, 0.1, 0.3, 0.6, 0.7], 
#                             "NanAll": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
# })
# metadata_practice_df.index = metadata_practice_df["HTB_ID"]
# print(metadata_practice_df)
# results = get_binary_variable(metadata_practice_df, "NanAll", ["Binary", "Cat_Three", "Cat_Four"], ["Binary_Num", "Numerical"])
# print(results )
# if results.empty:
#     print("EMPTY")

# results = prep_metadata_fortest(metadata_practice_df, "Cat_Three", ["Binary", "Cat_Three", "Cat_Four"], ["Binary_Num", "Numerical"])