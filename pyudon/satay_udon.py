import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


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


def fishers_clinical_feats(clinical_metadata_df, clinical_measure_key, udon_clusters=None, adata=None, udon_clusters_key='udon_clusters', p_val=0.1, n_samples=3):

    # Check input validity
    if udon_clusters is None and adata is None:
        raise ValueError("Either udon_clusters or adata must be provided.")
    if udon_clusters is not None and adata is not None:
        raise ValueError("Provide only one of udon_clusters or adata, not both.")
    if adata is not None and udon_clusters_key is None:
        raise ValueError("udon_cluster_key must be provided when adata is used.")

    # Get udon clusters depending on the input type
    if adata is not None:
        udon_clusters = adata.uns[udon_clusters_key]

    udon_clusters['cell_types'] = udon_clusters.index.str.split('__').str[0]
    cell_types = udon_clusters['cell_types'].unique()
    cell_types.sort()
    clusters = udon_clusters['cluster'].unique()
    udon_clusters['donors'] = udon_clusters.index.str.split('__').str[1]

    # Remove rows with NA in the specified column
    metadata_df = clinical_metadata_df.dropna(subset=clinical_measure_key, axis=0)

    # project the donor metadata to the udon clusters
    # create donor to clinical measure dictionary
    donor_to_clinical_measure = metadata_df[clinical_measure_key].to_dict()
    udon_clusters['clinical_measure'] = udon_clusters['donors'].map(donor_to_clinical_measure)

    # Initialize matrices for p-values and odds ratios
    p_val_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)
    odds_ratio_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)

    # Loop through each cluster and cell type
    for cluster in clusters:
        for cell_type in cell_types:

            # Subset data for the current cluster and other clusters
            pseudobulks_in_cluster = udon_clusters[udon_clusters['cluster'] == cluster]
            pseudobulks_in_other_clusters = udon_clusters[udon_clusters['cluster'] != cluster]

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

    return {'p_val': p_val_matrix, 'OR': odds_ratio_matrix}


def cmh_clinical_feats(clinical_metadata_df, clinical_measure_key, batch_key, udon_clusters=None, adata=None, udon_clusters_key='udon_clusters', n_samples=3):

    # Check input validity
    if udon_clusters is None and adata is None:
        raise ValueError("Either udon_clusters or adata must be provided.")
    if udon_clusters is not None and adata is not None:
        raise ValueError("Provide only one of udon_clusters or adata, not both.")
    if adata is not None and udon_clusters_key is None:
        raise ValueError("udon_cluster_key must be provided when adata is used.")

    # Get udon clusters depending on the input type
    if adata is not None:
        udon_clusters = adata.uns[udon_clusters_key]

    udon_clusters['cell_types'] = udon_clusters.index.str.split('__').str[0]
    cell_types = udon_clusters['cell_types'].unique()
    cell_types.sort()
    clusters = udon_clusters['cluster'].unique()
    udon_clusters['donors'] = udon_clusters.index.str.split('__').str[1]

    # Remove rows with NA in the specified column
    metadata_df = clinical_metadata_df.dropna(subset=clinical_measure_key, axis=0)

    # project the donor metadata to the udon clusters
    # create donor to clinical measure dictionary
    donor_to_clinical_measure = metadata_df[clinical_measure_key].to_dict()
    donor_to_batch = metadata_df[batch_key].to_dict()
    udon_clusters['clinical_measure'] = udon_clusters['donors'].map(donor_to_clinical_measure)
    udon_clusters['batches'] = udon_clusters['donors'].map(donor_to_batch)
    donor_to_batch = metadata_df[batch_key].to_dict()
    udon_clusters['batches'] = udon_clusters['donors'].map(donor_to_batch)

    batches = udon_clusters['batches'].unique()

    # Initialize matrices for p-values and odds ratios
    p_val_matrix = pd.DataFrame(np.nan, index=cell_types, columns=clusters)

    udon_clusters_og = udon_clusters.copy()

    for cluster in clusters:
        for cell_type in cell_types:
            print(cluster, cell_type)
            contigencies = np.ndarray(shape=(2,2,2), dtype=int)

            for batch_component in batches:
                udon_clusters = udon_clusters_og.copy()
                print(batch_component)

                # limit analysis to the individual batch
                udon_clusters = udon_clusters[udon_clusters['batches'] == batch_component]

                # Subset data
                pseudobulks_in_cluster = udon_clusters[udon_clusters['cluster'] == cluster]
                pseudobulks_in_other_clusters = udon_clusters[udon_clusters['cluster'] != cluster]

                # Create contingency table
                # Compute the counts for the Fisher's test
                n_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 1) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
                n_not_cm_ct_in_cluster = len(pseudobulks_in_cluster[(pseudobulks_in_cluster['clinical_measure'] == 0) & (pseudobulks_in_cluster['cell_types'] == cell_type)])
                n_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 1) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])
                n_not_cm_ct_in_other_clusters = len(pseudobulks_in_other_clusters[(pseudobulks_in_other_clusters['clinical_measure'] == 0) & (pseudobulks_in_other_clusters['cell_types'] == cell_type)])

                # Create a contingency table for Fisher's test
                contigency_mat = np.array([[n_cm_ct_in_cluster, n_not_cm_ct_in_cluster],
                                         [n_cm_ct_in_other_clusters, n_not_cm_ct_in_other_clusters]])

                contigencies[batch_component] = contigency_mat

            # Determine if CMH or Fisher's test should be used
            if np.shape(contigencies)[2] < 2 or np.any(contigencies[0, 0, :] < n_samples):
                p_value = np.nan
            else:
                udon_clusters['cm_ct_binary'] = 'no'
                udon_clusters.loc[(udon_clusters['clinical_measure'] == 1 & udon_clusters['cell_types'] == cell_type), 'cm_ct_binary'] = 'yes'
                udon_clusters['cluster_binary'] = 'no'
                udon_clusters.loc[udon_clusters['cluster'] == cluster, 'cluster_binary'] = 'yes'

                _, _, p_value = CMH(udon_clusters, 'cm_ct_binary', 'cluster_binary', stratifier='batches', raw=True)

            # Store results in matrices
            p_val_matrix.at[cell_type, cluster] = p_value

    return p_val_matrix


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


