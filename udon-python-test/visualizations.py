import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import to_hex
import umap

# Define a colormap to use for marker expression heatmap
N = 256
vals = np.ones((N*2, 4))
vals[:N, 0] = np.linspace(15/256, 0, N)
vals[:N, 1] = np.linspace(255/256, 0, N)
vals[:N, 2] = np.linspace(255/256, 0, N)
vals[N:, 0] = np.linspace(0, 255/256, N)
vals[N:, 1] = np.linspace(0, 243/256, N)
vals[N:, 2] = np.linspace(0, 15/256, N)
blue_black_yellow_cmap = ListedColormap(vals)


def assign_rainbow_colors_to_groups(groups):
    """ Creates a dictionary of cluster names to hexadecimal color strings

    This function takes a list of groups and assigns each unique item in the
    groups a color (using the matplotlib.cm.rainbow color-map) as a hexadecimal
    string value. This is useful for storing a single color scheme for clusters
    to be used with downstream visualizations.

    Arguments
    ---------
    groups : numpy.Array[str]
        List of cluster names. (Required)

    Returns
    -------
    dict {str:str}
        Dictionary of cluster-names to assigned colors (hexadecimal value)

    """
    unique_groups = np.unique(groups)
    groups_to_num = pd.Series(list(range(len(unique_groups))),
                              index=unique_groups)
    n = len(unique_groups)
    groups_to_color = pd.Series([to_hex(cm.rainbow(item / n)) for item in
                                 groups_to_num.values], index=groups_to_num.index.values).to_dict()
    return groups_to_color


def plot_markers_df(marker_heatmap, markers_df, clusters, path_to_save_figure):
    """ Plots a heatmap of the MarkerFinder results
    Arguments
    ---------
    marker_heatmap : pandas.DataFrame
        Data to plot. The rows must intersect with features in the \
        ordered_markers_df and the columns must intersect with the cells in \
        clusters.
    markers_df : pandas.DataFrame
        The markers_df where "marker" column are the genes and the "top_cluster" is the cluster to which the gene belongs to  (Required)
    clusters : pandas.DataFrame
        The cluster assignments of the samples, index of this data frame are the sample names (Required)
    path_to_save_figure : str
        The path to save the figure. (Required)

    Returns
    -------
    None
        Saves the heatmap to the file specified by path_to_save_figure

    """
    try:
        # get Dictionary of cluster-names to assigned colors (hexadecimal value)
        groups_to_colors = assign_rainbow_colors_to_groups(groups=np.array(clusters["cluster"]))

        # Filter the input data matrix to the markers and cells of interest
        input_df = marker_heatmap.loc[markers_df["marker"].values, clusters.index]
        input_df = input_df.astype(float)

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in clusters["cluster"].values]},
                                      index=clusters.index).T

        # Build row annotation dataframe
        row_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in markers_df["top_cluster"].values]},
                                     index=markers_df["marker"].values)

        # Get top markers for each cluster
        top_markers = markers_df.groupby('top_cluster')['marker'].first()

        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": clusters["cluster"].values,
            "position": list(range(clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build df to add marker label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": markers_df["top_cluster"].values,
            "position": list(range(markers_df.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Z-score normalize the expression matrix
        z_df = input_df.apply(lambda x: pd.Series(zscore(x.values),
                                                  index=x.index.values), axis=1)

        # Build the heatmap
        plt.close("all")
        fig = plt.figure(constrained_layout=True, figsize=(9, 5))
        ax = fig.add_gridspec(2, 2, width_ratios=(1, 20), height_ratios=(1, 20),
                      wspace=0.005, hspace=0.01)  # Adjusted to create space for row annotation and top markers
        ax2 = fig.add_subplot(ax[1, 1])
        ax1 = fig.add_subplot(ax[0, 1])
        ax3 = fig.add_subplot(ax[1, 0])  # Subplot for row annotation

        heat1 = sns.heatmap(cell_labels_df,
                            yticklabels=False,
                            xticklabels=False,
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(z_df,
                            vmin=-3,
                            vmax=3,
                            cmap=blue_black_yellow_cmap,  # Assuming you have blue_black_yellow_cmap defined
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here
                            cbar=True,
                            cbar_kws={"shrink": 0.5},
                            ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")

        heat3 = sns.heatmap(row_labels_df,
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here as well
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax3)
        ax3.set_yticks(label_to_position2.values)
        ax3.set_yticklabels(label_to_position.index.values)
        # ax3.yaxis.tick_left()
        fig.suptitle('Marker Finder Heatmap', fontsize=16)
        # fig.tight_layout()
        plt.savefig(path_to_save_figure, dpi=10, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_markers_df. See above Exception.")

# Usage example (this part should be replaced with actual data in your script)
# marker_heatmap = pd.DataFrame(...)
# markers_df = pd.DataFrame(...)
# clusters = pd.DataFrame(...)
# groups_to_colors = {'cluster1': '#1f77b4', 'cluster2': '#ff7f0e', ...}
# path_to_save_figure = 'heatmap.png'
# plot_markers_df(marker_heatmap, markers_df, clusters, groups_to_colors, path_to_save_figure)


def plot_metadata_heatmap(udon_clusters, udon_metadata, metadata_col, path_to_save_figure):
    try:
        # get Dictionary of cluster-names to assigned colors (hexadecimal value)
        groups_to_colors = assign_rainbow_colors_to_groups(groups=np.array(udon_clusters["cluster"]))

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in udon_clusters["cluster"].values]},
                                      index=udon_clusters.index).T


        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": udon_clusters["cluster"].values,
            "position": list(range(udon_clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Get unique clusters and samples
        unique_clusters = udon_metadata[metadata_col].unique()
        samples = udon_metadata.index

        # Initialize the binary matrix with zeros
        binary_matrix = pd.DataFrame(0, index=unique_clusters, columns=samples)

        # Populate the binary matrix
        for sample in samples:
            cluster = udon_metadata.loc[sample, metadata_col]
            binary_matrix.loc[cluster, sample] = 1

        binary_matrix = binary_matrix.sort_index()
        print("sorted index")

        # Build df to add cluster label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": binary_matrix.index,
            "position": list(range(binary_matrix.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build the heatmap
        plt.close("all")
        fig = plt.figure(constrained_layout=True, figsize=(9, 6))
        ax = fig.add_gridspec(2, 1, height_ratios=(1, 20),
                                  wspace=0.005,
                                  hspace=0.01)  # Adjusted to create space for row annotation and top markers
        ax2 = fig.add_subplot(ax[1, 0])
        ax1 = fig.add_subplot(ax[0, 0])

        heat1 = sns.heatmap(cell_labels_df,
                                yticklabels=False,
                                xticklabels=False,
                                cmap=sns.color_palette(groups_to_colors.values()),
                                cbar=False,
                                ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(binary_matrix,
                                vmin=0,
                                vmax=1,
                                cmap='YlGn',
                                xticklabels=False,
                                yticklabels=False,  # Disable yticklabels here
                                cbar=True,
                                cbar_kws={"shrink": 0.5},
                                ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")
        ax2.set_yticks(label_to_position2.values)
        ax2.set_yticklabels(label_to_position2.index.values, ha="right")
        ax2.yaxis.tick_left()

        fig.suptitle('Metadata Heatmap', fontsize=16)
        # fig.tight_layout()
        plt.savefig(path_to_save_figure, dpi=10, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_metadata_heatmap. See above Exception.")


def plot_markers_df_subplots(marker_heatmap, markers_df, clusters, groups_to_colors, path_to_save_figure):
    """
    Plots a heatmap of the MarkerFinder results
    Arguments
    ---------
    marker_heatmap : pandas.DataFrame
        Data to plot. The rows must intersect with features in the \
        ordered_markers_df and the columns must intersect with the cells in clusters.
    markers_df : pandas.DataFrame
        The markers_df where "marker" column are the genes and the "top_cluster" is the cluster to which the gene belongs to  (Required)
    clusters : pandas.DataFrame
        The cluster assignments of the samples, index of this data frame are the sample names (Required)
    groups_to_colors : dict {str:str}
        Dictionary of cluster-names to assigned colors (hexadecimal value) \
        The pyInfinityFlow.Plotting_Utilities.assign_rainbow_colors_to_groups \
        can be used to generate this dictionary from a list of clusters. \
        (Required)
    path_to_save_figure : str
        The path to save the figure. (Required)

    Returns
    -------
    None
        Saves the heatmap to the file specified by path_to_save_figure
    """
    try:
        # Filter the input data matrix to the markers and cells of interest
        input_df = marker_heatmap.loc[markers_df["marker"].values, clusters.index]
        input_df = input_df.astype(float)

        # Map the order of the groups_to_colors to make a cmap
        groups_to_order = pd.Series(list(range(len(groups_to_colors))),
                                    index=groups_to_colors.keys())

        # Build top heatmap to label clusters
        cell_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in clusters["cluster"].values]},
                                      index=clusters.index).T

        # Build row annotation dataframe
        row_labels_df = pd.DataFrame({"cluster": [groups_to_order[item] for item in markers_df["top_cluster"].values]},
                                     index=markers_df["marker"].values)

        # Get top markers for each cluster
        top_markers = markers_df.groupby('top_cluster')['marker'].first()

        # Build df to add cluster label ticks
        label_to_position = pd.pivot_table(pd.DataFrame({
            "cluster": clusters["cluster"].values,
            "position": list(range(clusters.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Build df to add marker label ticks
        label_to_position2 = pd.pivot_table(pd.DataFrame({
            "cluster": markers_df["top_cluster"].values,
            "position": list(range(markers_df.shape[0]))}),
            index="cluster",
            values="position",
            aggfunc=np.mean)["position"].sort_values()

        # Z-score normalize the expression matrix
        z_df = input_df.apply(lambda x: pd.Series(zscore(x.values),
                                                  index=x.index.values), axis=1)

        # Build the heatmap
        plt.close("all")
        fig, axs = plt.subplots(2, 2, figsize=(18, 10),
                                gridspec_kw={'width_ratios': [1, 20], 'height_ratios': [1, 20], 'wspace': 0.005,
                                             'hspace': 0.01})
        ax1, ax2, ax3 = axs[0, 1], axs[1, 1], axs[1, 0]  # ax3 remains the same

        heat1 = sns.heatmap(cell_labels_df,
                            yticklabels=False,
                            xticklabels=False,
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax1)
        ax1.set_xticks(label_to_position.values)
        ax1.set_xticklabels(label_to_position.index.values, rotation=45, ha="left")
        ax1.xaxis.tick_top()

        heat2 = sns.heatmap(z_df,
                            vmin=-3,
                            vmax=3,
                            cmap=blue_black_yellow_cmap,  # Assuming you have blue_black_yellow_cmap defined
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here
                            cbar=True,
                            cbar_kws={"shrink": 0.5},
                            ax=ax2)
        ax2.collections[0].colorbar.set_label("Z-Score Normalized Expression")

        heat3 = sns.heatmap(row_labels_df,
                            xticklabels=False,
                            yticklabels=False,  # Disable yticklabels here as well
                            cmap=sns.color_palette(groups_to_colors.values()),
                            cbar=False,
                            ax=ax3)
        ax3.set_yticks(label_to_position2.values)
        ax3.set_yticklabels(label_to_position.index.values)

        fig.suptitle('Marker Finder Heatmap', fontsize=16)

        # Save the figure with compression
        plt.savefig(path_to_save_figure, dpi=100, bbox_inches='tight', format="pdf", pad_inches=0.9)

    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_markers_df. See above Exception.")

def plot_NMFbased_UMAP(basis_matrix, clusters_clean, path_to_save_figure, n_U_neighbors=15):
    """NMF/
    Takes an NMF basis matrix (components x samples) and plots a UMAP with dots
        colored by cluster, shape for cell type (each dot refers to a Diseased 
        SampleID-CT pair for UDON)
    Arguments
    ---------
    basis_matrix: Matrix where rows = NMF components and columns = samples from NMF
        (e.g. basis_matrix output from clustering_wrapper)
        NOTE: The sample names have to follow trend of CELLTYPE_SAMPLEID
    clusters_clean: Dataframe where rownames = Sample names as in columns of basis_matrix
        and column called cluster indicates the UDON-based cluster
    path_to_save_figure: Path where figure will be saved as a pdf
    n_U_neighbors: Value used as n_neighbors input for umap.UMAP
    
    Returns
    ---------
    None (plot)
    """
    try:
        # Transpose basis matrix (components x samples) to make rows = samples
        X = basis_matrix.T  # Now shape: (n_samples, n_components)   

        # Extract cell types from sample names
        umap_df = pd.DataFrame(X.index, columns=['sample'])
        split_sample = umap_df['sample'].str.split('__')
        umap_df['cell_type'] = split_sample.str[0] 
        umap_df['sampleid'] = split_sample.str[1] 

        # Add the clusters to the umap_df
        clusters_clean.index.name = 'sample'
        clusters_clean_reset = clusters_clean.reset_index()
        umap_df = umap_df.merge(clusters_clean_reset[['sample', 'cluster']], on='sample', how='left')

        # Run UMAP
        reducer = umap.UMAP(n_neighbors=n_U_neighbors, min_dist=0.1, metric='euclidean', random_state=42)
        embedding = reducer.fit_transform(X)  
        # Add UMAP coordinates
        umap_df['UMAP1'] = embedding[:, 0]
        umap_df['UMAP2'] = embedding[:, 1]


        # Plot
        plt.figure(figsize=(10, 6))
        ax2 = sns.scatterplot(data=umap_df, x='UMAP1', y='UMAP2', style='cell_type', 
                              hue="cluster", palette='tab10', s=40, alpha=0.8)
        for i, row in umap_df.iterrows():
            ax2.text(row['UMAP1'] + 0.2, row['UMAP2'], row['sampleid'], fontsize=7, alpha=0.7)
        plt.title("NMF-based UMAP Colored by Cluster")
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
        plt.tight_layout()

        # Save the figure with compression
        plt.savefig(path_to_save_figure, dpi=200, bbox_inches='tight', format="pdf", pad_inches=0.9)
        plt.close()
    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_NMFbased_UMAP. See above Exception.")


## DEPRECATED VERSION WHERE I PLOT TWO DIFFERENT PLOTS
def plot_NMFbased_UMAP_OG(basis_matrix, path_to_save_ctfigure, path_to_save_sifigure, n_U_neighbors=15):
    """NMF/
    Takes an NMF basis matrix (components x samples) and plots a UMAP with dots
        colored by cell type and SampleID (each dot refers to a Diseased SampleID-CT pair for UDON)
    Arguments
    ---------
    basis_matrix: Matrix where rows = NMF components and columns = samples from NMF
        (e.g. basis_matrix output from clustering_wrapper)
        NOTE: The sample names have to follow trend of CELLTYPE_SAMPLEID
    path_to_save_figure: Path where figure will be saved as a pdf
    n_U_neighbors: Value used as n_neighbors input for umap.UMAP
    
    Returns
    ---------
    None (plot)
    """
    try:
        # Transpose basis matrix (components x samples) to make rows = samples
        X = basis_matrix.T  # Now shape: (n_samples, n_components)   

        # Extract cell types from sample names
        umap_df = pd.DataFrame(X.index, columns=['sample'])
        split_sample = umap_df['sample'].str.split('_')
        umap_df['cell_type'] = split_sample.str[0] 
        umap_df['sampleid'] = "_".join(split_sample.str[1], split_sample.str[2])

        # Run UMAP
        reducer = umap.UMAP(n_neighbors=n_U_neighbors, min_dist=0.1, metric='euclidean', random_state=42)
        embedding = reducer.fit_transform(X)  
        # Add UMAP coordinates
        umap_df['UMAP1'] = embedding[:, 0]
        umap_df['UMAP2'] = embedding[:, 1]

        # Plot
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=umap_df, x='UMAP1', y='UMAP2', hue='cell_type', palette='tab10', s=40, alpha=0.8)
        plt.title("NMF-based UMAP Colored by Cell Type")
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        # Save the figure with compression
        plt.savefig(path_to_save_ctfigure, dpi=200, bbox_inches='tight', format="pdf", pad_inches=0.9)
        plt.close()

        # now the patient
        # Plot
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=umap_df, x='UMAP1', y='UMAP2', hue='sampleid', palette='tab10', s=40, alpha=0.8)
        plt.title("NMF-based UMAP Colored by Cell Type")
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        # Save the figure with compression
        plt.savefig(path_to_save_sifigure, dpi=200, bbox_inches='tight', format="pdf", pad_inches=0.9)
        plt.close()
    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_NMFbased_UMAP. See above Exception.")

def plot_top_marker_heatmap(markers_df, corr_df, path_to_save_figure, 
top_n=10, figsize=(10, 12), cmap='bwr'):
    """
    Plot heatmap of top markers by lowest p-value per cluster, using Pearson r from corr_df.

    Arguments
    ----------
    markers_df : pd.DataFrame
        Must contain: ['marker', 'top_cluster', 'pearson_r', 'p_value']
    corr_df : pd.DataFrame
        Genes (rows) × clusters (columns), with Pearson r values
    top_n : int
        Number of top markers per cluster
    figsize : tuple
        Heatmap figure size
    cmap : str
        Colormap for heatmap (e.g. 'Spectral', 'coolwarm')
    """
    
    required_cols = {'marker', 'top_cluster', 'pearson_r', 'p_value'}
    if not required_cols.issubset(markers_df.columns):
        raise ValueError(f"markers_df must have columns: {required_cols}")
    
    try:
        # Step 1: Get top N markers per cluster
        top_markers_df = (
            markers_df
            .sort_values(['top_cluster', 'p_value'])
            .groupby('top_cluster')
            .head(top_n)
        )

        unique_markers = top_markers_df['marker'].unique()

        # Step 2: Use the cluster-sorted marker list as row order
        ordered_markers = top_markers_df['marker'].values

        # Step 3: Filter correlation data
        heatmap_df = corr_df.loc[corr_df.index.intersection(ordered_markers)]

        # Ensure markers are ordered as in ordered_markers
        heatmap_df = heatmap_df.loc[ordered_markers]

        # Drop rows with all NaNs
        heatmap_df = heatmap_df.dropna(how='all')
        if heatmap_df.empty:
            raise ValueError("Heatmap matrix is empty after filtering.")

        # Plot using sns.heatmap
        plt.figure(figsize=figsize)
        ax = sns.heatmap(
            heatmap_df,
            cmap=cmap,
            vmin=-1,
            vmax=1,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": "Pearson r"},
            linewidths=0.3,
            linecolor="black"
        )

        # Clean up axis labels and title
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Top Markers")
        ax.set_title(f"Top {top_n} Markers per Cluster by Pearson R^2", pad=20)
        ax.grid(False)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)

        # Shrink colorbar
        colorbar = ax.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=8)  # Smaller font on colorbar ticks
        colorbar.ax.set_ylabel("Pearson R^2", fontsize=10)
        # Resize colorbar dimensions
        colorbar.ax.set_aspect(30) 

        # Remove border spines (optional)
        for _, spine in ax.spines.items():
            spine.set_visible(False)

        plt.tight_layout()

        # Save the figure with compression
        plt.savefig(path_to_save_figure, dpi=200, bbox_inches='tight', format="pdf", pad_inches=0.9)
    except Exception as e:
        print(str(e))
        print("Warning! Failed to run plot_top_marker_heatmap. See above Exception.")
