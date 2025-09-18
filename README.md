# udon
UDON is an unsupervised framework for defining disease subtypes from patient scRNA-seq data. It uses control-normalized pseudobulk expression and sparse NMF to identify stable cell-state perturbations. UDON includes SATAY-UDON for metadata enrichment and SASHIMI-UDON for cohorts with limited controls.

# Steps for using this code
SASHIMI-UDON files coming soon! 

Note before starting: since it is in python, you will either need to convert your [Seurat object to h5ad](https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html) (or use the `sceasy` package) or create an h5ad/anndata object using the scanpy pipeline. Currently UDON supports only human data. Please reach out in the Github Issues if you have any other species data and I can help with setting up database for those species. 
 

1. Create a virtual environment using conda or pip.
2. Install the `py-udon` package. Package can be found on PyPi. 
   `pip install py-udon`
3. Ensure you have the following packages (preferably with the versions mentioned) available in a virtual environment:  
 `scanpy=1.9.5 ; numpy=1.24.4 ; scikit-learn=1.3.0; pandas=2.1.0; nimfa=1.4.0 ; statsmodels=0.14.0; tqdm=4.66.1; time; seaborn=0.12.2; matplotlib.pyplot=3.7.2`
4. Download the database files attached as zip file here. Please unzip the database folder and note the path to the folder and its contents. 
4 .You can edit the run_example_dataset.py to include your data and run udon on it. 
Note: if loading the anndata converted from Seurat object, it is important that adata.X have normalized data (achieved with `NormalizeData` function in Seurat or `sc.pp.normalize_total(adata, target_sum=1e4) ; sc.pp.log1p(adata)`) . The dataset I have shown in this script is for the SJIA cohort data in the paper whose h5 files are at: GSE207633
5. Please review the metadata file [donor_metadata.xlsx](https://github.com/kairaveet/py-udon/blob/main/donor_metadata.xlsx) required for the above-mentioned script to ensure your metadata follows the same format.
