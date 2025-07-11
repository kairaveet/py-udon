# udon
UDON is an unsupervised framework for defining disease subtypes from patient scRNA-seq data. It uses control-normalized pseudobulk expression and sparse NMF to identify stable cell-state perturbations. UDON includes SATAY-UDON for metadata enrichment and SASHIMI-UDON for cohorts with limited controls.

## Steps for using this code (Updated for Refined Science)
WARNING: This code is forked from a repo of code that is more organized that than from the original paper, but is still under development. The code was forked in June 2025 so be aware that any changes to the original should likely be integrated.

1. Create two virtual environment using the udon_env.yaml file and udon_satay_env.yaml: `conda env create -f udon_env.yaml` and `conda env create -f udon_satay_env.yaml`
    * I have two separate environments because the initial UDON README specified a specific combination of seaborn and numpy that didn't actually work for graphing a heatmap while running UDON Satay. To avoid any strange complications of environments, I made a new environment just for running UDON Satay.
2. Unzip the pathway_files.zip folder under udon-python-test. This should get you two files:
  1. PathwayCommons-Hs.txt
  2. ProteinCoding-Hs-Mm.txt
  
The source code is found in udon-python-test. 

## Running UDON
The most recent example for running this code can be found at [here](https://github.com/RefinedScience/AML_atlas_exploration/tree/BDS-30-Gene-signature-extraction/analyses/gene-sig-extraction/UDON/BoneMarrow/Running_UDON.py) as either a jupyter notebook or python script.
The README in that repository folder also indicates the input and output.

## Running UDON-Satay
The most recent example for running this udon-satay can be found at [here](https://github.com/RefinedScience/AML_atlas_exploration/tree/BDS-50-Survival-Response-analysis/analyses/gene-sig-extraction/UDON/BoneMarrow/Run_UDON_Satay.py) as either a jupyter notebook or python script.
The README in that repository folder also indicates the input and output.

## Trouble Shooting
### Environment
If the environment doesn't work for some reason, feel free to start from scratch. The original instructions for UDON are:
 Ensure you have the following packages (preferably with the versions mentioned) available in a virtual environment: scanpy=1.9.5 ; numpy=1.24.4 ; scikit-learn=1.3.0; pandas=2.1.0; nimfa=1.4.0 ; statsmodels=0.14.0; tqdm=4.66.1; time; seaborn=0.12.2; matplotlib.pyplot=3.7.2
 
