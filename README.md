# Network-based Shallow Shotgun Metagenomics Functional Profiling Mining pipeline

This repository contains scripts and necessary data for a **network-based shallow shotgun metagenomics (SSMS) functional profiling pipeline**. The pipeline aims to provide a workflow to impute the missing values of shallow shotgun metagenomic data using information from metagenomic reference profiles. 

## Repository overview

The repository is organized as follows:

1. **Data**  
   You can find the following documents in this file:
   (1) the network files, which are in the adjacency matrix format
   (2) the original unstratified relative abundance of level4 EC Number and the pathway profiles of both SSMS and MGS data
   (3) the imputed level4 EC Number of SSMS
  （4）the imputed pathway profiles of SSMS

   **Note** that the raw fastq data of the metagenomics data is stored under SRA (Sequence Read Archive). Please find their project ID in our paper. In this repo, we only show the Linux script for preprocessing of the raw fastq data under /scripts/preprocessing.

2. **scripts**  
   You can find the following scripts in this file:
   (1) scripts for curating SSMS from MGS, and the preprocessing of the raw SSMS and MGS fastq data (Linux)
   (2) scripts for comparison, evaluation and visualization of the results generated using the imputation workflow (R, Jupyter Notebook)
   (3) scripts for the imputation workflow of level4 EC Number profiles of SSMS (Python)
   (4) adopted scripts from HUMAnN (a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community), which is used to generate pathway profiles from the imputed level4 EC Number profiles (Python)
  （5）scripts for network construction （R）

## Recommended order of scripts use/reading

The scripts in this repository are organized according to the main analysis workflow. We recommend using them in the following order:

1. **Preprocess raw FASTQ files**  
   First, we randomly subsample the reads from raw MGS sequencing data to curate the raw SSMS sequencing data. **Note that if you have paired sequencing data (perform both MGS and SSMS on the same sample), then you can skip this step.** Then, we use the Linux shell scripts to preprocess the raw FASTQ sequencing data. This step prepares the sequencing reads for downstream analysis. After preprocessing, we will get both level4 EC Number and pathway profiles of samples. 

2. **Construct functional association networks**  
   Next, use the R scripts for network construction. These scripts build level 4 EC Number association networks that are later used as adjacency matrices in the imputation workflow.

3. **Run the level 4 EC Number imputation workflow**  
   Use the Python scripts to perform network-guided imputation of level 4 EC Number profiles from SSMS data. This step takes the SSMS EC profile, the reference MGS EC profile, and the EC adjacency matrix as input, and returns an imputed level 4 EC profile.

4. **Generate pathway profiles from imputed EC profiles**  
   After obtaining the imputed level 4 EC Number profiles, use the adopted HUMAnN-related Python scripts to infer pathway-level profiles. These scripts are adapted from HUMAnN, a pipeline for profiling the presence/absence and abundance of microbial pathways in a community.

5. **Compare, evaluate, and visualize the results**  
   Finally, use the R scripts and Jupyter notebooks for downstream comparison, evaluation, and visualization. This step is used to assess the performance of the imputation workflow and to compare different methods, networks, or parameter settings.

## evaluation of the imputation pipeline
randomized_network.ipynb

This notebook is used to analyze the effect of network structure. It helps test whether the observed performance depends on meaningful functional associations or could also arise from randomized network relationships.

Typical questions addressed by this notebook include:

* Does the real functional network provide useful information?
* How does prediction performance change when the network is randomized?
* Is the network-guided result better than a random-network baseline?

## Conceptual inputation workflow
Input matrices
The workflow uses three feature-by-sample matrices:

adj_matrix: EC-EC adjacency matrix describing which ECs are associated with each other.
ssms: shallow shotgun metagenomic sequencing profile, where rows are EC features and columns are samples.
mgs: reference metagenomic sequencing profile, also with EC features as rows and samples as columns.
Feature alignment
The EC features are first aligned across the adjacency matrix, SSMS matrix, and MGS matrix. Only ECs shared by all required inputs are retained, ensuring that the network, SSMS data, and MGS reference data use the same feature space.

Neighbor selection for each target EC
For each EC feature, the function identifies its connected neighbors from the adjacency matrix. These neighbors are treated as informative correlated ECs that can help predict the abundance of the target EC.

Model training using MGS reference data
For each target EC, a regression model is trained on the MGS dataset:

Predictors: abundances of the neighboring ECs in MGS samples.
Response: abundance of the target EC in MGS samples.
The model can be either:

Random Forest regression (RF)
Linear Regression (LR)
Preserve observed SSMS values
For each SSMS sample, if the target EC already has a positive observed abundance, the original SSMS value is kept unchanged. The model only attempts to impute values that are zero or absent.

Filtering before imputation
Before predicting a missing value, two filters are applied:

Target prevalence filter: the target EC must be present in at least 10% of MGS reference samples.
Neighbor support filter: enough neighboring ECs must be observed in the SSMS sample, requiring at least 20% of neighbors and at least one neighbor to be present.
If either filter fails, the target EC is assigned zero for that sample.

Imputation of missing EC values
If the target EC passes the filters, the observed neighboring EC abundances in the SSMS sample are used as model input. The trained MGS-based regression model then predicts the missing abundance of the target EC.



## Example use case
input:
* /path/your_ssms.tsv
* /path/your_mgs.tsv
* /path/your_adj_matrix.txt

output:
* /path/your_imputed_ssms.tsv

python pipeline.py --ssms /path/your_ssms.tsv --mgs /path/your_mgs.tsvc --adj /path/your_adj_matrix.txt --output /path/your_imputed_ssms.tsv

## Prerequisites
R Prerequisites
R >= 4.0
Required R packages:
dplyr
ggplot2
stringr
Matrix
igraph
pheatmap
scales
proxy
ggVennDiagram
reshape2
Rfast
infotheo
mclust
boot
ComplexHeatmap
NetCoMi
SpiecEasi
GENIE3

Python >= 3.7
Required Python packages:
numpy
pandas
scikit-learn
scipy
matplotlib
joblib
tqdm
networkx
powerlaw
torch
pytorch-lightning
jupyter

## Cite 
If you use or mention this method in your research, please cite:

## Contact
If you have any questions or issues, please feel free to contact us via email zihua.huang@tum.de
