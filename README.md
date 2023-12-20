# DMMO

-------------------------------------
Introduction
--------------------------------------
We developed DMMO, a method to identify cancer driver modules based on multi-omics data. This method integrates gene expression, gene mutation, gene transcription regulation and protein interaction network data to accurately identify cancer driver modules.

--------------------------------------
How to use
--------------------------------------
## Data processing
1.Data_PCA.py: The dimension of gene expression data is reduced by principal component analysis

input:RNA_exp.csv The input is the gene expression matrix 

output: exp_pca.csv

2.Data_GAE.py: GAE Feature learning for PPI networks

input: The input is the gene expression matrix and the corresponding PPI edge relationship (maximal connected graph).

output: data/ppi_gae.csv

3.Date_CoReg.py: The PPI network calculates the weights of genes and screens the regulatory data to calculate the value of co-regulation

input：PPI network and gene transcriptional regulatory relationship network

output：node_coReg.txt

## The four data types are merged and the score is calculated
4.DMMO.R:

input：

data/exp_pca.csv  Gene expression matrix after dimensionality reduction

data/mut.txt   A two-dimensional matrix of somatic mutations

data/node_coReg.txt   Co-regulation data for gene pairs

data/ppi_gae.csv  Data after dimensionality reduction of ppi network

output：

p_sort%K_index.txt && p_sort%K_name.txt Preliminary screening of modules according to K

p_union_name(%K).txt The module after going to the superset

p_mergeModul_%K.txt  Merged modules according to different K values

com.txt && score.txt  Modules of different sizes are merged, and the collection file of the module and the score value corresponding to the module are output

## Reordering of modules is performed based on module size and score
5.Module_analysis.py  Reorder the set of modules based on their size and score

input: score_p10.txt com_p10.txt

output: p10_4_DMMO.txt
