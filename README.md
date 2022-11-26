# AMLproject
Type of data: flow cytometry, scRNAseq, spectral flow cytometry


# scRNAseq
# Step1 (preprocessing)
This script is on preprocessing of scRNAseq data.

In summary the script does the followings:
1) rename the genes
2) Add metadata
3) Remove undetected genes
3) Infer and remove doublets
4) Calculate QC metrics with diagnostic plots
5) Remove outliers
# Step2 (integration + add TCR)
This script does the followings:
1) Normalize and find the most variable features 
2) integrate the different datasets (I have split them according to batch).
3) Runs PCA, TSNE, UMAP on the integrated dataset
4) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
5) Adds the TCR data to the Seurat object 
# Step3 (subsets)
This script does the followings:
1) Pick the resolution for the clustering
2) Using scGate creates a gating model to subset the different cell subsets. 
# CD8
This script does the followings:
1) Loads in the CD8 subset, find the most variable features, scale the data, run dimensionality reduction, run clustering
2) Plots single genes and scores to identify different subsets
3) Labels clusters and look at clusters abundancies
4) Performs trajectory inference (slingshot)

In progress
1) DGE along and between lineages and across conditions (Res, NonRes, HD)
2) clonotype analysis

# Velocity (CD8)
This script use scVelo to infer velocity
