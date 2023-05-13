## AMLproject
Type of data: flow cytometry, scRNAseq, spectral flow cytometry

#AllCells

This script includes 3 parts (Preprocessing, integration/batch correction and clonotype data addition, subclustering):

A) Preprocessing
1) Load in CellRanger output
2) rename rowData, colData and barcodes
3) define dimnames(sce)
4) Add metadata
5) Remove undetected genes
6) Infer and remove doublets
7) Calculate QC metrics and diagnostic plots
8) Find and remove outliers

B) Integration/batch correction and clonotype data addition
1) Create seurat object from sce
2) Split by batch them normalize, find the most variable features and scale
3) integrate the different datasets
4) Scale the integrated datased and run PCA, TSNE and UMAP
5) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
6) Adds the data on TCR (it's a VDJ scRNAseq) to the Seurat object 

C) Subclustering
1) Pick the resolution for clustering
2) Visualize not-annotated clusters' distribution onto 2D UMAP
3) Using scGate (https://github.com/carmonalab/scGate) to plot the subsets of interest
4) Label the clusters according to scGate results
5) Use gating models to subset and save into single rds objects
6) Use cluster annotation to subset CD4 cells (identified by scGate). This because CD4 gene is not very well detected in scRNA-seq hence the gating model 
can be not enough senesitive to "gate" on CD4 (this was discussed a little bit here: https://github.com/carmonalab/scGate/issues/15)


#CD8

This script does the followings:
1) Loads in the CD8 subset, find the most variable features, scale the data, run dimensionality reduction, run clustering
2) Plots single genes and scores to identify different subsets
3) Labels clusters and look at clusters abundancies
4) Performs trajectory inference (slingshot)
1) DGE along and between lineages and across conditions (Res, NonRes, HD)
2) clonotype analysis

# Velocity (CD8)
This script use scVelo to infer velocity
