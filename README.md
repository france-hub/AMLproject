# AMLproject
Type of data: flow cytometry, scRNAseq, spectral flow cytometry

scRNAseq
# scRNAseq_step1
This script is on preprocessing of scRNAseq data. Later on I'll use Seurat to analyze my data. Here though I use singlecellexperiment because it allows me to use the doubletfinder package (https://github.com/chris-mcginnis-ucsf/DoubletFinder). 
In summary the script does the followings:
1) rename the genes
2) Add metadata
3) Remove undetected genes
3) Infer and remove doublets
4) Calculate QC metrics with diagnostic plots
5) Remove outliers
# scRNAseq_step2
This script does the followings:
1) Normalize and find the most variable features 
2) integrate the different datasets (I have split them according to batch). I use Seurat to do this; there are several integration workflows
#here: https://satijalab.org/seurat/articles/integration_mapping.html. The developers don't really suggest one 
#over the other. I was interested in their SCTransform algorithm, however it's computationally heavy (my computer
#always crush). I have run it in the past through a server and did not get different results
3) Runs PCA, TSNE, UMAP on the integrated dataset
4) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
5) Adds the data on TCR (it's a VDJ scRNAseq) to the Seurat object, for this I used scRepertoire: http://127.0.0.1:22003/library/scRepertoire/doc/vignette.html 
# scRNAseq_step3
This script does the followings:
1) Pick the resolution for the clustering
2) Using scGate (https://github.com/carmonalab/scGate) creates a gating model to subset the different cell subsets. This can be done also in Seurat but I prefer scGate because you can also exclude negative margers to make the model more precise and subset the object on the output 
3) According to the expression and the DGE (FindAllMarkers) labels the clusters and then creates a seurat object for each subset
# step4_CD8
This script does the followings:
1) Loads in the CD8 subset, find the most variable features, scale the data, run dimensionality reduction, run clustering
2) Plots single genes and scores to identify different subsets
3) Labels clusters and look at clusters abundancies
4) Performs trajectory inference (slingshot)

In progress
1) DGE along and between lineages and across conditions (Res, NonRes, HD)
2) clonotype analysis

# Velocity
This script use scVelo to infer velocity
