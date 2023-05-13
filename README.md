# AMLproject
Type of data:  flow cytometry and bulkRNA-seq (Exploratory), scRNAseq, spectral flow cytometry

# Exploratory
## EploratoryFlow
CD8 This script can be divided into 3 parts

A)Define/create the directories where our files are located, merge 2 experiments, perform Arcsinh transformation and create flowSet
1) Define the directory where this script and the fcs files are located
NB: we are using pregated CD8+ fcs files obtained with flowJo software
2) Define workingDirectory and create the directory
3) Create flowsets for the two experiments and merge the two flowsets
4) Perform Arcsinh transformation (#Arcsinh transformation (this is necessary when you have flow citometry data in fcs format instead of csv) It is well explained here https://wiki.centenary.org.au/display/SPECTRE/Data+transformation, choose cofactors for each channel)
5) Define the directory where the transformed files are located and save fcs
6) Read the transformed files as flowset

B) Prepare the data to create a single cell experiment (sce) using CATALYST
1) Create panel dataframe
2) Create metadata dataframe
3) Create sce object

C) Perform QC, clustering and dimensionality reduction
1) Visualize CATALYST QC plots 
2) Run FlowSOM and ConsensusClusterPlus
3) Run dimensionality reduction (UMAP) and visualize
4) Add annotations and visualize

## DotPlot_GSEA

This script plots gene signature enrichments from Zheng et al.

# scRNAseq
## AllCells

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


## CD8.R

This script can be organized into 8 subsections (A-H)

A)Subclustering and use of custom and published markers along with bulk-RNA seq results and DGE to annotate clusters
1) Subclustering as suggested here https://github.com/satijalab/Seurat/issues/2087 by timoast
2) Use of custom single markers and visualize as DotPlot and FeaturePlot
3) Use published signatures and signature of senescence obtained by bulk-seq and plot as feature and violin plots
4) Correlation plot of dysfunction vs senescence signature
5) DGE (top10 markers)
5) Gene coexpression matrix and find genes with correlation > q99
6) Annotate clusters

B) Analyze annotated subsets
1) Plot and look at clusters proportions and densities across samples, between conditions (group_id)
2) Combine top10 and custom markers and plot them as heatmap and as scores with FeaturePlot
3) DA analysis using permutation test, speckle (Anova) and CNA

C) Compute trajectories to infer CD8+ cells differentiation path
1) Run Slingshot and plot the two trajectories
2) Use tradeSeq to fit a generalized additive model (GAM) and then use patternTest to study gene expression

D) Compute gene set variation analysis and infer pathway activity
1) Use GSVA and limma to study gene set enrichments between couples of subsets (SenLvsStL, SenLvsActEx, ActExvsStL)
2) Use SCpubr and decoupleR for pathway activity inference

E) TCR repertoire analysis
1) scRepertoire analysis
2) Putative tumor reactive

F) Reclustering with increased resolution and further analyses
1) Reclustering and clusters frequencies across conditions
2) STARTRAC method
3) DGE Int vs. SenL
4) Gating model for SLEC markers 
5) Trm markers
6) Gating model for Int markers (from spectral flow data)

G) Velocity inference
1) Load in loom files and prepare to analyze in Python via scVelo
2) import from python the csv with the info on latent time and different subsets and plot

H) SCENIC
1) Download the required databases
2) Run SCENIC workflow
3) Use exportsForArboreto to export matrix and scenicOptions and analyze in Python using GRNBoost 
4) Read in GRNBoost output from Python and proceed with the workflow

## GRNboost.py
This script uses the GRNBoost2 function to infer gene regulatory networks

## velocity.py
This script

1) loads the h5ad file obtained from the Seurat object in python as AnnData 
2) performs RNA velocity using the scvelo: it uses scv.tl.recover_dynamics to recover the full dynamics and scv.tl.velocity to infer the velocity for each gene
3) computes the velocity graph and the cell transitions predicted via the function scvelo.tl.velocity_graph 
4) computes the latent time using scvelo.tl.latent_time and setting the Naive subset as the key of root cells to be used

# SpectralFlow

## CD8_spec
This script can be organized into 4 parts

A)Define/create the directories where our files are located, transform the csv in fcs and create flowSet
1) Define the directory where this script is located
2) Define the directory where the csv files are located 
NB: we are using pregated CD8+ csv files obtained with the flowJo csv channel values option and then transform them here 
in fcs files so that we won't need to arcsinh transform the data as explained here https://wiki.centenary.org.au/display/SPECTRE/Data+transformation)
3) Define and create the WorkingDirectory (here we'll save the plots) and the directory that will contain the fcs files
4) Transform the csv in fcs, save and create the flowSet

B) Prepare the data to create a single cell experiment (sce) using CATALYST
1) Create panel dataframe
2) Create metadata dataframe
3) Visualize density plots and use warpSet to normalize 
4) Create sce object

C) Perform QC, clustering and dimensionality reduction
1) Visualize CATALYST QC plots and remove Ki67 channel because of the very low NRS
2) Run FlowSOM and ConsensusClusterPlus
3) Run dimensionality reduction - PCA, UMAP and visualize
4) Add annotations and visualize

D) Trajectory inference
1) Use PCA as dimensionality reduction algorithm and run Slingshot
2) Visualize trajectories on first 3 PCs
3) Visualize subsets distribution along PT (jitter plot)
4) Use UMAP as dimensionality reduction algorithm and run Slingshot
5) Visualize trajectories onto 2D UMAP

