# AMLproject
Type of data:  flow cytometry (Exploratory), scRNAseq, spectral flow cytometry

# Exploratory


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


## CD8

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
