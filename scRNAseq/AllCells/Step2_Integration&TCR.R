#This script does the followings:
#1) Normalize and find the most variable features 
#2) integrate the different datasets (I have split them according to batch). I use Seurat to do this; there are several integration workflows
#here: https://satijalab.org/seurat/articles/integration_mapping.html. The developers don't really suggest one 
#over the other. I was interested in their SCTransform algorithm, however it's computationally heavy (my computer
#always crush). I have run it in the past through a server and did not get different results
#3) Runs PCA, TSNE, UMAP on the integrated dataset
#4) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
#5) Adds the data on TCR (it's a VDJ scRNAseq) to the Seurat object, for this I used scRepertoire: http://127.0.0.1:22003/library/scRepertoire/doc/vignette.html 


rm(list = ls())
library(cowplot)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scRepertoire)

setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")
sce <- readRDS("scRNAseq_step1_test.rds")

sobj <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "AML")

# split by batch
cells_by_batch <- split(colnames(sce), sce$batch)
so.list <- lapply(cells_by_batch, function(i) subset(sobj, cells = i))

# normalize, find variable genes, and scale
so.list <- lapply(so.list, NormalizeData, verbose = FALSE)
so.list <- lapply(so.list, FindVariableFeatures, nfeatures = 2e3, 
                  selection.method = "vst", verbose = FALSE)
so.list <- lapply(so.list, ScaleData, verbose = FALSE)

#adjust the limit for allowable object sizes within R
options(future.globals.maxSize = 4000 * 1024^2)

# find anchors & integrate
as <- FindIntegrationAnchors(so.list, verbose = FALSE) #ignore warnings the function does not use random seeds
sobj <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

# scale integrated data
DefaultAssay(sobj) <- "integrated"
sobj <- ScaleData(sobj, verbose = FALSE)

#Run dimensionality reduction
sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
sobj <- RunTSNE(sobj, reduction = "pca", dims = seq_len(20),
                seed.use = 1, do.fast = TRUE, verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = seq_len(20),
                seed.use = 1, verbose = FALSE)

DimPlot(sobj, reduction = "pca", group.by = "batch")
DimPlot(sobj, reduction = "umap", group.by = "batch", split.by = "group_id")

# Plot the elbow plot
ElbowPlot(object = sobj, ndims = 30)

#Clustering
sobj <- FindNeighbors(sobj, reduction = "pca", dims = seq_len(20), verbose = FALSE)
sobj <- FindClusters(object = sobj, random.seed = 1, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 2))

#Add clonotype data
#Load data
contig_dir <- paste0(getwd(), "/../data_VDJ/")
files <- as.vector(list.files(contig_dir, pattern = "*.csv"))

# List contig files
contig_list <- list()
for (i in seq_along(files)){
  contig_list[[i]] <- read.csv(paste0(contig_dir, files[[i]]))
}

View(contig_list[[1]])

# Combine TCR data and add to Seurat object
combined <- combineTCR(contig_list, 
                       samples = c("p219", "p229", "p264", "pV03", "pV09", "p219", "p229", "p264", "pV03", "pV09"), 
                       ID = c(rep("base", 5), rep("post", 5)), 
                       cells = "T-AB")

combined <- addVariable(combined, name = "batch", 
                        variables = c("b2", "b1", "b2", "b2", "b3", "b2", "b1", "b2", "b2", "b3"))

sobj <- combineExpression(combined, sobj, 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

#Organizing the order of the factor cloneType
update_geom_defaults("point", list(stroke=0.5))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(sobj, reduction = "umap", group.by = "cloneType") + scale_color_manual(values = c(colorblind_vector(5)), na.value="grey")

x <- table(sobj[[]]$seurat_clusters, sobj[[]]$cloneType, useNA = "ifany")
for (i in 1:nrow(x)) {
  x[i,] <- x[i,]/sum(x[i,])
}
x <- data.frame(x)
ggplot(x, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_manual(values = colorblind_vector(5), na.value="grey") +
  theme_classic()

saveRDS(sobj, file = "scRNAseq_step2_test.rds")
