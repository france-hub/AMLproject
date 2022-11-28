###This script includes 3 parts (Preprocessing, integration/batch correction and clonotype data addition, subclustering):

##A) Preprocessing
#1) Load in CellRanger output
#2) rename rowData, colData and barcodes
#3) define dimnames(sce)
#4) Add metadata
#5) Remove undetected genes
#6) Infer and remove doublets
#7) Calculate QC metrics and diagnostic plots
#8) Find and remove outliers

##B) Integration/batch correction and clonotype data addition
#1) Create seurat object from sce
#2) Split by batch them normalize, find the most variable features and scale
#3) integrate the different datasets
#4) Scale the integrated datased and run PCA, TSNE and UMAP
#5) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
#6) Adds the data on TCR (it's a VDJ scRNAseq) to the Seurat object 

##C) Subclustering
#1) Pick the resolution for clustering
#2) Visualize not-annotated clusters' distribution onto 2D UMAP
#3) Using scGate (https://github.com/carmonalab/scGate) to plot the subsets of interest
#4) Label the clusters according to scGate results
#5) Use gating models to subset and save into single rds objects
#6) Use cluster annotation to subset CD4 cells (identified by scGate). This because CD4 gene is not very well detected in scRNA-seq hence the gating model 
#can be not enough senesitive to "gate" on CD4 (this was discussed a little bit here: https://github.com/carmonalab/scGate/issues/15)

rm(list = ls())

#Load packages
library(rstudioapi)
library(DropletUtils)
library(SingleCellExperiment)
library(readxl)
library(scds)
library(scater)
library(cowplot)
library(ggplot2)
library(LSD)
library(Matrix)
library(Seurat)
library(scRepertoire)
library(scCustomize)
library(scGate)
library(magrittr)

#Set working directory where the script is located
setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")

##################
#A) Preprocessing
##################

#1) Load in CellRanger output
dirs <- list.dirs("data_scRNAseq", recursive = FALSE, full.names = TRUE)
names(dirs) <- basename(dirs)
sce <- read10xCounts(dirs)
dim(sce) #check sce dimension

#2) rename rowData, colData and barcodes
rowData(sce) <- rowData(sce)[,1:2]
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))

#3) define dimnames(sce)
dimnames(sce) <- list(
  with(rowData(sce), SYMBOL), 
  with(colData(sce), paste(sample_id, barcode, sep = "_"))) #merge sample_ids with barcodes

#4) Add metadata
md <- file.path("../md_dir", "metadata.xlsx")
md <- read_excel(md)
m <- match(sce$sample_id, md$sample_id)
sce$group_id <- md$group_id[m]
sce$timepoint <- md$timepoint[m]
sce$RespTmp <- md$RespTmp[m]
sce$batch <- md$Batch[m]

#5) Remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#6) Infer and remove doublets (reference: https://github.com/chris-mcginnis-ucsf/DoubletFinder)
# split SCE by sample 
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

# run 'scds' for each sample
sce_by_s <- lapply(sce_by_s, function(u) 
  cxds_bcds_hybrid(bcds(cxds(u))))

# remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
  # compute expected nb. of doublets (10x)
  n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
  # remove 'n_dbl' cells w/ highest doublet score
  o <- order(u$hybrid_score, decreasing = TRUE)
  u[, -o[seq_len(n_dbl)]]
})

# merge back into single SCE
sce <- do.call("cbind", sce_by_s)

#7) Calculate QC metrics and diagnostic plots
(mito <- grep("MT-", rownames(sce), value = TRUE)) #find mitochondrial genes
sce <- addPerCellQC(sce, subsets = list(Mt = mito)) 
plotHighestExprs(sce, n = 20)

#8) Find and remove outliers
cols <- c("sum", "detected", "subsets_Mt_percent")
log <- c(TRUE, TRUE, FALSE)
type <- c("both", "both", "higher")

drop_cols <- paste0(cols, "_drop")
for (i in seq_along(cols))
  colData(sce)[[drop_cols[i]]] <- isOutlier(sce[[cols[i]]], 
                                            nmads = 2.5, type = type[i], log = log[i], batch = sce$sample_id)
sapply(drop_cols, function(i) 
  sapply(drop_cols, function(j)
    sum(sce[[i]] & sce[[j]])))

cd <- data.frame(colData(sce))
ps <- lapply(seq_along(cols), function (i) {
  p <- ggplot(cd, aes_string(x = cols[i], alpha = drop_cols[i])) +
    geom_histogram(bins = 100, show.legend = FALSE) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4)) +
    facet_wrap(~sample_id, ncol = 1, scales = "free") + 
    theme_classic() + theme(strip.background = element_blank())
  if (log[i]) 
    p <- p + scale_x_log10()
  return(p)
})

plot_grid(plotlist = ps, ncol = 3)

layout(matrix(1:2, nrow = 1))
out <- rowSums(as.matrix(colData(sce)[drop_cols])) != 0
x <- sce$sum
y <- sce$detected
heatscatter(x, y, log="xy", main = "unfiltered", 
            xlab = "Total counts", ylab = "Non-zero features")
heatscatter(x[!out], y[!out], log="xy", main = "filtered", 
            xlab = "Total counts", ylab = "Non-zero features")

# summary of cells kept
ns <- table(sce$sample_id)
ns_fil <- table(sce$sample_id[!out])
print(rbind(
  unfiltered = ns, filtered = ns_fil, 
  "%" = ns_fil / ns * 100), digits = 0)

# drop outlier cells
sce <- sce[, !out]
dim(sce)

# require count > 1 in at least 20 cells
sce <- sce[rowSums(counts(sce) > 1) >= 20, ]

dim(sce) #Dimension of sce after filtering

###############################################################
##B) Integration/batch correction and clonotype data addition
###############################################################

#1) Create seurat object from sce
sobj <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "AML")

#2) Split by batch then normalize, find the most variable features and scale
cells_by_batch <- split(colnames(sce), sce$batch) #split by batch
so.list <- lapply(cells_by_batch, function(i) subset(sobj, cells = i))
so.list <- lapply(so.list, NormalizeData, verbose = FALSE)
so.list <- lapply(so.list, FindVariableFeatures, nfeatures = 2e3, 
                  selection.method = "vst", verbose = FALSE)
so.list <- lapply(so.list, ScaleData, verbose = FALSE)

options(future.globals.maxSize = 4000 * 1024^2) #adjust the limit for allowable object sizes within R

#3) integrate the different datasets 
as <- FindIntegrationAnchors(so.list, verbose = FALSE) #ignore warnings the function does not use random seeds
sobj <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#4) Scale the integrated datased and run PCA, TSNE and UMAP
DefaultAssay(sobj) <- "integrated"
sobj <- ScaleData(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)
sobj <- RunTSNE(sobj, reduction = "pca", dims = seq_len(20),
                seed.use = 1, do.fast = TRUE, verbose = FALSE)
sobj <- RunUMAP(sobj, reduction = "pca", dims = seq_len(20),
                seed.use = 1, verbose = FALSE)

DimPlot(sobj, reduction = "pca", group.by = "batch")
DimPlot(sobj, reduction = "umap", group.by = "batch", split.by = "group_id")

#5) Computes the k.param nearest neighbors and then identifies clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm (Louvain). 
ElbowPlot(object = sobj, ndims = 30)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = seq_len(20), verbose = FALSE)
sobj <- FindClusters(object = sobj, random.seed = 1, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 2))

#6) Add the data on TCR (it's a VDJ scRNAseq) to the Seurat object 
contig_dir <- paste0(getwd(), "/../data_VDJ/") #data_VDJ is the folder where my VDJ csv files are located
files <- as.vector(list.files(contig_dir, pattern = "*.csv"))

contig_list <- list()
for (i in seq_along(files)){
  contig_list[[i]] <- read.csv(paste0(contig_dir, files[[i]]))
} #list contig files
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

#Organize the order of the factor cloneType
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

###################
##C) Subclustering
###################

#1) Pick the resolution for clustering
sobj <- SetIdent(sobj, value = "integrated_snn_res.0.8")
DefaultAssay(sobj) <- "RNA" #Put "RNA" as default assay

#2) Visualize not-annotated clusters' distribution onto 2D UMAP
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256) #define palette for expression to use later on
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq") #define palette for clusters to use later on
#Plot and save
tiff("../plots_CD8/moannot.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
noAnn <- DimPlot_scCustom(sobj, colors_use = pal_ident, label.size = 6) + theme_void() & theme(legend.position = "none")
noAnn
dev.off()

#3) Using scGate (https://github.com/carmonalab/scGate) to plot the subsets of interest
#CD3
CD3 <- gating_model(name = "CD3", signature = c("CD3D", "CD3G"))
CD3 <- scGate(sobj, model = CD3)
pcd3 <- DimPlot_scCustom(CD3, group.by = "is.pure", colors_use = c("red", "blue")) + 
  ggtitle(expression(paste("CD3"^{"+"}, " T cells")))

#CD8
CD8 <- gating_model(name = "CD8", signature = c("CD8A", "CD8B"))
CD8 <- scGate(sobj, model = CD8)
pcd8 <- DimPlot_scCustom(CD8, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle(expression(paste("CD8"^{"+"}, " T cells"))) + theme(legend.position="none")

#CD4
CD4 <- gating_model(name = "CD4", signature = "CD4")
CD4 <- scGate(sobj, model = CD4)
pcd4 <- DimPlot_scCustom(CD4, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle(expression(paste("CD4"^{"+"}, " T cells"))) + theme(legend.position="none")

#MAIT
MAIT <- c("TRAV1-2","SLC4A10")
MAIT <- gating_model(name = "MAIT", signature = MAIT)
MAIT <- scGate(sobj, model = MAIT)
pmait <- DimPlot_scCustom(MAIT, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle("MAIT") + theme(legend.position="none")

#NK
NK <- c("NCAM1","KLRD1","KLRG1")
CD3 <- c("CD3D", "CD3G")
NK <- gating_model(level=1, name = "NK", signature = NK)
NK <- gating_model(model=NK, level=1, name = "CD3", signature = CD3, negative=TRUE)
NK <- scGate(sobj, model = NK)
pnk <-DimPlot_scCustom(NK, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle("NK cells") + theme(legend.position="none")

#GD
GD <- c("CD3D", "TRGC1", "TRDV1", "TRDV2", "TRDV3")
CD4 <- c("CD4")
CD8 <- c("CD8A", "CD8B")
GD <- gating_model(level=1, name="GD", signature=GD)
GD <- gating_model(model=GD, level=1, name = "CD8", signature = CD8, negative=TRUE)
GD <- gating_model(model=GD, level=1, name = "CD4", signature = CD4, negative=TRUE)
GD <- scGate(sobj, model = GD)
pgd <- DimPlot_scCustom(GD, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle("GD cells") + theme(legend.position="none")

#NKT
#NKT
NKT <-c("CD3D", "CD3G", "NCAM1")
NKT <- gating_model(level=1, name = "NKT", signature = NKT)
NKT <- gating_model(model=NKT, level=1, name = "CD8", signature = CD8, negative=TRUE)
NKT <- gating_model(model=NKT, level=1, name = "CD4", signature = CD4, negative=TRUE)
NKT <- scGate(sobj, model = NKT)
pnkt <- DimPlot_scCustom(NKT, group.by = "is.pure", colors_use = c("red", "blue")) +
  ggtitle("NKT cells")

#Formatting for plotting
pcd3 <- pcd3 + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
pcd8 <- pcd8 + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
pcd4 <- pcd4 + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
pnk <- pnk + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
pmait <- pmait + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
pgd <- pgd + theme_void() + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))

# extract the legend from one of the plots
legend <- get_legend(
  pnkt + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p.row1 <- cowplot::plot_grid(NULL, noAnn, pcd3, pcd8, NULL, rel_widths = c(1/3,1.5,1,1,1/3), nrow =1)
p.row2 <- cowplot::plot_grid(NULL,  pcd4, pnk, pmait, NULL, rel_widths = c(1/3,1,1,1,1/3), nrow =1)
p.row3 <- cowplot::plot_grid(NULL, NULL,  pgd, legend, NULL, rel_widths = c(1, 1/3,1,1,1/3), nrow =1)

tiff("../plots_CD8/umap_allT.tiff", width = 5*700, height = 5*700, res = 300, pointsize = 5)     
cowplot::plot_grid(p.row1, p.row2, p.row3, nrow = 3) +
  draw_label("UMAP_1", x=0.5, y=  0, vjust=-0.5, angle= 0) +
  draw_label("UMAP_2", x=  0, y=0.5, vjust= 1.5, angle=90)
dev.off()

#4) Label the clusters according to scGate results
sobj <- RenameIdents(object = sobj,
                     "0" = "CD4+ T cells",
                     "1" = "CD4+ T cells", 
                     "4" = "CD4+ T cells",
                     "8" = "CD4+ T cells",
                     "10" = "CD4+ T cells",
                     "12" = "MAIT", 
                     "11" = "GD T cells",
                     "5" = "NK cells",
                     "3" = "CD8+ T cells",
                     "16" = "CD8+ T cells",
                     "9" = "CD8+ T cells",
                     "7" = "CD8+ T cells",
                     "6" = "NKT cells",
                     "2" = "CD8+ T cells")

sobj$clusters <- sobj@active.ident #define a variable "clusters" and assign active.ident to it

#Plot clusters proportion across all samples
tiff("../plots_CD8/barplot.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
ggplot(sobj@meta.data, aes(x=clusters, fill=sample_id)) + geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values= pal_ident) + xlab("") + ylab("") + coord_flip()
dev.off()

#Remove non-T, non-NK cells 
plot <- DimPlot(sobj, reduction = "umap")
T_NK <- CellSelector(plot=plot)
sobj <- subset(sobj, cells = T_NK)

p.subsets <- cowplot::plot_grid(pcd3, pcd8, pcd4, pnk, pgd, pclus, nrow = 2) 

tiff("../plots/all.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p.subsets #plot and save subsets
dev.off()

tiff("../plots_CD8/all_annot.tiff", width = 5*250, height = 5*250, res = 300, pointsize = 5)     
pclus2 <- DimPlot_scCustom(sobj, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001,
                           label.size = 7) & theme(legend.position = 'none')
pclus2 #plot and save annotated clusters
dev.off()

tiff("../plots_CD8/dotplot.tiff", width = 5*500, height = 5*180, res = 300, pointsize = 5)     
DotPlot_scCustom(sobj, features = c("CD3D","CD3G","CD8A","CD8B","CD4",
                                    "TRAV1-2","SLC4A10","NCAM1","KLRD1",
                                    "KLRG1","TRGC1", "TRDV1", "TRDV2", "TRDV3"),
                 colors_use = pal_exp, flip_axes = TRUE)
dev.off()

#5) Use gating models to subset and save into single rds objects

#CD8 (take all CD8+ excluding CD4, GD)
CD8 <- c("CD8A", "CD8B")
CD4 <- "CD4"
MAIT <- c("TRAV1-2","SLC4A10")
mmCD8 <- scGate::gating_model(level=1, name="CD8T", signature = CD8)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="CD4T", signature = CD4, negative=TRUE)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="MAIT", signature = MAIT, negative=TRUE)
CD8 <- scGate(sobj, model = mmCD8)
CD8sub <- subset(CD8, subset = `is.pure` == "Pure")
saveRDS(CD8sub, file = "CD8sub.rds")

#NK
NK <- c("NCAM1","KLRD1", "KLRG1")
CD3 <- c("CD3D","CD3G")
mmNK <- gating_model(level=1, name = "NK", signature = NK)
mmNK <- gating_model(model=mmNK, level=1, name = "CD3", signature = CD3, negative=TRUE)
NK <- scGate(sobj, model=mmNK)
NKsub <- NK
NKsub <- subset(NK, subset = `is.pure.level1` == "Pure")
saveRDS(NKsub, file = "NKsub.rds")

#GD
GD <- c("CD3D","TRGC1","TRDV1","TRDV2","TRDV3")
CD4 <- c("CD4")
CD8 <- c("CD8A","CD8B")
mmGD <- gating_model(level=1, name="GD", signature=GD)
mmGD <- gating_model(model=mmGD, level=1, name = "CD8", signature = CD8, negative=TRUE)
mmGD <- gating_model(model=mmGD, level=1, name = "CD4", signature = CD4, negative=TRUE)
GD <- scGate(sobj, model = mmGD)
GDsub <- GD
GDsub <- subset(GDsub, subset = `is.pure.level1` == "Pure")
saveRDS(GDsub, file = "GDsub.rds")

#Use cluster annotation to subset CD4 cells (identified by scGate)
CD4sub <- subset(sobj, clusters == "CD4+ T cells")
saveRDS(CD4sub, file = "CD4sub.rds")