#This script does the following:
#1) Pick the resolution for clustering
#2) Using scGate (https://github.com/carmonalab/scGate) creates a gating model to subset the different cell subsets
# This can be done also in Seurat but I prefer scGate because you can also exclude negative margers to make the model 
#more precise and subset the object on the output 
#3) According to the expression and the DGE (FindAllMarkers) labels the clusters and then creates a seurat object for 
#each subset

rm(list = ls())

#Load packages
library(Seurat)
library(scCustomize)
library(scGate)
library(ggplot2)
library(magrittr)

setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")
sobj <- readRDS("scRNAseq_step2_test.rds")

# set cluster IDs to resolution 1 clustering
sobj <- SetIdent(sobj, value = "integrated_snn_res.0.8")

DefaultAssay(sobj) <- "RNA" #Put "RNA" as default assay

pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")

#UMAP (Visualize clusters)
tiff("../plots_CD8/moannot.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
noAnn <- DimPlot_scCustom(sobj, colors_use = pal_ident, label.size = 6) + theme_void() & theme(legend.position = "none")
noAnn
dev.off()
#Plot different phenotypes
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

#Label subsets
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

sobj$clusters <- sobj@active.ident

tiff("../plots_CD8/barplot.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
ggplot(sobj@meta.data, aes(x=clusters, fill=sample_id)) + geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values= pal_ident) + xlab("") + ylab("") + coord_flip()
dev.off()

#Remove non-T, non-NK cells (cluster 13)
plot <- DimPlot(sobj, reduction = "umap")
T_NK <- CellSelector(plot=plot)
sobj <- subset(sobj, cells = T_NK)

DimPlot(sobj, split.by = "group_id")
pclus <- DimPlot_scCustom(sobj, label = TRUE, colors_use = pal_ident, figure_plot = TRUE)  
p.subsets <- cowplot::plot_grid(pcd3, pcd8, pcd4, pnk, pgd, pclus, nrow = 2)

tiff("../plots/all.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p.subsets
dev.off()

tiff("../plots_CD8/all_annot.tiff", width = 5*250, height = 5*250, res = 300, pointsize = 5)     
pclus2 <- DimPlot_scCustom(sobj, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001,
                           label.size = 7) & theme(legend.position = 'none')
pclus2
dev.off()

tiff("../plots_CD8/dotplot.tiff", width = 5*500, height = 5*180, res = 300, pointsize = 5)     
DotPlot_scCustom(sobj, features = c("CD3D","CD3G","CD8A","CD8B","CD4",
                                    "TRAV1-2","SLC4A10","NCAM1","KLRD1",
                                    "KLRG1","TRGC1", "TRDV1", "TRDV2", "TRDV3"),
                 colors_use = pal_exp, flip_axes = TRUE)
dev.off()

#Subsets
#CD8 (take all CD8+ excluding CD4, GD)
CD8 <- c("CD8A", "CD8B")
CD4 <- "CD4"
MAIT <- c("TRAV1-2","SLC4A10")
mmCD8 <- scGate::gating_model(level=1, name="CD8T", signature = CD8)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="CD4T", signature = CD4, negative=TRUE)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="MAIT", signature = MAIT, negative=TRUE)
CD8 <- scGate(sobj, model = mmCD8)
CD8sub <- subset(CD8, subset = `is.pure` == "Pure")
saveRDS(CD8sub, file = "CD8sub_test.rds")
DimPlot(CD8sub)

CD8sub <- subset(sobj, clusters == "CD8+ T cells")

#CD4
CD4sub <- subset(sobj, clusters == "CD4+ T cells")
saveRDS(CD4sub, file = "CD4sub.rds")

#NK
NK <- c("NCAM1","KLRD1", "KLRG1")
CD3 <- c("CD3D","CD3G")
mmNK <- gating_model(level=1, name = "NK", signature = NK)
mmNK <- gating_model(model=mmNK, level=1, name = "CD3", signature = CD3, negative=TRUE)
NK <- scGate(sobj, model=mmNK)
NKsub <- NK
NKsub <- subset(NK, subset = `is.pure.level1` == "Pure")

tiff("./plots/NKsub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(NKsub)
dev.off()

dim(NKsub)

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

tiff("./plots/GDsub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(GDsub)
dim(GDsub)

saveRDS(GDsub, file = "GDsub.rds")
