##CD8 This script can be organized into 7 parts

##A)Subclustering and use of custom and published markers along with bulk-RNA seq results and DGE to annotate clusters
#1) Subclustering as suggested here https://github.com/satijalab/Seurat/issues/2087 by timoast
#2) Use of custom single markers and visualize as DotPlot and FeaturePlot
#3) Use published signatures and signature of senescence obtained by bulk-seq and plot as feature and violin plots
#4) Correlation plot of dysfunction vs senescence signature
#5) DGE (top10 markers)
#5) Gene coexpression matrix and find genes with correlation > q99
#6) Annotate clusters

##B) Analyze annotated subsets
#1) Plot and look at clusters proportions and densities across samples, between conditions (group_id)
#2) DA analysis using permutation test, speckle (Anova) and CNA
#3) Combine top10 and custom markers and plot them as heatmap and as scores with FeaturePlot
#4) Plot scores obtained combining top10 and custom markers onto 2D UMAP
#5) Compute drivers of SL and ActEx via lasso regression

##C) Compute trajectories to infer CD8+ cells differentiation path
#1) Run Slingshot and plot the two trajectories
#2) Use tradeSeq to fit a generalized additive model (GAM) and then use patternTest to study gene expression

##D) Compute gene set variation analysis and infer pathway activity
#1) Use GSVA and limma to study gene set enrichments between couples of subsets (SLvsStL, SLvsActEx, ActExvsStL)
#2) Use SCpubr and decoupleR for pathway activity inference

##E) TCR repertoire analysis
#1) Export for scirpy
#2) scRepertoire analysis

##F) Reclustering with increased resolution and further analyses
#1) Reclustering and clusters frequencies across conditions
#2) STARTRAC method
#3) DGE including clusters SL1 and SL2

##G) Velocity inference
#1) Load in loom files and prepare to analyze in Python via scVelo
#2) import from python the csv with the info on latent time and different subsets and plot

##H) SCENIC
#1) Download the required databases
#2) Run SCENIC workflow
#3) Use exportsForArboreto to export matrix and scenicOptions and analyze in Python using GRNBoost 
#4) Read in GRNBoost output from Python and proceed with the workflow

rm(list = ls())

library(Seurat)
library(scGate)
library(slingshot)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(magrittr)
library(tidyr)
library(scRepertoire)
library(SeuratWrappers)
library(hdf5r)
library(circlize)
library(ComplexHeatmap)
library(tibble)
library(msigdbr)
library(scater)
library(gam)
library(tradeSeq)
library(stringr)
library(SeuratDisk)
library(viridis)
library(EnhancedVolcano)
library(scCustomize)
library(vip)
library(tidymodels)
library(scProportionTest)
library(decoupleR)
library(OmnipathR)
library(ggdist)
library(dplyr)
library(limma)
library(ggplotify)
library(ggthemes)
library(rcna)
library(glue)
library(ggnetwork)
library(ggforce)
library(speckle)
library(scry)
library(tidyverse)
library(AUCell)
library(RcisTarget)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
library(SCENIC)
library(arrow)

setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")
CD8 <-readRDS("CD8sub_test.rds")

#######################################################################################################################
#A) Subclustering and use of custom and published markers along with bulk-RNA seq results and DGE to annotate clusters
#######################################################################################################################

#1) Subclustering as suggested here https://github.com/satijalab/Seurat/issues/2087 by timoast
DefaultAssay(CD8) <- "RNA"
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", 
                            nfeatures = 2000, 
                            verbose = FALSE)
# Scale the counts in the integrated assay
DefaultAssay(CD8) <- "integrated"
CD8 <- ScaleData(CD8)
# Run PCA and UMAP
CD8 <- RunPCA(CD8)
CD8 <- RunUMAP(CD8, dims = 1:30, verbose = FALSE)
# Plot the elbow plot
ElbowPlot(object = CD8, ndims = 30)
# Determine the K-nearest neighbor graph

CD8 <- FindNeighbors(object = CD8, dims = 1:20)
# Determine the clusters for various resolutions                                
CD8 <- FindClusters(object = CD8, resolution = c(0.4, 0.6, 0.7, 0.8, 0.9, 1.0,1.2))
# set cluster IDs to resolution 0.4 clustering
CD8 <- SetIdent(CD8, value = "integrated_snn_res.0.4")
#Plot UMAP
DimPlot(CD8, label = T)
DefaultAssay(CD8) <- "RNA" #Put "RNA" as default assay
DimPlot(CD8, split.by = "group_id")
#Visualize the new subcluster and define custom palette
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
p.cd8 <- DimPlot_scCustom(CD8, label = TRUE, colors_use = pal_ident, pt.size = 0.00001, label.size = 6) + 
  theme_void() +
  theme(legend.position="none") 

#Save Fig. 2A 
tiff("../plots_CD8/cd8clus.tiff", width = 5*200, height = 5*200, res = 300, pointsize = 5)     
p.cd8
dev.off()

#2) Use of custom single markers and visualize as DotPlot and FeaturePlot

#DotPlot
features <- c("CCR7", "TCF7", "LEF1", "SELL", "IL7R", "SLAMF6", 
              "CXCR3", "GZMK", "CCL3", "CCL4", "CCL5", "XCL1", 
              "TIGIT", "PDCD1", "CD160", "LAG3", "CD69", "DUSP2", 
              "NR4A2", "ENTPD1", "ITGAE", "NKG7", "CX3CR1", "GNLY",
              "PRF1", "GZMB", "ZEB2", "CD226", "ZNF683",
              "TRAV1-2", "SLC4A10")   
#Define blue-white-red palette to be used consistently for plotting expression values
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)
tiff("../plots_CD8/cd8dot.tiff", width = 5*500, height = 5*250, res = 300, pointsize = 5)     
DotPlot_scCustom(CD8, features, colors_use = pal_exp, x_lab_rotate = TRUE) + 
  theme(axis.text.x = element_text(size = 8))
dev.off()

#FeaturePlot
p.CD8 <- FeaturePlot(CD8, features = "CD8A", pt.size = 0.00001, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")), breaks=c(0.3, 3.9), label = c("Min", "Max"))

#get legend to be used later on (same scale palette)
#get the legend from plot with same scale palette
legend <- get_legend(
  p.CD8 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "bottom",
                  legend.justification = "center") 
)

p <-FeaturePlot(CD8, features = c("TCF7", "IL7R", "GZMK", "CD69", "PDCD1", "GZMB", "PRF1", "GNLY"), combine=F, pt.size=0.00001, order=T) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position = "none")
}

tiff("../plots_CD8/p.markers.umap.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p.markers <- cowplot::plot_grid(plotlist = p, nrow =2)
p.mar_leg <- plot_grid(p.markers, legend, ncol = 1, rel_heights = c(1, .1))
p.mar_leg
dev.off()

#3) Use published signatures and signature of senescence obtained by bulk-seq and plot as feature and violin plots
#signatures
sig <- readxl::read_xlsx("../signatures/sig.xlsx", sheet = 1)
sig_naive <- list(sig$`Naive (from Szabo et al. 2019)`[!is.na(sig$`Naive (from Szabo et al. 2019)`)])
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive")
p.naive <- FeaturePlot(CD8, "sig_naive1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")), breaks=c(0.03, 0.18), label = c("Min", "Max")) +   
  ggtitle("Naive (from Szabo et al. 2019)") + theme(plot.title = element_text(size = 15, face = "bold")) + NoLegend() + NoAxes()
p.naive

sig_stem <- list(sig$`Stemness (from Pace et al. 2018)`[!is.na(sig$`Stemness (from Pace et al. 2018)`)])
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem")
p.stem <- FeaturePlot(CD8, "sig_stem1", pt.size = 0.000001, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("Stemness (from Pace et al. 2018)") + theme(plot.title = element_text(size = 15, face = "bold")) + NoLegend()+ NoAxes()

p.sig <- plot_grid(p.naive, p.stem, nrow = 1)
tiff("../plots_CD8/p.sig.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plot_grid(p.sig, legend, ncol = 1, rel_heights = c(1, .1)) 
dev.off()

#Bulk senescence signature
genes <- readxl::read_xlsx("../signatures/sig.xlsx", sheet = 2)
sen <- list(genes$gene_name[1:100])
CD8 <- AddModuleScore(CD8, features = sen, name = "sen")
genes$gene_name[1:100][(genes$gene_name[1:100] %in% rownames(CD8))]
p.sen <- FeaturePlot(CD8, features = "sen1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position="none") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("") + theme(plot.title = element_text(size = 15, face = "bold"))

tiff("../plots_CD8/score_sen.tiff", width = 5*180, height = 5*200, res = 300, pointsize = 5)     
p.sen <- p.sen + theme_void() + theme(legend.position = "none") + ggtitle("n = 43 genes") + theme(plot.title = element_text(hjust = 0.5))
plot_grid(p.sen, legend, ncol = 1, rel_heights = c(1, .1))
dev.off()

#Signature dysfunction June et al.
dys.genes <- readxl::read_xlsx("../signatures/sig.xlsx", sheet = 1)
dys <- list(dys.genes$`Dysfunction (Good et al. 2021)`)
CD8 <- AddModuleScore(CD8, features = dys, name = "dys")
p.dys <- FeaturePlot(CD8, features = "dys1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position="none") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) 

tiff("../plots_CD8/scores_dys.tiff", width = 5*180, height = 5*200, res = 300, pointsize = 5)     
p.dys <- p.dys + theme_void() + theme(legend.position = "none") + ggtitle("Dysfunction score (Good et al. 2021)") + 
  theme(plot.title = element_text(hjust = 0.5))
plot_grid(p.dys, legend, ncol = 1, rel_heights = c(1, .2))
dev.off()

#Distribution of senescence an disfunction across the 12 clusters
md <- CD8@meta.data %>% rename(sen1 = "Senescence", dys1 = "Dysfunction")
colnames(md)
CD8@meta.data <- md

#Stacked violin 
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab("") + ggtitle(feature) + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

# extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

#Function for stacked violin
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

tiff("../plots_CD8/stackvln.tiff", width = 5*180, height = 5*150, res = 300, pointsize = 5)     
StackedVlnPlot(obj = CD8, features = c("Senescence", "Dysfunction")) & scale_fill_manual(values = pal_ident)
dev.off()

#4) Correlation plot of dysfunction vs senescence signature
data<- CD8@assays$integrated@scale.data 
data<- t(data) %>% as.data.frame()
data$Senescence <- CD8$Senescence
data$Dysfunction <- CD8$Dysfunction
  
data$pc <- predict(prcomp(~Senescence+Dysfunction, data))[,1]

#Corr plot
tiff("../plots_CD8/corr.tiff", width = 5*200, height = 5*200, res = 300, pointsize = 5)     
ggplot(data, aes(x = Senescence, y = Dysfunction, color = pc)) +
  geom_point(shape = 16, size = 3, alpha = .4) +
  geom_smooth(method = "lm", color = "black", size=0.3, se = TRUE) +
  ggpubr::stat_cor() +
  theme_classic() +  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  ylab("Dysfunction score") + xlab("Senescence score") + 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))
dev.off()

#5) DGE (top10 markers)
mark <- FindAllMarkers(CD8)
mark %>% dplyr::filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

top10 <- top10[!duplicated(top10$gene),]
df <- data.frame(top10$cluster, top10$gene)
df  <-  reshape(transform(df, indx = ave(as.character(top10.cluster), top10.cluster, FUN = seq)), 
                idvar = "indx", timevar = "top10.cluster", direction = "wide") 
colnames(df) <- gsub("top10.gene.", "", colnames(df))
all_markers <- df %>%
  select(-indx) 
colnames(all_markers) <- as.character(colnames(all_markers))
all_m <- lapply(all_markers, function(x) x[!is.na(x)])

#Prepare Heatmap of mean marker-exprs. by cluster
un_cm <- unlist(all_m) 
num_mark <- vapply(all_m, length, numeric(1))
rep_clus <- rep.int(names(all_m), num_mark)
labs <- sprintf("%s(%s)", un_cm, rep_clus)

# split cells by cluster
cel_by_clus <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

m_by_clus <- lapply(all_m, function(un_cm)
  vapply(cel_by_clus, function(i)
    Matrix::rowMeans(log.ct[un_cm, i, drop = FALSE]), 
    numeric(length(un_cm))))

# prep. for plotting 
mat <- do.call("rbind", m_by_clus)

#Z-score
mat <- t(scale(t(mat)))
mat <- mat %>% as.data.frame() %>%  select(names(all_m)) %>% as.matrix()

lgd_aes <- list(direction = "horizontal", legend_width = unit(2.2, "cm"),
                title = "Expression")
h.heat.NoAnnot <- Heatmap(t(mat),
                       cluster_rows = TRUE,
                       cluster_columns = TRUE,
                       row_names_side = "left",
                       row_names_gp = grid::gpar(fontsize = 10),
                       heatmap_legend_param = lgd_aes,
                       column_names_gp = gpar(fontsize = 7))

tiff("../plots_CD8/heatNoAnnot.tiff", width = 5*500, height = 5*200, res = 300, pointsize = 5)     
p <- draw(h.heat.NoAnnot, heatmap_legend_side = "bottom", align_heatmap_legend = "heatmap_center", 
     show_annotation_legend = FALSE)
p
dev.off()

#5) Gene coexpression matrix and find genes with correlation > q99
mat <- cor(t(as.matrix(CD8@assays$RNA@data[VariableFeatures(CD8)[1:500],])),method = "spearman")
mat_corr <- mat
diag(mat_corr) <- 0 #give 0 to the diagonal
q <- quantile(mat_corr, probs = .99)
mat2 <-mat_corr %>% 
  as.data.frame() %>% 
  mutate(Res = ifelse(rowSums(mat_corr >= q) > 0, "Yes", "No"))
mat2 <- mat2 %>% filter(Res == "Yes")
mat2 <- mat2 %>% filter(!str_detect(rownames(mat2), "^RP[SL]"))

#169 genese above the 99th percentile
position1 <- which(rownames(mat) %in% rownames(mat2)[0:25])
position2 <- which(rownames(mat) %in% rownames(mat2)[26:50])
position3 <- which(rownames(mat) %in% rownames(mat2)[51:75])
position4 <- which(rownames(mat) %in% rownames(mat2)[76:100])
position5 <- which(rownames(mat) %in% rownames(mat2)[100:125])
position6 <- which(rownames(mat) %in% rownames(mat2)[126:150])
position7 <- which(rownames(mat) %in% rownames(mat2)[151:169])

#0-25
row_an1 <- rowAnnotation(Genes = anno_mark(at = position1,
                                          labels = rownames(mat)[position1],
                                          labels_gp = gpar(fontsize = 7),
                                          link_width = unit(2.5, "mm"),
                                          padding = unit(1, "mm"),
                                          link_gp = gpar(lwd = 0.5)))

ht.1 <- Heatmap(mat, name = "Spearman correlation",
              column_names_gp = grid::gpar(fontsize = 0),
              row_names_gp = grid::gpar(fontsize = 0),
              right_annotation = row_an1,
              heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.1, heatmap_legend_side = "bottom")

#26-50
row_an2 <- rowAnnotation(Genes = anno_mark(at = g1_pos,
                                           labels = rownames(mat)[g1_pos],
                                           labels_gp = gpar(fontsize = 7),
                                           link_width = unit(2.5, "mm"),
                                           padding = unit(1, "mm"),
                                           link_gp = gpar(lwd = 0.5)))

ht.2 <- Heatmap(mat, name = "Spearman correlation",
                column_names_gp = grid::gpar(fontsize = 0),
                row_names_gp = grid::gpar(fontsize = 0),
                right_annotation = row_an2,
                heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.2, heatmap_legend_side = "bottom")

#51-100
row_an3 <- rowAnnotation(Genes = anno_mark(at = position3,
                                           labels = rownames(mat)[position3],
                                           labels_gp = gpar(fontsize = 7),
                                           link_width = unit(2.5, "mm"),
                                           padding = unit(1, "mm"),
                                           link_gp = gpar(lwd = 0.5)))

ht.3 <- Heatmap(mat, name = "Spearman correlation",
                column_names_gp = grid::gpar(fontsize = 0),
                row_names_gp = grid::gpar(fontsize = 0),
                right_annotation = row_an3,
                heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.3, heatmap_legend_side = "bottom")

#101-125
row_an4 <- rowAnnotation(Genes = anno_mark(at = position4,
                                           labels = rownames(mat)[position4],
                                           labels_gp = gpar(fontsize = 7),
                                           link_width = unit(2.5, "mm"),
                                           padding = unit(1, "mm"),
                                           link_gp = gpar(lwd = 0.5)))

ht.4 <- Heatmap(mat, name = "Spearman correlation",
                column_names_gp = grid::gpar(fontsize = 0),
                row_names_gp = grid::gpar(fontsize = 0),
                right_annotation = row_an4,
                heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.4, heatmap_legend_side = "bottom")

#126-150
row_an5 <- rowAnnotation(Genes = anno_mark(at = position5,
                                           labels = rownames(mat)[position5],
                                           labels_gp = gpar(fontsize = 7),
                                           link_width = unit(2.5, "mm"),
                                           padding = unit(1, "mm"),
                                           link_gp = gpar(lwd = 0.5)))

ht.5 <- Heatmap(mat, name = "Spearman correlation",
                column_names_gp = grid::gpar(fontsize = 0),
                row_names_gp = grid::gpar(fontsize = 0),
                right_annotation = row_an5,
                heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.5, heatmap_legend_side = "bottom")

#151-169
row_an6 <- rowAnnotation(Genes = anno_mark(at = position6,
                                           labels = rownames(mat)[position6],
                                           labels_gp = gpar(fontsize = 7),
                                           link_width = unit(2.5, "mm"),
                                           padding = unit(1, "mm"),
                                           link_gp = gpar(lwd = 0.5)))

ht.6 <- Heatmap(mat, name = "Spearman correlation",
                column_names_gp = grid::gpar(fontsize = 0),
                row_names_gp = grid::gpar(fontsize = 0),
                right_annotation = row_an6,
                heatmap_legend_param = list(legend_direction = "horizontal")) 
draw(ht.6, heatmap_legend_side = "bottom")

DefaultAssay(CD8) <- "RNA" #Put "RNA" as default assay

#clean workspace
rm(list=setdiff(ls(), "CD8"))

#6) Annotate clusters
CD8 <- RenameIdents(CD8,
                    "0" = "Naive",
                    "1" = "StL",
                    "2" = "SL",
                    "3" = "ActEx",
                    "4" = "SL",
                    "5" = "Naive", 
                    "6" = "SL", 
                    "7" = "SL", 
                    "8" = "Naive",
                    "9" = "SL",
                    "10" = "StL",
                    "11" = "StL",
                    "12" = "ActEx")

#Add clusters variable to mddata
CD8$clusters<- CD8@active.ident

################################
##B) Analyze annotated subsets
################################

#1) Plot and look at clusters proportions and densities across samples, between conditions (group_id)
#Cluster distribution across all samples
tiff("../plots_CD8/barplotCD8.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
ggplot(CD8@meta.data, aes(x=clusters, fill=sample_id)) + geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values= pal_ident) + xlab("") + ylab("") + coord_flip()
dev.off()

tiff("../plots_CD8/UMAP_ann.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p.ann <- DimPlot_scCustom(CD8, label = TRUE, label.size = 4, colors_use = pal_ident[1:4], pt.size = 0.00001, figure_plot = T) + NoLegend()
p.ann 

tiff("../plots_CD8/UMAP_groupId.tiff", width = 5*300, height = 5*600, res = 300, pointsize = 5)     
p.abund <- DimPlot_scCustom(CD8, label = TRUE, split.by = "group_id", colors_use = pal_ident, split_seurat = TRUE, label.size = 6, num_columns = 1) + 
  theme_minimal(base_size = 35) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_discrete(breaks = NULL) +
  scale_x_discrete(breaks = NULL) + NoLegend()
p.abund 
dev.off()

tiff("../plots_CD8/UMAP_density.tiff", width = 5*500, height = 5*180, res = 300, pointsize = 5)    
CD8.dens <- subset(CD8, group_id %in% c("Res", "NonRes"))
p.dens <- DimPlot(CD8.dens, reduction = 'umap', split.by = "RespTmp")  + NoLegend() + NoAxes() 
p.dens +  geom_density_2d_filled(p.dens$data, mapping = aes(x = p.dens$data[,"UMAP_1"], y = p.dens$data[,"UMAP_2"]), contour_var = "ndensity") + 
  facet_wrap(vars(RespTmp), nrow = 1)
dev.off()

#2) DA analysis using permutation test, speckle (Anova) and CNA
#Permutation test 
prop_test <- sc_utils(CD8)

prop_test_ResVSNonRes <- permutation_test(
  prop_test, cluster_identity = "clusters",
  sample_1 = "Res", sample_2 = "NonRes",
  sample_identity = "group_id"
)

prop_test_HDVSRes <- permutation_test(
  prop_test, cluster_identity = "clusters",
  sample_1 = "HD", sample_2 = "Res",
  sample_identity = "group_id"
)

prop_test_HDVSNonRes <- permutation_test(
  prop_test, cluster_identity = "clusters",
  sample_1 = "HD", sample_2 = "NonRes",
  sample_identity = "group_id"
)

p.perm1 <- permutation_plot(prop_test_ResVSNonRes, log2FD_threshold = log2(2)) + 
  ggtitle("Res VS NonRes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold")) +
  scale_color_manual(values = c("red", "blue")) + theme(legend.position="none")
p.perm2 <- permutation_plot(prop_test_HDVSRes, log2FD_threshold = log2(2))+
  ggtitle("HD VS Res") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold")) +
  scale_color_manual(values = c("red", "blue")) +  theme(legend.position="none")
p.perm3 <- permutation_plot(prop_test_HDVSNonRes, log2FD_threshold = log2(2)) + 
  ggtitle("HD VS NonRes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face = "bold")) +
  scale_color_manual(values = c("red", "blue"))
?permutation_plot
p.perm <- plot_grid(p.perm1, p.perm2, p.perm3, rel_widths = c(1,1,1.9), nrow = 1)

tiff("../plots_CD8/sig.tiff", width = 5*300, height = 5*60, res = 150, pointsize = 5)     
p.perm
dev.off()

#Speckle
props <- getTransformedProps(CD8$clusters, CD8$sample_id,
                             transform="logit")

tiff("../plots_CD8/MeanVar.tiff", width = 5*150, height = 5*150, res = 300, pointsize = 5)     
plotCellTypePropsMeanVar(props$Counts)
dev.off()

condition <- c(rep("HD",2), rep("NonRes",2), rep("Res", 4), rep("NonRes", 2), 
               rep("Res", 2))
pair <- c(1,2, rep(c(3,4,5,6,7), each = 2))
batch <- c(rep(1,2), rep(2,2), rep(3,2), rep(2,2), rep(1,2), rep(2,2))
design <- model.matrix(~0 + condition + pair  + batch)
df.anova <- propeller.anova(prop.list=props, design=design, coef = c(1,2,3), 
                            robust=TRUE, trend=FALSE, sort=TRUE) %>% as.data.frame() 
df.anova <- df.anova %>% set_colnames(c("HD", "NonRes", "Res", "Fstat", "p.value", "FDR"))
fdr = 0.05
s <- factor(
  ifelse(df.anova$FDR < fdr, "yes", "no"), 
  levels = c("no", "yes"))
fdr_pal <- c("blue", "red3")
names(fdr_pal) <- levels(s)
F_pal <- c("blue",  "red3")
F_lims <- range(df.anova$Fstat)
F_brks <- c(0, F_lims[2])
Fstat = df.anova$Fstat
anno_cols <- list(Fstat = colorRamp2(F_brks, F_pal))
anno_cols$significant <- fdr_pal
right_anno <- rowAnnotation(
  Fstat = Fstat,
  significant = s,
  "foo" = row_anno_text(
    gp = gpar(fontsize = 8),
    scientific(df.anova$FDR, 2)),
  col = anno_cols,
  gp = gpar(col = "white"),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm"))

mat <- as.matrix(df.anova[,1:3])
png("../plots/HeatAnova.png", width = 5*250, height = 5*250, res = 300, pointsize = 5)     
Heatmap(mat, name = "Transformed proportions",
        right_annotation = right_anno,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        heatmap_legend_param = list(title_gp = gpar(
          fontsize = 10, fontface = "bold", lineheight = 0.8)))
dev.off()

#as boxplots
fqs <- props$Proportions
df.plot <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df.plot$group_id <- CD8$group_id[match(df.plot$sample_id, CD8$sample_id)]
df.plot <- df.plot %>% dplyr::mutate(positions = case_when(cluster_id == "StL" ~ "Anova: FDR ***, F 43.271719", TRUE ~ ""))
tiff("../plots_CD8/boxAnova.tiff", width = 5*400, height = 5*120, res = 300, pointsize = 5)     
pal_box <- pal_ident[c(10, 14, 35)]
p.box <- ggplot(df.plot, aes(x = group_id, y = frequency, fill = group_id)) +
  labs(x = NULL, y = "Proportion") + 
  geom_boxplot() + scale_fill_manual(values= pal_box) +
  geom_point(aes_string(x = "group_id"), position = position_jitter(width = 0.2)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.key.height  =  unit(0.8, "lines")) +
  facet_wrap(~ cluster_id, scales = "free_y", nrow = 1) + theme(legend.position = "none")
p.box + geom_text(aes(x=1.9,y=0.45,label=positions), size = 1.6)
dev.off()

#CNA (Res vs NonRes)
CD8.cna <- subset(CD8, group_id %in% c("Res", "NonRes"))
CD8.cna <- FindVariableFeatures(CD8.cna,selection.method = "vst", nfeatures = 2000) 
DefaultAssay(CD8.cna) <- "integrated"
CD8.cna <- ScaleData(CD8.cna) 
CD8.cna <- RunPCA(CD8.cna, npcs = 20, verbose = FALSE) 
CD8.cna <- FindNeighbors(CD8.cna, dims = 1:20, verbose = FALSE) 
CD8.cna <- FindClusters(CD8.cna, resolution = 0.4, verbose = FALSE)
CD8.cna <- RunUMAP(CD8.cna, dims = 1:10, verbose = FALSE)

CD8.cna@meta.data$status_val <- as.numeric(factor(CD8.cna@meta.data$group_id, c('Res', 'NonRes')))
CD8.cna@meta.data$batch <- as.numeric(factor(CD8.cna@meta.data$batch, c('b1', 'b2', 'b3')))
obj <- association.Seurat(
  seurat_object = CD8.cna, 
  test_var = 'status_val', 
  samplem_key = 'sample_id', 
  graph_use = 'integrated_snn',
  batches = "batch")
options(repr.plot.width=14, repr.plot.height=4)
p1 <- FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'Disease association', color = 'Correlation',
    subtitle = sprintf('Global p=%0.3f', obj@reductions$cna@misc$p)) 
p2 <- FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(title = '', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') 
p3 <- DimPlot(CD8.cna, group.by = "group_id") + ggtitle("Condition")

tiff("../plots_CD8/cna.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p2 + theme_void(base_size = 18)
dev.off()

DefaultAssay(CD8.cna) <- "RNA"
sig <- readxl::read_xlsx("../signatures/sig.xlsx", sheet = 1)
sig_naive <- list(sig$`Naive (from Szabo et al. 2019)`[!is.na(sig$`Naive (from Szabo et al. 2019)`)])
CD8.cna <- AddModuleScore(CD8.cna, features = sig_naive, name = "sig_naive")
p.naive.cna <- FeaturePlot(CD8.cna, "sig_naive1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("Naive (from Szabo et al. 2019)") + theme(plot.title = element_text(size = 15, face = "bold"))

sig_stem <- list(sig$`Stemness (from Pace et al. 2018)`[!is.na(sig$`Stemness (from Pace et al. 2018)`)])
CD8.cna <- AddModuleScore(CD8.cna, features = sig_stem, name = "sig_stem")
p.stem.cna <- FeaturePlot(CD8.cna, "sig_stem1", pt.size = 0.000001, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("Stemness (from Pace et al. 2018)") + theme(plot.title = element_text(size = 15, face = "bold"))

tiff("../plots_CD8/stem_cna.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p.stem.cna + theme_void(base_size = 18) + theme(legend.position = "none")
dev.off()

#3) Combine top10 and custom markers and plot them as heatmap and as scores with FeaturePlot
#Custom Markers Heatmap
DefaultAssay(CD8) <- "RNA"
cust_mark <- list(
  Naive = c("CCR7", "MAL", "LEF1", "SELL", "TCF7"),
  StL = c("BCL2", "BACH2", "CD27","IL7R", "SLAMF6", "CXCR3", "GZMK"),
  ActEx = c("EOMES", "CCL4", "XCL2", "CCL3", "XCL1", "KLF6", "TIGIT", "CD69", "CD160", "PDCD1", "TOX", "NR4A2",  "DUSP2"),
  SL = c("NKG7", "TBX21",  "CD38","FCRL6", "FCGR3A", "C1orf21", "PRF1", "ENO1",
         "GNLY", "ZEB2", "CX3CR1", "FGFBP2", "KLRD1", "GZMB","ZNF683", "CD226")
)

#Plot custom markers as scores
sig_naive <- list(cust_mark$Naive)
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive")
p.naive <- FeaturePlot(CD8, "sig_naive1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("Naive") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_stem <- list(cust_mark$StL)
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem")
p.stl <- FeaturePlot(CD8, "sig_stem1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("StL") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_actex <- list(cust_mark$ActEx)
CD8 <- AddModuleScore(CD8, features = sig_actex, name = "sig_actex")
p.actex <- FeaturePlot(CD8, "sig_actex1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("ActEx") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_sl <- list(cust_mark$SL)
CD8 <- AddModuleScore(CD8, features = sig_sl, name = "sig_sl")
p.sl <- FeaturePlot(CD8, "sig_sl1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("SL") + theme(plot.title = element_text(size = 10, face = "bold"))

tiff("../plots_CD8/p.featcust.tiff", width = 5*210, height = 5*210, res = 300, pointsize = 5)     
p.sig <- plot_grid(p.naive, p.stl, p.actex, p.sl, nrow = 2)
plot_grid(p.sig, legend, ncol = 1, rel_heights = c(1, .2)) 
dev.off()

#Prepare Heatmap of mean marker-exprs. by cluster
un_cm <- unlist(cust_mark)
num_mark <- vapply(cust_mark, length, numeric(1))
rep_clus <- rep.int(names(cust_mark), num_mark)
labs <- sprintf("%s(%s)", un_cm, rep_clus)

# split cells by cluster
cel_by_clus <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

m_by_clus <- lapply(cust_mark, function(un_cm)
  vapply(cel_by_clus, function(i)
    Matrix::rowMeans(log.ct[un_cm, i, drop = FALSE]), 
    numeric(length(un_cm))))

# prep. for plotting 
mat <- do.call("rbind", m_by_clus)

#Z-score
mat <- t(scale(t(mat)))

mat <- mat %>% as.data.frame() %>%  select(names(cust_mark)) %>% as.matrix()

cols <- pal_ident[seq_along(levels(CD8$clusters))]
cols <- setNames(cols, levels(CD8$clusters))
col_anno <- HeatmapAnnotation(
  df = data.frame(cluster_id = c("Naive", "StL", "ActEx", "SL")),
  col = list(cluster_id = cols, gp = gpar(col = "white"))) 
lgd_aes <- list(direction = "horizontal", legend_width = unit(2.2, "cm"),
                title = "Expression")
h.heat.cust <- Heatmap(mat,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       row_names_side = "left",
                       bottom_annotation = col_anno,
                       row_names_gp = grid::gpar(fontsize = 10),
                       heatmap_legend_param = lgd_aes,
                       column_names_gp = gpar(fontsize = 15))

tiff("../plots_CD8/p.heat_custom.tiff", width = 5*200, height = 5*400, res = 300, pointsize = 5)     
draw(h.heat.cust, heatmap_legend_side = "bottom", align_heatmap_legend = "heatmap_center", 
     show_annotation_legend = FALSE)
dev.off()

##Combine custom with top10 markers
#DGE all markers
mark <- FindAllMarkers(CD8)

#top10
mark %>% dplyr::filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

df <- data.frame(top10$cluster, top10$gene)
df <- df %>% distinct(top10.gene, .keep_all = TRUE) 
df  <-  reshape(transform(df, indx = ave(as.character(top10.cluster), top10.cluster, FUN = seq)), 
              idvar = "indx", timevar = "top10.cluster", direction = "wide") 
colnames(df) <- gsub("top10.gene.", "", colnames(df))
df <- df %>% relocate(SL, .after = ActEx)

all_markers <- df %>%
  select(-indx) 
all_markers <- lapply(all_markers, function(x) x[!is.na(x)])
comb_mark <- cust_mark
for(i in seq_along(cust_mark)){
  comb_mark[[i]] <- c(cust_mark[[i]], all_markers[[i]][!all_markers[[i]]%in% cust_mark[[i]]])
}

#check for duplicates in subsets
names(which(table(unlist(comb_mark)) > 1)) #CX3CR1
lapply(comb_mark, function(x) grep("CX3CR1", x))
comb_mark$StL <- comb_mark$StL[-11] #remove duplicate

idx <- which(unlist(comb_mark) %in% unlist(all_markers)) #to use after to color genes

#Prepare Heatmap of mean marker-exprs. by cluster
un_cm <- unlist(comb_mark)
num_mark <- vapply(comb_mark, length, numeric(1))
rep_clus <- rep.int(names(comb_mark), num_mark)
labs <- sprintf("%s(%s)", un_cm, rep_clus)

# split cells by cluster
cel_by_clus <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

m_by_clus <- lapply(comb_mark, function(un_cm)
  vapply(cel_by_clus, function(i)
    Matrix::rowMeans(log.ct[un_cm, i, drop = FALSE]), 
    numeric(length(un_cm))))

# prep. for plotting 
mat <- do.call("rbind", m_by_clus)

#Z-score
mat <- t(scale(t(mat)))

mat <- mat %>% as.data.frame() %>%  select(names(comb_mark)) %>% as.matrix()

cols <- pal_ident[seq_along(levels(CD8$clusters))]
cols <- setNames(cols, levels(CD8$clusters))
col_anno <- HeatmapAnnotation(
  df = data.frame(cluster_id = c("Naive", "StL", "ActEx", "SL")),
  col = list(cluster_id = cols, gp = gpar(col = "white")), 
  show_legend = c(FALSE, TRUE)) 
graphics = list(
  "Top10 genes" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "darkgreen"), pch = 16)
  },
  "Custom genes" = function(x, y, w, h) {
    grid.points(x, y, gp = gpar(col = "darkviolet"), pch = 16)
  }
)
lgd = Legend(title = "", at = names(graphics), graphics = graphics)

lgd_aes <- list(direction = "vertical", legend_width = unit(2.2, "cm"),
                title = "Expression")

colLab <- rep("darkviolet", nrow(mat)) 
colLab[idx] <- "darkgreen"

h.dge <- Heatmap(mat,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 row_names_side = "left",
                 bottom_annotation = col_anno,
                 row_names_gp = grid::gpar(col = colLab, fontsize = 8),
                 heatmap_legend_param = lgd_aes,
                 column_names_gp = gpar(fontsize = 14))


tiff("../plots_CD8/p.heat_combo.tiff", width = 5*250, height = 5*400, res = 300, pointsize = 5)     
draw(h.dge,annotation_legend_list = lgd)
dev.off()

#plot scores
sig_naive <- list(comb_mark$Naive)
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive")
p.naive <- FeaturePlot(CD8, "sig_naive1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("Naive") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_stem <- list(comb_mark$StL)
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem")
p.stl <- FeaturePlot(CD8, "sig_stem1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("StL") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_actex <- list(comb_mark$ActEx)
CD8 <- AddModuleScore(CD8, features = sig_actex, name = "sig_actex")
p.actex <- FeaturePlot(CD8, "sig_actex1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("ActEx") + theme(plot.title = element_text(size = 10, face = "bold"))

sig_sl <- list(comb_mark$SL)
CD8 <- AddModuleScore(CD8, features = sig_sl, name = "sig_sl")
p.sl <- FeaturePlot(CD8, "sig_sl1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  NoAxes() + NoLegend() + ggtitle("SL") + theme(plot.title = element_text(size = 10, face = "bold"))

tiff("../plots_CD8/p.featcust.tiff", width = 5*210, height = 5*210, res = 300, pointsize = 5)     
p.sig <- plot_grid(p.naive, p.stl, p.actex, p.sl, nrow = 2)
plot_grid(p.sig, legend, ncol = 1, rel_heights = c(1, .2)) 
dev.off()

#Drivers of SL and ActEx
CD8.new <- CD8
SL<- subset(CD8.new, clusters %in% c("SL", "ActEx"))
Idents(SL)<- SL$clusters

SL<- SL %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.1, verbose = FALSE) %>%
  RunUMAP(dims = 1:10, verbose = FALSE)

data <- SL@assays$RNA@scale.data

# let's transpose the matrix and make it to a dataframe
dim(data)
data <- t(data) %>% as.data.frame()

## add the cell type/the outcome/y to the dataframe
data$cell_type<- SL$clusters
data$cell_barcode<- rownames(data)
## it is important to turn it to a factor for classification
data$cell_type<- factor(data$cell_type)
View(data)
set.seed(123)
data_split <- initial_split(data, strata = "cell_type")
data_train <- training(data_split)
data_test <- testing(data_split)

# 10 fold cross validation
data_fold <- vfold_cv(data_train, v = 10)

#Lasso regression
lasso_spec<-
  logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

lasso_recipe <- 
  recipe(formula = cell_type ~ ., data = data_train) %>% 
  step_zv(all_predictors())

# step_normalize(all_predictors())
## the expression is already scaled, no need to do step_normalize

lasso_workflow <- workflow() %>% 
  add_recipe(lasso_recipe) %>% 
  add_model(lasso_spec)

penalty_grid <- grid_regular(penalty(range = c(-5, 5)), levels = 50)
#penalty_grid
tune_res <- tune_grid(
  lasso_workflow,
  resamples = data_fold, 
  grid = penalty_grid
)

best_penalty <- select_best(tune_res, metric = "accuracy")
best_penalty

lasso_final <- finalize_workflow(lasso_workflow, best_penalty)
lasso_final_fit <- fit(lasso_final, data = data_train)

## confusion matrix
predict(lasso_final_fit, new_data = data_test) %>%
  bind_cols(data_test %>% select(cell_type)) %>%
  conf_mat(truth = cell_type, estimate = .pred_class)

lasso_features<- tidy(lasso_final_fit) %>% 
  arrange(desc(abs(estimate))) %>%
  filter(estimate != 0) 

Idents(SL) <- SL$clusters

tiff("../plots/lasso_genes.tiff", width = 5*150, height = 5*400, res = 300, pointsize = 5)     
scCustomize::Stacked_VlnPlot(SL, features = lasso_features %>% pull(term) %>% head(n = 20),
                             colors_use = c("blue", "red")) 
dev.off()

##################################################################
##C) Compute trajectories to infer CD8+ cells differentiation path
##################################################################

#1) Run Slingshot and plot the two trajectories
clusterLabels <- CD8$clusters
sceCD8 <- as.SingleCellExperiment(CD8, assay = "RNA")
sds <- slingshot(sceCD8, clusterLabels = clusterLabels, 
                  allow.breaks = TRUE, stretch = 2, reducedDim = "UMAP", start.clus = "Naive") #Calcualting the trajectory
sds <- SlingshotDataSet(sds)

df <- bind_cols(
  as.data.frame(reducedDim(sds, "UMAP")),
  slingPseudotime(sds) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(sds) %>% as.data.frame(),
) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst))

curves <- slingCurves(sds, as.df = TRUE)

#Look at single lineages
#Lineage 1 
p1 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes_string(color = df$Lineage1_pst),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  theme_minimal() + labs(colour = "Pseudotime") 

#Lineage 2
p2 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes_string(color = df$Lineage2_pst),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  theme_minimal() + labs(colour = "Pseudotime") 
plot_grid(p1,p2)

#Together
p.traj <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = -7.7, y = 3.7, label = "ActEx", size = 5) +
  annotate("text", x = -8, y = -5.6, label = "SL", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +  theme_minimal()  

tiff("../plots_CD8/umap_traj.tiff", width = 5*350, height = 5*250, res = 300, pointsize = 5)     
p.traj +theme_void()
dev.off()

#2) Use tradeSeq to fit a generalized additive model (GAM) and then use patternTest to study gene expression
set.seed(5)
par(mar=c(1,1,1,1))
icMat <- evaluateK(counts = counts(sceCD8), sds = sds, k = 3:7, 
                   nGenes = 100, verbose = T, plot = TRUE)
print(icMat[1:2,])
set.seed(7)

###Parallel computing
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 2

#Use deviance for feature selection
sce.dev <-devianceFeatureSelection(sceCD8, assay="counts", sorted=TRUE)
plot(rowData(sce.dev)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
#fitGAM
sce.dev <- sce.dev[1:1000,] #one thousand genes with highest deviance
genes <- which(rownames(sceCD8) %in% rownames(sce.dev))
sce.gam <- fitGAM(counts = counts(sceCD8), sds = sds, nknots = 5, verbose = TRUE, 
                  genes = genes, BPPARAM = BPPARAM)
#Save
#saveRDS(sce.gam, "scegam.rds")
#sce.gam <- readRDS("scegam.rds")

# plot our Slingshot lineage trajectories, this time illustrating the new tradeSeq knots
tiff("./plots/traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plotGeneCount(curve = sds, counts = counts,
              clusters = CD8@active.ident,
              models = sce.gam)
dev.off()

### Discovering differentiated cell type markers
# discover marker genes for the differentiated cell types
#Genes with different expression patterns (most interesting part)
patternRes <- patternTest(sce.gam)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
rownames(patternRes)[oPat][1:20]

p.gzmk <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][1], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("GZMK") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

legend <- cowplot::get_legend(p.leg + theme(legend.position = "right"))

p.prf1 <-plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][5], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("PRF1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p.gzmb <-plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][6], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("GZMB") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p.fgfbp2 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][7], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("FGFBP2") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p.gnly <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][9], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("GNLY") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 


p.nkg7 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][15], xlab = "", ylab = "") + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "SL", "2"="ActEx")) +
  theme(legend.position = "none") +
  ggtitle ("NKG7") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

p.patt <- plot_grid(p.gzmk, p.gzmb, p.gnly, p.prf1, p.nkg7, p.fgfbp2) +
  draw_label("Pseudotime", x=0.5, y=  0, vjust=-0.5, angle= 0) +
  draw_label("Log(Expression +1)", x=  0, y=0.5, vjust= 1.5, angle=90)
tiff("../plots_CD8/genes_traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plot_grid(p.patt, legend, rel_widths = c(4,1))
dev.off()

#Heatmap
genes <- rownames(patternRes)[oPat]
df.tr <- predictSmooth(sce.gam, gene = genes, tidy = FALSE)
df <- as.data.frame(df.tr)

#sl = senescent-like; ae = activated-exhausted
df.ls <- map(set_names(c("lineage1", "lineage2")),~select(df,starts_with(.x)))
mat.sl <- df.ls$lineage1 
mat.sl <- mat.sl %>% filter(!str_detect(rownames(mat.sl), "^RP[SL]|^MT-")) %>% 
  as.data.frame.matrix()
mat.sl <- t(scale(t(mat.sl)))
genes_int <- c("GZMK", "PRF1", "GZMB", "FGFBP2", "GNLY", "CTSW", "CD69", "NKG7", "KLRG1", "KLRD1", "DUSP2",
               "NR4A2", "CD27","GZMA", "CCL5", "SELL", "LGALS1", "NFKBIA", "LTB", "TCF7",
               "TNFAIP3", "NELL2", "GZMH", "CST7","MALAT1", "IL7R", "MAL", "ZNF683",
               "CCR7", "LDHA", "KLRB1", "CCL4", "CX3CR1", "LEF1", "ADGRG1", "CCL3",
               "TIGIT", "FCRL6", "XCL1", "FCGR3A", "GZMM", "ZEB2")

heat.sl <- Heatmap(mat.sl[1:100,],
                   show_column_names = FALSE,
                   row_names_gp = grid::gpar(fontsize = 0),
                   cluster_columns = FALSE,
                   cluster_rows = TRUE,
                   show_heatmap_legend = FALSE,
                   column_title = "SL")

mat.ae <- df.ls$lineage2 
mat.ae <- mat.ae %>% filter(!str_detect(rownames(mat.ae), "^RP[SL]|^MT-")) %>% 
  as.data.frame.matrix()
mat.ae <- t(scale(t(mat.ae)))

position <- which(rownames(mat.ae) %in% genes_int)
row_an.ae <- rowAnnotation(Genes = anno_mark(at = which(rownames(mat.ae) %in% genes_int),
                                             labels = rownames(mat.ae)[position],
                                             labels_gp = gpar(fontsize = 7),
                                             link_width = unit(2.5, "mm"),
                                             padding = unit(1, "mm"),
                                             link_gp = gpar(lwd = 0.5)))


lgd_aes <- list(legend_width = unit(2.2, "cm"),
                title = "Expression")

heat.ae <- Heatmap(mat.ae[1:100,],
                   show_column_names = FALSE,
                   row_names_gp = grid::gpar(fontsize = 0),
                   right_annotation = row_an.ae,
                   cluster_columns = FALSE,
                   cluster_rows = T,
                   column_title = "ActEx",
                   heatmap_legend_param = lgd_aes)

tiff("../plots_CD8/heatTraj.tiff", width = 5*350, height = 5*400, res = 300, pointsize = 5)     
heat.sl + heat.ae 
dev.off()

###################################################################
##D) Compute gene set variation analysis and infer pathway activity
###################################################################

#1) Use GSVA and limma to study gene set enrichments between couples of subsets (SLvsStL, SLvsActEx, ActExvsStL)
# function to read GMT file
read_GMT_file <- function(file) {
  gmt <- readr::read_delim(
    file,
    delim = ';',
    col_names = c('X1'),
    col_types = readr::cols()
  )
  
  gene_set_genes <- list()
  for ( i in seq_len(nrow(gmt)) )
  {
    temp_genes <- strsplit(gmt$X1[i], split = '\t')[[1]] %>% unlist()
    temp_genes <- temp_genes[3:length(temp_genes)]
    gene_set_genes[[i]] <- temp_genes
  }
  gene_set_loaded <- list(
    genesets = gene_set_genes,
    geneset.names = lapply(strsplit(gmt$X1, split = '\t'), '[', 1) %>% unlist(),
    geneset.description = lapply(
      strsplit(gmt$X1, split = '\t'), '[', 2
    ) %>% unlist()
  )
  
  return(gene_set_loaded)
}

# load gene sets from GMT file
gene_sets <- read_GMT_file('../h.all.v7.2.symbols.gmt')

# set gene set names
names(gene_sets$genesets) <- gene_sets$geneset.names

##Senescence-like VS Stem-like
# get indices of cells which are either SL or StL
c_SLvStL <- CD8@meta.data %>%
  mutate(row_number = row_number()) %>%
  dplyr::filter(grepl(clusters, pattern = 'SL|StL')) %>%
  arrange(clusters) %>%
  pull(row_number)

# get list of genes unique genes across all gene sets
genes_to_analyze <- gene_sets$genesets %>% unlist() %>% unique()
# filter gene list for those which are present in the data set
genes_to_analyze <- genes_to_analyze[which(genes_to_analyze %in% rownames(CD8@assays$RNA@counts))]

# get expression matrix and reduce it to cells and genes of interest
mat_SLvStL <- CD8@assays$RNA@counts[ genes_to_analyze , c_SLvStL] %>% as.matrix()

# perform GSVA
g_SLvStL <- GSVA::gsva(
  mat_SLvStL,
  gset.idx.list = gene_sets$genesets,
  parallel.sz = 1
)

# generate design matrix
dm_SLvStL <- tibble(
  control = 1,
  test = c(
    rep(0, CD8@meta.data %>% dplyr::filter(clusters == 'StL') %>% nrow()),
    rep(1, CD8@meta.data %>% dplyr::filter(clusters == 'SL') %>% nrow())
  )
)

# fit linear model, followed by empirical Bayes statistics for differential
# enrichment analysis
fit_SLvStL <- lmFit(g_SLvStL, dm_SLvStL)
fit_SLvStL <- eBayes(fit_SLvStL)

# prepare data for plotting
data_SLvStL <- topTable(fit_SLvStL, coef = 'test', number = 50) %>%
  mutate(gene_set = rownames(fit_SLvStL$t)) %>%
  arrange(t) %>%
  mutate(
    gene_set = factor(gene_set, levels = gene_set),
    just = ifelse(t < 0, 0, 1),
    nudge_y = ifelse(t < 0, 1, -1),
    color = ifelse(t < -5 | t > 5, 'black', 'grey')
  )

# plot t-value
slVSstl <- ggplot(data = data_SLvStL, aes(x = gene_set, y = t, fill = t)) +
  geom_col() +
  geom_hline(yintercept = c(-5,5), linetype = 'dashed', color = 'grey80') +
  geom_text(
    aes(
      x = gene_set,
      y = 0,
      label = gene_set,
      hjust = just,
      color = color
    ),
    nudge_y = data_SLvStL$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data_SLvStL$t), max(data_SLvStL$t))) +
  scale_color_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y =  element_blank(),
    legend.position = 'none'
  ) + ggtitle("Senescent-like VS Stem-like") +
  theme(plot.title = element_text(hjust = 0.5)) 

##Senescence-like VS Activated Exhausted
# get indices of cells which are either SL or StL
c_SLvActEx <- CD8@meta.data %>%
  mutate(row_number = row_number()) %>%
  dplyr::filter(grepl(clusters, pattern = 'SL|ActEx')) %>%
  arrange(clusters) %>%
  pull(row_number)

# get expression matrix and reduce it to cells and genes of interest
mat_SLvActEx <- CD8@assays$RNA@counts[ genes_to_analyze , c_SLvActEx] %>% as.matrix()

# perform GSVA
g_SLvActEx <- GSVA::gsva(
  mat_SLvActEx,
  gset.idx.list = gene_sets$genesets,
  parallel.sz = 1
)

# generate design matrix
dm_SLvActEx <- tibble(
  control = 1,
  test = c(
    rep(0, CD8@meta.data %>% dplyr::filter(clusters == 'ActEx') %>% nrow()),
    rep(1, CD8@meta.data %>% dplyr::filter(clusters == 'SL') %>% nrow())
  )
)

# fit linear model, followed by empirical Bayes statistics for differential
# enrichment analysis
fit_SLvActEx <- lmFit(g_SLvActEx, dm_SLvActEx)
fit_SLvActEx <- eBayes(fit_SLvActEx)

# prepare data for plotting
data_SLvActEx <- topTable(fit_SLvActEx, coef = 'test', number = 50) %>%
  mutate(gene_set = rownames(fit_SLvActEx$t)) %>%
  arrange(t) %>%
  mutate(
    gene_set = factor(gene_set, levels = gene_set),
    just = ifelse(t < 0, 0, 1),
    nudge_y = ifelse(t < 0, 1, -1),
    color = ifelse(t < -5 | t > 5, 'black', 'grey')
  )

# plot t-value
SlVSActEx <- ggplot(data = data_SLvActEx, aes(x = gene_set, y = t, fill = t)) +
  geom_col() +
  geom_hline(yintercept = c(-5,5), linetype = 'dashed', color = 'grey80') +
  geom_text(
    aes(
      x = gene_set,
      y = 0,
      label = gene_set,
      hjust = just,
      color = color
    ),
    nudge_y = data_SLvActEx$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data_SLvActEx$t), max(data_SLvActEx$t))) +
  scale_color_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y =  element_blank(),
    legend.position = 'none'
  ) + ggtitle("Senescent-like VS Activated-exhausted") +
  theme(plot.title = element_text(hjust = 0.5)) 

##Activated-Exhausted VS Stem-like
# get indices of cells which are either SL or StL
c_ActExvStL <- CD8@meta.data %>%
  mutate(row_number = row_number()) %>%
  dplyr::filter(grepl(clusters, pattern = 'ActEx|StL')) %>%
  arrange(clusters) %>%
  pull(row_number)

# get expression matrix and reduce it to cells and genes of interest
mat_ActExvStL <- CD8@assays$RNA@counts[ genes_to_analyze , c_ActExvStL] %>% as.matrix()

# perform GSVA
g_ActExvStL <- GSVA::gsva(
  mat_ActExvStL,
  gset.idx.list = gene_sets$genesets,
  parallel.sz = 1
)

# generate design matrix
dm_ActExvStL <- tibble(
  control = 1,
  test = c(
    rep(0, CD8@meta.data %>% dplyr::filter(clusters == 'StL') %>% nrow()),
    rep(1, CD8@meta.data %>% dplyr::filter(clusters == 'ActEx') %>% nrow())
  )
)

# fit linear model, followed by empirical Bayes statistics for differential
# enrichment analysis
fit_ActExvStL<- lmFit(g_ActExvStL, dm_ActExvStL)
fit_ActExvStL<- eBayes(fit_ActExvStL)

# prepare data for plotting
data_ActExvStL <- topTable(fit_ActExvStL, coef = 'test', number = 50) %>%
  mutate(gene_set = rownames(fit_ActExvStL$t)) %>%
  arrange(t) %>%
  mutate(
    gene_set = factor(gene_set, levels = gene_set),
    just = ifelse(t < 0, 0, 1),
    nudge_y = ifelse(t < 0, 1, -1),
    color = ifelse(t < -5 | t > 5, 'black', 'grey')
  )

# plot t-value
ActExvStL <- ggplot(data = data_ActExvStL, aes(x = gene_set, y = t, fill = t)) +
  geom_col() +
  geom_hline(yintercept = c(-5,5), linetype = 'dashed', color = 'grey80') +
  geom_text(
    aes(
      x = gene_set,
      y = 0,
      label = gene_set,
      hjust = just,
      color = color
    ),
    nudge_y = data_ActExvStL$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data_ActExvStL$t), max(data_ActExvStL$t))) +
  scale_color_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y =  element_blank(),
    legend.position = 'none'
  ) + ggtitle("Activated-exhausted VS Stem-like") +
  theme(plot.title = element_text(hjust = 0.5)) 

tiff("../plots_CD8/GSVA_plots.tiff", width = 5*830, height = 5*300, res = 300, pointsize = 5)     
plot_grid(slVSstl, SlVSActEx, nrow = 1)
dev.off()

tiff("../plots_CD8/GSVA_plots2.tiff", width = 5*415, height = 5*300, res = 300, pointsize = 5)     
ActExvStL
dev.off()

#2) Use SCpubr and decoupleR for pathway activity inference
# Retrieve prior knowledge network.
network <- decoupleR::get_progeny(organism = "human")

# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(CD8@assays[["RNA"]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)

# General heatmap.
out <- SCpubr::do_PathwayActivityPlot(sample = CD8,
                                      activities = activities,
                                      split.by = "group_id")
p.patHeat <- out$heatmaps$average_scores

tiff("../plots_CD8/heat_hyp.tiff", width = 5*400, height = 5*500, res = 300, pointsize = 5)     
p.patHeat
dev.off()

#Since hypoxia is up in GSVA and pathway analysis LDHA might be increased in SL
CD8@active.ident <- CD8$clusters
tiff("../plots_CD8/DP_LDHA.tiff", width = 5*200, height = 5*180, res = 300, pointsize = 5)     
DotPlot_scCustom(CD8, features = "LDHA",colors_use = pal_exp, x_lab_rotate = TRUE) + 
  theme(axis.text.x = element_text(size = 8))
dev.off()

############################
##E) TCR repertoire analysis
############################

#1) Export for scirpy
CD8.clono <- CD8
CD8.clono <- subset(CD8.clono, subset = group_id == "HD", invert = TRUE) #we do not have clonotype data for the HDs
CD8.clono$batch <- NULL
md <- CD8.clono@meta.data
md$clusters <- as.vector(md$clusters)
CD8.clono@meta.data <- md
SaveH5Seurat(CD8.clono, filename = "CD8_clonoscv.h5Seurat")
Convert("CD8_clonoscv.h5Seurat", dest = "h5ad") 

SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat")
Convert("pbmc3k.h5Seurat", dest = "h5ad")

#2) scRepertoire analysis

#Clonotype analysis
#Delete HD (we have VDJ only for responders and nonresponders)
CD8.clono <- subset(CD8, group_id %in% c("Res", "NonRes"))
md <- CD8.clono@meta.data

#Define hyper and large as expanded, all the others as not expanded
md <- md %>% mutate(`Clonal expansion` = case_when(cloneType %in% c("Hyperexpanded (100 < X <= 500)","Large (20 < X <= 100)") ~ 'Expanded',
                                           TRUE ~ "Not Expanded"))
CD8.clono@meta.data <- md
CD8.clono$`Clonal expansion` %>% as.factor()
tiff("../plots_CD8/clonalExp.tiff", width = 5*300, height = 5*250, res = 300, pointsize = 5)     
DimPlot_scCustom(CD8.clono, group.by = "Clonal expansion", colors_use =  c("lightgrey", "darkred"), order = c("Expanded", "Not Expanded"),
                 pt.size = 0.00001) + ggtitle("") + theme_void()
dev.off()

tiff("../plots_CD8/clonoContour.tiff", width = 5*350, height = 5*250, res = 300, pointsize = 5)     
p <- clonalOverlay(CD8.clono, reduction = "umap", 
                   freq.cutpoint = 20, bins = 25, facet = "group_id") + 
  guides(color = FALSE) 
p + theme_void(base_size = 20) + scale_color_manual(values=pal_ident)
dev.off()

CD8.clono_list <- SplitObject(CD8.clono, split.by = "RespTmp")
occrep <- lapply(CD8.clono_list, function(x) occupiedscRepertoire(x, x.axis = "clusters", label = FALSE))
occrep_res.bas <- occrep[[1]] + ggtitle ("Res_bas") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                                                            legend.title=element_blank(),
                                                            axis.text.x = element_text(size = 20),
                                                            axis.text.y = element_text(size = 20),
                                                            legend.text=element_text(size=15)) 
occrep_res.post <- occrep[[2]] + ggtitle ("Res_post")+ theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                                                             legend.title=element_blank(),
                                                             axis.text.x = element_text(size = 15),
                                                             axis.text.y = element_text(size = 15),
                                                             legend.text=element_text(size=15))
occrep_Nres.bas <- occrep[[3]] + ggtitle ("NonRes_bas") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                                                                legend.title=element_blank(),
                                                                axis.text.x = element_text(size = 15),
                                                                axis.text.y = element_text(size = 15),
                                                                legend.text=element_text(size=15)) 
occrep_Nres.post <- occrep[[4]] + ggtitle ("NonRes_post")+  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                                                                  legend.title=element_blank(),
                                                                  axis.text.x = element_text(size = 20),
                                                                  axis.text.y = element_text(size = 20),
                                                                  legend.text=element_text(size=15))

tiff("../plots_CD8/occrep.tiff", width = 5*850, height = 5*400, res = 300, pointsize = 5)     
cowplot::plot_grid(occrep_res.bas, occrep_res.post, occrep_Nres.bas, occrep_Nres.post, nrow =2) 
dev.off()

#################
##F) Reclustering
#################

#1) Reclustering and clusters frequencies across conditions
#Different SL like clusters
CD8 <- SetIdent(CD8, value = "integrated_snn_res.1.2")
p.umap2 <- DimPlot(CD8, label = T, split.by = "group_id")

#Better look at senescent cells
CD8 <- RenameIdents(CD8, 
                    "0" = "Naive",
                    "1" = "StL",
                    "2" = "SL1",
                    "3" = "Naive",
                    "4" = "SL2",
                    "5" = "ActEx",
                    "6" = "Naive",
                    "7" = "StL",
                    "8" = "SL1",
                    "9" = "ActEx",
                    "10" = "SL1",
                    "11" = "SL2",
                    "12" = "Naive",
                    "13" = "SL2",
                    "14" = "StL",
                    "15" = "SL2",
                    "16" = "StL",
                    "17" = "StL",
                    "18" = "SL2",
                    "19" = "ActEx",
                    "20" = "ActEx")

CD8$clusters2 <- CD8@active.ident
CD8$clusters2 <- factor(CD8$clusters2, levels = c("Naive", "StL",   "SL1", "ActEx", "SL2"))
CD8@active.ident <- CD8$clusters2
tiff("../plots_CD8/umap_sen.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
p <- DimPlot_scCustom(CD8, label = TRUE, label.size = 8, colors_use = pal_ident, pt.size = 0.00001) + NoAxes()+
  theme(legend.position="none") 
p
#df <- data.frame(x1 = -5.4, x2 = -0.4, y1 = -4.3, y2 = -0.2)
#p +  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), face = "bold", linetype = 2, data = df)
dev.off()

cluster_stats <- Cluster_Stats_All_Samples(CD8, group_by_var = "group_id") %>%  
  mutate(freq_groupID = log10(Res/NonRes)) %>%  
  dplyr::filter(!str_detect("Total", Cluster))

cluster_stats$Cluster <-  factor(cluster_stats$Cluster, levels = c("Naive", "StL",  "SL1","ActEx", "SL2"))

colScale <- scale_colour_manual(name = "Cluster",values = pal_ident)

tiff("../plots_CD8/dot_group.tiff", width = 5*80, height = 5*60, res = 150, pointsize = 5)     
cluster_stats %>% ggplot(aes(x=Cluster, y = freq_groupID)) + 
  geom_point(aes(size = abs(freq_groupID), colour = Cluster)) + theme_classic() + ylim(-1.5, 1.5) + geom_hline(yintercept=0,linetype=2) + ylab("log10(Res/NonRes freq)") +
  theme(legend.position = "none")+ colScale
dev.off()

#2) STARTRAC method 
CD8.clono <- subset(CD8, group_id %in% c("Res", "NonRes")) #subset again to include clusters2 variable

table.div <- StartracDiversity(CD8.clono, 
                       type = "group_id", 
                       sample = "sample_id", 
                       by = "overall", exportTable = TRUE)

table.div <- melt(table.div)
tran <- table.div %>% filter(variable == "tran")

#Check for normality 
shapiro.test(tran$value)#not normal
kruskal.test(value ~ majorCluster, data = tran)

p.tran <- ggplot(tran, aes(x=majorCluster, y=value)) +
  geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) + scale_fill_manual(values= pal_ident[c(4,1,3,5,2)]) +
  theme_classic() +
  ylab("Tran score") +
  guides(fill="none") +
  theme(axis.title.x = element_blank()) + ggtitle("Kruskal-Wallis chi-squared,\n p = 3.02e-05") +
  theme(plot.title = element_text(size = 5, face = "bold")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Check for normality
exp <- table.div %>% filter(variable == "expa")
shapiro.test(exp$value) #not normal
kruskal.test(value ~ majorCluster, data = exp)

p.exp <- ggplot(exp, aes(x=majorCluster, y=value)) +
  geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) + scale_fill_manual(values= pal_ident[c(4,1,3,5,2)]) +
  theme_classic() +
  ylab("Exp score") +
  guides(fill="none") +
  theme(axis.title.x = element_blank()) + ggtitle("Kruskal-Wallis chi-squared, \n p = 4.194e-05") +
  theme(plot.title = element_text(size = 5, face = "bold")) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tiff("../plots_CD8/startrac.tiff", width = 5*83, height = 5*50, res = 150, pointsize = 5)     
plot_grid(p.tran, p.exp, ncol = 2) 
dev.off()

DefaultAssay(CD8) <- "RNA"

#3) DGE including clusters SL1 and SL2
mark2 <- FindAllMarkers(CD8)
mark2 %>% dplyr::filter(!str_detect(rownames(mark2), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10 <- top10[!duplicated(top10$gene),]
df <- data.frame(top10$cluster, top10$gene)
df  <-  reshape(transform(df, indx = ave(as.character(top10.cluster), top10.cluster, FUN = seq)), 
                idvar = "indx", timevar = "top10.cluster", direction = "wide") 
colnames(df) <- gsub("top10.gene.", "", colnames(df))
all_markers <- df %>%
  select(-indx) 
colnames(all_markers) <- as.character(colnames(all_markers))
all_m <- lapply(all_markers, function(x) x[!is.na(x)])

#Prepare Heatmap of mean marker-exprs. by cluster
un_cm <- unlist(all_m) 
num_mark <- vapply(all_m, length, numeric(1))
rep_clus <- rep.int(names(all_m), num_mark)
labs <- sprintf("%s(%s)", un_cm, rep_clus)

# split cells by cluster
cel_by_clus <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

m_by_clus <- lapply(all_m, function(un_cm)
  vapply(cel_by_clus, function(i)
    Matrix::rowMeans(log.ct[un_cm, i, drop = FALSE]), 
    numeric(length(un_cm))))

# prep. for plotting 
mat <- do.call("rbind", m_by_clus)

#Z-score
mat <- t(scale(t(mat)))

mat <- mat %>% as.data.frame() %>%  select(c("Naive", "StL", "ActEx", "SL1", "SL2")) %>% 
  as.matrix()

cols <- pal_ident[seq_along(levels(CD8$clusters2))]
cols <- setNames(cols, levels(CD8$clusters2))
row_anno <- rowAnnotation(
  df = data.frame(cluster_id = c("Naive", "StL", "ActEx", "SL1", "SL2")),
  col = list(cluster_id = cols, gp = gpar(col = "white")), 
  show_legend = c(FALSE, TRUE)) 

idx <- which(rownames(mat) %in% c("KLRB1", "ZNF683", "GNLY", "FGFBP2"))
colLab <- rep("black", nrow(mat))
colLab[idx] <- "red"

lgd_aes <- list(direction = "horizontal", legend_width = unit(2.2, "cm"),
                title = "Expression")
h.heat.sl1sl2 <- Heatmap(t(mat),
                         cluster_rows = FALSE,
                         cluster_columns = TRUE,
                         row_names_side = "left",
                         left_annotation = row_anno,
                         row_names_gp = grid::gpar(fontsize = 10),
                         heatmap_legend_param = lgd_aes,
                         column_names_gp = gpar(col = colLab,fontsize = 7))

tiff("../plots_CD8/heatSL.tiff", width = 5*330, height = 5*150, res = 300, pointsize = 5)     
p <- draw(h.heat.sl1sl2, heatmap_legend_side = "bottom", align_heatmap_legend = "heatmap_center", 
          show_annotation_legend = FALSE)
p
dev.off()

tiff("../plots_CD8/vlnSL.tiff", width = 5*250, height = 5*100, res = 300, pointsize = 5)     
Stacked_VlnPlot(CD8, features = c("CX3CR1", "ZNF683"), colors_use = pal_ident)
dev.off()

#Trm markers
p <-FeaturePlot(CD8, features = c("CXCR6", "ITGAE"), combine=F, pt.size=0.00001, order=T) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position = "none")
}

tiff("../plots_CD8/p.Trm.umap.tiff", width = 5*300, height = 5*500, res = 300, pointsize = 5)     
p.markers <- cowplot::plot_grid(plotlist = p, nrow =2)
p.mar_leg <- plot_grid(p.markers, legend, ncol = 1, rel_heights = c(1, .1))
p.mar_leg
dev.off()

tiff("../plots_CD8/p.Trm.vln.tiff", width = 5*300, height = 5*500, res = 300, pointsize = 5)     
Stacked_VlnPlot(CD8, features = c("CD69", "CXCR6", "ITGAE"), colors_use = pal_ident)
dev.off()

#######################
##G) Velocity inference
#######################

#1) Load in loom files and prepare to analyze in Python via scVelo
ld.list <- list.files("../velocity", pattern = "loom", all.files = TRUE)
ld <- lapply(ld.list, function(x){
  ReadVelocity(paste0("../velocity/", x))
})
ld.list <- gsub(".loom", "_", ld.list)

#Rename colnames for each object according to the colnames of the integrated CD8 object
so <- lapply(ld, as.Seurat)
so <- lapply(so, function(x) RenameCells(x, new.names = gsub("x", "-1", colnames(x)), for.merge = T))
so <- lapply(so, function(x) RenameCells(x, new.names = gsub(".*:", "", colnames(x)), for.merge = T))
so <-  lapply(seq_along(ld.list), function(y) RenameCells(so[[y]], new.names = paste0(ld.list[[y]], colnames(so[[y]]))))

#Merge the objects
som <- merge(so[[1]], so[-(1)], merge.data = TRUE)

# Extract only the intersection
som <- som[, intersect(colnames(CD8), colnames(som[, colnames(CD8)]))] # This works

#Create the specific assays
spliced <- CreateAssayObject(GetAssayData(som, assay = "spliced"))
unspliced <- CreateAssayObject(GetAssayData(som, assay = "unspliced"))
ambiguous <- CreateAssayObject(GetAssayData(som, assay = "ambiguous"))

#Import the assays in the original CD8 object
CD8[["spliced"]] <- spliced
CD8[["unspliced"]] <- unspliced
CD8[["ambiguous"]] <- ambiguous

#Assign the spliced assay to the RNA assay and put it as default
CD8[['RNA']] <- CD8[["spliced"]]
DefaultAssay(CD8) <- "RNA"
CD8@active.ident <- CD8$clusters
DimPlot(CD8)
md <- CD8@meta.data
md$clusters <- as.vector(md$clusters)
md$clusters2 <- as.vector(md$clusters2)
CD8@meta.data <- md
CD8.velo <- CD8
SaveH5Seurat(CD8, filename = "CD8_veloscv.h5Seurat")
Convert("CD8_veloscv.h5Seurat", dest = "h5ad") ## --> continue in python

#2) import from python the csv with the info on latent time and different subsets and plot
df <- read.csv("l_time.csv")
df <- df %>% arrange(factor(clusters2, levels = c("Naive", "StL", "ActEx", "SL1", "SL2")))
df$clusters2 <- factor(df$clusters2, levels = c("Naive", "StL", "ActEx", "SL1", "SL2")) 

tiff("../plots_CD8/jitt.tiff", width = 5*100, height = 5*60, res = 150, pointsize = 5)     
ggplot(df, 
       aes(x = latent_time, y = clusters2, colour = clusters2)) + geom_jitter(size = 0.0000001) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_colour_manual(values = pal_ident[c(1,2,4,3,5)]) +
  theme_classic() +
  xlab("Latent Time") +
  ylab("") + labs(col = "") +theme(legend.text=element_text(size=10))
dev.off()

############
##H) SCENIC
############

#1) Download the required databases:
#I used the command-line interface downloaded (wget) the files in a folder called cisTarget_databases
# The files are:
#a)https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather;
#b) https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
# I also renamed the files as hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather and hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

#2) Run SCENIC workflow
#Proceed as here: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
exprMat <- CD8@assays$RNA@data #expresion matrix
cellInfo <- CD8@meta.data #metadata

levels(cellInfo$clusters2)

colVars <- list(CellType=c("Naive"="#F0E442",
                           "StL"="#E69F00",
                           "SL1"="#0072B2",
                           "SL2"='#56B4E9',
                           "ActEx"="#009E73"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$clusters2)]
myDatasetTitle <- "SCENIC AML" 
dbVersion <- 'v9'
dbDir <- 'cisTarget_databases'

#Scenic expects names(db@rankings)[1] to be "features" instead of motifs. To avoid errors and use hg38 instead of hg19 we need to add some more steps to SCENIC tutorial 
#After having modified the files we can write them and then read again in inizializeScenic
db <- importRankings("cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",  indexCol = "motifs")
names(db@rankings)[1] <- "features" 
db@org <- "hgnc"
db@genome <- "hg38"
arrow::write_feather(db@rankings, "cisTarget_databases/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

db <- importRankings("cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(db@rankings)[1] <- "features"
db@org <- "hgnc"
db@genome <- "hg38"
arrow::write_feather(db@rankings,"cisTarget_databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")

#assign the two files we are going to read to hg38Dbs
hg38Dbs <- c('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
             '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')

#Proceed as in the tutorial
scenicOptions <- initializeScenic(org=org, dbs = hg38Dbs, dbDir = dbDir,
                                   datasetTitle=myDatasetTitle, nCores=10)

scenicOptions@inputDatasetInfo$colVars <- colVars
scenicOptions@inputDatasetInfo$colVars <- cellInfo


scenicOptions@settings$dbs <- hg38Dbs
scenicOptions@settings$dbDir <- dbDir
scenicOptions@settings$db_mcVersion <- dbVersion

exprMat <- as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions2,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#Check if the interesting genes in my dataset are in genesKept
interestingGenes <- c("TCF7", "IL7R", "GZMK", "GNLY", "NR4A2", "ZNF683", "CX3CR1")
interestingGenes[which(!interestingGenes %in% genesKept)]

#Filter the matrix
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)

#3) Use exportsForArboreto to export matrix and scenicOptions and analyze in Python using GRNBoost 
#(this will markedly reduce the time to analyze the data)
exportsForArboreto(exprMat = exprMat_filtered, scenicOptions = scenicOptions, dir = ".")

#4) Read in GRNBoost output from Python and proceed with the workflow
GRNBoost_output <- read.delim("./adjacencies.tsv", header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight") #rename colnames to follow SCENIC R tutorial

saveRDS(GRNBoost_output, file="int/1.4_GENIE3_linkList.Rds") #save

exprMat_log <- log2(exprMat_filtered+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/SO.Rds") # save

nPcs <- c(5,15,50)

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")

saveRDS(cellInfo, file="int/cellInfo.Rds")
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))

par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="clusters2", cex=.5)

#Obtain matrix of regulon activity and plot as heatmap
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$clusters2),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
mat <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rownames(mat) <- gsub("_extended", "", rownames(mat))
mat <- mat %>% as.data.frame() %>%  select(c("Naive", "StL", "ActEx", "SL1", "SL2")) %>% as.matrix()
cols <- pal_ident[seq_along(levels(CD8$clusters2))]
cols <- setNames(cols, c("Naive", "StL", "SL1", "ActEx", "SL2"))
col_anno <- HeatmapAnnotation(
  df = data.frame(cluster_id = c("Naive", "StL", "ActEx", "SL1", "SL2")),
  col = list(cluster_id = cols, gp = gpar(col = "white"))) 
lgd_aes <- list(direction = "vertical", legend_width = unit(2.2, "cm"),
                title = "Regulon activity")
h.reg <- Heatmap(mat,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 row_names_side = "right",
                 bottom_annotation = col_anno,
                 row_names_gp = grid::gpar(fontsize = 7),
                 heatmap_legend_param = lgd_aes,
                 column_names_gp = gpar(fontsize = 15))

tiff("../heat.tiff", width = 5*300, height = 5*500, res = 300, pointsize = 5)     
h.reg
dev.off()

#Add the info obtained with SCENIC to Seurat metadata for additional plotting options
AUC.df <- t(AUCell::getAUC(regulonAUC)) %>% as.data.frame()
colnames(AUC.df) <- gsub("_extended|\\s", "", colnames(AUC.df))
md <- CD8@meta.data
AUC.ord <- AUC.df %>% slice(match(rownames(md), rownames(AUC.df)))
all(rownames(md) == rownames(AUC.ord))
md <- cbind(md, AUC.ord)
CD8@meta.data <- md
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, "HIF1A"), plots = "AUC")

FeaturePlot_scCustom(CD8, features = "HIF1A(50g)", colors_use =  viridis_plasma_light_high, order = TRUE, 
                     split.by = "group_id")



