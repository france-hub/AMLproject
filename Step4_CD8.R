##STEP4: CD8
#This script does the followings:
#Loads in the CD8 subset, find the most variable features, scale the data, run dimensionality reduction, run clustering
#Plots single genes and scores to identify different subsets
#Labels clusters and look at clusters abundancies

###In progress
#Trajectory inference
#single cell pathway analysis
#DGE along and across trajectories (TradeSeq)
#DGE between conditions (condiments)
#velocity inference
#clonotype analysis

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
library(CD8Disk)
library(viridis)
library(EnhancedVolcano)
library(scCustomize)
library(vip)
library(tidymodels)
library(scProportionTest)
library(fgsea)
library(dplyr)
library(limma)
library(ggplotify)
library(miloR)
library(rcna)
library(glue)
library(ggnetwork)
library(ggforce)
library(speckle)

setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")
getwd()
CD8 <-readRDS("CD8sub_test")

#Subclustering the integrated assay as suggested here https://github.com/satijalab/CD8/issues/2087 by timoast
# Identify the 2000 most variable genes in the RNA assay 
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

# Look at the batches
DimPlot(CD8, reduction = "pca", group.by = "batch")
DimPlot(CD8, reduction = "umap", group.by = "batch")

# Plot the elbow plot
ElbowPlot(object = CD8, ndims = 30)

# Determine the K-nearest neighbor graph
CD8 <- FindNeighbors(object = CD8, dims = 1:20, graph.name = c("RNA_nn", "integrated_nn"))

# Determine the clusters for various resolutions                                
CD8 <- FindClusters(object = CD8, resolution = c(0.4, 0.6, 0.8, 1.0))

# set cluster IDs to resolution 0.4 clustering
CD8 <- SetIdent(CD8, value = "integrated_snn_res.0.4")
DimPlot(CD8, label = T)
DefaultAssay(CD8) <- "RNA" #Put "RNA" as default assay
DimPlot(CD8, split.by = "group_id")

#Visualize the new subcluster
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
p.cd8 <- DimPlot_scCustom(CD8, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")


tiff("../plots/cd8clus.tiff", width = 5*400, height = 5*400, res = 300, pointsize = 5)     
p.cd8
dev.off()

###Visualize single markers and scores distributions
#Check they are CD8
FeaturePlot(CD8, "CD8A") 

#Check there are no CD4
FeaturePlot(CD8, "CD4") 

#DotPlot
features <- c("CCR7", "TCF7", "LEF1", "SELL", "IL7R", "SLAMF6", 
              "CXCR3", "GZMK", "CCL3", "CCL4", "CCL5", "XCL1", 
              "TIGIT", "PDCD1", "CD160", "LAG3", "CD69", "DUSP2", 
              "NR4A2", "ENTPD1", "ITGAE", "NKG7", "CX3CR1", "GNLY",
              "PRF1", "GZMB", "ZEB2", "CD226", "ZNF683",
              "TRAV1-2", "SLC4A10")   

pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

tiff("../plots/cd8dot.tiff", width = 5*400, height = 5*350, res = 300, pointsize = 5)     
Clustered_DotPlot(CD8, features, colors_use_exp = pal_exp, 
                  colors_use_idents = pal_ident) 
dev.off()

#Scores
naive <- c("TCF7", "CCR7", "SELL", "LEF1")
CD8 <- AddModuleScore(CD8, features = naive, name = "naive")

GI7 <- c("GZMK", "IL7R")
CD8 <- AddModuleScore(CD8, features = GI7, name = "GI7")

MAIT <- c("TRAV1-2", "SLC4A10")
MAIT <- AddModuleScore(CD8, features = MAIT, name = "MAIT")

#Scores
#signatures
sig <- readxl::read_xlsx("../signatures/sig.xlsx", sheet = 1)
sig_naive <- list(sig$`Naive (from Szabo et al. 2019)`[!is.na(sig$`Naive (from Szabo et al. 2019)`)])
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive")
p.naive <- FeaturePlot(CD8, "sig_naive1", pt.size = 0.00001, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("Naive (from Szabo et al. 2019)") + theme(plot.title = element_text(size = 15, face = "bold"))

sig_stem <- list(sig$`Stemness (from Pace et al. 2018)`[!is.na(sig$`Stemness (from Pace et al. 2018)`)])
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem")
p.stem <- FeaturePlot(CD8, "sig_stem1", pt.size = 0.000001, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("Stemness (from Pace et al. 2018)") + theme(plot.title = element_text(size = 15, face = "bold"))

#Bulk senescence signature
genes <- readxl::read_xlsx("../signatures/Model_5_results.xlsx", sheet = 2)
sen <- list(genes$gene_name[1:100])
CD8 <- AddModuleScore(CD8, features = sen, name = "sen")
p.sen <- FeaturePlot(CD8, features = "sen1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position="none") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("Senescence score") + theme(plot.title = element_text(size = 15, face = "bold"))

#Signature dysfunction June et al.
dys.genes <- readxl::read_xlsx("../signatures/Dysfunction.xlsx")
dys <- list(dys.genes$Dysfunction)
CD8 <- AddModuleScore(CD8, features = dys, name = "dys")
tiff("./plots/sen_group.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
p.dys <- FeaturePlot(CD8, features = "dys1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position="none") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("Dysfunction score (Good et al. 2021)")

tiff("../plots/scores.tiff", width = 5*550, height = 5*450, res = 300, pointsize = 5)     
cowplot::plot_grid(p.naive, p.stem, p.sen, p.dys)
dev.off()

md <- CD8@meta.data %>% rename(sen1 = "Senescence", dys1 = "Dysfunction")
colnames(md)
CD8@meta.data <- md

tiff("../plots/stackvln.tiff", width = 5*400, height = 5*200, res = 300, pointsize = 5)     
Stacked_VlnPlot(CD8, features = c("Senescence", "Dysfunction"), plot_legend = TRUE) 
dev.off()

p.sc <- FeatureScatter(CD8, feature1 = "Dysfunction", feature2 = "Senescence", pt.size = 0.1)

#Calculate principal component to scale colors
data <- p.sc$data
data$pc <- predict(prcomp(~Senescence+Dysfunction, data))[,1]

#Nicer picture
tiff("../plots/corr.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
ggplot(data, aes(x = Senescence, y = Dysfunction, color = pc)) +
  geom_point(shape = 16, size = 3, alpha = .4) +
  geom_smooth(method = "lm", color = "black", size=0.3, se = TRUE) +
  ggpubr::stat_cor() +
  theme_classic() +  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  ylab("Dysfunction score") + xlab("Senescence score") + theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Plots for figure 2
gene.list <- c("CCR7", "TCF7", "IL7R", "GZMK", "CD69", "GZMB", "GNLY", "CX3CR1", "CD226")

p <-FeaturePlot(CD8, features = c("CCR7", "TCF7", "IL7R", "GZMK", "CD69", "GZMB", "GNLY", "CX3CR1", "CD226"), combine=F, pt.size=0.001, order=T) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))
  
}
tiff("../plots/p.markers.umap.tiff", width = 5*400, height = 5*450, res = 300, pointsize = 5)     
cowplot::plot_grid(plotlist = p, nrow =3)
dev.off()

#DGE all markers
mark <- FindAllMarkers(CD8)

mark %>% filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
df <- data.frame(top10$cluster, top10$gene)
df  <-  reshape(transform(df, indx = ave(as.character(top10.cluster), top10.cluster, FUN = seq)), 
                idvar = "indx", timevar = "top10.cluster", direction = "wide") 
colnames(df) <- gsub("top10.gene.", "", colnames(df))

all_markers <- df %>%
  select(-indx) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

#Use same colors used for UMAP
tiff("../plots/p.dotDGE_noAnnot.tiff", width = 5*400, height = 5*600, res = 300, pointsize = 5)     
Clustered_DotPlot(CD8, all_markers, colors_use_exp = pal_exp,
                  exp_color_min = -1, exp_color_max = 1, colors_use_idents = pal_ident,
                  x_lab_rotate = TRUE, k =4)
dev.off()

DimPlot(CD8)
#Clean workspace
rm(list=setdiff(ls(), "CD8"))
p.umap <- DimPlot(CD8, label = T)
cowplot::plot_grid(p.umap, p.vln)

CD8$NotAnnClus <- CD8@active.ident

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
                    "12" = "SL")

#Add clusters variable to metadata
CD8$clusters<- CD8@active.ident

tiff("../plots/UMAP_ann.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p.ann <- DimPlot_scCustom(CD8, label = TRUE, figure_plot = TRUE, colors_use = pal_ident[1:4], pt.size = 0.00001)
p.ann
dev.off()

tiff("../plots/UMAP_groupId.tiff", width = 5*600, height = 5*250, res = 300, pointsize = 5)     
p1 <- DimPlot_scCustom(CD8, label = TRUE, split.by = "group_id", colors_use = pal_ident)
p1
dev.off()

tiff("../plots/UMAP_density.tiff", width = 5*600, height = 5*400, res = 300, pointsize = 5)     
p.dens <- DimPlot(CD8, reduction = 'umap', split.by = "RespTmp")  + NoLegend() + NoAxes() 
p.dens +  geom_density_2d_filled(p.dens$data, mapping = aes(x = p.dens$data[,"UMAP_1"], y = p.dens$data[,"UMAP_2"]), contour_var = "ndensity") + 
  facet_wrap(vars(RespTmp))
dev.off()

#Plot network
snn <- CD8@graphs$integrated_snn %>%
  as.matrix() %>%
  ggnetwork() %>%
  left_join(CD8@meta.data %>% mutate(vertex.names = rownames(.)), by = 'vertex.names')

p1 <- ggplot(snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = group_id), size = 0.5) +
  scale_color_manual(
    name = 'group_id', values = pal_ident,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'right') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(CD8@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )

p2 <- ggplot(snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = clusters), size = 0.5) +
  scale_color_manual(
    name = 'clusters', values = pal_ident,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'right') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(CD8@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  )
tiff("../plots/umap_snn.tiff", width = 5*400, height = 5*200, res = 150, pointsize = 5)     
p1 + p2 + plot_layout(ncol = 2)
dev.off()

#Look for antigen specific
#neoTCR
DefaultAssay(CD8) <- "RNA"
neoT <- readxl::read_xlsx("../signatures/neoTCR.xlsx", sheet = 2)
neoT <- list(neoT$Gene)
CD8 <- AddModuleScore(CD8, features = neoT, name = "neoT", random.seed = 1)
neoT <- FeaturePlot(CD8, features = "neoT1", pt.size = 0.0000001, order = T) 
tiff("../plots/neoT.tiff", width = 5*250, height = 5*200, res = 300, pointsize = 5)     
neoT + scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + theme(legend.position="none") +
  ggtitle("neoTCR8")
dev.off()

v <- VlnPlot(CD8, features = "neoT1")

tiff("../plots/boxNeoTCR8.tiff", width = 5*70, height = 5*80, res = 150, pointsize = 5)     
ggplot(v$data, aes(x=ident, y=neoT1, fill=ident)) +
  labs(x = NULL, y = NULL) +
  geom_boxplot() +  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height  =  unit(0.8, "lines")) + theme_classic() + ggtitle("NeoTCR8") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Relative abundance of the each clusters per condition (How to check for "significance"? What test should we apply?)
freq_table <- table(CD8@active.ident, CD8$group_id)
fqs <- prop.table(table(CD8@active.ident, CD8$sample_id), 2)
df <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df$group_id <- CD8$group_id[match(df$sample_id, CD8$sample_id)]
View(df)
tiff("../plots/box.tiff", width = 5*250, height = 5*40, res = 150, pointsize = 5)     
p.box <- ggplot(df, aes(x = group_id, y = frequency, color = group_id)) +
  labs(x = NULL, y = "Proportion [%]") +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height  =  unit(0.8, "lines")) +
  geom_boxplot(aes_string(color = "group_id", fill = "group_id"), position = position_dodge(), alpha = 0.2, 
               outlier.color = NA, show.legend = FALSE) + 
  geom_point(aes_string(x = "group_id", col = "group_id"), position = position_jitter(width = 0.2)) +
  facet_wrap(~ cluster_id, scales = "free_y", ncol = 6) +
  theme_classic()
p.box + ggtitle("Differential abundance")  + theme(plot.title = element_text(hjust = 0.5))

#Differential abundance testing
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

p.perm <- plot_grid(p.perm1, p.perm2, p.perm3, rel_widths = c(1,1,1.9), nrow = 1)

tiff("../plots/sig.tiff", width = 5*300, height = 5*60, res = 150, pointsize = 5)     
p.perm
dev.off()

#Speckle
props <- getTransformedProps(CD8$clusters, CD8$sample_id,
                             transform="logit")

tiff("../plots/MeanVar.tiff", width = 5*150, height = 5*150, res = 300, pointsize = 5)     
plotCellTypePropsMeanVar(props$Counts)
dev.off()

condition <- c(rep("HD",2), rep("NonRes",2), rep("Res", 4), rep("NonRes", 2), 
               rep("Res", 2))
pair <- rep(c(1,2,3,4,5,6), each = 2)
batch <- c(rep(1,2), rep(2,2), rep(3,2), rep(2,2), rep(1,2), rep(2,2))
design <- model.matrix(~0 + condition + pair + batch)
df.anova <- propeller.anova(prop.list=props, design=design, coef = c(1,2,3), 
                            robust=TRUE, trend=FALSE, sort=TRUE) %>% as.data.frame() 

df.anova <- df.anova %>% set_colnames(c("HD", "NonRes", "Res", "Fstat", "FDR"))
fdr = 0.05
s <- factor(
  ifelse(df.anova$FDR < fdr, "yes", "no"), 
  levels = c("no", "yes"))
fdr_pal <- c("blue", "red3")
names(fdr_pal) <- levels(s)
F_pal <- c("blue",  "red3")
F_lims <- range(df.anova$Fstat)
F_brks <- c(0, f_lims[2])
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

tiff("../plots/boxAnova.tiff", width = 5*200, height = 5*350, res = 300, pointsize = 5)     
p.box <- ggplot(df.plot, aes(x = group_id, y = frequency, color = group_id)) +
  labs(x = NULL, y = "Logit Transformed Proportions") +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height  =  unit(0.8, "lines")) +
  geom_boxplot(aes_string(color = "group_id", fill = "group_id"), position = position_dodge(), alpha = 0.2, 
               outlier.color = NA, show.legend = FALSE) + 
  geom_point(aes_string(x = "group_id", col = "group_id"), position = position_jitter(width = 0.2)) +
  facet_wrap(~ cluster_id, scales = "free_y", ncol = 1) +
  theme_classic()
p.box + ggtitle("Differential abundance")  + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#CNA
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
  labs(title = 'Disease association (filtered)', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') 
p3 <- DimPlot(CD8.cna, group.by = "group_id") + ggtitle("Condition")

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

genes <- readxl::read_xlsx("../signatures/Model_5_results.xlsx", sheet = 2)
sen <- list(genes$gene_name[1:100])
CD8.cna <- AddModuleScore(CD8.cna, features = sen, name = "sen")
p.sen.cna <- FeaturePlot(CD8.cna, features = "sen1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position="none") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("Senescence score") + theme(plot.title = element_text(size = 15, face = "bold"))

tiff("../plots/cna.tiff", width = 5*800, height = 5*400, res = 300, pointsize = 5)     
plot_grid(p1,p2,p3, p.naive.cna, p.stem.cna, p.sen.cna)
dev.off()

#Look at local densities
use_condaenv("~/miniconda3/")
# make sure reticulate package can find correct version of "umap" model of python
umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
umap_import$pkg_resources$get_distribution("umap-learn")$version # here should print 0.5.2 in your case
CD8 <- RunUMAP(CD8, dims = 1:30, umap.method = "umap-learn", densmap=TRUE)
tiff("../plots/densMAP.tiff", width = 5*300, height = 5*150, res = 150, pointsize = 5)     
DimPlot(CD8, split.by = "group_id")
dev.off()

#DGE all markers
mark <- FindAllMarkers(CD8)

#signatures for TCGA
#top20
mark %>% filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
df20 <- data.frame(top20$cluster, top20$gene)
df20  <-  reshape(transform(df20, indx = ave(as.character(top20.cluster), top20.cluster, FUN = seq)), 
                idvar = "indx", timevar = "top20.cluster", direction = "wide") %>% select(-indx)
colnames(df20) <- gsub("top20.gene.", "", colnames(df))

#top50
mark %>% filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
df50 <- data.frame(top50$cluster, top50$gene)
df50  <-  reshape(transform(df50, indx = ave(as.character(top50.cluster), top50.cluster, FUN = seq)), 
                idvar = "indx", timevar = "top50.cluster", direction = "wide") %>% select(-indx)
colnames(df50) <- gsub("top50.gene.", "", colnames(df))
writexl::write_xlsx(df20, path="top20.xlsx")
writexl::write_xlsx(df50, path="top50.xlsx")

mark %>% filter(!str_detect(rownames(mark), "^RP[SL]")) %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
df <- data.frame(top10$cluster, top10$gene)
df  <-  reshape(transform(df, indx = ave(as.character(top10.cluster), top10.cluster, FUN = seq)), 
              idvar = "indx", timevar = "top10.cluster", direction = "wide") 
colnames(df) <- gsub("top10.gene.", "", colnames(df))

all_markers <- df %>%
  select(-indx) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

#Use same colors used for UMAP
tiff("../plots/p.dotDGE.tiff", width = 5*280, height = 5*300, res = 300, pointsize = 5)     
Clustered_DotPlot(CD8, all_markers, colors_use_exp = pal_exp,
                  exp_color_min = -1, exp_color_max = 1, colors_use_idents = pal_ident,
                  x_lab_rotate = TRUE, k =4)
dev.off()

#by condition
tiff("../plots/p.dot_cond.tiff", width = 5*280, height = 5*300, res = 300, pointsize = 5)     
Clustered_DotPlot(CD8, all_markers, colors_use_exp = pal_exp, group.by = "group_id",
                  exp_color_min = -1, exp_color_max = 1, colors_use_idents = pal_ident,
                  x_lab_rotate = TRUE, k =2)
dev.off()
DefaultAssay(CD8) <- "RNA"
#Custom Markers
clust_mark <- list(
  Naive = c("CCR7", "MAL", "LEF1", "SELL", "TCF7"),
  StL = c("BCL2", "BACH2", "CD27","IL7R", "SLAMF6", "CXCR3", "GZMK",  "ITGAE","ENTPD1"),
  ActEx = c("EOMES", "CCL4", "XCL2", "CCL3", "XCL1", "KLF6", "TIGIT", "CD69", "CD160", "PDCD1", "TOX", "NR4A2",  "DUSP2"),
  SL = c("NKG7", "TBX21",  "CD38","FCRL6", "FCGR3A", "C1orf21", "PRF1", "ENO1",
                     "GNLY", "ZEB2", "CX3CR1", "FGFBP2", "KLRD1", "GZMB","ZNF683", "CD226")
)
#Prepare Heatmap of mean marker-exprs. by cluster
un_cm <- unlist(clust_mark)
num_mark <- vapply(clust_mark, length, numeric(1))
rep_clus <- rep.int(names(clust_mark), num_mark)
labs <- sprintf("%s(%s)", un_cm, rep_clus)

# split cells by cluster
cel_by_clus <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

m_by_clus <- lapply(clust_mark, function(un_cm)
  vapply(cel_by_clus, function(i)
    Matrix::rowMeans(log.ct[un_cm, i, drop = FALSE]), 
    numeric(length(un_cm))))

# prep. for plotting 
mat <- do.call("rbind", m_by_clus)

#Z-score
mat <- t(scale(t(mat)))

mat <- mat %>% as.data.frame() %>%  select(names(clust_mark)) %>% as.matrix()

h.heat <- Heatmap(mat,
             name = "Z-score",
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_names_side = "left",
             row_names_gp = grid::gpar(fontsize = 6))

tiff("../plots/p.heat_custom.tiff", width = 5*150, height = 5*250, res = 300, pointsize = 5)     
h.heat
dev.off()

#Correlation Heatmap
custom_colors <- c(pal_exp)

# calculate average expression value for all variable genes for each cluster
average_expression_profiles_by_cluster <- CD8@assays$RNA@data[CD8@assays$RNA@var.features,] %>%
  t() %>%
  as.matrix() %>%
  as_tibble() %>%
  mutate(cluster = CD8@meta.data$clusters) %>%
  select(cluster, everything()) %>%
  group_by(cluster) %>%
  summarize_all(~mean(.))

# calculate Spearman correlation matrix
correlation_matrix <- average_expression_profiles_by_cluster %>%
  select(-1) %>%
  as.matrix() %>%
  t() %>%
  cor(method = 'spearman')

# assign row and column names
rownames(correlation_matrix) <- levels(CD8@meta.data$clusters)
colnames(correlation_matrix) <- levels(CD8@meta.data$clusters)

# save cluster names for later
cluster <- rownames(correlation_matrix)

# assign a color to each cluster
colors_for_clusters <- c(pal_ident[1:length(cluster)])
names(colors_for_clusters) <- cluster

# create annotation function
func_cell_cluster <- function(i, j, x, y, width, height, fill) {
  grid.text(cluster[j], x = x, y = y, gp = gpar(fontsize = 8))
}

# create main heatmap
ht_matrix <- Heatmap(
  correlation_matrix,
  name = 'Spearman\ncorrelation',
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    title = 'Spearman correlation',
    legend_height = unit(2, 'cm'),
    legend_width = unit(1, 'cm')
  )
)

tiff("../plots/p.corr.tiff", width = 5*250, height = 5*250, res = 300, pointsize = 5)     
ht_matrix 
dev.off()

#Coexp gene matrix
mat_corr <- cor(t(as.matrix(CD8@assays$RNA@data[VariableFeatures(CD8)[1:500],])),
                method = "spearman")
tiff("../plots/h.coexp.tiff", width = 5*750, height = 5*750, res = 300, pointsize = 5)     
Heatmap(mat_corr, name = "Spearman correlation",
        column_names_gp = grid::gpar(fontsize = 2),
        row_names_gp = grid::gpar(fontsize = 2))
dev.off()

#Identify cell cycling cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

CD8cc <- CellCycleScoring(CD8, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

tiff("../plots/ccUMAP.tiff", width = 5*200, height = 5*150, res = 300, pointsize = 5)     
DimPlot(CD8cc,
        reduction = "umap",
        group.by= "Phase")
dev.off()

freq_table <- CD8cc[[]]
colnames(freq_table)
freq_table <- freq_table[,c("group_id", "clusters", "Phase")]
freq_table <- subset(freq_table, Phase != "Undecided") #removing undecided phases
freq_table <- freq_table %>%
  group_by(group_id, clusters, Phase) %>%
  summarise(n = n())
freq_table$Phase <- factor(freq_table$Phase, levels = c("G1", "G2M", "S")) #ordering phases

tiff("../plots/ccAbund.tiff", width = 5*150, height = 5*100, res = 300, pointsize = 5)     
ggplot(freq_table, aes(x=clusters, y=n, fill=Phase)) + 
  labs(x = NULL, y = NULL)+
  stat_summary(geom="bar", position="fill", color="black", lwd=0.25) + 
  theme(axis.title.x = element_blank()) + 
  theme_classic() + coord_flip()
dev.off()

#Condition DE within cluster(Try also https://github.com/egarren/scTfh/blob/main/code/01_gex.R)
CD8$celltype.group <- paste(CD8$clusters, CD8$group_id, sep = "_")
Idents(CD8) <- "celltype.group"
clusters <- levels(CD8$clusters)
df <- list()
for (i in clusters){
  Idents(CD8) <- "celltype.group" #setting idents to new metadata column
  df[[i]] <- FindMarkers(CD8, ident.1 = paste0(i,"_Res"), ident.2 = paste0(i,"_NonRes"), 
                    min.pct=0,logfc.threshold = -Inf, base =2)
}

p.volc <- lapply(df, function(x)
  EnhancedVolcano(x,
                  lab = rownames(x),
                  x = "avg_log2FC",
                  y = 'p_val_adj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  labSize = 6.0,
                  colAlpha = 1,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75))


#Drivers of SL and Tex
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

set.seed(123)
View(data)
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
  update_role(cell_barcode, new_role = "ID")  %>%
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


#TRAJECTORY INFERENCE using Slingshot
clusterLabels <- CD8@active.ident
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
ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes_string(color = df$Lineage1_pst),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  theme_minimal() + labs(colour = "Pseudotime") 

#Lineage 2
ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes_string(color = df$Lineage2_pst),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  theme_minimal() + labs(colour = "Pseudotime") 

#Together
p.traj <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = -5, y = 5.8, label = "Activation-Exhaustion", size = 5) +
  annotate("text", x = -7, y = -5.6, label = "Senescence", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +  theme_minimal()  

tiff("../plots/umap_traj.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
p.traj
dev.off()

## Identifying differentially expressed genes along a trajectory
#Fit negative binomial model
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

#fitGAM
sce.gam <- fitGAM(counts = counts(sceCD8), sds = sds, nknots = 5, verbose = TRUE, BPPARAM = BPPARAM)

# plot our Slingshot lineage trajectories, this time illustrating the new tradeSeq knots
tiff("./plots/traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plotGeneCount(curve = sds, counts = counts,
              clusters = CD8@active.ident,
              models = sce.gam)
dev.off()

#Association test
assoRes <- associationTest(sce.gam)
head(assoRes)

### Discovering differentiated cell type markers
# discover marker genes for the differentiated cell types
endRes <- diffEndTest(sce.gam) #Nothing interesting with this analysis
head(endRes)

o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce.gam)[o[2]]
plotSmoothers(sceCD8, counts(sceCD8), sigGene) 

plotGeneCount(sds, counts(sce.gam), gene = sigGene)

# Marker genes between specific roots   
earlyDERes <- earlyDETest(sce.gam, knots = c(3, 4))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(earlyDERes)[oEarly][1])
plotGeneCount(sds, counts(sce.gam), gene = rownames(earlyDERes)[oEarly][3])

#Genes with different expression patterns (most interesting part)
patternRes <- patternTest(sce.gam)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
rownames(patternRes)[oPat][1:10]

p.leg <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][1]) + 
  scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.text = element_text(size=15),
        legend.title = element_text(size=16))

legend <- cowplot::get_legend(p.leg + theme(legend.position = "right"))

p2 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][1], 
                    xlab = "",
                    ylab = "") + 
  ggtitle ("NKG7") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p3 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][2],
                    xlab = "",
                    ylab = "") + 
  ggtitle ("GZMK") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p4 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][6],
                    xlab = '', ylab = '') + 
  ggtitle ("PRF1") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p5 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][7],
                    xlab = '',
                    ylab = '') + 
  ggtitle ("GNLY") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())


p6 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][8],
                    xlab = '',
                    ylab = '') + 
  ggtitle ("GZMB") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

#NR4A2
p7 <- plotSmoothers(sce.gam, counts(sce.gam), gene = rownames(patternRes)[oPat][68],
                    xlab = '',
                    ylab = '') + 
  ggtitle ("NR4A2") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_color_viridis_d(name = "Lineages", labels=c("1" = "Senescence", "2"="Exhaustion")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p.patt <- plot_grid(p2, p4, p5, p6, p3, p7) + theme_void()
legend <- cowplot::get_legend(p.leg)
p.patt.stack <- plot_grid(p.traj, NULL, legend, rel_widths = c(1, -0.2, 1), nrow = 1)

ptrajcoord <- ggplot(data.frame(x = 100, y = 100), aes(x = x, y = y)) +
  geom_point() +
  xlim(c(0, 10)) + ylim(c(0,10)) +
  theme_classic() +
  ylab("Log(expression + 1") + xlab("Pseudotime") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(
          arrow = arrow(angle = 15, length = unit(0.5, "cm"), type = "closed")))

layout <- c(
  area(t = 1, l = 3, b = 11, r = 11),
  area(t = 10, l = 2, b = 12, r = 2)
)
plot(layout)

p.patt <- p.patt.stack + ptrajcoord + plot_layout(design = layout)

p.umapcoord <- ggplot(data.frame(x = 100, y = 100), aes(x = x, y = y)) +
  geom_point() +
  xlim(c(0, 10)) + ylim(c(0,10)) +
  theme_classic() +
  ylab("UMAP_2") + xlab("UMAP_1") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(
          arrow = arrow(angle = 15, length = unit(0.5, "cm"), type = "closed")))

layout <- c(
  area(t = 1, l = 3, b = 11, r = 11),
  area(t = 10, l = 2, b = 12, r = 2)
)

p.traj <- p1 + p.umapcoord + plot_layout(design = layout)

tiff("./plots/traj.tiff", width = 5*1400, height = 5*500, res = 300, pointsize = 5)     
plot_grid(p.traj, p.patt, rel_widths = c(1,2), nrow = 1)
dev.off()
ls()

#Different SL like clusters
CD8 <- SetIdent(CD8, value = "integrated_snn_res.0.4")
p.umap2 <- DimPlot(CD8, label = T)

freq_table <- table(CD8@active.ident, CD8$group_id)
fqs <- prop.table(table(CD8@active.ident, CD8$sample_id), 2)
df <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df$group_id <- CD8$group_id[match(df$sample_id, CD8$sample_id)]

p.box_2 <- ggplot(df, aes(x = group_id, y = frequency, color = group_id)) +
  labs(x = NULL, y = "Proportion [%]") +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height  =  unit(0.8, "lines")) +
  geom_boxplot(aes_string(color = "group_id", fill = "group_id"), position = position_dodge(), alpha = 0.2, 
               outlier.color = NA, show.legend = FALSE) + 
  geom_point(aes_string(x = "group_id", col = "group_id"), position = position_jitter(width = 0.2)) +
  facet_wrap(~ cluster_id, scales = "free_y", ncol = 6) +
  theme_classic()

tiff("../plots/box_sen.tiff", width = 5*600, height = 5*250, res = 300, pointsize = 5)     
p.box_2 + ggtitle("")  + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Better look at senescent cells
CD8.sen <- RenameIdents(CD8, "0" = "Naive",
                        "1" = "StL",
                        "2" = "SL2",
                        "3" = "ActEx",
                        "4" = "SL1",
                        "5" = "Naive",
                        "6" = "SL1",
                        "7" = "SL2",
                        "8" = "Naive",
                        "9" = "SL1",
                        "10" = "StL",
                        "11" = "StL",
                        "12" = "SL2")
tiff("../plots/umap_sen.tiff", width = 5*250, height = 5*300, res = 300, pointsize = 5)     
DimPlot_scCustom(CD8.sen, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")
dev.off()

levels(CD8.sen@active.ident) <- c("Naive", "StL", "ActEx",  "SL1", "SL2")
CD8.sen$clusters <- CD8.sen@active.ident
cluster_stats <- Cluster_Stats_All_Samples(CD8.sen, group_by_var = "group_id") %>%  
  mutate(freq_groupID = log10(Res/NonRes)) %>%  
  filter(!str_detect("Total", Cluster))

cluster_stats$Cluster <-  factor(cluster_stats$Cluster, levels = c("Naive", "StL", "ActEx", "SL1", "SL2"))

colScale <- scale_colour_manual(name = "Cluster",values = pal_ident[1:5])

tiff("../plots/dot_group.tiff", width = 5*80, height = 5*60, res = 150, pointsize = 5)     
cluster_stats %>% ggplot(aes(x=Cluster, y = freq_groupID, size = abs(freq_groupID), colour = Cluster)) + 
  geom_point() + theme_classic() + ylim(-1.5, 1.5) + geom_hline(yintercept=0,linetype=2) +
  theme(legend.position="none") + ylab("log10(Res/NonRes)") + colScale
dev.off()

#Find markers of SL1 and SL2
SL <- subset(CD8.sen, group_id %in% c("Res", "NonRes") & clusters %in%
               c("SL1", "SL2"))
SL<- SL %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.1, verbose = FALSE) %>%
  RunUMAP(dims = 1:10, verbose = FALSE)

gene.list <- c("CX3CR1", "ZNF683", "TNFAIP3")

p <-FeaturePlot(CD8, features = gene.list, combine=F, pt.size=0.001, order=T) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))
  
}
tiff("../plots/p.markers_sen.umap.tiff", width = 5*500, height = 5*150, res = 300, pointsize = 5)     
cowplot::plot_grid(plotlist = p, nrow =1)
dev.off()

tiff("../plots/p.markers_sen.vln.tiff", width = 5*500, height = 5*150, res = 300, pointsize = 5)     
VlnPlot_scCustom(CD8.sen, features = gene.list, colors_use = pal_ident, pt.size = 0)
dev.off()

rm(list = ls(pattern = "CD8"))
data <- SL@assays$RNA@scale.data

# let's transpose the matrix and make it to a dataframe
dim(data)
data <- t(data) %>% as.data.frame()

## add the cell type/the outcome/y to the dataframe
data$cell_type<- SL$clusters
data$cell_barcode<- rownames(data)
## it is important to turn it to a factor for classification
data$cell_type<- factor(data$cell_type)

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
  update_role(cell_barcode, new_role = "ID")  %>%
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

tiff("../plots/lasso_genes.tiff", width = 5*150, height = 5*400, res =300, pointsize = 5)
scCustomize::Stacked_VlnPlot(SL, features = lasso_features %>%
                               pull(term) %>% head(n = 20),
                             colors_use = c("blue", "red"))
dev.off()

elastic_recipe <-
  recipe(formula = cell_type ~ ., data = data_train) %>%
  update_role(cell_barcode, new_role = "ID") %>%
  step_zv(all_predictors())

# we will tune both penalty and the mixture
elastic_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

elastic_workflow <- workflow() %>%
  add_recipe(elastic_recipe) %>%
  add_model(elastic_spec)

## tune both the penalty and the mixture from 0-1
penalty_grid <- grid_regular(penalty(range = c(-2, 2)), mixture(), levels = 50)

doParallel::registerDoParallel()
tune_res <- tune_grid(
  elastic_workflow,
  resamples = data_fold,
  grid = penalty_grid
)

best_penalty <- select_best(tune_res, metric = "accuracy")
elastic_final <- finalize_workflow(elastic_workflow, best_penalty)
elastic_final_fit <- fit(elastic_final, data = data_train)

## confusion matrix
predict(elastic_final_fit, new_data = data_test) %>%
  bind_cols(data_test %>% select(cell_type)) %>%
  conf_mat(truth = cell_type, estimate = .pred_class)

elastic_features<- tidy(elastic_final_fit) %>%
  arrange(desc(abs(estimate))) %>%
  filter(estimate != 0) 

merged_markers<- left_join(elastic_features, lasso_features, by = c("term" = "term")) %>%
  dplyr::rename(estimate.elastic = estimate.x, estimate.lasso= estimate.y) %>%
  select(-penalty.x, - penalty.y) 

Idents(pbmc_subset)<- pbmc_subset$cell_type
tiff("lasso_elasto.tiff", width = 5*250, height = 5*300, res = 300, pointsize = 5)     
scCustomize::Stacked_VlnPlot(SL, features = merged_markers %>% slice(-1) %>% pull(term) %>% head(n = 10),
                             colors_use = c("blue", "red") )

dev.off()
#Random Forest
rf_recipe <- 
  recipe(formula = cell_type ~ ., data = data_train) %>%
  update_role(cell_barcode, new_role = "ID") %>%
  step_zv(all_predictors())

## feature importance sore to TRUE
rf_spec <- rand_forest() %>%
  set_engine("randomForest", importance = TRUE) %>%
  set_mode("classification")

rf_workflow <- workflow() %>% 
  add_recipe(rf_recipe) %>% 
  add_model(rf_spec)

rf_fit <- fit(rf_workflow, data = data_train)

## confusion matrix, perfect classification! 
predict(rf_fit, new_data = data_test) %>%
  bind_cols(data_test %>% select(cell_type)) %>%
  conf_mat(truth = cell_type, estimate = .pred_class)

rf_fit %>%
  extract_fit_parsnip() %>%
  vip::vip(geom = "col", num_features = 25) + 
  theme_bw(base_size = 14)+
  labs(title = "Random forest variable importance") 

rf_fit %>%
  extract_fit_parsnip() %>%
  vip::vi_model() %>%
  arrange(desc(abs(Importance))) %>%
  head(n = 20)

rf_features<- rf_fit %>%
  extract_fit_parsnip() %>%
  vip::vi_model() %>%
  arrange(desc(abs(Importance))) %>%
  head(n = 20) %>%
  pull(Variable)

tiff("rf.tiff", width = 5*250, height = 5*500, res = 300, pointsize = 5)     
scCustomize::Stacked_VlnPlot(SL, features = rf_features,
                             colors_use = c("blue", "red")) 
dev.off()

#GSVA
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
  filter(grepl(clusters, pattern = 'SL|StL')) %>%
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
    rep(0, CD8@meta.data %>% filter(clusters == 'StL') %>% nrow()),
    rep(1, CD8@meta.data %>% filter(clusters == 'SL') %>% nrow())
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
    nudge_y = data$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data$t), max(data$t))) +
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
  filter(grepl(clusters, pattern = 'SL|ActEx')) %>%
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
    rep(0, CD8@meta.data %>% filter(clusters == 'ActEx') %>% nrow()),
    rep(1, CD8@meta.data %>% filter(clusters == 'SL') %>% nrow())
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
    nudge_y = data$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data$t), max(data$t))) +
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
  filter(grepl(clusters, pattern = 'ActEx|StL')) %>%
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
    rep(0, CD8@meta.data %>% filter(clusters == 'StL') %>% nrow()),
    rep(1, CD8@meta.data %>% filter(clusters == 'ActEx') %>% nrow())
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
    nudge_y = data$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data$t), max(data$t))) +
  scale_color_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y =  element_blank(),
    legend.position = 'none'
  ) + ggtitle("Activated-exhausted VS Stem-like") +
  theme(plot.title = element_text(hjust = 0.5)) 

##Senescence-like2 VS Senescence-like1
# get indices of cells which are either SL2 or SL1
# get indices of cells which are either SL or StL
c_SL2vSL1 <- CD8.sen@meta.data %>%
  mutate(row_number = row_number()) %>%
  filter(grepl(clusters, pattern = 'SL2|SL1')) %>%
  arrange(clusters) %>%
  pull(row_number)

# get expression matrix and reduce it to cells and genes of interest
mat_SL2vSL1<- CD8.sen@assays$RNA@counts[ genes_to_analyze , c_SL2vSL1] %>% as.matrix()

# perform GSVA
g_SL2vSL1 <- GSVA::gsva(
  mat_SL2vSL1,
  gset.idx.list = gene_sets$genesets,
  parallel.sz = 1
)

# generate design matrix
dm_SL2vSL1 <- tibble(
  control = 1,
  test = c(
    rep(0, CD8.sen@meta.data %>% filter(clusters == 'SL1') %>% nrow()),
    rep(1, CD8.sen@meta.data %>% filter(clusters == 'SL2') %>% nrow())
  )
)

# fit linear model, followed by empirical Bayes statistics for differential
# enrichment analysis
fit_SL2vSL1<- lmFit(g_SL2vSL1, dm_SL2vSL1)
fit_SL2vSL1<- eBayes(fit_SL2vSL1)

# prepare data for plotting
data_SL2vSL1 <- topTable(fit_SL2vSL1, coef = 'test', number = 50) %>%
  mutate(gene_set = rownames(fit_SL2vSL1$t)) %>%
  arrange(t) %>%
  mutate(
    gene_set = factor(gene_set, levels = gene_set),
    just = ifelse(t < 0, 0, 1),
    nudge_y = ifelse(t < 0, 1, -1),
    color = ifelse(t < -5 | t > 5, 'black', 'grey')
  )

# plot t-value
SL2vSL1 <- ggplot(data = data_SL2vSL1, aes(x = gene_set, y = t, fill = t)) +
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
    nudge_y = data$nudge_y, size = 3
  ) +
  scale_x_discrete(name = '', labels = NULL) +
  scale_y_continuous(name = 't-value', limits = c(-55,55)) +
  scale_fill_distiller(palette = 'Spectral', limits = c(-max(data$t), max(data$t))) +
  scale_color_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y =  element_blank(),
    legend.position = 'none'
  ) + ggtitle("Senescence-like2 VS Senescence-like1") +
  theme(plot.title = element_text(hjust = 0.5)) 

tiff("../plots/GSVA_plots.tiff", width = 5*830, height = 5*650, res = 300, pointsize = 5)     
plot_grid(slVSstl, SlVSActEx, ActExvStL, SL2vSL1, nrow = 2)
dev.off()

#GSEA
##SLvsStL
SLvStL_ranked <- FindMarkers(CD8, ident.1 = "SL", ident.2 = "StL", min.pct = 0.1, logfc.threshold = 0)

# order list, pull out gene name and log2fc, and convert genes to uppercase
SLvStL_ranked <- SLvStL_ranked[order(SLvStL_ranked$avg_log2FC, decreasing = T),]
SLvStL_ranked$Gene.name <- str_to_upper(rownames(SLvStL_ranked))
SLvStL_ranked <- SLvStL_ranked[,c("Gene.name", "avg_log2FC")]
rownames(SLvStL_ranked) <- NULL

# read in file containing lists of genes for each pathway
getwd()
hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

SLvStL_ranked <- prepare_ranked_list(SLvStL_ranked)
head(SLvStL_ranked)

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = SLvStL_ranked,
                       minSize = 15,
                       maxSize = 500)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}

wf_SLvStL <- waterfall_plot(fgsea_results, "Pathways enriched in SL vs StL")

ActExvStL_ranked <- FindMarkers(CD8, ident.1 = "ActEx", ident.2 = "StL", min.pct = 0.1, logfc.threshold = 0)

# order list, pull out gene name and log2fc, and convert genes to uppercase
ActExvStL_ranked <- ActExvStL_ranked[order(ActExvStL_ranked$avg_log2FC, decreasing = T),]
ActExvStL_ranked$Gene.name <- str_to_upper(rownames(ActExvStL_ranked))
ActExvStL_ranked <- ActExvStL_ranked[,c("Gene.name", "avg_log2FC")]
rownames(ActExvStL_ranked) <- NULL

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

ActExvStL_ranked <- prepare_ranked_list(ActExvStL_ranked)
head(ActExvStL_ranked)

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = ActExvStL_ranked,
                       minSize = 15,
                       maxSize = 500)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}

wf_ActExvStL <- waterfall_plot(fgsea_results, "Pathways enriched in SL vs StL")

#SLvActEx
SLvActEx_ranked <- FindMarkers(CD8, ident.1 = "SL", ident.2 = "ActEx", min.pct = 0.1, logfc.threshold = 0)

# order list, pull out gene name and log2fc, and convert genes to uppercase
SLvActEx_ranked <- SLvActEx_ranked[order(SLvActEx_ranked$avg_log2FC, decreasing = T),]
SLvActEx_ranked$Gene.name <- str_to_upper(rownames(SLvActEx_ranked))
SLvActEx_ranked <- SLvActEx_ranked[,c("Gene.name", "avg_log2FC")]
rownames(SLvActEx_ranked) <- NULL

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

SLvActEx_ranked <- prepare_ranked_list(SLvActEx_ranked)
head(SLvActEx_ranked)

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = SLvActEx_ranked,
                       minSize = 15,
                       maxSize = 500)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}

wf_SLvActEx <- waterfall_plot(fgsea_results, "Pathways enriched in SL vs ActEx")

#SL2 vs SL1
SL2vSL1_ranked <- FindMarkers(CD8.sen, ident.1 = "SL2", ident.2 = "SL1", min.pct = 0.1, logfc.threshold = 0)

# order list, pull out gene name and log2fc, and convert genes to uppercase
SL2vSL1_ranked <- SL2vSL1_ranked[order(SL2vSL1_ranked$avg_log2FC, decreasing = T),]
SL2vSL1_ranked$Gene.name <- str_to_upper(rownames(SL2vSL1_ranked))
SL2vSL1_ranked <- SL2vSL1_ranked[,c("Gene.name", "avg_log2FC")]
rownames(SL2vSL1_ranked) <- NULL

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt.txt")
head(names(hallmark_pathway))

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

SL2vSL1_ranked <- prepare_ranked_list(SL2vSL1_ranked)
head(SL2vSL1_ranked)

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = SL2vSL1_ranked,
                       minSize = 15,
                       maxSize = 500)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}

wf_SL2vSL1 <- waterfall_plot(fgsea_results, "Pathways enriched in SL2 vs SL1")

tiff("../plots/wf_plots.tiff", width = 5*900, height = 5*650, res = 300, pointsize = 5)     
plot_grid(wf_SLvStL, wf_ActExvStL, wf_SLvActEx, wf_SL2vSL1, nrow = 2)
dev.off()

#SCPA
library(SCPA)
pathways <- msigdbr("Homo sapiens", "C5") %>%
  format_pathways()


# The populations here just need to be your normalized expression matrices
scpa_out1 <- compare_seurat(CD8,
                           group1 = "clusters", 
                           group1_population = c("SL", "StL"),
                           pathways = pathways)
head(scpa_out,10)

p1 <- plot_rank(scpa_out1, "lact", 
                highlight_point_size = 3.5, highlight_point_color = "#60c5f7") + ggtitle("SL vs StL")

p2 <- plot_rank(scpa_out, "TNF",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")

SLvsStL <- patchwork::wrap_plots(p1, p2)

scpa_out2 <- compare_seurat(CD8,
                           group1 = "clusters", 
                           group1_population = c("SL", "ActEx"),
                           pathways = pathways)

p1 <- plot_rank(scpa_out2, "hypoxia", 
                highlight_point_size = 3.5, highlight_point_color = "#60c5f7") + ggtitle("SL vs ActEx")

p2 <- plot_rank(scpa_out2, "TNF",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")

SLvsActEx <-patchwork::wrap_plots(p1, p2)

scpa_out3 <- compare_seurat(CD8,
                            group1 = "clusters", 
                            group1_population = c("ActEx", "StL"),
                            pathways = pathways)

p1 <- plot_rank(scpa_out3, "hypoxia", 
                highlight_point_size = 3.5, highlight_point_color = "#60c5f7") + ggtitle("ActEx vs StL")

p2 <- plot_rank(scpa_out3, "TNF",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")

ActExvsStL <-patchwork::wrap_plots(p1, p2)

scpa_out4 <- compare_seurat(CD8.sen,
                            group1 = "clusters", 
                            group1_population = c("SL2", "SL1"),
                            pathways = pathways)

p1 <- plot_rank(scpa_out4, "hypoxia", 
                highlight_point_size = 3.5, highlight_point_color = "#60c5f7") + ggtitle("SL2 vs SL1")

p2 <- plot_rank(scpa_out4, "TNF",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")

SL2vsSL1 <-patchwork::wrap_plots(p1, p2)

tiff("../plots/rank_plots.tiff", width = 5*400, height = 5*500, res = 300, pointsize = 5)     
plot_grid(SLvsStL, SLvsActEx, ActExvsStL, SL2vsSL1, nrow = 4)
dev.off()

pathways <- "h_k_r_go_pid_reg_wik.csv"



#Velocity
#Load in the files
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
CD8$clusters <- CD8@active.ident
md <- CD8@meta.data
md$clusters <- as.vector(md$clusters)
CD8@meta.data <- md

CD8.velo <- CD8
SaveH5Seurat(CD8, filename = "CD8_veloscv.h5Seurat")
Convert("CD8_veloscv.h5Seurat", dest = "h5ad") 

#Clonotyep for scirpy
CD8.clono <- CD8
CD8.clono <- subset(CD8.clono, subset = group_id == "HD", invert = TRUE)
CD8.clono$batch <- NULL
SaveH5CD8(CD8, filename = "CD8_clonoscv.h5CD8")
Convert("CD8_clonoscv.h5CD8", dest = "h5ad") 

SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat")
Convert("pbmc3k.h5Seurat", dest = "h5ad")

#Clonotype analysis
#Delete HD (we have VDJ only for responders and nonresponders)
CD8 <- subset(CD8, group_id == c("Res", "NonRes"))
CD8@meta.data$cloneType <- factor(CD8@meta.data$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(CD8, group.by = "cloneType", split.by = "group_id") + scale_color_manual(values = c(colorblind_vector(5)), na.value="grey")

alluvialClonotypes(CD8, cloneCall = "gene", 
                   y.axes = c("cluster", "group_id"), 
                   color = "cluster") 

CD8@meta.data$cluster_id <- CD8@active.ident

tiff("../plots/clonoUMAP.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p <- clonalOverlay(CD8, reduction = "umap", 
                   freq.cutpoint = 30, bins = 10, facet = "group_id") + 
  guides(color = FALSE) 
p
dev.off()

CD8_list <- SplitObject(CD8, split.by = "group_id")
occrep <- lapply(CD8_list, function(x) occupiedscRepertoire(x, x.axis = "clusters"))

tiff("../plots/occrep.tiff", width = 5*600, height = 5*300, res = 300, pointsize = 5)     
cowplot::plot_grid(plotlist =occrep, nrow =1, labels = names(occrep), align = 'v', hjust= -2)
dev.off()

tiff("../plots/alluv.tiff", width = 5*600, height = 5*300, res = 300, pointsize = 5)     
alluvialClonotypes(CD8, cloneCall = "gene", 
                   y.axes = c("sample_id", "clusters", "group_id"), 
                   color = "clusters") 
dev.off()


circles <- getCirclize(CD8, 
                       group.by = "clusters")

#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(CD8@active.ident)))
names(grid.cols) <- levels(CD8@active.ident)

#Graphing the chord diagram
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)

saveRDS(CD8, file = "CD8_final.rds")
