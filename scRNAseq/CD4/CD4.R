##CD4 This script can be organized into n parts

##A)Subclustering and use of DGE to annotate clusters
#1) Subclustering as suggested here https://github.com/satijalab/Seurat/issues/2087 by timoast
#2) DGE (top10 markers) 
#3) Annotate clusters
#4) Remove patient-specific cluster (#) and recluster
#5) DGE (top10 markers) and reannotate


##B) Analyze annotated subsets
#1) DGE (top10 markers)
#2) DA analysis using permutation test and speckle (Anova) 

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
CD4 <- readRDS("CD4sub.rds")


######################################################
#A) Subclustering and use of DGE to annotate clusters
######################################################

#1) Subclustering as suggested here https://github.com/satijalab/Seurat/issues/2087 by timoast
DefaultAssay(CD4) <- "RNA"
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", 
                            nfeatures = 2000, 
                            verbose = FALSE)

# Scale the counts in the integrated assay
DefaultAssay(CD4) <- "integrated"
CD4 <- ScaleData(CD4)

# Run PCA and UMAP
CD4 <- RunPCA(CD4)
CD4 <- RunUMAP(CD4, dims = 1:30, verbose = FALSE)

# Look at the batches
DimPlot(CD4, reduction = "pca", group.by = "batch")
DimPlot(CD4, reduction = "umap", group.by = "batch")

# Plot the elbow plot
ElbowPlot(object = CD4, ndims = 30)

# Determine the K-nearest neighbor graph
CD4 <- FindNeighbors(object = CD4, dims = 1:20)

# Determine the clusters for various resolutions                                
CD4 <- FindClusters(object = CD4, resolution = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2))

# set cluster IDs to resolution 1 clustering
CD4 <- SetIdent(CD4, value = "integrated_snn_res.0.4")

DefaultAssay(CD4) <- "RNA" #Put "RNA" as default assay

#Custom palettes for visualization purpose
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

#Visualize the CD4 subset
p.cd4 <- DimPlot_scCustom(CD4, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")

tiff("../plots_CD4/cd4clus.tiff", width = 5*400, height = 5*400, res = 300, pointsize = 5)     
p.cd4
dev.off()

#2) DGE (top10 markers) 
mark <- FindAllMarkers(CD4)

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
tiff("../plots_CD4/p.dotDGE_noAnnot.tiff", width = 5*400, height = 5*600, res = 300, pointsize = 5)     
Clustered_DotPlot(CD4, all_markers, colors_use_exp = pal_exp,
                  exp_color_min = -1, exp_color_max = 1, colors_use_idents = pal_ident,
                  x_lab_rotate = TRUE, k =4)
dev.off()

#3) Annotate clusters
CD4 <- RenameIdents(object = CD4,
                    "0" = "Naive",
                    "2" = "Naive", 
                    "6" = "Naive",
                    "10" = "Naive",
                    "1" = "EM",
                    "3" = "EM",
                    "4" = "ActEff",
                    "5" = "EM NR4A2+",
                    "7" = "Tregs",
                    "8" = "Tregs",
                    "9" = "Naive NR4A2+")

CD4$clusters <- CD4@active.ident
tiff("../plots_CD4/cd4ann.tiff", width = 5*400, height = 5*400, res = 300, pointsize = 5)     
p.cd4 <- DimPlot_scCustom(CD4, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")
p.cd4
dev.off()

tiff("../plots_CD4/barplotCD4.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
ggplot(CD4@meta.data, aes(x=clusters, fill=sample_id)) + geom_bar(position = "fill") + theme_bw() + 
  scale_fill_manual(values= pal_ident) + xlab("") + ylab("") + coord_flip()
dev.off()

#4) Remove patient-specific cluster and recluster
CD4 <- subset(CD4, subset = clusters == "Naive NR4A2+", invert = TRUE )

DefaultAssay(CD4) <- "RNA"
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", 
                            nfeatures = 2000, 
                            verbose = FALSE)

# Scale the counts in the integrated assay
DefaultAssay(CD4) <- "integrated"
CD4 <- ScaleData(CD4)

# Run PCA and UMAP
CD4 <- RunPCA(CD4)
CD4 <- RunUMAP(CD4, dims = 1:30, verbose = FALSE)

# Look at the batches
DimPlot(CD4, reduction = "pca", group.by = "batch")
DimPlot(CD4, reduction = "umap", group.by = "batch")

# Plot the elbow plot
ElbowPlot(object = CD4, ndims = 30)

# Determine the K-nearest neighbor graph
CD4 <- FindNeighbors(object = CD4, dims = 1:20)

# Determine the clusters for various resolutions                                
CD4 <- FindClusters(object = CD4, resolution = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2))

# set cluster IDs to resolution 1 clustering
CD4 <- SetIdent(CD4, value = "integrated_snn_res.0.4")

DefaultAssay(CD4) <- "RNA" #Put "RNA" as default assay

pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

#Visualize the CD4 subset
p.cd4 <- DimPlot_scCustom(CD4, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")

tiff("../plots_CD4/cd4clus.tiff", width = 5*400, height = 5*400, res = 300, pointsize = 5)     
p.cd4
dev.off()

#5) DGE (top10 markers) and reannotate
mark <- FindAllMarkers(CD4)
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
cel_by_clus <- split(colnames(CD4), CD4@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD4, slot = "counts")
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
DimPlot(CD4, label = TRUE)

CD4 <- RenameIdents(object = CD4,
                    "0" = "Naive",
                    "2" = "Naive", 
                    "6" = "Naive",
                    "9" = "Naive",
                    "5" = "EM NR4A2+",
                    "1" = "EM",
                    "3" = "EM",
                    "4" = "ActEff",
                    "7" = "Tregs",
                    "8" = "Tregs")

CD4$clusters <- CD4@active.ident
tiff("../plots_CD4/cd4ann.tiff", width = 5*400, height = 5*400, res = 300, pointsize = 5)     
p.cd4 <- DimPlot_scCustom(CD4, label = TRUE, colors_use = pal_ident, figure_plot = TRUE, pt.size = 0.00001) + 
  theme(legend.position="none")
p.cd4
dev.off()

##############################
#2) Analyze annotated subsets
##############################

#1) DGE (top10 markers) 
mark <- FindAllMarkers(CD4)
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
cel_by_clus <- split(colnames(CD4), CD4@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD4, slot = "counts")
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
heat.Annot <- Heatmap(t(mat),
                      cluster_rows = TRUE,
                      cluster_columns = TRUE,
                      row_names_side = "left",
                      row_names_gp = grid::gpar(fontsize = 10),
                      heatmap_legend_param = lgd_aes,
                      column_names_gp = gpar(fontsize = 7))

#2) DA analysis using permutation test and speckle (Anova) 
#permutation test
prop_test <- sc_utils(CD4)
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

tiff("../plots_CD4/perm.tiff", width = 5*300, height = 5*60, res = 150, pointsize = 5)     
p.perm
dev.off()

#Speckle
props <- getTransformedProps(CD4$clusters, CD4$sample_id,
                             transform="logit")

tiff("../plots_CD4/MeanVar.tiff", width = 5*150, height = 5*150, res = 300, pointsize = 5)     
plotCellTypePropsMeanVar(props$Counts)
dev.off()

condition <- c(rep("HD",2), rep("NonRes",2), rep("Res", 4), rep("NonRes", 2), 
               rep("Res", 2))
pair <- rep(c(1,2,3,4,5,6), each = 2)
batch <- c(rep(1,2), rep(2,2), rep(3,2), rep(2,2), rep(1,2), rep(2,2))
design <- model.matrix(~0 + condition + pair + batch)
df.anova <- propeller.anova(prop.list=props, design=design, coef = c(1,2,3), 
                            robust=TRUE, trend=FALSE, sort=TRUE) %>% as.data.frame() 

#as boxplots
fqs <- props$Proportions
df.plot <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df.plot$group_id <- CD4$group_id[match(df.plot$sample_id, CD4$sample_id)]
df.plot <- df.plot %>% dplyr::mutate(positions = case_when(cluster_id == "EM NR4A2+" ~ "Anova: FDR ***, F 18.976471", TRUE ~ ""))
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

