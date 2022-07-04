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
library(mclust)
library(SCPA)
library(dyno)
library(scRepertoire)
library(SeuratWrappers)
library(monocle3)
library(hdf5r)
library(circlize)
library(ComplexHeatmap)
library(tibble)
library(msigdbr)
library(statmod)
library(scater)
library(mgcv)
library(gam)
library(tradeSeq)
library(condiments)
library(stringr)
library(scuttle)
library(destiny)
library(SeuratDisk)
library(cmapR)
library(escape)
library(dittoSeq)
library(snow)
library(UCell)

setwd("~/Documents/AML_project/scRNA_AMLproj")


#Subclustering the integrated assay as suggested here https://github.com/satijalab/seurat/issues/2087 by timoast
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
#CD8 <- RunTSNE(CD8, dims = 1:30, verbose = FALSE)

# Look at the batches
DimPlot(CD8, reduction = "pca", group.by = "batch")
DimPlot(CD8, reduction = "umap", group.by = "batch")

# Plot the elbow plot
ElbowPlot(object = CD8, ndims = 30)

# Determine the K-nearest neighbor graph
CD8 <- FindNeighbors(object = CD8, dims = 1:20)

# Determine the clusters for various resolutions                                
CD8 <- FindClusters(object = CD8, resolution = 1.2)

# set cluster IDs to resolution 1 clustering
CD8 <- SetIdent(CD8, value = "integrated_snn_res.1.2")

DefaultAssay(CD8) <- "RNA" #Put "RNA" as default assay

#Visualize the new subcluster
p <- DimPlot(object = CD8, 
             reduction = "umap", 
             label = TRUE,
             label.size = 3,
             repel = TRUE)
p

###Visualize single markers and scores distributions
#Check they are CD8
T_CD8 <- paste0("CD8A$|CD8B$")
T_CD8 <- rownames(CD8)[grep(T_CD8, rownames(CD8))]
FeaturePlot(CD8, T_CD8) 

#Check there are no CD4
T_CD4 <- paste0("[.]CD4$")
T_CD4 <- rownames(CD8)[grep(T_CD4, rownames(CD8))]
FeaturePlot(CD8, T_CD4) 

#Single Markers
TCF7 <- "TCF7$"
TCF7 <- rownames(CD8)[grep(TCF7, rownames(CD8))]
FeaturePlot(CD8, TCF7,min.cutoff = "q10", max.cutoff = "q90") + ggtitle("TCF7")

CX3CR1 <- "CX3CR1$"
CX3CR1 <- rownames(CD8)[grep(CX3CR1, rownames(CD8))]
FeaturePlot(CD8, CX3CR1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("CX3CR1")

tiff("./plots/CX3CR1.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, CX3CR1, pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.8, 1.9), label = c("0", "Maximum")) +
  ggtitle("Effectors CX3CR1+")
dev.off()

LEF1 <- "LEF1$"
LEF1 <- rownames(CD8)[grep(LEF1, rownames(CD8))]
FeaturePlot(CD8, LEF1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("LEF1")

CD69 <- "CD69$"
CD69 <- rownames(CD8)[grep(CD69, rownames(CD8))]

tiff("./plots/CD69.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, CD69, min.cutoff = "q10", max.cutoff = "q90") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("CD69")
dev.off()

PD1 <- "PDCD1$"
PD1 <- rownames(CD8)[grep(PD1, rownames(CD8))]
FeaturePlot(CD8, PD1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("PD1")

IL7R <- "IL7R$"
IL7R <- rownames(CD8)[grep(IL7R, rownames(CD8))]
tiff("./plots/IL7R.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, IL7R, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("IL7R") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("IL7R")
dev.off()

GZMK <- "GZMK$"
GZMK <- rownames(CD8)[grep(GZMK, rownames(CD8))]
tiff("./plots/GZMK.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GZMK, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GZMK")
dev.off()

EOMES <- "EOMES$"
EOMES <- rownames(CD8)[grep(EOMES, rownames(CD8))]
FeaturePlot(CD8, EOMES, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("EOMES")

TBX21 <- "TBX21$"
TBX21 <- rownames(CD8)[grep(TBX21, rownames(CD8))]
FeaturePlot(CD8, TBX21, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("TBX21")

GZMB <- "GZMB$"
GZMB <- rownames(CD8)[grep(GZMB, rownames(CD8))]
tiff("./plots/GZMB.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GZMB, min.cutoff = "q10", max.cutoff = "q90") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GZMB")
dev.off()

GZMA <- "GZMA$"
GZMA <- rownames(CD8)[grep(GZMA, rownames(CD8))]
FeaturePlot(CD8, GZMA, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("GZMA")

GNLY <- "GNLY$"
GNLY <- rownames(CD8)[grep(GNLY, rownames(CD8))]
tiff("./plots/GNLY.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GNLY, min.cutoff = "q10", max.cutoff = "q90") +   
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")), breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GNLY")
dev.off()

FCGR3A <- "FCGR3A$"
FCGR3A <- rownames(CD8)[grep(FCGR3A, rownames(CD8))]
FeaturePlot(CD8, FCGR3A, min.cutoff = "q10", max.cutoff = "q90") 

CRTAM <- "CRTAM$"
CRTAM <- rownames(CD8)[grep(CRTAM, rownames(CD8))]
FeaturePlot(CD8, CRTAM, min.cutoff = "q10", max.cutoff = "q90") 

KLF2 <- "KLF2$"
KLF2 <- rownames(CD8)[grep(KLF2, rownames(CD8))]
FeaturePlot(CD8, KLF2, min.cutoff = "q10", max.cutoff = "q90") 

CXCR6 <- "CXCR6$"
CXCR6 <- rownames(CD8)[grep(CXCR6, rownames(CD8))]
FeaturePlot(CD8, CXCR6, min.cutoff = "q10", max.cutoff = "q90") 

#Scores
naive_genes <- "TCF7$|CCR7|SELL|LEF1|GZMK"
naive <- list(rownames(CD8)[grep(naive_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = naive, name = "naive", random.seed = 1)
FeaturePlot(CD8, "naive1", min.cutoff = "q10", max.cutoff = "q90") + ggtitle("Tnaive")

Tcm_genes <- "CCR7|SELL|IL7R"
Tcm <- list(rownames(CD8)[grep(Tcm_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tcm, name = "Tcm", random.seed = 1)
FeaturePlot(CD8, "Tcm1", min.cutoff = "q10", max.cutoff = "q90") + ggtitle("Tcm")

Teff_genes <- "CX3CR1|FGFBP2|FCGR3A"
Teff <- list(rownames(CD8)[grep(Teff_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Teff, name = "Teff", random.seed = 1)
FeaturePlot(CD8, "Teff1") + ggtitle("Teff")

Tpex <- "GZMK|IL7R"
Tpex <- list(rownames(CD8)[grep(Tpex, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tpex, name = "Tpex", random.seed = 1)
FeaturePlot(CD8, "Tpex1") + ggtitle("Tpex")

tiff("./plots/stem-like.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "Tpex1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + NoLegend() +
  ggtitle("stem-like IL7R+GZMK+")
dev.off()

# signatures
#Scores
#####Change --> all AddModuleScore or all AddModuleScore_UCell
#signatures
sig <- readxl::read_xlsx("signatures/sig.xlsx", sheet = 1)
sig_naive <- as.list(sig$`Naive (from Szabo et al. 2019)`[!is.na(sig$`Naive (from Szabo et al. 2019)`)])
sig_naive <- lapply(sig_naive, function(x){paste0("[.]", x, "$")})
sig_naive <- paste0(sig_naive, collapse = "|")
sig_naive <- list(rownames(CD8)[grep(sig_naive, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive", random.seed = 1)

tiff("./plots/naive.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "sig_naive1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.02, 0.16), label = c("0", "Maximum")) +
  ggtitle("Naive (from Szabo et al. 2019)")
dev.off()

sig_stem <- as.list(sig$`Stemness (from Pace et al. 2018)`[!is.na(sig$`Stemness (from Pace et al. 2018)`)])
sig_stem <- lapply(sig_stem, function(x){paste0("[.]", x, "$")})
sig_stem <- paste0(sig_stem, collapse = "|")
sig_stem <- list(rownames(CD8)[grep(sig_stem, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem", random.seed = 1)

tiff("./plots/stem.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "sig_stem1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.015, 0.125), label = c("0", "Maximum")) +
  ggtitle("Stemness (from Pace et al. 2018)")
dev.off()

Tex <- as.list(sig$`Terminal exhausted (from Miller et al. 2019)`[!is.na(sig$`Terminal exhausted (from Miller et al. 2019)`)])
Tex <- lapply(Tex, function(x){paste0("[.]",x, "$")})
Tex <- paste0(Tex, collapse = "|")
Tex <- list(rownames(CD8)[grep(Tex, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tex, name = "Tex", random.seed = 1)


tiff("./plots/Tex.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "Tex1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.06, 0.32), label = c("0", "Maximum")) +
  ggtitle("Terminally exhausted (from Miller et al. 2019)")
dev.off()

#Bulk senescence signature
genes <- readxl::read_xlsx("signatures/Model_5_results.xlsx", sheet = 2)
sen <- as.list(genes$gene_name[1:100])
sen <- lapply(sen, function(x){paste0("[.]",x, "$")})
sen <- paste0(sen, collapse = "|")
sen <- list(rownames(CD8)[grep(sen, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sen, name = "sen", random.seed = 1)
tiff("./plots/sen_group.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, features = "sen1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") + labs(col = "EXPRESSION") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("")
dev.off()

tiff("./plots/UMAP.tiff", width = 8*250, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD8, label = TRUE)
dev.off()

#Remove patient dependent subset
CD8 <- subset(CD8, idents = '10', invert = TRUE)

#Clean workspace
rm(list=setdiff(ls(), "CD8"))

CD8 <- RenameIdents(CD8,
                    "0" = "Naive",
                    "1" = "Naive",
                    "2" = "Int",
                    "3" = "Stem-like",
                    "4" = "Senescent-like",
                    "5" = "Tex", 
                    "6" = "Stem-like", 
                    "7" = "Naive", 
                    "8" = "Int", 
                    "9" = "Int-GZMK+",
                    "11" = "Senescent-like", 
                    "12" = "Senescent-like",
                    "13" = "Naive",
                    "14" = "Senescent-like",
                    "15" = "Stem-like",
                    "16" = "Stem-like",
                    "17" = "Stem-like",
                    "18" = "Stem-like",
                    "19" = "Tex",
                    "20" = "Int",
                    "21" = "Int")

tiff("./plots/UMAP_ann.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p <- DimPlot(CD8, label = TRUE)
p
dev.off()

tiff("./plots/UMAP_groupId.tiff", width = 5*600, height = 5*250, res = 300, pointsize = 5)     
p <- DimPlot(CD8, label = TRUE, split.by = "group_id")
p
dev.off()

tiff("./plots/UMAP_density.tiff", width = 5*600, height = 5*400, res = 300, pointsize = 5)     
a <- DimPlot(CD8, reduction = 'umap', split.by = "RespTmp")  + NoLegend() + NoAxes() 
a +  geom_density_2d_filled(a$data, mapping = aes(x = a$data[,"UMAP_1"], y = a$data[,"UMAP_2"]), contour_var = "ndensity") + 
  facet_wrap(vars(RespTmp))
dev.off()

mark <- FindAllMarkers(CD8)
mark %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
df <- data.frame(top20$cluster, top20$gene)
df <- reshape(transform(df, indx = ave(as.character(top20.cluster), top20.cluster, FUN = seq)), 
              idvar = "indx", timevar = "top20.cluster", direction = "wide")
View(top20)
writexl::write_xlsx(df, "./genes.xlsx")

tiff("./plots/DoHeat.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8, features = top20$gene, size = 3, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 6)) + ggtitle("Top 20")
dev.off()

#Relative abundance of the each clusters per condition (How to check for "significance"? What test should we apply?)
freq_table <- table(CD8@active.ident, CD8$group_id)
fqs <- prop.table(table(CD8@active.ident, CD8$sample_id), 2)
df <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df$group_id <- CD8$group_id[match(df$sample_id, CD8$sample_id)]

tiff("./plots/box.tiff", width = 5*400, height = 5*60, res = 150, pointsize = 5)     
ggplot(df, aes(x = group_id, y = frequency, color = group_id)) +
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
dev.off()

#Custom Markers
fs <- list(
  Naive = c("CCR7", "MAL", "LEF1", "SELL", "TCF7"),
  GZMK_IL7R = c("BCL2", "BACH2", "CD27","IL7R", "SLAMF6", "CXCR3", "GZMK",  "ITGAE","ENTPD1"),
  Int_GZMK = c("EOMES", "CCL4", "XCL2", "CCL3", "XCL1", "KLF6", "TIGIT"),
  Tex = c("CD69", "CD160", "PDCD1", "NR4A2",  "DUSP2"),
  Int = c("NKG7", "TBX21"),
  Senescent_like = c("CD38","FCRL6", "FCGR3A", "C1orf21", "PRF1", 
                     "GNLY", "ZEB2", "CX3CR1", "FGFBP2", "KLRD1", "GZMB","ZNF683", "CD226")
)
fs <- lapply(fs, sapply, function(g) 
  grep(paste0("\\.",g, "$"), rownames(CD8), value = TRUE))

#Prepare Heatmap of mean marker-exprs. by cluster
gs <- gsub(".*\\.", "", unlist(fs))
ns <- vapply(fs, length, numeric(1))
ks <- rep.int(names(fs), ns)
labs <- sprintf("%s(%s)", gs, ks)

# split cells by cluster
cbk <- split(colnames(CD8), CD8@active.ident)

# compute cluster-marker means
ct <- GetAssayData(CD8, slot = "counts")
libsizes <- colSums(ct)
sf <- libsizes/mean(libsizes)
log.ct <- log2(t(t(ct)/sf) + 1)

mbk <- lapply(fs, function(gs)
  vapply(cbk, function(i)
    Matrix::rowMeans(log.ct[gs, i, drop = FALSE]), 
    numeric(length(gs))))

# prep. for plotting 
mat <- do.call("rbind", mbk)
rownames(mat) <- gsub(".*\\.", "", rownames(mat))

#Z-score(?)
mat <- t(scale(t(mat)))

#0-1 scaling
qs <- rowQuantiles(mat, probs = c(0.01, 0.99), na.rm = TRUE)
mat <- (mat - qs[, 1])/(qs[, 2] - qs[, 1])
mat[mat < 0] <- 0
mat[mat > 1] <- 1

mat <- mat %>% as.data.frame() %>%  select(names(fs)) %>% as.matrix()

h <- Heatmap(mat,
             name = "Z-score",
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_names_side = "left",
             column_title = "cluster_id",
             column_title_side = "bottom",
             #rect_gp = gpar(col = "white"),
             row_names_gp = grid::gpar(fontsize = 8))
#left_annotation = row_anno)

h 
#Gene set enrichment analysis
#Download geneset acidosis from harmonizome
df.acidosis <- readxl::read_xlsx("harmonizome.xlsx")

acidosis <- df.acidosis$acidosis

gene.sets <- getGeneSets(library = "H")
gene.sets$ACIDOSIS <- acidosis

sceCD8 <- as.SingleCellExperiment(CD8, assay = "RNA")
ss1 <- strsplit(rownames(sceCD8), ".", fixed = TRUE)
rownames(sceCD8) <- sapply(ss1, .subset, 2)

ES <- enrichIt(obj = sceCD8, gene.sets = gene.sets, groups = 1000, cores = 6, method = "UCell")

CD8 <- AddMetaData(CD8, ES)

CD8@meta.data$active.idents <- CD8@active.ident
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

dittoHeatmap(CD8, genes = NULL, metas = names(ES),
             heatmap.colors = rev(colorblind_vector(50)),
             annot.by = c("active.idents", "group_id"),
             cluster_cols = TRUE,
             fontsize = 7)

tiff("./plots/hypoxia.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "HALLMARK_HYPOXIA", group.by = "group_id", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

tiff("./plots/hypoxia_RspTmp.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "HALLMARK_HYPOXIA", group.by = "RespTmp", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

tiff("./plots/hypoxia_idents.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "HALLMARK_HYPOXIA", group.by = "active.idents", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

tiff("./plots/acidosis.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "Acidosis", group.by = "group_id", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

tiff("./plots/acidosis_RspTmp.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "Acidosis", group.by = "RespTmp", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

tiff("./plots/acidosis_idents.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
dittoBoxPlot(CD8, "Acidosis", group.by = "active.idents", jitter.size = 0) + 
  scale_fill_manual(values = colorblind_vector(6)) + theme_bw()
dev.off()

#TRAJECTORY INFERENCE using Slingshot
clusterLabels <- CD8@active.ident
sceCD8 <- as.SingleCellExperiment(CD8, assay = "RNA")
sds <- slingshot(sceCD8, clusterLabels = clusterLabels, 
                 allow.breaks = TRUE, stretch = 2, reducedDim = "UMAP", start.clus = "Naive") #Calcualting the trajectory
sds <- SlingshotDataSet(sds)

#Create dataframe pseudotimes
pt <- as.data.frame(slingPseudotime(sds))
colnames(pt) <- c("pseudoT1", "pseudoT2")
names <- rownames(pt)
   
#Extract the curves coordinates
curve1 <- slingCurves(sds)[[1]]
curve2 <- slingCurves(sds)[[2]]
   
#Create dataframe with UMAP coordinates
umap_df <- reducedDim(sceCD8, "UMAP")
   
#Create dataframe with UMAP coordinates and pseudotimes
df <- cbind(umap_df, pt)
   
#Plot according to pseudotime values  (pseudoT1)
p1 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
   geom_point(aes_string(color = df$pseudoT1),
              alpha = 0.5) +
   scale_colour_viridis_c() +
   theme_minimal() + labs(colour = "Pseudotime") 
   
   
#Plot 1st lineage
tiff("./plots/traj1.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p1 <- p1 + geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve1$s[curve1$ord, ] %>% as.data.frame(),
                     col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime") + ggtitle("Senescence")
dev.off()
   
#Plot according to pseudotime (pseudoT2)
p2 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
   geom_point(aes_string(color = df$pseudoT2),
              alpha = 0.5) +
   scale_colour_viridis_c() +
   theme_minimal() + labs(colour = "Pseudotime")
   
tiff("./plots/traj2.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p2 <- p2 +   geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve2$s[curve2$ord, ] %>% as.data.frame(),
                       col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime") + ggtitle("Exhaustion")
dev.off()
   
tiff("./plots/traj.tiff", width = 8*130, height = 5*300, res = 150, pointsize = 5)  
plot_grid(p1, p2, ncol = 1, align = "h")
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
sceCD8 <- fitGAM(counts = counts(sceCD8), sds = sds, nknots = 5, verbose = TRUE, BPPARAM = BPPARAM)

# plot our Slingshot lineage trajectories, this time illustrating the new tradeSeq knots
tiff("./plots/traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plotGeneCount(curve = sds, counts = counts,
             clusters = CD8@active.ident,
             models = sceCD8)
dev.off()
   
#Association test
assoRes <- associationTest(sceCD8)
head(assoRes)
   
### Discovering differentiated cell type markers
# discover marker genes for the differentiated cell types
endRes <- diffEndTest(sceCD8) #Nothing interesting with this analysis
head(endRes)
   
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sceCD8)[o[2]]
plotSmoothers(sceCD8, counts(sceCD8), sigGene) 
   
plotGeneCount(sds, counts(sceCD8), gene = sigGene)

# Marker genes between specific roots   
earlyDERes <- earlyDETest(sceCD8, knots = c(3, 4))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(earlyDERes)[oEarly][1])
plotGeneCount(sds, counts(sceCD8), gene = rownames(earlyDERes)[oEarly][3])

#Genes with different expression patterns (most interesting part)
patternRes <- patternTest(sceCD8)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])

tiff("./plots/DEtrajGZMK.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(patternRes)[oPat][1]) + ggtitle ("GZMK")
dev.off()

tiff("./plots/DEtrajNKG7.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(patternRes)[oPat][2]) + ggtitle ("NKG7")
dev.off()

tiff("./plots/DEtrajGNLY.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(patternRes)[oPat][7]) + ggtitle ("GNLY")
dev.off()

tiff("./plots/DEtrajPRF1.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(patternRes)[oPat][8]) + ggtitle ("PRF1")
dev.off()

tiff("./plots/DEtrajGZMB.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sceCD8, counts(sceCD8), gene = rownames(patternRes)[oPat][9]) + ggtitle ("GZMB")
dev.off()

#DGE between lineages conditions according to pseudotime
#First consider only Res and NonRes 
sce <- as.SingleCellExperiment(CD8, assay = "RNA") #restart before running Slingshot and fitGAM
sce.cond <- subset(sce,,group_id %in% c("Res", "NonRes")) #filter out HD
sds <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$ident,
                    start.clus = 'Tnaive', approx_points = 150)
sds <- SlingshotDataSet(sds)

#Trajectory by condition (compare Res and NonRes)
df <- bind_cols(
 as.data.frame(reducedDim(sce.cond, "UMAP")),
 as.data.frame(colData(sce.cond))
) %>%
 sample_frac(1)

scores <- condiments::imbalance_score(
 Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
 conditions = df$group_id)
   
df$scores <- scores$scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
 geom_point(size = .7) +
 scale_color_viridis_c(option = "C") +
 labs(col = "Scores")
p3

df$scaled_scores <- scores$scaled_scores
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scaled_scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores")
p4

topologyTest(SlingshotDataSet(sds), sce.cond$group_id, rep = 100,
             methods = "KS_mean", threshs = .01)
knitr::kable(top_res)

psts <- slingPseudotime(sds) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$group_id) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

ggplot(psts, aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")

prog_res <- progressionTest(sds, conditions = df$group_id, global = TRUE, lineages = TRUE)
knitr::kable(prog_res)

df$weight_1 <- slingCurveWeights(sds, as.probs = TRUE)[, 1]
ggplot(df, aes(x = weight_1, fill = group_id)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for the first lineage")

dif_res <- differentiationTest(sds, conditions = df$group_id, global = FALSE, pairwise = TRUE)
knitr::kable(dif_res)

#Run fitGAM according to conditions and plot heatmaps 
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 10

set.seed(5)
icMat <- evaluateK(counts(sce), sds = sds, k = 3:7, nGenes = 100, verbose = T, conditions = factor(colData(sce)$group_id), plot = T)

set.seed(3)
sce <-fitGAM(counts = counts(sce), sds = sds, nknots = 5, conditions = factor(colData(sce)$group_id), verbose = TRUE, BPPARAM = BPPARAM)

#Differential expression between lineages
patRes <- patternTest(sce, l2fc = log2(2))
patRes$padj <- p.adjust(patRes$pvalue, "fdr")
mean(patRes$padj <= 0.05, na.rm = TRUE)

sum(patRes$padj <= 0.05, na.rm = TRUE)

patternGenes <- rownames(patRes)[patRes$padj <= 0.05]
patternGenes <- patternGenes[!is.na(patternGenes)]

### based on mean smoother
#DGE NonRes lineage1 vs lineage2
yhatSmooth <- predictSmooth(sce, gene = patternGenes, nPoints = 50, tidy = TRUE) %>%
  mutate(yhat = log1p(yhat)) %>%
  group_by(gene) %>%
  mutate(yhat = scales::rescale(yhat)) %>%
  filter(condition == "NonRes") %>%
  select(-condition) %>%
  ungroup()

heatSmooth_Lineage1 <- pheatmap(
  yhatSmooth %>%
    filter(lineage == 1) %>%
    select(-lineage) %>%
    arrange(-time) %>%
    pivot_wider(names_from = time, values_from = yhat) %>%
    select(-gene),
  cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
  main = "Lineage 1", legend = FALSE, silent = TRUE
)

heatSmooth_Lineage2 <- pheatmap(
  (yhatSmooth %>%
     filter(lineage == 2) %>%
     select(-lineage) %>%
     arrange(-time) %>%
     pivot_wider(names_from = time, values_from = yhat) %>%
     select(-gene))[heatSmooth_Lineage1@row_order, ],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "Lineage 2",
  legend = FALSE, silent = TRUE 
)

heatSmooth_Lineage1
heatSmooth_Lineage2

#DGE Responders libeage1 vs lineage2
yhatSmooth <- predictSmooth(sce, gene = patternGenes, nPoints = 50, tidy = TRUE) %>%
  mutate(yhat = log1p(yhat)) %>%
  group_by(gene) %>%
  mutate(yhat = scales::rescale(yhat)) %>%
  filter(condition == "NonRes") %>%
  select(-condition) %>%
  ungroup()

heatSmooth_Lineage1 <- pheatmap(
  yhatSmooth %>%
    filter(lineage == 1) %>%
    select(-lineage) %>%
    arrange(-time) %>%
    pivot_wider(names_from = time, values_from = yhat) %>%
    select(-gene),
  cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
  main = "Lineage 1", legend = FALSE, silent = TRUE
)

heatSmooth_Lineage2 <- pheatmap(
  (yhatSmooth %>%
     filter(lineage == 2) %>%
     select(-lineage) %>%
     arrange(-time) %>%
     pivot_wider(names_from = time, values_from = yhat) %>%
     select(-gene))[heatSmooth_Lineage1@row_order, ],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "Lineage 2",
  legend = FALSE, silent = TRUE 
)

heatSmooth_Lineage1
heatSmooth_Lineage2

#Differential expression between conditions
condRes <- conditionTest(sce, l2fc = log2(2), lineages = TRUE)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)

sum(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

scales <- brewer.pal(3, "Accent")[1:4]

# plot genes
oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
p1 <- plotSmoothers(sce, assays(sce)$counts,
                    gene = rownames(assays(sce)$counts)[oo[1]],
                    alpha = 1, border = TRUE, curvesCols = scales) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(sce)$counts)[oo[1]])

p1

# Second most significant gene
p2 <- plotSmoothers(sce, assays(sce)$counts,
                    gene = rownames(assays(sce)$counts)[oo[2]],
                    alpha = 1, border = TRUE, curvesCols = scales) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(sce)$counts)[oo[2]])

p2

### based on mean smoother
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
conditionGenes_lineage1 <- rownames(condRes)[condRes$padj_lineage1 <= 0.05]
conditionGenes_lineage1 <- conditionGenes_lineage1[!is.na(conditionGenes_lineage1)]

yhatSmooth <- predictSmooth(sce, gene = conditionGenes_lineage1, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth[, c(1:50, 101:150)], 1, scales::rescale))
heatSmooth_NonRes <- pheatmap(yhatSmoothScaled[, 1:50],
                              cluster_cols = FALSE,
                              show_rownames = FALSE, show_colnames = FALSE, main = "NonRes", legend = FALSE,
                              silent = TRUE
)

matchingHeatmap_Res <- pheatmap(yhatSmoothScaled[heatSmooth_NonRes@row_order, 51:100],
                                cluster_cols = FALSE, cluster_rows = FALSE,
                                show_rownames = FALSE, show_colnames = FALSE, main = "Res",
                                legend = FALSE, silent = TRUE
)

heatSmooth_NonRes
matchingHeatmap_Res

condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
conditionGenes_lineage2 <- rownames(condRes)[condRes$padj_lineage2 <= 0.05]
conditionGenes_lineage2 <- conditionGenes_lineage2[!is.na(conditionGenes_lineage2)]

yhatSmooth <- predictSmooth(sce, gene = conditionGenes_lineage2, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth[, c(51:100, 151:200)], 1, scales::rescale))
colnames(yhatSmoothScaled)
heatSmooth_NonRes <- pheatmap(yhatSmoothScaled[, 1:50],
                              cluster_cols = FALSE,
                              show_rownames = FALSE, show_colnames = FALSE, main = "NonRes", legend = FALSE,
                              silent = TRUE
)

matchingHeatmap_Res <- pheatmap(yhatSmoothScaled[heatSmooth_NonRes$tree_row$order, 51:100],
                                cluster_cols = FALSE, cluster_rows = FALSE,
                                show_rownames = FALSE, show_colnames = FALSE, main = "Res",
                                legend = FALSE, silent = TRUE
)

p10 <- plot_grid(heatSmooth_NonRes[[4]], matchingHeatmap_Res[[4]],
                 NULL, NULL, NULL, ncol = 2, rel_widths = c(1.4, 1, 1), rel_heights = c(10, 1)) +
  draw_text("Lineage 2", x = .5, y = .5)
p10

plot_grid(p9, p10, ncol = 2)

##Fix the error (GSEA does not work)
## C7 category is according to gene ontology grouping: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
geneSets <- msigdbr(species = "Homo sapiens", category = "C7")
### filter background to only include genes that we assessed.
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
ss1 <- strsplit(rownames(sce), ".", fixed=TRUE)
genenames <- sapply(ss1, .subset, 2)
geneSets <- geneSets[geneSets$gene_symbol %in% genenames,]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
condRes <- conditionTest(sce, l2fc = log2(2))
statsCond <- condRes$waldStat
ss2 <- strsplit(rownames(condRes), ".", fixed=TRUE)
names(statsCond) <- sapply(ss1, .subset, 2)
eaRes <- fgsea(pathways = m_list, stats = statsCond, nperm = 5e4, minSize = 10)
ooEA <- order(eaRes$pval, decreasing = FALSE)
kable(head(eaRes[ooEA, 1:3], n = 20))
##????

#Velocity
#Load in the files
ld.list <- list.files("velocity", pattern = "loom", all.files = TRUE)
ld <- lapply(ld.list, function(x){
  ReadVelocity(paste0("velocity/", x))
})

ld.list <- gsub(".loom", "_", ld.list)

#Rename colnames for each object according to the colnames of the integrated Seurat object
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

#Import the assays in the original Seurat object
CD8[["spliced"]] <- spliced
CD8[["unspliced"]] <- unspliced
CD8[["ambiguous"]] <- ambiguous

#Assign the spliced assay to the RNA assay and put it as default
CD8[['RNA']] <- CD8[["spliced"]]
DefaultAssay(CD8) <- "RNA"
CD8$celltype <- CD8@active.ident

#save seurat and convert into h5ad
SaveH5Seurat(CD8, filename = "CD8scv.h5Seurat")
Convert("CD8scv.h5Seurat", dest = "h5ad") #--> continue this using Python (script velocity.py in the Github repo)

#Clonotype analysis
#Delete HD (we have VDJ only for responders and nonresponders)
CD8 <- subset(CD8, group_id == c("Res", "NonRes"))

CD8@meta.data$cloneType <- factor(CD8@meta.data$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(CD8, group.by = "cloneType", split.by = "group_id") + scale_color_manual(values = c(colorblind_vector(5)), na.value="grey")
TCD8 <- expression2List(CD8, group = "cluster")
Samples <- expression2List(CD8, group = "RespTmp")

p229 <- subset(CD8, orig.ident == "p229")
p229 <- expression2List(p229, group = "cluster")

p219 <- subset(CD8, orig.ident == "p219")
p219 <- expression2List(p219, group = "cluster")

p264 <- subset(CD8, orig.ident == "p264")
p264 <- expression2List(p264, group = "cluster")

pV03 <- subset(CD8, orig.ident == "pV03")
pV03 <- expression2List(pV03, group = "cluster")

pV09 <- subset(CD8, orig.ident == "pV09")
pV09 <- expression2List(pV09, group = "cluster")

clonalProportion(p229, cloneCall = "gene+nt")
clonalProportion(p219, cloneCall = "gene+nt")
clonalProportion(p264, cloneCall = "gene+nt")
clonalProportion(pV03, cloneCall = "gene+nt")
clonalProportion(pV09, cloneCall = "gene+nt")

tiff("./pipelines/plots/clono.tiff", width = 5*10, height = 5*10, res = 300, pointsize = 5)     
compareClonotypes(p229, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(p219, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(p264, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(pV03, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(pV09, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")

quantContig(TCD8, cloneCall="gene+nt", scale = T) + guides(fill = F)
quantContig(Samples, cloneCall="gene+nt", scale = TRUE) + guides(fill = F)

clonalHomeostasis(TCD8, cloneCall = "gene+nt")
clonalProportion(TCD8, cloneCall = "gene+nt")

clonalOverlap(TCD8, cloneCall = "gene+nt", method = "morisita") + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11)), na.value = "white")

clonalOverlap(TCD8, cloneCall = "gene+nt", method = "overlap") + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11)), na.value = "white")

clonalDiversity(CD8, cloneCall = "gene+nt", group = "cluster")

clonalDiversity(CD8, cloneCall = "gene+nt", group = "RespTmp")

tiff("./plots/clonoDivClus.tiff", width = 5*1500, height = 5*800, res = 300, pointsize = 5)     
clonalDiversity(CD8, cloneCall = "gene+nt", group.by = "RespTmp",
                x.axis = "cluster", n.boots = 100)
dev.off()

alluvialClonotypes(CD8, cloneCall = "gene", 
                   y.axes = c("cluster", "group_id"), 
                   color = "cluster") 

CD8@meta.data$cluster_id <- CD8@active.ident



tiff("./plots/clonoUMAP.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p <- clonalOverlay(CD8, reduction = "umap", 
                   freq.cutpoint = 30, bins = 10, facet = "group_id") + 
  guides(color = FALSE) 
p
dev.off()


library(RColorBrewer)
cols<- brewer.pal(5, "Reds")
pal <- colorRampPalette(cols)

slot(CD8, "meta.data")$cloneType <- factor(slot(CD8, "meta.data")$cloneType, 
                                           levels = c(NA, "Single (0 < X <= 1)", 
                                                      "Small (1 < X <= 5)",
                                                      "Medium (5 < X <= 20)",
                                                      "Large (20 < X <= 100)",
                                                      "Hyperexpanded (100 < X <= 500)"))

table(CD8@meta.data$cloneType, CD8@active.ident)
meta <- melt(table(CD8@meta.data[!is.na(CD8@meta.data$Frequency),
                                 c("functional.cluster", "cloneType")]), varnames = c("functional.cluster", "cloneType"))
ggplot(data = meta, aes(x = CD8@active.ident, y = value, fill = cloneType)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_manual(values = c(palette(5)),
                                                                                            na.value = "grey")


CD8_split <- SplitObject(CD8, split.by = "group_id")
CD8_R <- CD8_split$Res
CD8_NR <- CD8_split$NonRes



mark_res <- FindAllMarkers(CD8_R)
mark_nr <- FindAllMarkers(CD8_NR)
mark_res %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

mark %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

tiff("./plots/DoHeat1.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8_R, features = top20$gene, size = 8, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 14)) + ggtitle("Top 20")
dev.off()

tiff("./plots/DoHeatNR.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8_NR,  group.by = "cloneType", features = top20$gene, size = 8, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 14)) + ggtitle("Top 20")
dev.off()

# get log normalized values for top 50 variable genes
to_plot <- GetAssayData(CD8_R, slot = "data")[top20$gene, ]

cluster_anno<- CD8@meta.data$cloneType


Heatmap(as.matrix(to_plot),
        show_column_names = TRUE)


comb2 <- expression2List(CD8, split.by = "cluster")

length(comb2) 

tiff("./plots/diversity.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
clonalDiversity(comb2, cloneCall = "nt")
dev.off()
library(scRepertoire)

CD8_list <- SplitObject(CD8, split.by = "group_id")

tiff("./plots/occupiedResp.tiff", width = 5*450, height = 5*150, res = 300, pointsize = 5)     
occupiedscRepertoire(CD8_list$Res, x.axis = "cluster") 
dev.off()

tiff("./plots/occupiedNonResp.tiff", width = 5*450, height = 5*150, res = 300, pointsize = 5)     
occupiedscRepertoire(CD8_list$NonRes, x.axis = "cluster") 
dev.off()

#ggraph needs to be loaded due to issues with ggplot
library(ggraph)

clonalHomeostasis(comb2, cloneCall = "nt")
clonalProportion(comb2, cloneCall = "nt")

tiff("./plots/overlap.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
clonalOverlap(comb2, cloneCall="aa", method="overlap")
dev.off()

#Look for epitope specific TCRs
library(Trex)

##Setting up Tensor Flow
library(reticulate)
conda_create("r-reticulate") ##If first time using reticulate
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(tensorflow)

install_tensorflow()

Trex_vectors <- maTrex(CD8, AA.properties = "AF")
qplot(data = as.data.frame(Trex_vectors), Trex_2, Trex_3) + theme_classic()

CD8 <- runTrex(CD8, chains = "TRB", AA.properties = "KF", reduction.name = "Trex.KF")

#Generating UMAP from Trex Neighbors
CD8 <- RunUMAP(CD8, reduction = "Trex.KF", dims = 1:30, reduction.name = 'Trex.umap', 
               reduction.key = 'trexUMAP_')

#Trex UMAP
plot1 <- DimPlot(CD8, reduction = "Trex.umap") + NoLegend()
plot2 <- DimPlot(CD8, group.by = "CTaa", reduction = "Trex.umap") + 
  scale_color_viridis(discrete = TRUE, option = "B") + 
  NoLegend()

plot1 + plot2

#Cluster clones
CD8<- clonalCommunity(CD8, reduction.name = "Trex.KF", cluster.parameter = KNNGraphParam(k=20))

#Run Harmony and new UMAP for RNA
CD8 <- RunHarmony(CD8, "trex.clusters", verbose = FALSE, assay.use = "integrated")
CD8 <- RunUMAP(CD8, reduction = "harmony", dims = 1:20, reduction.key = "RFS_UMAP_", reduction.name = "regressClone")

plot3 <- DimPlot(CD8, 
                 reduction = "umap", 
                 group.by = "CTaa") + 
  NoLegend() + 
  scale_color_viridis(option="B", discrete = TRUE)

plot4 <- DimPlot(CD8, 
                 reduction = "regressClone", 
                 group.by = "CTaa") + 
  NoLegend() + 
  scale_color_viridis(option="B", discrete = TRUE)

plot5 + plot6

CoNGA.seurat <- CoNGAfy(CD8, 
                        method = "dist")

CoNGA.seurat <- runTrex(CoNGA.seurat, 
                        AA.properties = "KF", 
                        reduction.name = "Trex.KF")

CoNGA.seurat <- CoNGA.seurat %>%
  FindNeighbors(reduction = "Trex.KF") %>%
  FindClusters(algorithm = 3)

CoNGA.seurat <- RunUMAP(CoNGA.seurat, 
                        reduction = "Trex.KF", 
                        dims = 1:20, 
                        reduction.name = 'Trex.umap', 
                        reduction.key = 'trexUMAP_')

DimPlot(CoNGA.seurat, reduction = "Trex.umap") + NoLegend()

CoNGA.seurat <- annotateDB(CoNGA.seurat, 
                           chains = "TRB")

tiff("./plots/ep_species.tiff", width = 8*1000, height = 5*600, res = 150, pointsize = 5)  
DimPlot(CoNGA.seurat, 
        reduction = "Trex.umap", 
        group.by = "TRB_Epitope.species")
dev.off()

DimPlot(CoNGA.seurat, 
        reduction = "Trex.umap", 
        group.by = "TRB_Epitope.sequence")

CD8 <- annotateDB(CD8, chains = "TRB")

tiff("./plots/ep_species.tiff", width = 8*1000, height = 5*600, res = 150, pointsize = 5)  
DimPlot(CD8, reduction = "umap", group.by = "TRB_Epitope.species") +theme(element_text(size = 8))
dev.off()

tiff("./plots/UMAP.tiff", width = 8*150, height = 5*150, res = 150, pointsize = 5)  
DimPlot(CD8, reduction = "umap", label = TRUE)
dev.off()