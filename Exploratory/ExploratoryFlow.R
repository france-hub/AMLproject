##Exploratory flow: CD8 This script can be divided into 4 parts

##A)Define/create the directories where our files are located, merge 2 experiments, perform Arcsinh transformation and create flowSet
#1) Define the directory where this script is located
#2) Define the directory where the fcs files, for the two experiments, are located 
#NB: we are using pregated CD8+ fcs files obtained with flowJo software
#3) Define workingDirectory and create the directory
#4) Create flowsets for the two experiments and merge the two flowsets
#5) Perform Arcsinh transformation (#Arcsinh transformation (this is necessary when you have flow citometry data in fcs format instead of csv)
# It is well explained here https://wiki.centenary.org.au/display/SPECTRE/Data+transformation. 
#choose cofactors for each channel)
#6) Define the directory where the transformed files are located and save fcs
#7) Read the transformed files as flowset

##B) Prepare the data to create a single cell experiment (sce) using CATALYST
#1) Create panel dataframe
#2) Create metadata dataframe
#3) Create sce object

##C) Perform QC, clustering and dimensionality reduction
#1) Visualize CATALYST QC plots and remove Ki67 channel because of the very low NRS
#2) Run FlowSOM and ConsensusClusterPlus
#3) Run dimensionality reduction - PCA, UMAP and visualize
#4) Add annotations and visualize

##D) DA analysis using diffcyt
#1) setup model formulas and contrast matrix
#2) Run GLMM and plot

##E) Trajectory inference
#1) Use PCA as dimensionality reduction algorithm and run Slingshot
#2) Visualize trajectories on first 3 PCs
#3) Visualize subsets distribution along PT (jitter plot)
#4) Use UMAP as dimensionality reduction algorithm and run Slingshot
#5) Visualize trajectories onto 2D UMAP
rm(list = ls())

library(rstudioapi)
library(flowCore)
library(FlowSOM)
library(dplyr)
library(flowDensity) 
library(flowStats) 
library(flowVS)
library(ggplot2)
library(dplyr)
library(CATALYST)
library(diffcyt)
library(stringr)
library(SingleCellExperiment)
library(scCustomize)
library(ComplexHeatmap)
library(cowplot)

#####################################################################################################################################
#A)Define/create the directories where our files are located, merge 2 experiments, perform Arcsinh transformation and create flowSet
#####################################################################################################################################

#1) Define the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
PrimaryDirectory <- getwd()
PrimaryDirectory

fcs_1 <- "fcs_1"
FCS1Directory <- paste(PrimaryDirectory, fcs_1, sep = "/")
dir.create(FCS1Directory)

fcs_2 <- "fcs_2"
FCS2Directory <- paste(PrimaryDirectory, fcs_2, sep = "/")
dir.create(FCS2Directory)

#3) Define workingDirectory and create the directory
wdName <- "Working_Directory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#4) Create flowsets for the two experiments and merge the two flowsets
#List FCSfiles_1
FCSfiles_1 <- list.files(FCS1Directory, pattern = ".fcs$", full = FALSE)
## flowSet_1 
fs_1 <- read.flowSet(files = FCSfiles_1, path = FCS1Directory, truncate_max_range = FALSE)
fs_1 <- fs_1[, -c(13,14)] #delete channels 13 and 14 (not included in second flow panel)

#List FCSfiles_2
FCSfiles_2 <- list.files(FCS2Directory, pattern = ".fcs$", full = FALSE)
## flowSet_2 
fs_2 <- read.flowSet(files = FCSfiles_2, path = FCS2Directory, truncate_max_range = FALSE)
## flowSet merged 
fs <- rbind2(fs_1, fs_2)

#5) Perform Arcsinh transformation (#Arcsinh transformation (this is necessary when you have flow citometry data in fcs format instead of csv)
cfs <- c(750, 2200, 2300, 800, 1000, 3200, 1500, 700, 1100, 1100, 600)
#Rename the colnames with the markers corresponding to the fluorochrome (some of this channels, like FSC, SSC,.. are not of interest)
#they do not need to be renamed
colnames(fs)[c(7:9, 11:18)] <- c("KLRG1", "CD45RA","CD27","CD160","PD1","CD28","CD56",
                                 "CD57","CD8","CD3","CCR7")
channels <- colnames(fs)[c(7:9, 11:18)]
#Transform data
fs_t <- transFlowVS(fs, channels = channels, cofactor = cfs)
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)
#Look at the density plots for each channel after transformation
densityplot(~KLRG1, fs_t)
densityplot(~CD45RA, fs_t)
densityplot(~CD27, fs_t)
densityplot(~CD160, fs_t)
densityplot(~PD1, fs_t)
densityplot(~CD28, fs_t)
densityplot(~CD56, fs_t)
densityplot(~CD57, fs_t)
densityplot(~CD8, fs_t)
densityplot(~CD3, fs_t)
densityplot(~CCR7, fs_t)

## Try to correct the signal in the channels with the biggest technical issues
#Warpset from fdaNorm
fs_fda <- warpSet(fs_t, stains = channels[-c(3,4)])
densityplot(~KLRG1, fs_fda)
densityplot(~CD45RA, fs_fda)
densityplot(~CD27, fs_fda)
densityplot(~CD160, fs_fda)
densityplot(~PD1, fs_fda)
densityplot(~CD28, fs_fda)
densityplot(~CD56, fs_fda)
densityplot(~CD57, fs_fda)
densityplot(~CD8, fs_fda)
densityplot(~CD3, fs_fda)
densityplot(~CCR7, fs_fda)

#6) Define the directory where the transformed files are located and save fcs
if(!dir.exists('fcs_t')){dir.create("fcs_t", showWarnings = FALSE)}
setwd("fcs_t")
write.flowSet(fs_fda, outdir='fcs_t', filename = sampleNames(fs_fda))
fcs_t <- "fcs_t"
FCSDirectory <- paste(PrimaryDirectory, fcs_t, sep = "/")

#7) Read the transformed files as flowset
FCSfiles <- list.files(path = FCSDirectory, pattern = ".fcs", full.names= FALSE)
fs <- read.flowSet(files = FCSfiles, path = FCSDirectory, transformation = FALSE, truncate_max_range = FALSE)

# Keyword ($CYT = "FACS"): CATALYST is a mass spectrometry package. 
# We need to change the keyword in order to adapt it to flow cytometry data
#fs[[1]]@description$`$CYT` <- "FACS" this step should not be anymore necessary for recent CATALYST releases
#read here: https://github.com/HelenaLC/CATALYST/issues/103

##########################################################################################################
#B)Prepare the data to create a single cell experiment (sce) using CATALYST
##########################################################################################################

#1) Create panel dataframe
fcs_colname <- colnames(fs)[c(7:9, 11:18)]
marker_class <- c(rep("type", 8), rep("state",2), "type")
antigen <- fcs_colname
length(marker_class) == length(fcs_colname)
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(fs))

#2) Create metadata dataframe
condition <- FCSfiles
condition <- word(condition, 2,4, sep = "_")
condition[grepl("HC", condition)] <- "HC"
condition <- gsub("DG", "base", condition)
patient_id <- FCSfiles
patient_id <- word(patient_id, 1, sep = "_")
sample_id <- paste(patient_id, condition, sep = "_")
file_name <- FCSfiles
md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

#3) Create sce object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE)
sce@assays@data$exprs <- sce@assays@data$counts

########################################################
#C) Perform QC, clustering and dimensionality reduction
########################################################

#1) Visualize CATALYST QC plots and remove Ki67 channel because of the very low NRS
# Density
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Counts
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Multidimensional scaling (MDS)
CATALYST::pbMDS(sce, color_by = "condition", label = NULL)

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

#Use same colors of scRNAseq
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

#2) Run FlowSOM and ConsensusClusterPlus
## Clustering
#FlowSOM
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta7", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last")

#Filter out subset <1%
sce <- filterSCE(sce, cluster_id %in% c(1:4,6), k = "meta7")

#3) Run dimensionality reduction - PCA, UMAP and visualize
#Set parameters for UMAP
set.seed(1234)
n_cells <- 1000
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

#4) Add annotations and visualize
annotation_table <- readxl::read_excel("annotation_1.xlsx")
# convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster,
                                       levels = c("Naive", "StL", "SenL",
                                                  "Act"))
#apply manual annotation
sce <- mergeClusters(sce, k = "meta7", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
#Change condition labels
sce$condition <- ifelse(sce$condition == "CR_BM_base", "Res_bas", 
                        ifelse(sce$condition == "NR_BM_base", "NonRes_bas", 
                               ifelse(sce$condition == "CR_BM_post", "Res_post",
                                      ifelse(sce$condition == "NR_BM_post", "NonRes_post", "HD"))))

#Heatmap with annotations 
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta7", m = "cluster_annotation", scale = "last")
p <- plotExprHeatmap(sce, features = "type",  hm_pal = pal_exp, k_pal = pal_ident,
                by = "cluster_id", k = "cluster_annotation", scale = "last") 

#make the heatmap vertical
cols <- pal_ident[seq_along(levels(annotation_table$new_cluster))]
cols <- setNames(cols, levels(annotation_table$new_cluster))
col_anno <- HeatmapAnnotation(
  df = data.frame(cluster_id = c("Naive", "StL", "SenL", "Act")),
  col = list(cluster_id = cols, gp = gpar(col = "white")),
  annotation_name_gp= gpar(fontsize = 0),
  show_legend = FALSE) 
lgd_aes <- list(direction = "horizontal", legend_width = unit(2.2, "cm"),
                title = "scaled median \n expression")
heat <- Heatmap(t(p@matrix), 
                rect_gp = gpar(col = "white"),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                bottom_annotation = col_anno,
                row_names_gp = grid::gpar(fontsize = 12),
                heatmap_legend_param = lgd_aes,
                column_names_gp = gpar(fontsize = 12))

#Save Fig. 1B
tiff("./plots/heat.tiff", width = 5*150, height = 5*300, res = 300, pointsize = 5)     
draw(heat, heatmap_legend_side = "bottom") 
dev.off()

#Reorder UMAP and change labeling
p <- plotDR(sce, "UMAP", facet_by = "condition", color_by = "cluster_annotation", k_pal = pal_ident, ncol = 2) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) + 
  geom_density2d(binwidth = 0.006, colour = "black") +
  theme(axis.title=element_text(size=0),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank())

splitFacet <- function(x){
  facet_vars <- names(x$facet$params$facets)         
  x$facet    <- ggplot2::ggplot()$facet 
  datasets   <- split(x$data, x$data[facet_vars])    
  new_plots  <- lapply(datasets,function(new_data) { 
    x$data <- new_data
    x})
}    

df <- p$data
df <- aggregate(
  x = df[c("x", "y")], 
  by = list(z = df$cluster_annotation, w = df$condition), 
  FUN = mean)
df.split <- split(df, df$w)
titles <- as.list(levels(as.factor(df$w)))
p.list <- splitFacet(p)
p.list <- lapply(seq_along(p.list), function(x) p.list[[x]] + 
         geom_label(
           aes(x, y, label = z),
           data = df.split[[x]], show.legend = FALSE,
           inherit.aes = FALSE,
           fill = NA,
           label.size = NA,
           size = 10) + ggtitle(titles[[x]]) +
           theme(plot.title = element_text(size = 27, hjust = 0.5, face = "bold"))) 
p.listHD <- p.list[[1]]
p.listAML <- p.list
p.listAML[[1]] <- NULL

#Save Fig. 1C
tiff("./plots/umap_AML.tiff", width = 5*700, height = 5*700, res = 300, pointsize = 5)     
plot_grid(p.listAML[[1]],p.listAML[[2]], p.listAML[[3]],p.listAML[[4]], nrow = 2)
dev.off()

#Save Supplementary Fig. 1C
tiff("./plots/umap_HD.tiff", width = 5*350, height = 5*350, res = 300, pointsize = 5)     
p.listHD 
dev.off()

p <- plotAbundances(sce, k = "cluster_annotation")
df <- p$data

barplot <- ggplot(df, aes(x = condition, y = Freq, fill = cluster_id)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ylab("") + xlab("") + scale_fill_manual(values = pal_ident) + theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm")) 

#Save Fig. 1D
tiff("./barplot.tiff", width = 5*200, height = 5*200, res = 300, pointsize = 5)     
barplot
dev.off()

saveRDS(sce, "sce.rds") #FROM HERE
sce <- readRDS("sce.rds")
