##Spectral flow: CD8 This script can be organized into n parts

##A)Define/create the directories where our files are located, transform the csv in fcs and create flowSet
#1) Define the directory where this script is located
#2) Define the directory where the csv files are located 
#NB: we are using pregated CD8+ csv files obtained with the flowJo csv channel values option and then transform them here 
#in fcs files so that we won't need to arcsinh transform the data as explained here https://wiki.centenary.org.au/display/SPECTRE/Data+transformation)
#3) Define and create the WorkingDirectory (here we'll save the plots) and the directory that will contain the fcs files
#4) Transform the csv in fcs, save and create the flowSet

##B) Prepare the data to create a single cell experiment (sce) using CATALYST
#1) Create panel dataframe
#2) Create metadata dataframe
#3) Visualize density plots and use warpSet to normalize 
#4) Create sce object

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

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2) 
  library(flowCore) 
  library(cytofCore)
  library(dplyr)
  library(devtools)
  library(FlowSOM)
  library(cluster)
  library(dplyr)
  library(ggthemes)
  library(RColorBrewer)
  library(uwot)
  library(CATALYST)
  library(diffcyt)
  library(stringr)
  library(scran)
  library(scater)
  library(SingleCellExperiment)
  library(flowWorkspace)
  library(reshape2)
  library(ggrepel)
  library(slingshot)
  library(knn.covertree)
  library(readxl)
  library(flowVS)
  library(flowStats)
  library(ggpubr)
  library(scCustomize)
  library(patchwork)
})

##########################################################################################################
#A)Define/create the directories where our files are located, transform the csv in fcs and create flowSet
##########################################################################################################

#1) Define the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

#2) Define the directory where the csv files are located 
csv <- "csv" #name of the directory
csvDirectory <- paste(PrimaryDirectory, csv, sep = "/")

#3) Define and create the WorkingDirectory (here we'll save the plots) and the directory that will contain the fcs files
# Define working directory
wdName <- "WorkingDirectory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#Define csv to fcs directory
csv2fcsDirectory <- "csv2fcsDirectory"
csv2fcsDirectory <- paste(PrimaryDirectory, csv2fcsDirectory, sep = "/")
dir.create(csv2fcsDirectory)

#4) Transform the csv in fcs, save and create the flowSet
#List CSV files and add .fcs extension
CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)
#Obtain fcs files from csv files
for(i in c(1:length(CSVfiles))){
  data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
  print(CSVfiles[i])
  print(fileName[i])
  cytofCore.write.FCS(as.matrix(data), 
                      filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                      what = "numeric")
}
# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
fs <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)

#fs[[1]]@description$`$CYT` <- "FACS" this step should not be anymore necessary for recent CATALYST releases
#read here: https://github.com/HelenaLC/CATALYST/issues/103

#######################################################################
#B) Prepare the data to create a single cell experiment using CATALYST
#######################################################################
#1) Create panel dataframe
fcs_colname <- colnames(fs)[10:33] # Define channels of interest and marker_class
antigen <- fcs_colname
marker_class <- c("state", "type", "state", rep("type", 20), "state") #define types and states as per CATALYST workflow
length(marker_class) == length(fcs_colname)
panel <- as.data.frame(cbind(fcs_colname, antigen, marker_class))
all(panel$fcs_colname %in% colnames(fs)) #check

#Set conditions
condition <- FCSfiles
condition <- word(condition, 2,3, sep = "_") 

#Set patient_id
patient_id <- FCSfiles
patient_id <- word(patient_id, 1, sep = "_")

#Set sample_id
sample_id <- paste(patient_id, condition, sep = "_")

#Set file_name
file_name <- FCSfiles

#2) Create metadata dataframe
md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)
#Check if ids and md$file_name are the same
ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

#3) Visualize density plots and use warpSet to normalize 
tiff("./WorkingDirectory/density1.tiff", width = 5*700, height = 5*1200, res = 300, pointsize = 5)     
densityplot(~., fs, channels = antigen[1:6]) #CD8, GZMK, GZMB
dev.off()

tiff("./WorkingDirectory/density2.tiff", width = 5*700, height = 5*1200, res = 300, pointsize = 5)     
densityplot(~., fs, channels = antigen[7:12])
dev.off()

tiff("./WorkingDirectory/density3.tiff", width = 5*700, height = 5*1200, res = 300, pointsize = 5)     
densityplot(~., fs, channels = antigen[13:18])
dev.off()

tiff("./WorkingDirectory/density4.tiff", width = 5*700, height = 5*1200, res = 300, pointsize = 5)     
densityplot(~., fs, channels = antigen[19:24])
dev.off()

#Normalize 
fs <- warpSet(fs, stains = antigen[c(3,4,5,17,18)])

#4) Create sce object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, 
                          features = panel$fcs_colname)

sce@assays@data$exprs <- sce@assays@data$counts
type_markers(sce) #visualize type markers
#save(sce, file = "Spectral_test_CD8.rds")
#load("Spectral_test_CD8.rds")

#Rename conditions to make them consistent with scRNAseq labeling
sce$condition <- ifelse(sce$condition == "CR_base", "Res_bas", 
                        ifelse(sce$condition == "NR_base", "NonRes_bas", 
                               ifelse(sce$condition == "CR_post", "Res_post", "NonRes_post")))
##############################
#C) Perform QC and clustering
##############################
#1) Visualize CATALYST QC plots and remove Ki67 channel because of the very low NRS
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")
# Pseudo-bulk multidimensional scaling (pbMDS)
CATALYST::pbMDS(sce, color_by = "patient_id")
# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")
# Filter out Ki67 (very low NRS)
sce <- filterSCE(sce, rownames(sce) != "Ki67")
rownames(sce)
#Use same colors of scRNAseq
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

#2) Run FlowSOM and ConsensusClusterPlus
set.seed(1234)
sce <- cluster(sce, xdim = 10, ydim = 10, maxK = 20, seed = 1234)
tiff("./WorkingDirectory/delta.tiff", width = 5*400, height = 5*200, res = 300, pointsize = 5)     
delta_area(sce)
dev.off()
tiff("./WorkingDirectory/heat10.tiff", width = 5*400, height = 5*200, res = 300, pointsize = 5)     
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta8", k_pal = pal_ident,  hm_pal = pal_exp,
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last")
dev.off()

#3) Run dimensionality reduction - PCA, UMAP and visualize
n_cells <- 3000
n_events <- min(n_cells(sce))
if(!(n_cells < n_events))
  n_cells <- n_events
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")
sce <- runDR(sce, dr = "PCA", cells = n_cells, features = "type")

# Plot UMAP
#metaclusters distribution
umap.meta8 <- plotDR(sce, "UMAP", color_by = "meta20", k_pal = pal_ident) + theme_void()
p.coord <- ggplot(data.frame(x = 100, y = 100), aes(x = x, y = y)) +
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
umap.meta8 <- umap.meta8 + p.coord + plot_layout(design = layout)
tiff("./WorkingDirectory/umap_8.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
umap.meta8
dev.off()

indx.mark <- match(c("CCR7", "TCF1", "CD127", "GZMK", "CD69", "CD57", 
                     "GZMB", "CX3CR1", "DNAM1"), type_markers(sce))

tiff("./WorkingDirectory/umap_markers.tiff", width = 5*500, height = 5*500, res = 300, pointsize = 5)     
plotDR(sce, "UMAP", color_by = type_markers(sce)[indx.mark], a_pal = pal_exp) +
  theme(strip.text = element_text(size=25)) +  xlab("UMAP_1") + ylab("UMAP_2") +
  theme(strip.text = element_text(size=25), axis.text=element_text(size=12),
        axis.title=element_text(size=14)) 
dev.off()

png("./WorkingDirectory/umap_clus.png", width = 5*900, height = 5*900, res = 300, pointsize = 5)     
plotDR(sce, "UMAP", color_by = "meta8", facet_by = "group_id", k_pal = pal_ident) +
  geom_point(size = 0.0000001)+
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) + 
  geom_density2d(binwidth = 0.006, colour = "black") +  
  xlab("UMAP_1") + ylab("UMAP_2") +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25)) + 
  theme(legend.text=element_text(size=30)) 
dev.off()

plotDR(sce, "UMAP", color_by = "meta8", facet_by = "sample_id") #by sample_id

#4) Add annotations and visualize
#Read annotation file
annotation_table <- readxl::read_excel("annotation.xlsx")
annotation_table$new_cluster <- factor(annotation_table$new_cluster, levels = c("Naive",
                                                                                "StL","SL", "ActEx",  "Int"))

#Apply manual annotation
sce <- mergeClusters(sce, k = "meta8", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

#Plot heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta8", m = "cluster_annotation", scale = "last", bin_anno = TRUE, perc = TRUE)

tiff("./WorkingDirectory/heat_def.tiff", width = 5*500, height = 5*250, res = 300, pointsize = 5)     
plotExprHeatmap(sce, features = "type",  hm_pal = pal_exp, k_pal = pal_ident,
                by = "cluster_id", k = "cluster_annotation", scale = "last") 
dev.off()

#Plot UMAP
tiff("./WorkingDirectory/umap_cond.tiff", width = 5*900, height = 5*900, res = 300, pointsize = 5)     
plotDR(sce, "UMAP", color_by = "cluster_annotation", facet_by = "condition", k_pal = pal_ident) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) + 
  geom_density2d(binwidth = 0.006, colour = "black") +  
  xlab("UMAP_1") + ylab("UMAP_2") +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25),
        legend.text=element_text(size=30),
        legend.title = element_text(size=0)) 
dev.off()

#Plot boxplot
tiff("./WorkingDirectory/box20.tiff", width = 5*300, height = 5*180, res = 300, pointsize = 5)     
p <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id")
pal_box <- pal_ident[c(6,19, 14, 35)]
p.box <- ggplot(p$data, aes(x = condition, y = Freq, fill = condition)) +
  labs(x = NULL, y = "Proportion") + 
  geom_boxplot() + scale_fill_manual(values= pal_box) +
  geom_point(aes_string(x = "condition"), position = position_jitter(width = 0.2)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    legend.key.height  =  unit(0.8, "lines")) +
  facet_wrap(~ cluster_id, scales = "free_y", nrow = 1) + theme(legend.position = "none")
p.box
dev.off()

##############################
#D) DA analysis using diffcyt
##############################

#1) setup model formulas and contrast matrix
ei <- sce@metadata$experiment_info
ei <- ei %>%  mutate(group_id = case_when(grepl("NR", sample_id) ~ "NonRes",
                                          grepl("CR", sample_id) ~ "Res")) 

ei$group_id <- factor(ei$group_id)

#Change condition names
ei$condition <- factor(ifelse(ei$condition == "CR_base", "Res_bas", 
                              ifelse(ei$condition == "NR_base", "NonRes_bas", 
                                     ifelse(ei$condition == "CR_post", "Res_post", "NonRes_post"))))

(da_formula <- createFormula(
  ei, cols_fixed = "group_id", cols_random = "patient_id"))

#Contrast Res vs NonRes
contrast <- createContrast(c(0,1))

#2) Run GLMM and plot
out_DA1 <- diffcyt(sce, formula = da_formula, contrast = contrast, method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "cluster_annotation")
da_1 <- topTable(out_DA1, format_vals = TRUE)

tiff("./WorkingDirectory/heatmapSig.tiff", width = 5*900, height = 5*300, res = 300, pointsize = 5)     
plotDiffHeatmap(sce, rowData(out_DA1$res), all = TRUE)
dev.off()
pal_box <- pal_ident[c(14, 35)]
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "group_id")
df.plot <- bxp$data %>% dplyr::mutate(positions = case_when(cluster_id == "StL" ~ "*", TRUE ~ ""))
bxp <- ggplot(df.plot, aes(x = group_id, y = Freq, fill = group_id)) +
  labs(x = NULL, y = "Proportion") + 
  geom_boxplot() + scale_fill_manual(values= pal_box) +
  geom_point(aes_string(x = "group_id"), position = position_jitter(width = 0.2)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
    legend.key.height  =  unit(0.8, "lines")) +
  facet_wrap(~ cluster_id, scales = "free_y", nrow = 1) + theme(legend.position = "none")
bxp <- bxp + geom_text(aes(x=1.5,y=35,label=positions), size = 5)
tiff("./WorkingDirectory/BoxplotSig.tiff", width = 5*600, height = 5*150, res = 300, pointsize = 5)     
bxp
dev.off()

########################
#E Trajectory inference
########################

#1) Use PCA as dimensionality reduction algorithm and run Slingshot
sce_sling <- filterSCE(sce, complete.cases(reducedDim(sce, "PCA"))) #remove NA
#Prepare clusterLabels
clusters <- sce_sling$cluster_id
levels(clusters) <- cluster_codes(sce_sling)$cluster_annotation
#Run slingshot
sds <- slingshot(sce_sling, reducedDim = "PCA", clusterLabels = clusters, 
                 start.clus = "Naive", stretch = 0)
df <- bind_cols(
  setNames(as.data.frame(reducedDim(sds, "PCA")), c("PCA_1", "PCA_2", "PCA_3")),
  slingPseudotime(sds) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(sds) %>% as.data.frame(),
) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst))

curves <- slingCurves(sds, as.df = TRUE)
curves <- setnames(curves, old = c('Dim.1','Dim.2', 'Dim.3'), new = c('PCA_1','PCA_2', 'PCA_3'))

#2) Visualize trajectories on first 3 PCs
pca.pseudo <- ggplot(df, aes(x = PCA_1, y = PCA_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 680, y = 600, label = "SL", size = 5) +
  annotate("text", x = 30, y = 20, label = "ActEx", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) + theme_classic()

df$clusters <- cluster_ids(sce_sling, "cluster_annotation")
pal_ident <- DiscretePalette_scCustomize(num_colors = 40, palette = "ditto_seq")

p1.pca <- ggplot(df, aes(x = PCA_1, y = PCA_2)) +
  geom_point(size = 0.0000001, aes(col = clusters)) +
  scale_colour_manual(values = pal_ident) +
  labs(col = "") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 680, y = 600, label = "SL", size = 5) +
  annotate("text", x = 30, y = 20, label = "ActEx", size = 5) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25)) + 
  theme(legend.text=element_text(size=30)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  theme_classic()

p2.pca <- ggplot(df, aes(x = PCA_1, y = PCA_3)) +
  geom_point(size = 0.0000001, aes(col = clusters)) +
  scale_colour_manual(values = pal_ident) +
  labs(col = "") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 550, y = 500, label = "SL", size = 5) +
  annotate("text", x = 30, y = -600, label = "ActEx", size = 5) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25)) + 
  theme(legend.text=element_text(size=30)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  theme_classic()

p3.pca <- ggplot(df, aes(x = PCA_2, y = PCA_3)) +
  geom_point(size = 0.0000001, aes(col = clusters)) +
  scale_colour_manual(values = pal_ident) +
  labs(col = "") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 500, y = 450, label = "SL", size = 5) +
  annotate("text", x = 30, y =-550, label = "ActEx", size = 5) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25)) + 
  theme(legend.text=element_text(size=30)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  theme_classic()

p.pca.traj <- cowplot::plot_grid(pca.pseudo, p1.pca, p2.pca, p3.pca, ncol = 2)

tiff("./WorkingDirectory/p.pca.traj.tiff", width = 5*300, height = 5*200, res = 150, pointsize = 5)     
p.pca.traj
dev.off()

#3) Visualize subsets distribution along PT (jitter plot)
p.jitt1 <- ggplot(df, 
                  aes(x = Lineage1_pst, 
                      y = clusters, colour = clusters)) + geom_jitter(size = 0.0000001) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_colour_manual(values = pal_ident) +
  theme_classic() +
  xlab("Pseudotime (SL)") +
  ylab("") + labs(col = "") +theme(legend.text=element_text(size=15))


p.jitt2 <- ggplot(df, 
                  aes(x = Lineage2_pst, 
                      y = clusters, colour = clusters)) + geom_jitter(size = 0.0000001) +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_colour_manual(values = pal_ident) +
  theme_classic() +
  xlab("Pseudotime (ActEx)") +
  ylab("") + labs(col = "") + theme(legend.text=element_text(size=15))

tiff("./WorkingDirectory/jitt.tiff", width = 5*300, height = 5*100, res = 150, pointsize = 5)     
cowplot::plot_grid(p.jitt1, p.jitt2)
dev.off()

#4) Use UMAP as dimensionality reduction algorithm and run Slingshot
sce_sling <- filterSCE(sce, complete.cases(reducedDim(sce, "UMAP"))) #remove NA
#Prepare clusterLabels
clusters <- sce_sling$cluster_annotation
levels(clusters) <- cluster_codes(sce_sling)$cluster_annotation
#Run slingshot
sds <- slingshot(sce_sling, reducedDim = "UMAP", clusterLabels = clusters, 
                 start.clus = "Naive", stretch = 0)
df.umap <- bind_cols(
  setNames(as.data.frame(reducedDim(sds, "UMAP")), c("UMAP_1", "UMAP_2")),
  slingPseudotime(sds) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(sds) %>% as.data.frame(),
) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst))
curves <- setNames(slingCurves(sds, as.df = TRUE), c("UMAP_1", "UMAP_2", "Order", "Lineage"))

#5) Visualize trajectories onto 2D UMAP
umap.traj <- ggplot(df.umap, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 6, y = 8, label = "ActEx", size = 5) +
  annotate("text", x = 6, y = -7, label = "SL", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +  theme_void()  

umap.traj.ann <- ggplot(df.umap, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.0000001, aes(col = clusters)) +
  scale_colour_manual(values = pal_ident) +
  labs(col = "") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black",  arrow = arrow(), lineend = "round", size = 1.5) +
  annotate("text", x = 6, y = 8, label = "ActEx", size = 5) +
  annotate("text", x = 6, y = -7, label = "SL", size = 5) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) +
  theme(strip.text = element_text(size=40), axis.text=element_text(size=12),
        axis.title=element_text(size=25)) + 
  theme(legend.text=element_text(size=30)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  theme_void()

p.traj <- cowplot::plot_grid(umap.traj, umap.traj.ann, ncol = 1)

tiff("./WorkingDirectory/umapTraj.tiff", width = 5*100, height = 5*150, res = 150, pointsize = 5)     
p.traj + labs(x = "UMAP_1", y = "UMAP_2")
dev.off()