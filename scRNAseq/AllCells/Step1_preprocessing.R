#This script is on preprocessing of scRNAseq data. Later on I'll use Seurat to analyze my data. Here though I use singlecellexperiment because
#it allows me to use the doubletfinder package (https://github.com/chris-mcginnis-ucsf/DoubletFinder). There is now an interface of the same package
#that works with Seurat, but when I did this I do not think was available. 
# In summary the script does the followings:
#1) rename the genes
#2) Add metadata
#3) Remove undetected genes
#3) Infer and remove doublets
#4) Calculate QC metrics with diagnostic plots
#5) Remove outliers

#Nice guide: https://bioconductor.org/books/release/OSCA/

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

#Set working directory where the script is located
setwd("~/Documents/AML_project/scRNA_AMLproj/scripts")

#Load raw counts for all the experiments (each dir contains the files barcodes.tsv, features.tsv, matrix.mtx for each patient)
dirs <- list.dirs("data_scRNAseq", recursive = FALSE, full.names = TRUE)
names(dirs) <- basename(dirs)
sce <- read10xCounts(dirs)

#Check sce dimension
dim(sce)

#rename row/colData colnames & SCE dimnames the genes will be in the format ensembl.symbol
#this will allow to keep as rows the info of ENSEMBLs and SYMBOLs
rowData(sce) <- rowData(sce)[,1:2]
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), SYMBOL), 
  with(colData(sce), paste(sample_id, barcode, sep = "_")))

#load metadata and add to sce object
md <- file.path("../md_dir", "metadata.xlsx")
md <- read_excel(md)
m <- match(sce$sample_id, md$sample_id)
sce$group_id <- md$group_id[m]
sce$timepoint <- md$timepoint[m]
sce$RespTmp <- md$RespTmp[m]
sce$ELN <- md$ELN[m]
sce$batch <- md$Batch[m]
colData(sce)

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#Doublets removal (reference: https://github.com/chris-mcginnis-ucsf/DoubletFinder)
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

#Calculate QC metrics QC metrics like number of unique genes detected per cell and total number of reads
#Passing in the list of mitochondrial genes, we will also calculate percent mitochondrial reads.
(mito <- grep("MT-", rownames(sce), value = TRUE)) #find mitochondrial genes
sce <- addPerCellQC(sce, subsets = list(Mt = mito)) 

##Visualize colData(sce) structure
plotHighestExprs(sce, n = 20)

##Filtering
#1)sum contains the total count for each cell 
#2)detected column contains the number of detected genes 
#3)subsets_Mt_percent column contains the percentage of reads mapped to mitochondrial transcripts

# get sample-specific outliers and plot
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

# Look at the filtered and unfiltered heatscatter
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

#Dimension sce after filtering
dim(sce)

saveRDS(sce, file = "scRNAseq_step1_test.rds")
