#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# Analyze single cell RNA-seq: TNBC_scRNA_GSE169246
# Aim: Cell typing based spatial data
# 10/02/2023
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# preprocess data ##############################################################
# Clear the environment
rm(list = ls())

# Load the libraries
library(Seurat)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# The input directory
# Please download the data from 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246
# And store the data in the following directory
myinf1 <- "./data/TNBC_scRNA_GSE169246/RNA"

# THe output directory
myoutd <- "./Fig1/result/TNBC_scRNA_GSE169246/"
dir.create(myoutd)

myoutf1 <- paste0(myoutd, "QC.pdf")
myoutf2 <- paste0(myoutd, "PC.pdf")
myoutf3 <- paste0(myoutd, "UMAP_cluster.pdf")
myoutf4 <- paste0(myoutd, "UMAP_patient.pdf")
myoutf5 <- paste0(myoutd, "UMAP_timepoint.pdf")
myoutf6 <- paste0(myoutd, "Spatial_marker_snn0.8_dotplot.pdf")
myoutf7 <- paste0(myoutd, "Spatial_marker_snn0.8_featureplot.pdf")
myoutf8 <- paste0(myoutd, "NKT_marker_snn1_dotplot.pdf")
myoutf9 <- paste0(myoutd, "NKT_marker_snn1_featureplot.pdf")
myoutf10 <- paste0(myoutd, "All_NKT_UMAP.pdf")
myoutf11 <- paste0("./data/TNBC_scRNA_GSE169246/", "celltyping_all.rds")
myoutf12 <- paste0("./data/TNBC_scRNA_GSE169246/", "celltyping_NKT.rds")

# Load data --------------------------------------------------------------------
tenx.raw <- Read10X(data.dir = myinf1, gene.column = 1,cell.column = 1)
tenx.data <- CreateSeuratObject(counts = tenx.raw,
                                min.cells=3,
                                min.features=200,
                                project="TNBC_aPDL1")

#  filter & QC  --------------------------------------------------------------
mydata <- tenx.data
# Clear the large object
rm(tenx.raw, tenx.data)


mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-|^mt-")
mydata@meta.data$Sample <- rownames(mydata@meta.data) %>%
  strsplit(split = ".",fixed = T) %>%
  lapply("[[", 2) %>%
  unlist()
unique(mydata@meta.data$Sample )


mydata@meta.data$TimePoint <-  mydata@meta.data$Sample %>%
  strsplit(split = "_",fixed = T) %>%
  lapply("[[", 1) %>%
  unlist()

unique(mydata@meta.data$TimePoint)

mydata@meta.data$Patient <- mydata@meta.data$Sample %>%
  strsplit(split = "_",fixed = T) %>%
  lapply("[[", 2) %>%
  unlist()


unique(mydata@meta.data$Patient)

mydata@meta.data$Tissue <- mydata@meta.data$Sample %>%
  strsplit(split = "_",fixed = T) %>%
  lapply("[[", 3) %>%
  unlist()

unique(mydata@meta.data$Tissue)


# Only keep TNBC tumor samples
comxx <- mydata@meta.data %>%
  filter(TimePoint == "Pre" | TimePoint == "Post") %>%
  filter(Tissue == "t") %>%
  pull(Sample)

Idents(mydata) <- "Sample"
filtered.data <- subset(x = mydata, idents = comxx)

# QC figures
plot.nF<- VlnPlot(filtered.data, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)
plot.nRead <- VlnPlot(filtered.data, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
plot.pmito <- VlnPlot(filtered.data, features = "percent.mt", pt.size = 0.1, log = TRUE)


plot.viability <- FeatureScatter(object = filtered.data,
                                 feature1 = "nCount_RNA",
                                 feature2 = "percent.mt")
plot.doublets <- FeatureScatter(object = filtered.data,
                                feature1 = "nCount_RNA",
                                feature2 = "nFeature_RNA")
gp.list <- list(plot.nF,
                plot.nRead,
                plot.pmito,
                plot.viability,
                plot.doublets)
combined.gp <- do.call(ggarrange, c(gp.list, ncol = 1, nrow = 1))

print('...Print Figures...')
pdf(myoutf1, width = 20, height = 10)
print(combined.gp,height = 10,
      width = 20)
dev.off()


range(filtered.data@meta.data$nFeature_RNA) # 398 7847
range(filtered.data@meta.data$percent.mt) #  0.000000 9.996214
hist(filtered.data@meta.data$percent.mt)


# This step is redundant, the author has filtered the data based above QC
# I just to double check to confirm the data is clean
filtered.data <- subset(x = filtered.data,
                        subset = nFeature_RNA > 200 &
                        nFeature_RNA < 8000 &
                        percent.mt <= 10)
                        
# remove the large object
rm(mydata)
# Remove P10 and P16 with poor quality
Idents(filtered.data) <- "Patient"
filtered.data <- subset(x = filtered.data, subset = Patient != "P010" & Patient != "P016")
unique(filtered.data@meta.data$Patient)

# Cell typing ##############################################################
#  normalizing -----------------------------------------------------------------

# normalize
normalized.data <- NormalizeData(object = filtered.data,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 100000)
dim(normalized.data@meta.data)
# remove the large object
rm(filtered.data)

# No batch correction Clustering -----------------------------------------------
# determine number of PCs to use
seed =1234
normalized.data <- FindVariableFeatures(object = normalized.data,
                                        selection.method = 'vst',
                                        nfeatures = 2000)
length(x = VariableFeatures(object = normalized.data))
hvg <- VariableFeatures(object = normalized.data)
normalized.data <- ScaleData(object = normalized.data, features = hvg,
                             vars.to.regress = c("nCount_RNA", "percent.mt"))
normalized.data <- RunPCA(object = normalized.data, features = hvg,seed.use=seed,
                          verbose = FALSE)

pdf(myoutf2)
VizDimLoadings(normalized.data, dims = 1:2)
DimPlot(normalized.data)
normalized.data <- ProjectDim(object = normalized.data)
DimHeatmap(normalized.data, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(object =normalized.data,ndims = 50)
dev.off()
# snn: select PC1:30

normalized.data <- FindNeighbors(object = normalized.data, dims = 1:30,k.param = 20)
for (i in seq(0.1,0.8,0.1)) {
  normalized.data <- FindClusters(object = normalized.data,resolution = i)
}


# uMAP visualization

normalized.data  <- RunUMAP(object = normalized.data,
                            reduction = "pca",
                            dims = 1:30,
                            min.dist = 0.8,
                            n.neighbors = 30,
                            seed.use = seed)
pdf(myoutf3, width = 18,height = 10)
DimPlot(normalized.data , reduction = "umap",
        group.by = c(paste0("RNA_snn_res.",seq(0.1,0.8,0.1)),"TimePoint","Patient"))
dev.off()

pdf(myoutf4, width = 50,height = 5)
DimPlot(normalized.data , reduction = "umap",
        group.by = c("RNA_snn_res.0.5"), label = T,
        split.by = "Patient")
dev.off()

# Patient has difference but considering biology > batch

pdf(myoutf5)
DimPlot(normalized.data , reduction = "umap",group.by = c("RNA_snn_res.0.5","TimePoint"),
        split.by = "TimePoint")
dev.off()
# TimePoint has a little difference, considering biology > batch

# Cell Typing 1: all -----------------------------------------------------------
print("-------Cell Typing----------")
warnings("---------Cell Typing----------")

normalized.data <- readRDS(myinf1)

# spatial markers:
# Figure 1g: https://www.nature.com/articles/s41586-023-06498-3/figures/1
# the single cell is immune cell atlas so we only focused on TME
# Stroma cell markers should have not expr but included for plotting,
# In order to consistent with the spatial cell typing
ref <- c( "CD8A", "CD8B", # "CD8"
          "TCF7", # "TCF1"
          "TBX21", # "T-bet"
          "CD3E",  "CD3D", "CD3G" , # "CD3"
          "GATA3", # "GATA3"
          "GZMK", "GZMA", "GZMH", "GZMB", "GZMM", # "GZMB"
          # GZMK NK + (CD8+ T); GZMA CD8+ CTL + NK; GZMH NK + CD8+ T
          # GZMB NK + (CD8+ T); GZMM NK
          "ITGAX", # "CD11c", DC
          "CD68", # "CD68"
          "CD163", #  "CD163"
          "CD4"  , # "CD4"
          "TOX",  "TOX2",  "TOX3", "TOX4" , # "TOX"
          "PDCD1", # "PD-1"
          "TNFRSF4", #"OX40
          "HLA-A"   , "HLA-B" , "HLA-C"   ,# "HLA-ABC"
          "PTPRC", #"CD45"
          "HLA-DRA" ,  "HLA-DRB5" ,    "HLA-DRB1"  ,# "HLA-DR"
          "CD274", # "PD-L1"
          "IDO1", #  "IDO1"
          "ICOS", # "ICOS"
          "FOXP3", # "FOXP3"
          "IKZF2", # "helios"
          "NCAM1", # "CD56"
          "CA9" , # "CA9"
          "PDPN", # "PDPN"
          "CNN1", # "Calponin"
          "MS4A1", #"CD20"
          "CD79A", # "CD79a"
          "MPO", # MPO"
           "FUT4", # "CD15"
          "CAV1", # "Caveolin-1"
          "SMN1", # "SMA"
          "VIM", # "Vimentin"
          "PDGFRB" , # "PDGFRB"
          "PECAM1", # "CD31",
          "NCR1", "GNLY", "NKG7", # Additional NK Cell markers
          "CLEC4C", # Additional pDC
          "KIT","MS4A2" # Additional Mast )
          )
# Check the gene symbol
# features <- rownames(normalized.data@assays$RNA)
# features[grep("FCGR3A" ,features)]
normalized.data@meta.data$RNA_snn_res.0.8 <- factor(normalized.data@meta.data$RNA_snn_res.0.8, levels = 0:30)
Idents(normalized.data) <- "RNA_snn_res.0.8"

pdf(myoutf6, width = 15, height = 10)
DotPlot(object =normalized.data , cols = c("blue","red"),
        feature= ref) + RotatedAxis()
dev.off()

pdf(myoutf7, width = 24, height = 30)
FeaturePlot(normalized.data, features = ref,ncol = 6,
            max.cutoff = 4, label = T)
dev.off()


# Find markers + enrichr to assist cell typeing
# plan("multisession", workers = 20)
# options(future.globals.maxSize = 8000 * 1024^2)
# cluster.markers <- FindMarkers(normalized.data, ident.1 = 27,
#                                #ident.2 = c(12,17,21,22,23,24,25),
#                                min.pct = 0.25)
# head(cluster.markers, n = 20)

# Level 1 cell typing
# The cells' name is consistent with the spatial data
celltype1 <- c("NK & T",  # 0
               "NK & T",  # 1
               "CD20^+B",  # 2
               "CD20^+B",  # 3
               "NK & T",  # 4
               "NK & T",  # 5
               "NK & T",  # 6
               "NK & T",  # 7
               "NK & T",  # 8
               "NK & T",  # 9
               "CD79a^+Plasma",  # 10
               "NK & T",  # 11
               "M1 Mac",  # 12
               "NK & T",  # 13
               "NK & T",  # 14
               "NK & T",  # 15
               "CD20^+B",  # 16
               "M2 Mac",  # 17
               "CD79a^+Plasma",  # 18
               "NK & T",  # 19
               "CD20^+B",  # 20
               "M2 Mac",  # 21
               "NK & T",  # 22
               "M2 Mac",  # 23
               "DCs",  # 24
               "pDC",  # 25
               "CD20^+B",  # 26
               "Mast",  # 27
               "CD79a^+Plasma",  # 28
               "NK & T",  # 29
               "CD20^+B"  # 30
               )
names(celltype1) <- levels(normalized.data)
normalized.data <- RenameIdents(normalized.data, celltype1)
normalized.data@meta.data$CellType1 <- normalized.data@active.ident
DimPlot(normalized.data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Cell Typing 2: NK T -----------------------------------------------------------
seed <-  1234
NKT <- subset(normalized.data, idents = "NK & T")

NKT <- FindNeighbors(object = NKT, dims = 1:30,k.param = 20)
NKT <- FindClusters(object = NKT,resolution = 1)


# uMAP visualization

NKT  <- RunUMAP(object = NKT,
                reduction = "pca",
                dims = 1:30,
                min.dist = 0.8,
                future.seed = T,
                n.neighbors = 30,
                seed.use = seed)

DimPlot(NKT , reduction = "umap", label = T,
        group.by = c("RNA_snn_res.1","TimePoint","Patient"))

NKT_markers <- c( "CD8A", "CD8B", # "CD8"
                  "TCF7", # "TCF1"
                  "TBX21", # "T-bet"
                  "CD3E",  "CD3D", "CD3G" , # "CD3"
                  "GATA3", # "GATA3"
                  "GZMK", "GZMA", "GZMH", "GZMB", "GZMM", # "GZMB"
                  "CD4"  , # "CD4"
                  "TOX",
                  "PDCD1", # "PD-1"
                  "TNFRSF4", #"OX40
                  "CD274", # "PD-L1"
                  "IDO1", #  "IDO1"
                  "ICOS", # "ICOS"
                  "FOXP3", # "FOXP3"
                  "IKZF2", # "helios"
                  "NCAM1", # "CD56"
                  "NCR1", "GNLY", "NKG7", # Additional NK Cell markers
                  "CTLA4", "LAG3", "TIGIT" # additional exhausted markers
)

NKT@meta.data$RNA_snn_res.1 <- factor(NKT@meta.data$RNA_snn_res.1, levels = 0:21)
Idents(NKT) <- "RNA_snn_res.1"

pdf(myoutf8, width = 15, height = 10)
DotPlot(object =NKT , cols = c("blue","red"),
        feature= NKT_markers) + RotatedAxis()
dev.off()

pdf(myoutf9, width = 24, height = 30)
FeaturePlot(NKT, features = NKT_markers,ncol = 6,
            max.cutoff = 4, label = T)
dev.off()

celltype2 <- c("CD8^+TCF1^+T",  # 0
               "CD4^+TCF1^+T",  # 1
               "Treg",  # 2
               "CD8^+PD1^+T_{Ex}",  # 3
               "CD4^+TCF1^+T",  # 4
               "CD4^+PD1^+T",  # 5
               "CD4^+TCF1^+T",  # 6
               "CD8^+TCF1^+T",  # 7
               "CD4^+TCF1^+T", # 8
               "CD4^+TCF1^+T",  # 9
               "CD8^+TCF1^+T",  # 10
               "CD56^+NK",  # 11
               "Uncharacterized",  # 12
               "CD8^+TCF1^+T",  # 13
               "CD56^+NK",  # 14
               "CD8^+GZMB^+T",  # 15
               "CD8^+TCF1^+T",  # 16
               "CD8^+GZMB^+T",  # 17
               "CD8^+PD1^+T_{Ex}",  # 18
               "CD4^+TCF1^+T",  # 19
               "Uncharacterized",  # 20
               "CD8^+PD1^+T_{Ex}"  # 21
)
names(celltype2) <- levels(NKT)
NKT <- RenameIdents(NKT, celltype2)
NKT@meta.data$CellType2 <- NKT@active.ident
DimPlot(NKT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Map T cell subtype ID to the normalized data
# Sample data
xx <- normalized.data@meta.data
xx$CellType1 <- as.character(xx$CellType1)
yy <- NKT@meta.data
yy$CellType2 <- as.character(yy$CellType2)
# Convert the row names of xx to numeric
xx_row_names <- row.names(xx)

# Match the row names of xx with yy
matching_rows <- match(xx_row_names, row.names(yy))

# Replace values in CellType1 of xx with values from yy
xx$CellType1[!is.na(matching_rows)] <- yy$CellType2[matching_rows[!is.na(matching_rows)]]

# Add the new column CellType2 to xx
xx$CellType2 <- NULL  # Remove CellType2 if it exists (to ensure no duplicate columns)
xx$CellType2 <- ifelse(!is.na(matching_rows), yy$CellType2[matching_rows], xx$CellType1)

#  updated normalized data
normalized.data@meta.data$CellType2 <- xx$CellType2

# Export res
pdf(myoutf10)
DimPlot(normalized.data, reduction = "umap",
        group.by = "CellType1", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(NKT, reduction = "umap",
        group.by = "CellType2", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(normalized.data, reduction = "umap",
        group.by = "CellType2", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

saveRDS(normalized.data,file = myoutf11)
saveRDS(NKT,file = myoutf12)
# These two files were stored through zenodo
# in the "data" folders
# Please download the subfolder "TNBC_scRNA_GSE169246"
# and save them in the following path
#  "./data/TNBC_scRNA_GSE169246/"

sessionInfo()