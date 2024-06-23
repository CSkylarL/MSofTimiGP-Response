#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# Analyze single cell RNA-seq: TNBC_scRNA_GSE169246
# Aim: Generate markers for TimiGP
# 10/02/2023
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Generate markers for TimiGP ##################################################
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# Find markers -----------------------------------------------------------------

# Load libraries
rm(list = ls())
library(Seurat)
library(dplyr)

print("-------Find markers for TimiGP----------")
warnings("---------Find markers for TimiGP----------")


# The input data is the Seurat object normalized.data from the previous script
# Fig1.script1_Preprocess_scRNAseq.R
# Please run the previous script before running this script or 
# Download the preprocessed data from zenodo
# https://doi.org/10.5281/zenodo.12209783
# Please decompress "TNBC_scRNA_GSE169246.zip"
# and move them into the `data` folder
# The path should be
#  "./data/TNBC_scRNA_GSE169246/*.rda"

myind <- "./data/TNBC_scRNA_GSE169246/"
myinf1 <- paste0(myind, "celltyping_all.rds")

# The output directory
myoutd <- "./Fig1/result/customized_markers_for_TimiGP/"
dir.create(myoutd, showWarnings = F)
myoutf1 <- paste0(myoutd, "unprocessed_markers_for_TimiGP.rda")

# Load the data
normalized.data <- readRDS(myinf1)

plan("multisession", workers = 40)
options(future.globals.maxSize = 150 * 1024^3)
Idents(normalized.data) <- "CellType2"

# Find markers
print("-------markers_wilcox----------")
warnings("---------markers_wilcox----------")
markers_wilcox <- FindAllMarkers(normalized.data, test.use = "wilcox",
                          only.pos = F, min.pct = 0.25, 
                          logfc.threshold = 0.25)

markers_roc<- FindAllMarkers(normalized.data, test.use = "roc",
                                  only.pos = F, min.pct = 0.25,
                                  logfc.threshold = 0.25)


# Save the markers
save(markers_wilcox,markers_roc, file = myoutf1)

# Get the spatial markers from IMC data ------------------------------------------------------
# Load libraries
rm(list = ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(fst)
library(parallel)
library(reshape)

# Please download the IMC data from https://zenodo.org/record/7990870
# Save the data in the folder "data"
# And rename the folder to "TNBC_IMC_zenodo.7990870"

myind <- "./data/TNBC_IMC_zenodo.7990870/"
myfun1 <- "./Fig1/code/Fig1_function.R"

myoutd <- "./Fig1/result/customized_markers_for_TimiGP/"
myoutf1 <- paste0(myoutd, "spatial_markers.rda")

# functions to analyze IMC data
source(myfun1)

## The code is modified from Fig 2 of 
# https://zenodo.org/record/7990870 (AUTHOR: H R Ali)
cells <- getCells(dir = myind,curated = TRUE, allCells = TRUE)
panel <- getPanel(dir = myind)
thresholds <- getPositiveThresholds()

cells[, eval(panel$panel) := lapply(.SD, function(x) scale(clip99(x))), .SDcols = panel$panel]
cells[, isKi67Pos := Ki67 > quantile(Ki67, probs = thresholds['Ki67'])]

setkeyv(cells, c('ImageNumber', 'ObjectNumber'))

# TME Heatmap

SpatialMarkers <- function(dat, clusterColumns = NULL, ByVar = NULL){
  
  Zmedians <- dat[!is.na(get(ByVar)), lapply(.SD, median), by = ByVar, .SDcols = clusterColumns]
  Zmedians <- Zmedians[!is.na(get(ByVar)), eval(clusterColumns) := 
                         lapply(.SD, scale), .SDcols = clusterColumns][order(get(ByVar))]
  
  cellAnnot <- dat[, .N, by = .(PrintOrder, Colour, Label)][order(PrintOrder)]    
  
  toHM <- as.matrix(Zmedians[,-1, with = F])
  rownames(toHM) <- cellAnnot[, Label]
  
  return(toHM)
}

cellClustersAnnot <- getCellClusters(dir = myind)
TME <- grep('H3|DNA|Carboplatin|Vimentin|Calponin|Ki67|c-PARP|CK|AR|pH2AX', 
            panel$panel, invert = TRUE, value = TRUE)
protein_ref <- SpatialMarkers(cells[!(isEpithelial)], clusterColumns = c(TME, 'Vimentin', 'Calponin'),
                              ByVar = 'PrintOrder')
# select the markers
cell <- c("CD20^+B" ,"CD79a^+Plasma",
          "CD4^+TCF1^+T", "CD4^+PD1^+T", "Treg",
          "CD8^+TCF1^+T", "CD8^+GZMB^+T","CD8^+PD1^+T_{Ex}",
          "CD56^+NK" ,
          "DCs"  ,"pDC" ,       
          "M1 Mac" ,"M2 Mac",
          "Mast")

spatial_markers <- melt(protein_ref) %>% data.frame() %>% 
  filter(X1 %in% cell)  %>%
  arrange(X1, -value)
colnames(spatial_markers) <- c("CellType","Protein","Value")
save(spatial_markers, file = myoutf1)


# Genearate final markers --------------------------------------------------------------

# Load libraries
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# The input and the output directories
myoutd <- "./Fig1/result/customized_markers_for_TimiGP/"
myinf1 <- paste0(myoutd, "unprocessed_markers_for_TimiGP.rda")
myinf2 <- paste0(myoutd, "spatial_markers.rda")
myoutf1 <- paste0(myoutd, "TNBC_scRNA_marker_for_spatial.rda")

# Load the markers from scRNA-seq
load(myinf1)
rownames(markers_roc) <- NULL
rownames(markers_wilcox) <- NULL

# Check cutoff AUC >0.7 
tmp <- markers_roc %>%
  filter(myAUC > 0.7)
table(tmp$cluster)


# CD4^+TCF1^+T     CD8^+GZMB^+T           CD8^+TCF1^+T             Treg 
# 2  (more)       19 (keep)               8    (more)           16 (keep) 
# pDC                     CD4^+PD1^+T                 CD56^+NK         CD20^+B 
# 287  (less)               9  (more)             35 (keep)               38 (keep)
# M1 Mac                Mast                 CD79a^+Plasma  CD8^+PD1^+T_{Ex} 
# 44  (less)             74 (less)              70  (less)          24 (keep) 
# M2 Mac              DCs  Uncharacterized 
# 359 (less)     253 (less)          30 (drop)
keep <- c("CD8^+GZMB^+T",
          "Treg",
          "CD56^+NK",
          "CD20^+B",
          "CD8^+PD1^+T_{Ex}")

more <- c("CD4^+TCF1^+T" ,      
          "CD8^+TCF1^+T",    
          "CD4^+PD1^+T"   )

geneset <- markers_roc %>%
  filter(myAUC > 0.7) %>% 
  filter(cluster %in% c(keep,more)) %>%
  select(cluster,  gene)

table(geneset$cluster)


# Select top 30 markers following AUC for cell types with too many markers (less):
less <- c("pDC",
          "M1 Mac",
          "Mast",
          "CD79a^+Plasma",
          "M2 Mac", 
          "DCs" )
tmp <- markers_roc %>%
  filter(cluster %in% less) %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = myAUC) %>% 
  ungroup() %>%
  data.frame() %>%
  select(cluster,  gene)

table(tmp$cluster)
geneset <- rbind(geneset,tmp)
table(geneset$cluster)

# Select top 30 markers following Pvalue and avg_log2FC for cell types with too less markers (more):

tmp <- markers_wilcox %>%
  filter(cluster %in% more) %>%
  filter(p_val_adj == 0) %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC,) %>% 
  ungroup() %>%
  data.frame() %>%
  select(cluster,  gene)

table(tmp$cluster)
geneset <- rbind(geneset,tmp)

# Load the spatial markers
load(myinf2) # spatial ref
# Add spatial markers used to define the cluster
tmp <- data.frame(
  cluster =  c(rep("CD20^+B" , 8       )    ,
               rep("CD79a^+Plasma",2   )    ,   
               rep("CD4^+TCF1^+T",  6   )   ,    
               rep("CD4^+PD1^+T",  9    )   ,  
               rep("Treg",      10      ),
               rep("CD8^+TCF1^+T",  9   )   ,    
               rep("CD8^+GZMB^+T",  7   )   ,  
               rep("CD8^+PD1^+T_{Ex}", 9)   ,      
               rep("CD56^+NK" ,  1)      ,
               rep("DCs"  ,    4        )  ,
               rep("M2 Mac",    4       )  ),
  gene =  c("MS4A1","CD79A",  "HLA-A"   , "HLA-B" , "HLA-C" , "HLA-DRA" ,  "HLA-DRB5" ,    "HLA-DRB1", # 8
            "CD79A", "PECAM1", # 2
            "TCF7", "CD4"  ,"CD3E",  "CD3D", "CD3G" ,  "GATA3", # 6
            "CD4"  ,"CD3E",  "CD3D", "CD3G" ,  "GATA3",  "PDCD1",  "TOX", "ICOS", "TBX21", # 9
            "CD4"  ,"CD3E",  "CD3D", "CD3G" ,  "GATA3",  "ICOS", "TBX21", "FOXP3","IKZF2","TNFRSF4", # 10
            "CD8A", "CD8B",  "CD3E",  "CD3D", "CD3G" , "TCF7","TBX21","GATA3", "TOX", # 9
            "CD8A", "CD8B",  "CD3E",  "CD3D", "CD3G" , "GZMB","TBX21", # 7
            "CD8A", "CD8B",  "CD3E",  "CD3D", "CD3G" , "TBX21", "GATA3",  "PDCD1",  "TOX", # 9
            "NCAM1", # 1
            "ITGAX", "CD68", "HLA-DRB5" ,    "HLA-DRB1", # 4
            "CD68", "CD163", "HLA-DRB5" ,    "HLA-DRB1"# 4 
            ))

# remove duplicates (spatial marker has a lot duplicates with scRNA markers)
geneset <- rbind(geneset,tmp)
geneset <- unique(geneset)
table(geneset$cluster)


# CD4^+TCF1^+T     CD8^+GZMB^+T     CD8^+TCF1^+T             Treg 
# 36               25               37               23 
# pDC      CD4^+PD1^+T         CD56^+NK          CD20^+B 
# 30               37               36               41 
# M1 Mac             Mast    CD79a^+Plasma CD8^+PD1^+T_{Ex} 
# 31               30               31               28 
# M2 Mac              DCs  Uncharacterized 
# 32               32                0 

# Finalize format for TimiGP
geneset <- geneset %>% 
  mutate(cluster = as.character(cluster)) %>%
  mutate(gene = as.character(gene))
colnames(geneset) <- c("CellType","Gene")
geneset$Dataset <- "TNBC_scRNA_marker_for_spatial"
geneset$CellType <- geneset$CellType %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T)

# rename the variable and save the data
TNBC_scRNA_marker_for_spatial <- geneset
save(TNBC_scRNA_marker_for_spatial, file = myoutf1)


sessionInfo()