#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Skylar Li
# 10/06/2023
# Generate Figures in Figure 1
# and the corresponding Supplementary tables and figures
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# scRNA UMAP ###################################################################

# Clear all variables
rm(list = ls())

# Load libraries
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(tibble)
library(scales)
library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape)

# The input files

myconfig1 <- "./Fig1/code/Fig1_config.R"


# These two files were preprocessed in Fig1.script1_Preprocess_scRNAseq.R
# And stored through zenodo
# in the "data" folders
# Please download the subfolder "TNBC_scRNA_GSE169246"
# and save them in the following path
#  "./data/TNBC_scRNA_GSE169246/"

myinf1 <- "./data/TNBC_scRNA_GSE169246/celltyping_all.rds"
myinf2 <- "./data/TNBC_scRNA_GSE169246/celltyping_NKT.rds"

# The clinical info was downloaded from Table S1 of the orginal paper:
# Zhang, Yuanyuan, et al. "Single-cell analyses reveal key immune cell subsets associated with response to PD-L1 blockade in triple-negative breast cancer." Cancer cell 39.12 (2021): 1578-1593.
# Here is the link to download the Table S1
# https://www.cell.com/cms/10.1016/j.ccell.2021.09.010/attachment/91cf6d86-3d31-41d4-8597-04aaff65572c/mmc2.xlsx
# This table was reformatted and saved as a csv file
# which is strode through zenodo in the subfolder "TNBC_scRNA_GSE169246"
myinf3 <- "./data/TNBC_scRNA_GSE169246/clinical_info.csv"

# The output files
myoutd <- "./Fig1/result/Figures_Tables"
dir.create(myoutd)
myoutf1 <- paste0(myoutd,"/Fig1b_UMAP_all_cell.pdf")
myoutf2 <- paste0(myoutd,"/FigS2a_UMAP_clinical_info.pdf")
myoutf3 <- paste0(myoutd,"/Fig1e_scRNA_marker_dotplot.pdf")
myoutf4 <- paste0(myoutd,"/FigS3_scRNA_marker_featureplot.pdf")
myoutf5 <- paste0(myoutd,"/FigS2b_scRNA_cell_proportion.pdf")
# settings for visualization  --------------------------------------------------
source(myconfig1)
cell <- c(cell,
          "NK & T",
          "Uncharacterized") 

color <- c(color,
           "#8856a7",
           "#d9d9d9"
           )
color <- alpha(color,0.8)
names(color) <- cell

# load data --------------------------------------------------------------------
normalized.data <- readRDS(myinf1)

# Fig1b All cell type in one figure --------------------------------------------
normalized.data@meta.data$CellType <- normalized.data@meta.data$CellType2 %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell[-15])
pdf(myoutf1, width = 8,height = 6)
DimPlot(normalized.data, reduction = "umap",
        group.by = "CellType", label = F, pt.size = 0.1,
        cols = color[-15], raster = F)
dev.off()


# ExtendedDataFig.2A Other info ------------------------------------------------
normalized.data@meta.data$CellType <- normalized.data@meta.data$CellType2 %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell[-15])

# combine clinical info
info <- read.csv(myinf3) 
yy <- normalized.data@meta.data %>% rownames_to_column("cellID")
xx <- merge(yy,info,by.x = "Patient", by.y = "Patient.ID") %>%
  column_to_rownames("cellID")

xx <- xx[rownames(normalized.data@meta.data),] %>%
  select(Treatment, Biopsied.lesion, Clinical.efficacy.)
all(rownames(normalized.data@meta.data) == xx$cellID)  # T
normalized.data@meta.data  <- cbind(normalized.data@meta.data , xx)

table(normalized.data$Clinical.efficacy.,exclude = F)
normalized.data@meta.data$Response <- ifelse(normalized.data@meta.data$Clinical.efficacy. == "PR", "Responder","Non-Responder") %>%
  factor(levels = c("Responder","Non-Responder"))
table(normalized.data$Response,exclude = F)

pdf(myoutf2, width = 16,height = 12)
mycols <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Dark2"), brewer.pal(8,"YlGnBu")[7])
mycols <- alpha(mycols,0.3)
DimPlot(normalized.data, reduction = "umap", ncol = 2,cols = mycols,
        group.by = c("Patient","TimePoint","Treatment","Response"), label = F, 
        pt.size = 0.1, raster = T)
dev.off()

# scRNA Dotplot ================================================================
# Fig1e dot plot for main figure -----------------------------------------------

normalized.data <- normalized.data %>% 
  subset(subset = CellType2 != "Uncharacterized")

normalized.data@meta.data$CellType <- normalized.data@meta.data$CellType2 %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell[length(cell):1])

Idents(normalized.data) <- "CellType"

pdf(myoutf3, width = 10,height = 5)
DotPlot(object =normalized.data , scale = T,
        feature= gene) + 
  theme(axis.text.x = element_text(angle = 90,hjust=1)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  guides(size=guide_legend(override.aes=list(shape=21, 
                                             colour="black", 
                                             fill="white"),
                           title = "Percent Expressed" ))
dev.off()


pdf(myoutf4, width = 24, height = 35)
FeaturePlot(normalized.data, features = gene,ncol = 5,
            pt.size = 0.1,
            raster=T, max.cutoff = 4, label = T)
dev.off()

# ExtendeData 2b scRNA stackplot ===============================================
# Reset the cell type
source(myconfig1)

normalized.data <- normalized.data %>% 
  subset(subset = CellType2 != "Uncharacterized")

normalized.data@meta.data$CellType <- normalized.data@meta.data$CellType2 %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell) 

info <- normalized.data@meta.data[,c("CellType","Patient","Treatment","Response","TimePoint")] %>%
  mutate(Treatment = ifelse(Treatment == "Chemo", "C", "C&I")) %>% 
  data.frame()
dim(info)

mylist <- data.frame(Treatment = c(rep("C&I",2)) ,
                     TimePoint = c("Pre","Post"))

color <- alpha(color,0.8)
p.list <- list()
for ( ii in 1:nrow(mylist)) {
  
  ref <- info %>%
    filter(Treatment == mylist$Treatment[ii] & 
             TimePoint == mylist$TimePoint[ii] ) %>%
    select(Patient, Response) %>%
    unique()
  rownames(ref) <- NULL
  
  
  p.data  <- info %>%
    filter(Treatment == mylist$Treatment[ii] & 
             TimePoint == mylist$TimePoint[ii] ) %>%
    group_by(Patient,CellType,.drop = F) %>%
    count(name = "Count") %>%
    data.frame() %>%
    merge(y=ref, by ="Patient")
  
  p.list[[ii]] <- p.data %>%
    ggplot(aes(fill=CellType, y=Count, x=Patient)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(name="Cell Type", 
                      values= color) +
    theme_bw(base_size = 12, base_family = "serif") + 
    facet_grid( ~ Response, space="free", scales="free",drop = TRUE) +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Patient ID",y="Percentage of Immune Cells",
         title=paste0( "Cell Proportion - Treatment: ",
                       mylist$Treatment[ii],
                       "(",
                       mylist$TimePoint[ii],
                       ")")) +
    theme(legend.position="right",
          legend.direction ="vertical",
          legend.box = "vertical",
          panel.grid=element_blank(),
          legend.title = element_text(face="bold", color="black",family = "serif", size=12),
          legend.text= element_text(face="bold", color="black",family = "serif", size=12),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="black",family = "serif", size=10),
          axis.text.y = element_text(face="bold", color="black", family = "serif",size=10),
          axis.title.x = element_text(face="bold", color="black",family = "serif", size=12),
          axis.title.y = element_text(face="bold",color="black",family = "serif", size=12)) +
    scale_y_continuous(expand = c(0,0)) 
  
}


pdf(myoutf5, width = 12,height = 4.5)
do.call("grid.arrange", c(plotlist = p.list, ncol=2))
dev.off()

# Spatial Heatmap #######################################################################
rm(list = ls())
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(reshape)
myoutd <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/Fig1"
myconfig1 <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/Fig1_config.R"
source(myconfig1)


myinf1 <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/TNBC_scRNA_GSE169246/spatial_markers.rda"

myoutf1 <- paste0(myoutd,"/Spatial_marker_HM.pdf")


# Spatial markers  -------------------------------------------------------------

load(myinf1)

# Choose 11 immune cell markers for visualization
markers_to_print <- spatial_markers %>% 
  filter(Value > 0.5) %>%
  filter(! Protein %in% c("Vimentin", "PDGFRB", "CD31")) %>%
  pull(Protein) %>%
  unique()


spatial_markers$CellType <- spatial_markers$CellType %>%
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell)

p.data <- spatial_markers %>% 
  filter(CellType %in% cell) %>%
  filter(Protein %in% markers_to_print) %>%
  mutate(Protein =   factor(Protein, levels = protein)) %>%
  cast(CellType ~ Protein, value = "Value") %>%
  column_to_rownames(var = "CellType")

# Heatmap -------------------------------------------------------------
color <- alpha(color,0.8)
p1 <-  Heatmap(as.matrix(p.data),  
               cluster_rows = F, 
               cluster_columns = F,
               heatmap_legend_param = list(title = "z-score"), 
               col = colorRamp2(c(-1.5, 0, 4.5), c("blue", "white", "red")),
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               row_split = factor(group[rownames(p.data)],levels = unique(group[rownames(p.data)])),
               left_annotation = 
                 rowAnnotation(
                   CellType = rownames(p.data),
                   col = list(CellType = color[rownames(p.data)]),
                   show_annotation_name = c(CellType = F),
                   show_legend = F),
               
               row_title=NULL,
               show_row_names = T,
               row_names_side = 'left',
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(p.data)*unit(5, "mm"),
               height = nrow(p.data)*unit(5, "mm"))
p1
pdf(myoutf1,
    width = 8,  
    height = 4)
print(p1)
dev.off()





# TimiGP IMC dotplot ###################################################################
# conclusion: only use the C&I and GSE194040 set
rm(list = ls())
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(tibble)

myPackage <- "~/Mypackage/TimiGP/"
library(devtools)
load_all(myPackage)

myoutd <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/Fig1"
myconfig1 <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/Fig1_config.R"
source(myconfig1)
myind <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/TNBC_scRNA_marker_for_spatial/"
myinf1 <- list.files(myind,pattern = "_TNBC_ENRICH_responder.rda", 
                     recursive = T, full.names = T) %>%.[c(6,4,1,2)]
myinf2 <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/TNBC_spatial_ratio_responder/Spatial_cell_ratio_response.rda"
myinf3 <- "~/ChaoCheng/A_TimiGP_Therapy/A05_TimiFisher_TNBC/TNBC_scRNA_GSE169246/TNBC_scRNA_marker_for_spatial.rda"


myoutf1 <- paste0(myoutd,"/TimiGP_score_dotplot.pdf")

basename(myinf1) 
treatment <- rep(c("C&I","C"),2)
# TimiGP interaction & IMC ratio  ---------------------------------------------
for ( dd in 1: length(myinf1)) {
  load(myinf1[dd])
  load(myinf2)
  rna_ds_id <- basename(myinf1[dd]) %>% strsplit("_") %>% 
    lapply("[[", 2) %>% unlist()
  
  myoutf1 <- paste0(myoutd,"/cell.interaction_",rna_ds_id,"_",treatment[dd],"_dotplot.pdf")
  cat("\n",basename(myoutf1), ":", basename(myinf1[dd]))
  
  # TimiGP CI
  res_ci <- enrich_res %>% 
    filter(Adjust.P.Value < 0.05) %>%
    pull(Cell.Interaction)
  length(res_ci) #49
  # extract spatial result

  tp <- names(spatial_ref[[treatment[dd]]])
  
  # select cell interaction to print
  print_ci <- NULL
  for ( ii in tp) {
    
    tmp <- spatial_ref[[treatment[dd]]][[ii]]$log_res_ratio %>%
      filter(PV < 0.05) %>%
      rownames(.)  %>%
      gsub(pattern = ".",replacement = "+",fixed = T)
    cat("\n",ii,
        "; IMC PV < 0.05: ", length(tmp),
        "; TimiGP QV < 0.05:", length(res_ci),
        "; intersect: ", length(intersect(tmp,res_ci)))
    print_ci <- unique(c(print_ci,tmp))
    
    
  }
  
  # dd=1
  # Post-treatment ; IMC PV < 0.05:  20 ; TimiGP QV < 0.05: 49 ; intersect:  8
  # Baseline ; IMC PV < 0.05:  14 ; TimiGP QV < 0.05: 49 ; intersect:  4
  # On-treatment ; IMC PV < 0.05:  19 ; TimiGP QV < 0.05: 49 ; intersect:  10
  
  
  length(print_ci) # 42
  # extract cell interaction from TimiGP result
  p.data <- enrich_res %>% 
    filter(Cell.Interaction %in% print_ci) %>%
    filter(P.Value < 0.05) %>%
    select(Cell.Interaction, Enrichment.Ratio, P.Value, Adjust.P.Value) %>%
    setnames(old = "Enrichment.Ratio",new = "ER_OR") %>%
    setnames(old = "P.Value",new = "PV") %>%
    setnames(old = "Adjust.P.Value",new = "QV") %>%
    mutate(Group = paste0("TimiGP_",rna_ds_id)) 
  
  print_ci <- p.data$Cell.Interaction
  length(print_ci) #16
  
  if(dd == 1) {
    load(myinf3)
    net.dir <- paste0(myoutd,"/netowrk_",rna_ds_id,"_",treatment[dd])
    dir.create(net.dir)
    sel <- enrich_res %>% 
      filter(Cell.Interaction %in% print_ci) %>%
      pull(Index)
    TimiCellNetwork(resdata = resdata,
                    dataset = "Other",
                    geneset = TNBC_scRNA_marker_for_spatial,
                    group = group,
                    select = sel,
                    export = T,
                    path = net.dir)
  }
  
  
  # cobime IMC and TimiGP res
  for ( ii in tp) {
    print(ii)
    p.data <-  spatial_ref[[treatment[dd]]][[ii]]$log_res_ratio %>%
      rownames_to_column("Cell.Interaction") %>% 
      mutate(Cell.Interaction = gsub(pattern = ".",replacement = "+", 
                                     x = Cell.Interaction,fixed = T)  ) %>%
      filter(Cell.Interaction %in% print_ci) %>%
      setnames(old = "OR",new = "ER_OR")  %>%
      mutate(Group = paste0("IMC_",ii))  %>%
      rbind(p.data)
    
  }
  
  # dot plot
  ER_max <- max(p.data$ER_OR)
  ds <- c(paste0("TimiGP_",rna_ds_id),
          "IMC_Baseline",
          "IMC_On-treatment",
          "IMC_Post-treatment")
  p.data <- p.data %>%
    mutate(Group = factor(Group,levels = ds[length(ds):1])) %>%
    mutate(Cell.Interaction = factor(Cell.Interaction, 
                                     levels = print_ci)) 
  ER_max <- ceiling(max(p.data$ER_OR))
  p1 <- ggplot() +
    geom_point(data = p.data,
               mapping = aes( y=Group, x=Cell.Interaction,color=PV, 
                              size=ER_OR),shape = 16) +
    geom_point(data = p.data[which(p.data$PV < 0.05),],
               mapping = aes( y=Group, x=Cell.Interaction,
                              size=ER_OR),shape = 1,fill = "white",color="black", stroke = 1) +
    geom_point(data = p.data[which(p.data$PV > 0.05),],
               mapping = aes( y=Group, x=Cell.Interaction,
                              size=ER_OR),shape = 1,fill = "white",color="black", stroke = 0.2) +
    scale_colour_gradientn(colours=c("#7fcdbb","#edf8b1","white","white"),
                           breaks = c(0,0.01,0.05,0.1),
                           labels = c(0,0.01,0.05,1),
                           limits = c( 0,0.1),
                           guide = guide_colorbar(barwidth = 1, barheight =8,
                                                  reverse=TRUE,
                                                  frame.colour = "black", ticks.colour = "black"),
                           oob = scales::squish,
                           name = "P-Value") +
    
    scale_size(range = c(0,5), 
               limits = c(1,ER_max),name="Enrichment/\n Odds Ratio") +
    geom_text(data = p.data[which(p.data$QV < 0.1),],
              aes(y=Group, x=Cell.Interaction),hjust = 0.5,
              vjust = 0, label = "*", size = 5,color = "black") +
    theme_bw(base_size = 8, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(y="Method_Dataset",x="Cell Interaction",
         title= paste0(treatment[dd]," (TimiGP:",rna_ds_id,")"),
         base_size = 8, base_family = "serif",face="bold") +
    theme(legend.background = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          legend.position="right",
          legend.box = "horizontal",
          legend.direction= "vertical",
          panel.grid.major = element_line(color = "#d9d9d9",
                                          linewidth = 0.5,
                                          linetype = 2),
          panel.grid=element_blank(),
          legend.key.width = unit(0.5,"cm"),
          legend.title = element_text(face="bold", color="black",
                                      family = "serif", size=8),
          legend.text= element_text(face="bold", color="black",
                                    family = "serif", size=8),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text( color="black", size=8, 
                                      angle = 40,  vjust = 1, hjust=1),
          axis.text.y = element_text( color="black", size=8),
          axis.title.x = element_text(face="bold", color="black", size=8.5),
          axis.title.y = element_text(face="bold",color="black", size=8.5))
  pdf(myoutf1, width = 8,height = 4)
  print(p1)
  dev.off()
}


# conclusion: only use the C&I and GSE194040 set

# TimiGP score & IMC count  ---------------------------------------------
# Describe CD8+GZMB+T is the only one significant after correction
# No figures (the analysis has been done by the IMC group, just describe our result and their result)
dd <- 1
load(myinf1[dd])
load(myinf2)
rna_ds_id <- basename(myinf1[dd]) %>% strsplit("_") %>% 
  lapply("[[", 2) %>% unlist()

myoutf1 <- paste0(myoutd,"/Cell.Type_",rna_ds_id,"_",treatment[dd],"_dotplot.pdf")
cat("\n",basename(myoutf1), ":", basename(myinf1[dd]))

# extract spatial result
imc_res <- NULL
tp <- names(spatial_ref[[treatment[dd]]])
for ( ii in tp) {
  print(ii)
  imc_res[[ii]] <- spatial_ref[[treatment[dd]]][[ii]]$log_res_cell %>%
    rownames_to_column("Cell.Type") %>% 
    mutate(Cell.Type = gsub(pattern = ".",replacement = "+", 
                            x = Cell.Type,fixed = T)  ) %>%
    mutate(Group = ii) %>%
    filter(PV < 0.05) 
  
}

# select cell type to print
print_ci <-  union(imc_res[[1]]$Cell.Type,imc_res[[2]]$Cell.Type) %>%
  union(imc_res[[3]]$Cell.Type)

#
# TimiGP bar plot

score <- enrich_res %>% TimiFS() %>%
  filter(Cell.Type %in% print_ci) 
p1 <- TimiFSBar(score)


# Spatial Heatmap 
# plot data
p.data <- matrix(ncol = length(print_ci),nrow = 3)
rownames(p.data) <- factor(names(spatial_ref[[treatment[dd]]]),
                           levels = c("Baseline","On-treatment",   "Post-treatment"))

colnames(p.data) <- factor(score$Cell.Type, 
                           levels = score$Cell.Type)
pv <- qv <- p.data
for ( ii in rownames(p.data) ) {
  print(ii)
  tmp <- spatial_ref[[treatment[dd]]][[ii]]$log_res_cell %>%
    rownames_to_column("Cell.Type") %>% 
    mutate(Cell.Type = gsub(pattern = ".",replacement = "+", 
                            x = Cell.Type,fixed = T)  ) %>%
    filter(Cell.Type %in% colnames(p.data)) %>%
    column_to_rownames("Cell.Type") 
  tmp <- tmp[colnames(p.data),]
  p.data[ii,] <- log(tmp$OR)
  pv[ii,] <- tmp$PV
  qv[ii,] <- tmp$QV
  
}
# heatmap 
p2 <-  Heatmap(as.matrix(p.data),  
               cluster_rows = F, 
               cluster_columns = F,
               heatmap_legend_param = list(title = "log(Odds Ratio)"), 
               col = colorRamp2(c(-0.07, 0, 0.02), c("blue", "white", "red")),
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               #  row_split = factor(group[colnames(p.data)],levels = unique(group[colnames(p.data)])),
               bottom_annotation = 
                 columnAnnotation(
                   CellType = colnames(p.data),
                   col = list(CellType = color[colnames(p.data)]),
                   show_annotation_name = c(CellType = F),
                   show_legend = F),
               
               row_title=NULL,
               show_row_names = T,
               row_names_side = 'left',
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(p.data)*unit(5, "mm"),
               height = nrow(p.data)*unit(5, "mm"))
p2
p.data
pv
qv

# CD8+GZMB+T  CD8+PD1+Tex       CD20+B            M2      CD56+NK   CD8+TCF1+T CD79a+Plasma         Treg    CD4+PD1+T   CD4+TCF1+T
# Post-treatment -0.029290460 -0.001563908 -0.002875893 -1.871094e-05 -0.002384876 -0.011816500 -0.008707292 -0.132506342 -0.016410058 -0.007185415
# Baseline        0.007895981  0.001161754  0.002319497 -3.975515e-04  0.018636400  0.001616396  0.001278921  0.003866025  0.003536268  0.001086224
# On-treatment    0.006488873  0.002500123  0.003054536  1.703499e-03  0.036664323  0.008911918  0.000710594  0.001887121  0.004125066  0.001050639
# > pv
# CD8+GZMB+T CD8+PD1+Tex     CD20+B          M2     CD56+NK CD8+TCF1+T CD79a+Plasma        Treg   CD4+PD1+T CD4+TCF1+T
# Post-treatment 0.015061052 0.346928758 0.12878899 0.985913676 0.687385346 0.10981340   0.10085621 0.001389132 0.028469552 0.04675083
# Baseline       0.003564723 0.453506067 0.07777441 0.731470808 0.027802624 0.42552383   0.04302342 0.014691672 0.024438706 0.14850767
# On-treatment   0.000513281 0.003060139 0.03968195 0.009283761 0.007955441 0.00201136   0.29363936 0.083492107 0.007682898 0.19494394
# > qv
# CD8+GZMB+T CD8+PD1+Tex     CD20+B         M2    CD56+NK CD8+TCF1+T CD79a+Plasma       Treg  CD4+PD1+T CD4+TCF1+T
# Post-treatment 0.082835784  0.46555436 0.20238269 0.98591368 0.75612388 0.20132457   0.20132457 0.01528046 0.10438836  0.1285648
# Baseline       0.039211955  0.55428519 0.14258641 0.75741171 0.07645721 0.55428519   0.09465153 0.07645721 0.07645721  0.2333692
# On-treatment   0.005646091  0.01122051 0.06235735 0.01702023 0.01702023 0.01106248   0.32300329 0.11480165 0.01702023  0.2382648

# Describe CD8+GZMB+T is the only one significant after correction

# Tables ######################################################################
# TODO
# Gene information used to generate TimiGP markers
  # ROC Gene
  # p.adj gene
  # Spatial value gene

# TODO IMC_TimiGP_dotplot
# Using adobe illustrator to 
# add "* FDR < 0.1" and" border circle, P-Value < 0.05" 