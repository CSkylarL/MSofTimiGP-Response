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
# Figures ######################################################################
# scRNA UMAP ===================================================================

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
myoutd <- "./Fig1/result/Figures_Tables/"
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

# Spatial Heatmap ==============================================================
# Clear all variables
rm(list = ls())

# Libraries
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(reshape)

# The input files

myconfig1 <- "./Fig1/code/Fig1_config.R"
source(myconfig1)
myinf1 <- "./Fig1/result/customized_markers_for_TimiGP/spatial_markers.rda"

# The output files
myoutd <- "./Fig1/result/Figures_Tables/"
myoutf1 <- paste0(myoutd,"/Fig1e_Spatial_marker_HM.pdf")

# Spatial markers  -------------------------------------------------------------

load(myinf1)

# Choose 11 immune cell markers for visualization
markers_to_print <- spatial_markers %>% 
  filter(Value > 0.5) %>%
  filter(! Protein %in% c("Vimentin", "PDGFRB", "CD31")) %>%
  pull(Protein) %>%
  unique()

# Update the cell type name
spatial_markers$CellType <- spatial_markers$CellType %>%
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>%
  factor(levels = cell)

# Heatmap data 
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

# TimiGP results Fig 1c-d ======================================================
# These figures were automatically generated by the TimiGP R package
# Which were generated with the following script
# Fig1.script3_customized_TimiGP_response_for_TNBC.R
# and stored in the following path
# "./Fig1/result/TimiGP_TNBC_scRNA_marker_for_spatial/Wolf_GSE194040_ISPY2_ptx_pem"
# Wolf_GSE194040_ISPY2_ptx_pem_TNBC_circle.pdf -> cell-cell interaction network Fig1c
# Wolf_GSE194040_ISPY2_ptx_pem_TNBC_score.pdf -> favorability score Fig 1d


# TimiGP IMC dotplot ===========================================================
# conclusion: only use the C&I and GSE194040 set
# clear all variables
rm(list = ls())

# load libraries
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(TimiGP)

# The input files

myconfig1 <- "./Fig1/code/Fig1_config.R"
source(myconfig1)
myind <- "./Fig1/result/TimiGP_TNBC_scRNA_marker_for_spatial/"
myinf1 <- list.files(myind,pattern = "_TNBC_ENRICH_responder.rda", 
                     recursive = T, full.names = T) 

myinf2 <- "./Fig1/result/TNBC_spatial_ratio_responder/Spatial_cell_ratio_response.rda"
myinf3 <- "./Fig1/result/customized_markers_for_TimiGP/TNBC_scRNA_marker_for_spatial.rda"

# The output files
myoutd <- "./Fig1/result/Figures_Tables/"

# The dataset name
basename(myinf1) 
treatment <- "C&I"

# TimiGP interaction & IMC ratio  ---------------------------------------------
for ( dd in 1: length(myinf1)) {
  load(myinf1[dd])
  load(myinf2)
  rna_ds_id <- basename(myinf1[dd]) %>% strsplit("_") %>% 
    lapply("[[", 2) %>% unlist()
  
  myoutf1 <- paste0(myoutd,"/Fig1g_cell.interaction_",rna_ds_id,"_",treatment[dd],"_dotplot.pdf")
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
  # Fig 1f --------------------------------------------------------------------
  if(dd == 1) {
    load(myinf3)
    net.dir <- paste0(myoutd,"/network_",rna_ds_id,"_",treatment[dd])
    dir.create(net.dir)
    sel <- enrich_res %>% 
      filter(Cell.Interaction %in% print_ci) %>%
      pull(Index)
    TimiCellNetwork(resdata = enrich_res,
                    dataset = "Other",
                    geneset = TNBC_scRNA_marker_for_spatial,
                    group = group,
                    select = sel,
                    export = T,
                    path = net.dir)
  }
  
  # Fig 1g --------------------------------------------------------------------
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
  # Using adobe illustrator to 
  # add "* FDR < 0.1" and" border circle, P-Value < 0.05" to the dot plot
}


# Tables ######################################################################

# Table S1 customized markers for TimiGP ======================================
# Clear all variables
rm(list = ls())
# The input files
myind <- "./Fig1/result/customized_markers_for_TimiGP/"
myinf1 <- paste0(myind, "TNBC_scRNA_marker_for_spatial.rda")
load(myinf1)

# The output files
myoutd <- "./Fig1/result/Figures_Tables/"
myoutf1 <- paste0(myoutd,"TableS1_scRNA_marker_for_TimiGP.csv")
write.csv(TNBC_scRNA_marker_for_spatial[1:2], quote=F,
          file = myoutf1, row.names = F)

# Table S2 TimiGP cell-cell interactions ======================================
# Clear all variables
rm(list = ls())

# The input files
myinf1 <-  "./Fig1/result/TimiGP_TNBC_scRNA_marker_for_spatial/Wolf_GSE194040_ISPY2_ptx_pem/Wolf_GSE194040_ISPY2_ptx_pem_TNBC_ENRICH_responder.rda"
load(myinf1)

# The output files
myoutd <- "./Fig1/result/Figures_Tables/"
myoutf1 <- paste0(myoutd,"TableS2_TimiGP_cell_interaction.csv")
write.csv(enrich_res, quote=F,
          file = myoutf1, row.names = F)

# Table S3 IMC cell-cell interactions =========================================
# Clear all variables
rm(list = ls())

# libraries
library(dplyr)

# The input files
myinf1 <-  "./Fig1/result/TNBC_spatial_ratio_responder/Spatial_cell_ratio_response.rda"
load(myinf1)

# The output files
myoutd <- "./Fig1/result/Figures_Tables/"
myoutf1 <- paste0(myoutd,"TableS3_IMC_cell_interaction.csv")

# combine list to make a table

tp <- c("Baseline","On-treatment","Post-treatment")
res <- data.frame()

for (ii in tp) {
  tmp <- spatial_ref$`C&I`[[ii]]$log_res_ratio %>%
    mutate(TimePoint = ii) %>%
    rownames_to_column("Cell.Interaction")
  colnames(tmp) <- c("Cell.Interaction","Odds.Ratio","P-Value","BH-Adjusted.P-Value","Time.Point")
  res <- rbind(res,tmp)
}

write.csv(res, quote = F,
          file = myoutf1, row.names = F)

