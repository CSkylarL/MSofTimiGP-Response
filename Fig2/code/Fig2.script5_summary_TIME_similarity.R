#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
#  Calculate TIME cell-cell interactions similarity between datasets
# and visualize using heatmap
# Generate Extended Data Fig.5
# 02/28/2024
# Compared to A06, this version add EC in A04
# The CI definition of few datasets were changed from default, see config.txt

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# clear workspace
rm(list = ls())

# load libraries
library(ggplot2)
library(TimiGP)
library(dplyr)
library(RColorBrewer)

library(ComplexHeatmap)
library(circlize)


# Input & Output ===============================================================
# load config file
myinf1 <- "./Fig2/code/config.rda"
load(myinf1)

# The output directory
myoutd <- "./Fig2/result/summary"


kk <- "Newman2015"
  cat("\n", kk)
  # output file
  myoutf1 <-  paste0(myoutd,"/immunotherapy_all_hm_similarity_",kk,".pdf")
  myinf <- config %>% filter(CI_geneset == kk)
  
  # load data ----
  resdata <- data.frame()
  
  for (i in myinf$order) {
    load(myinf$location[i])
    tmp <-enrich_res %>% filter (
      get(myinf$CI_condition[i]) < myinf$CI_cutoff[i])
    tmp$Dataset <- myinf$oldname[i]
    resdata <- rbind(resdata,tmp)
    rm(enrich_res,tmp)
  } 
  
  # Tversky index  --------------------------------------------------------
  DS <- myinf$oldname
  p.data <- matrix(nrow = length(DS),ncol = length(DS))
  colnames(p.data) <- rownames(p.data) <- DS
  p.value <- p.data
  
  for ( i in 1:nrow(p.data)){
    tmpi <- resdata %>% 
      filter(Dataset == rownames(p.data)[i]) %>% 
      pull(Cell.Interaction)
    for (j in 1:ncol(p.data)){
      tmpj <- resdata %>% 
        filter(Dataset == rownames(p.data)[j]) %>% 
        pull(Cell.Interaction)
      
      shared <- intersect(tmpi, tmpj) %>% length()
      union <- union(tmpi, tmpj) %>% length()
      p.data[i,j] <-  shared / min(length(tmpi),length(tmpj)) 
      
    }
  }
  
  round(p.data,2)
  
  # p.value  --------------------------------------------------------
  if(kk == "Bindea2013_Cancer" ) {
    data("CellType_Bindea2013_cancer")
    geneset <- CellType_Bindea2013_cancer
    
  }
  
  if(kk == "Zheng2021") {
    data("CellType_Zheng2021_Tcell")
    geneset <- CellType_Zheng2021_Tcell
    
    
  }
  
  if(kk == "Newman2015") {
    data("CellType_Newman2015_LM22")
    geneset <- CellType_Newman2015_LM22
  }
  
  cell_pair <- TimiCellPair(geneset = geneset,core = 20)
  nn <- unique(cell_pair$Cell.Pair) %>% length()# total number of interaction
  for ( i in 1:nrow(p.value)){
    tmpi <- resdata %>% 
      filter(Dataset == rownames(p.value)[i]) %>% 
      pull(Cell.Interaction)
    for (j in 1:ncol(p.value)){
      tmpj <- resdata %>% 
        filter(Dataset == rownames(p.value)[j]) %>% 
        pull(Cell.Interaction)
      
      shared <- intersect(tmpi, tmpj) %>% length()
      union <- union(tmpi, tmpj) %>% length()
      aa <- length(tmpi)
      bb <- length(tmpj)
      p.value[i,j] <-  sum(dhyper(shared:min(aa,bb), aa, nn - aa, bb)) # p-value
    }
  }
  #round(p.value, 4) 
  
  
  
  # heatmap - order  --------------------------------------------------------
  p1 <-  Heatmap(as.matrix(p.data),  cluster_rows = F, 
                 heatmap_legend_param = list(title = "Tversky index"), 
                 cluster_columns = F,
                 col =  colorRamp2(c( 0, 1), c("white", "red")),
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   
                   if(p.value[i, j] < 0.05) {
                     # select cell with *
                     gb = textGrob("*")
                     gb_w = convertWidth(grobWidth(gb), "mm")
                     gb_h = convertHeight(grobHeight(gb), "mm")
                     grid.text("*", x, y - gb_h*0.6 + gb_w*0.4,
                               gp = gpar(fontsize = 15))
                     # # select cell with border
                     # border <- "black" # Set the border color to black
                     # border_width <- unit(2, "mm") # Set the border width to 2 mm
                     # border_rect <- rectGrob(x, y, w, h, 
                     #                         gp = gpar(col = border, 
                     #                                   fill = NA,
                     #                                   lwd = border_width)) # Create a rectangular grob with black border
                     # grid.draw(border_rect) # Draw the rectangular grob
                   } 
                 },
                 #  border_gp = gpar(col = "black", lty = 1),
                 # rect_gp = gpar(col = "white", lwd = 0.5),
                 row_split = as.integer(myinf$cancer),
                 column_split = as.integer(myinf$cancer),
                 left_annotation =
                   rowAnnotation(
                     "Cancer Type" = anno_block(gp =
                                                  gpar(col = unique(myinf$color)),
                                                # fill = unique(myinf$color)),
                                                labels = unique(myinf$cancer),
                                                labels_gp = gpar(col =
                                                                   unique(myinf$color),
                                                                 fontsize = 20))),
                 row_title=NULL,
                 
                 top_annotation =
                   columnAnnotation(
                     "Cancer Type" = anno_block(gp =
                                                  gpar(col = unique(myinf$color)),
                                                # fill = unique(myinf$color)),
                                                labels = unique(myinf$cancer),
                                                labels_gp = gpar(col =
                                                                   unique(myinf$color),
                                                                 fontsize = 20))),
                 column_title=NULL,
                 
                 right_annotation = 
                   rowAnnotation(
                     Dataset = anno_text(myinf$newname,
                                         just = "left", 
                                         location = unit(0, "npc"), 
                                         show_name = F)),
                 show_row_names =F,
                 
                 bottom_annotation = 
                   columnAnnotation(
                     Dataset = anno_text(myinf$newname,
                                         just = "right", 
                                         rot=45,
                                         location = 1, 
                                         show_name = F)),
                 show_column_names =F,
                 
                 row_gap = unit(1, "mm"),
                 column_title_gp = gpar(fontsize = 15),
                 width = ncol(p.data)*unit(10, "mm"),
                 height = nrow(p.data)*unit(10, "mm"))

  # export plot
  pdf(myoutf1,
      width = ncol(p.data)*0.6,  
      height = nrow(p.data)*0.6)
  print(p1)
  dev.off()
  













