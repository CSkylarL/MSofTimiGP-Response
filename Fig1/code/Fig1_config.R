#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Skylar Li
# 10/06/2023
# Figure 1 visualization setting
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# settings for visualization  --------------------------------------------------
# Cell names
cell <- c("CD20+B" , "CD79a+Plasma",
          "CD4+TCF1+T" ,"CD4+PD1+T" ,"Treg" , 
          "CD8+TCF1+T" ,"CD8+GZMB+T" ,"CD8+PD1+Tex" ,
          "CD56+NK" ,"cDC" , "pDC" , 
          "M1" , "M2" , 
          "Mast" ) 

# Cell group
group <- c(rep("B",2),
           rep("CD4T", 3),
           rep("CD8T", 3),
           "NK",
           rep("DC", 2),
           rep("Mac", 2),
           "Mast")
# Cell color
color <- c("#66c2a4","#238b45",
           "#bcbddc","#807dba","#54278f",
           "#fa9fb5","#dd3497","#7a0177",
           "#2b8cbe",
           "#225ea8","#081d58",
           "#feb24c","#ec7014",
           "#e31a1c")
names(group) <- cell
names(color) <- cell

# Protein names
protein <- c("CD45", 
  "CD20","CD79a", "HLA-ABC", "HLA-DR",
  "CD3", "CD4", "CD8", "GZMB", "TCF1","Helios", "FOXP3", "GATA3", "T-bet","ICOS", "PD-1", "TOX","OX40",
  "CD56",
  "CD11c","CD68", "CD163")

# Gene names
gene <-  c( "PTPRC", #"CD45"
            "MS4A1", #"CD20"
            "CD79A", # "CD79a"
            "HLA-A"   , "HLA-B" , "HLA-C"   ,# "HLA-ABC"
            "HLA-DRA" ,  "HLA-DRB5" ,    # "HLA-DR"
            "CD3E", # "CD3"
            "CD4"  , # "CD4"
            "CD8A",  # "CD8"
            "GZMB",  # "GZMB"
            "TCF7", # "TCF1"
            "IKZF2", # "helios"
            "FOXP3", # "FOXP3"
            "GATA3", # "GATA3"
            "TBX21", # "T-bet"
            "ICOS", # "ICOS"
            "PDCD1", # "PD-1"
            "TOX",  # "TOX"
            "TNFRSF4", #"OX40
            "CTLA4", "LAG3", "TIGIT", # additional exhausted markers
            "NCAM1", # "CD56"
            "NCR1", "GNLY", "NKG7", # Additional NK Cell markers
            "ITGAX", # "CD11c", DC
            "CD68", # "CD68"
            "CD163", #  "CD163"
            "CLEC4C", # Additional pDC
            "KIT","MS4A2" # Additional Mast 
            )