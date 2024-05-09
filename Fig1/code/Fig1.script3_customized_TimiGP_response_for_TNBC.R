#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Skylar Li
# 10/06/2023
# Run_TimiFisher for TNBC cohort
# Wolf_GSE194040_ISPY2_ptx_pem (chemo + immuno) 
# Celltype: TNBC_scRNA_marker_for_spatial
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# TimiGP TNBC_scRNA_Spatial #######################################################################

# clear all variables
rm(list = ls())

# load the TimiGP package
library(TimiGP)

# The input files
myinf1 <-  "./data/bulk_transcriptomics/Wolf_GSE194040_ISPY2_ptx_pem_TNBC.rda"
myinf2 <- "./Fig1/result/customized_markers_for_TimiGP/TNBC_scRNA_marker_for_spatial.rda"

myconfig1 <- "./Fig1/code/Fig1_config.R"
# The output directory
myoutd <- "./Fig1/result/TimiGP_TNBC_scRNA_marker_for_spatial/"
dir.create(myoutd)

dataset <-  gsub(x = basename(myinf1),pattern = ".rda",replacement = "")

# the scRNA cell-type annotaion
load(myinf2)
geneset <- TNBC_scRNA_marker_for_spatial
unique(geneset$CellType)
marker <- unique(geneset$Gene)
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

# settings for visualization 
source(myconfig1)

# Perform TimiGP analysis
for (i in 1:length(myinf1)) {
  # input & output ------------------------------------------------------------
  load(myinf1[i])
  myoutd1 <- paste0(myoutd,dataset[i])
  dir.create(myoutd1)
  
  myoutf1 <- paste0(myoutd1,"/",dataset[i],"_TNBC_MPS.rda")
  myoutf2 <- paste0(myoutd1,"/",dataset[i],"_TNBC_FISHER_MP.rda")
  myoutf3 <- paste0(myoutd1,"/",dataset[i],"_TNBC_ENRICH_responder.rda")
  myoutf4 <- paste0(myoutd1,"/",dataset[i],"_TNBC_dotplot.pdf")
  myoutf5 <- paste0(myoutd1,"/",dataset[i],"_TNBC_circle.pdf")
  myoutf6 <- paste0(myoutd1,"/",dataset[i],"_TNBC_score.pdf")
  
  cat("\nResponder: ",i, dataset[i], "TNBC info",
      "R/NR" ,"rna", "id_pair","range(rna)\n")
  print(dim(info))
  print(table(info$Response))
  print(dim(rna))
  print(all(colnames(rna) == rownames(info)))
  print(range(rna,na.rm = T))
  
  if (max(rna,na.rm = T) < 100){
    log = FALSE
  } else {
    log = TRUE
  }
  cat("\n",dataset[i], max(rna,na.rm = T)," log = ", log,"\n")
  # association between GP and responder ---------------------------------------
  rna <- TimiPrePropress(marker = marker, cohort = rownames(info),
                         log = log, GMNorm = T, rna = rna)
  mps <- TimiGenePair(rna = rna)
  
  se <- which(colnames(info) == "Response")
  info <- info[,se,drop=F]
  tmp_res <- TimiFisher(mps = mps,info = info,p.adj = "BH")
  
  mps <- tmp_res$mps
  fisher_res <- tmp_res$fisher_res
  save(mps, file = myoutf1)
  save(fisher_res, file = myoutf2)
  
  #select GP for enrichment analysis -------------------------------------------
  
  cat("GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
      sum(fisher_res$PV < 0.05)/nrow(fisher_res) ,")")
  GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()
  
  
  # enrichment analysis --------------------------------------------------------
  background <- TimiBG(marker.pair = row.names(fisher_res))
  
  enrich_res <- TimiEnrich(gene = GP, background = background,
                           geneset = cell_pair, p.adj = "BH",core=20)
  
  save(enrich_res, file = myoutf3)
  
  #define cell-cell interaction -------------------------------------------
  CI_condition <- "Adjust.P.Value"
  CI_cutoff <- 0.05
  cat("CI_cutoff: ", CI_condition, " < ",
      CI_cutoff, " : #",sum(enrich_res[,CI_condition] < CI_cutoff))
  
  # figures --------------------------------------------------------------------
  
  NET <- TimiCellNetwork(resdata = enrich_res,
                         condition = CI_condition, 
                         cutoff = CI_cutoff,
                         geneset = geneset,
                         dataset = "Other",
                         group = group,
                         export =TRUE, path = myoutd1)
  
  nn <- min(sum(enrich_res[,CI_condition] < CI_cutoff),10)
  p1 <- TimiDotplot(resdata = enrich_res,
                    condition = CI_condition,
                    cutoff = CI_cutoff,
                    select = c(1:nn))
  pdf(myoutf4,width = 8,height = 6)
  print(p1)
  dev.off()
  
  
  pdf(myoutf5,width = 15,height = 15)
  TimiCellChord(resdata = enrich_res,
                dataset = "Other",group = group,color = color,
                condition = CI_condition, cutoff = CI_cutoff,)
  dev.off()
  
  score <- TimiFS(enrich_res,condition = CI_condition,cutoff = CI_cutoff)
  head(score)
  p2 <- TimiFSBar(score)
  pdf(myoutf6,width = 8,height = 5)
  print(p2)
  dev.off()
  
  cat("\nvvv-----Analysis-----Done----", dataset[i],"\n")
  rm(myoutf1,myoutf2,myoutf3,myoutf4,myoutf5,myoutf6,
     rna,info,tmp_res,mps,fisher_res,GP,
     background,enrich_res,nn,p1,p2)
  
}

sessionInfo()

sessionInfo()
