#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# Pan-cancer immune landscape
# Celltype: CellType_Zheng2021_Tcell
# pan_cancer_T_Cell_landscape
# 11/19/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# TimiGP #######################################################################
# Clear the environment
rm(list = ls())

# load package
library(TimiGP)
# config -----------------------------------------------------------------------
myconfig <- "./Fig2/code/config.rda"
load(myconfig)
config <- config [which(config$CI_geneset == "Zheng2021"),]

# input & output dir ------------------------------------------------------------
# please following the instruction to create the input directory
# path to instruction: ./data/README.md
myind <- "./data/bulk_transcriptomics/"
myoutd <- "./Fig2/result/pan_cancer_T_Cell_landscape/"
dir.create(myoutd)

myind <- list.files(myind,full.names = T,pattern = "2")
myinf1 <- paste0(myind,"/",config$oldname, ".rda")

dataset <-  config$oldname

data("CellType_Zheng2021_Tcell")
geneset <- CellType_Zheng2021_Tcell
marker <- unique(geneset$Gene)
cell_pair <- TimiCellPair(geneset = geneset,core = 20)

for (i in 1:nrow(config)) {
  # load data ------------------------------------------------------------
  load(myinf1[i])
  #  output ------------------------------------------------------------
  myoutd1 <- paste0(myoutd,dataset[i])
  dir.create(myoutd1)

  myoutf1 <- paste0(myoutd1,"/",dataset[i],"_MPS.rda")
  myoutf2 <- paste0(myoutd1,"/",dataset[i],"_FISHER_MP.rda")
  myoutf3 <- paste0(myoutd1,"/",dataset[i],"_ENRICH_responder.rda")
  myoutf4 <- paste0(myoutd1,"/",dataset[i],"_dotplot.pdf")
  myoutf5 <- paste0(myoutd1,"/",dataset[i],"_circle.pdf")
  myoutf6 <- paste0(myoutd1,"/",dataset[i],"_score.pdf")

  cat("\nResponder: ",i, dataset[i], "info",
      "R/NR" ,"rna", "id_pair","range(rna)\n")
  print(dim(info))
  print(table(info$Response))
  print(dim(rna))
  print(all(colnames(rna) == rownames(info)))
  print(range(rna))

  if (max(rna) < 100){
    log = FALSE
  } else {
    log = TRUE
  }
  cat("\n",dataset[i], max(rna)," log = ", log,"\n")
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
   if (dataset[i] == "Liu_phs000452.v3.p1") {
    cat("mannually set  GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
        sum(fisher_res$PV < 0.05)/nrow(fisher_res),")")
    GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()
  } else  if (dataset[i] == "Wolf_GSE194040_ISPY2_ptx_pem_non-TNBC") {
    cat("mannually set  GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
        sum(fisher_res$PV < 0.05)/nrow(fisher_res) ,")")
    GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()
  } else if (dataset[i] == "Wolf_GSE194040_ISPY2_paclitaxel_non-TNBC") {
    cat("mannually set  GP_cutoff: QV < 0.01(",sum(fisher_res$QV < 0.01), ":",
        sum(fisher_res$QV < 0.01)/nrow(fisher_res) ,")")
    GP <- fisher_res %>% filter(QV < 0.01)  %>% rownames()
  } else if (dataset[i] == "Wolf_GSE194040_ISPY2_ptx_pem_TNBC") {
    cat("mannually set GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
        sum(fisher_res$PV < 0.05)/nrow(fisher_res) ,")")
    GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()
  } else if (dataset[i] == "Wolf_GSE194040_ISPY2_paclitaxel_TNBC") {
     cat("mannually set GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
        sum(fisher_res$PV < 0.05)/nrow(fisher_res) ,")")
    GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()
  } else if(sum(fisher_res$QV < 0.05)/nrow(fisher_res) > 0.01) {

    cat("GP_cutoff: QV < 0.05(",sum(fisher_res$QV < 0.05), ":",
        sum(fisher_res$QV < 0.05)/nrow(fisher_res) ,")")
    GP <- fisher_res %>% filter(QV < 0.05)  %>% rownames()

  } else if(sum(fisher_res$PV < 0.01)/nrow(fisher_res) > 0.01) {

    cat("GP_cutoff: PV < 0.01(",sum(fisher_res$PV < 0.01), ":",
        sum(fisher_res$PV < 0.01)/nrow(fisher_res),")")
    GP <- fisher_res %>% filter(PV < 0.01)  %>% rownames()

  } else if (sum(fisher_res$PV < 0.05)/nrow(fisher_res) > 0) {

    cat("GP_cutoff: PV < 0.05(",sum(fisher_res$PV < 0.05), ":",
        sum(fisher_res$PV < 0.05)/nrow(fisher_res),")")
    GP <- fisher_res %>% filter(PV < 0.05)  %>% rownames()

  } else {

    cat("!!! No IMGP associated responder identified")
    next
  }
  # enrichment analysis --------------------------------------------------------
  background <- TimiBG(marker.pair = row.names(fisher_res))

  enrich_res <- TimiEnrich(gene = GP, background = background,
                     geneset = cell_pair, p.adj = "BH",core=20)

  save(enrich_res, file = myoutf3)
  # figures --------------------------------------------------------------------
    NET <- TimiCellNetwork(resdata = enrich_res,
                         dataset = config$CI_geneset[i],
                         condition = config$CI_condition[i],
                         cutoff = config$CI_cutoff[i],
                         export =TRUE, path = myoutd1)
  
  nn <- min(enrich_res %>% filter(Adjust.P.Value < 0.05) %>% nrow(),10)
  p1 <- TimiDotplot(resdata = enrich_res,select = c(1:nn))
  pdf(myoutf4,width = 12,height = 8)
  print(p1)
  dev.off()


  pdf(myoutf5,width = 20,height = 21)
    TimiCellChord(resdata = enrich_res,
                condition = config$CI_condition[i],
                cutoff = config$CI_cutoff[i],
                dataset = config$CI_geneset[i])
  dev.off()

  score <- TimiFS(enrich_res,                         
                  condition = config$CI_condition[i],
                  cutoff = config$CI_cutoff[i])
  head(score)
  p2 <- TimiFSBar(score)
  pdf(myoutf6,width = 14,height = 10)
  print(p2)
  dev.off()

  cat("\nvvv-----Analysis-----Done----", dataset[i],"\n")
  rm(myoutf1,myoutf2,myoutf3,myoutf4,myoutf5,myoutf6,
    rna,info,tmp_res,mps,fisher_res,GP,
     background,enrich_res,nn,p1,p2)

}
#############
sessionInfo()
