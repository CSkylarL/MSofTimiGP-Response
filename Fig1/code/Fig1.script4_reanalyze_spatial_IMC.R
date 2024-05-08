#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# 10/06/2023
# Analyze spatial data
# calculate the association between cell ratio and responder 
# using logistic regression
# from IMC data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up directories
# Please set the working directory to the git repository
# Using
# setwd("~/path/to/git/MSofTimiGP-Response")
# Re-analyze spatial ###########################################################

# Clear environment 
rm(list = ls())

# load libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(fst)
library(parallel)
library(glmnet)
library(reshape)

# input data
# Please download the IMC data from https://zenodo.org/record/7990870
# Save the data in the folder "data"
# And rename the folder to "TNBC_IMC_zenodo.7990870"

myind <- "./data/TNBC_IMC_zenodo.7990870/"
myfun1 <- "./Fig1/code/Fig1_function.R"

# The output directory
myoutd <- "./Fig1/result/TNBC_spatial_ratio_responder/"
dir.create(myoutd)
myoutf1 <- paste0(myoutd,"/Spatial_cell_ratio_response.rda")
source(myfun1)

# Get spatial data ======================================================

clinical <- getClinical(dir = myind)[(isPerProtocol)]
clinical$Response <- ifelse(clinical$pCR == "pCR",1,0)
table(clinical$pCR, exclude = F)

# RD pCR 
# 310 314 

clinical$Response <- ifelse(clinical$pCR == "pCR",1,0)
table(clinical$Response)
# 0   1 
# 310 314 


# cellClusters <- getCellClusters(dir = myind)[['All']]
# cellCounts <- getCellCounts(dir = myind)

# Mean cell phenotype density
Densities <- getDensities(dir = myind)
Densities <- merge(x = Densities,
                   y = clinical[, .(PatientID, BiopsyPhase, pCR, Arm)],
                   by = c('PatientID', 'BiopsyPhase'))
Densities[, predictor := sqrt(Density)]

# Preprocess spatial data ======================================================
cell <- c("CD20+B" , "CD79a+Plasma",
          "CD4+TCF1+T" ,"CD4+PD1+T" ,"Treg" , 
          "CD8+TCF1+T" ,"CD8+GZMB+T" ,"CD8+PD1+Tex" ,
          "CD56+NK" ,"cDC" , "pDC" , 
          "M1" , "M2" , 
          "Mast" ) %>% 
  gsub(pattern = "+",replacement = ".", fixed = T) 

# Update cell name
Densities$Label <- Densities$Label  %>% 
  gsub(pattern = "^",replacement = "", fixed = T) %>% 
  gsub(pattern = "_{E",replacement = "e", fixed = T) %>% 
  gsub(pattern = " Mac",replacement = "", fixed = T) %>% 
  gsub(pattern = "DCs",replacement = "cDC", fixed = T) %>% 
  gsub(pattern = "}",replacement = "", fixed = T) %>% 
  gsub(pattern = "+",replacement = ".", fixed = T) 

se <- cell %in% unique(Densities$Label)
cell <- cell[se]

# separate data by Arm and BiopsyPhase
arm <- unique(clinical$Arm)
timepoint <- unique(clinical$BiopsyPhase)
mydata <- list()
for (ii in arm) {
  for (jj in timepoint) {
    cat("\n", as.character(ii), "&", as.character(jj))
    info <- clinical %>% 
      filter(Arm == ii) %>% 
      filter(BiopsyPhase == jj) %>%
      data.frame()

    rownames(info) <- as.character(info$PatientID)
    
    cat("\t", nrow(mydata[[as.character(ii)]][[as.character(jj)]]$info))
    tmp <- Densities %>%
      filter(Label %in% cell) %>%
      filter(PatientID %in%  info$PatientID) %>%
      filter(BiopsyPhase == jj)
    
     count <-  tmp %>%
      select(Label, PatientID, CellCountPerPatient) %>% 
      cast(PatientID ~ Label, value = "CellCountPerPatient")  
     
     rownames(count) <- as.character(count$PatientID)
    
     density <- tmp %>%
       select(Label, PatientID, Density) %>% 
       cast(PatientID ~ Label, value = "Density")  
     
     rownames(density) <- as.character(density$PatientID)
     
     comxx <- intersect(rownames(info),rownames(count))
    mydata[[as.character(ii)]][[as.character(jj)]]$info <- info[comxx,-1]
    mydata[[as.character(ii)]][[as.character(jj)]]$count <- count[comxx,-1]
    mydata[[as.character(ii)]][[as.character(jj)]]$density <- density[comxx,-1]
    rm(count,info,tmp,density)

  }
  
  
}

# calculate Ratio ==============================================================

# Generate cell pair
cell_pair <- outer(cell,cell, paste, sep="_") 
cell_pair <- cell_pair[lower.tri(cell_pair)]

xx <- cell_pair
tmp <- unlist(strsplit(xx, "_"))
nn <- length(xx)
tmp1 <- tmp[(1:nn)*2-1]
tmp2 <- tmp[(1:nn)*2-0]
B_A <- paste(tmp2, tmp1, sep="_")
 
# Calculate ratio
for (ii in arm) {
  for (jj in timepoint) {
    cat("\n", as.character(ii), "&", as.character(jj))
    ratioA_B <- ratioB_A <- matrix(
      nrow = nrow(mydata[[as.character(ii)]][[as.character(jj)]]$count),
                    ncol=length(cell_pair))

    for (kk in 1:length(cell_pair)) {
      cellA <- strsplit(cell_pair[kk], "_") %>% unlist() %>% .[1]
      cellB <- strsplit(cell_pair[kk], "_") %>% unlist() %>% .[2]
      logA <- log2(as.numeric(mydata[[as.character(ii)]][[as.character(jj)]]$count[,cellA])+1)
      logB <- log2(as.numeric(mydata[[as.character(ii)]][[as.character(jj)]]$count[,cellB])+1)
      ratioA_B[,kk] <- logA-logB

    }
    rownames(ratioA_B) <- rownames(mydata[[as.character(ii)]][[as.character(jj)]]$count)
    colnames(ratioA_B) <- cell_pair
    
    ratioB_A <- -ratioA_B

    colnames(ratioB_A) <- B_A
    mydata[[as.character(ii)]][[as.character(jj)]]$ratioA_B <- ratioA_B
    mydata[[as.character(ii)]][[as.character(jj)]]$ratioB_A <- ratioB_A
    
  }
}

# Logistic regression =============================================================
# cell -------------------------------------------------------------------------
for (ii in arm) {
  for (jj in timepoint) {
    cat("\n", as.character(ii), "&", as.character(jj))
    name <- OR <- PV <- numeric(
      length = length(cell))
    name <- as.character(name)
    data <- cbind(mydata[[as.character(ii)]][[as.character(jj)]]$density,
                  mydata[[as.character(ii)]][[as.character(jj)]]$info) %>%
      as.data.frame()
    
    
    for (kk in 1:length(cell)) {
      
      myformular <- as.formula(paste0("Response ~ ", cell[kk]))
      logistic <- glm(myformular, data=data,family = "binomial")
      logistic <- summary(logistic)
      # Odds Ratio
        name[kk] <- cell[kk]
        OR[kk] <- exp(logistic$coefficients[2,1])
        PV[kk] <- logistic$coefficients[2,4]
       
      
    }
    log_res <- data.frame(OR,PV)
    log_res$QV <- p.adjust(log_res$PV, method="BH")
    rownames(log_res) <- name
    print(log_res)
    mydata[[as.character(ii)]][[as.character(jj)]]$log_res_cell <- log_res %>% arrange(PV)
    rm(log_res)
    
  }
}
# ratio ------------------------------------------------------------------------
for (ii in arm) {
  for (jj in timepoint) {
    cat("\n", as.character(ii), "&", as.character(jj))
    name <- OR <- PV <- numeric(
      length = length(cell_pair))
    name <- as.character(name)
    data <- cbind(mydata[[as.character(ii)]][[as.character(jj)]]$ratioA_B,
                  mydata[[as.character(ii)]][[as.character(jj)]]$info) %>%
      as.data.frame()
      
    
    for (kk in 1:length(cell_pair)) {
     
      myformular <- as.formula(paste0("Response ~ ", cell_pair[kk]))
      logistic <- glm(myformular, data=data,family = "binomial")
      logistic <- summary(logistic)
      # Odds Ratio
      or_tmp <- logistic$coefficients[2,1]
      if(or_tmp >= 0) {
        name[kk] <- cell_pair[kk]
        OR[kk] <- exp(or_tmp)
        PV[kk] <- logistic$coefficients[2,4]
      } else {
        name[kk] <- B_A[kk]
        OR[kk] <- exp(-or_tmp)
        PV[kk] <- logistic$coefficients[2,4]
      }
     
      
    }
    log_res <- data.frame(OR,PV)
    log_res$QV <- p.adjust(log_res$PV, method="BH")
    rownames(log_res) <- name
    mydata[[as.character(ii)]][[as.character(jj)]]$log_res_ratio <- log_res %>% arrange(PV)
    rm(log_res)
    
  }
}
spatial_ref <- mydata
save(spatial_ref, file = myoutf1)