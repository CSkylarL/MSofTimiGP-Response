#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Chenyang Li
# Generate tiny pie chart as summary of favorability score
# # Generate Extended Data Fig.4, Fig 2b and c
# 02/28/2024
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# clear workspace
rm(list = ls())

# load libraries
library(ggplot2)
library(TimiGP)
library(dplyr)
library(RColorBrewer)
library(VennDiagram)
library(scatterpie)

# Input & Output ===============================================================
# load config file
myinf1 <- "./Fig2/code/config.rda"
load(myinf1)

# The output directory
myoutd <- "./Fig2/result/summary"

# get the geneset ID which determines the analysis resolution
geneset <- unique(config$CI_geneset)

for (kk in geneset) {
  cat("\n", kk)
  # output file
  myoutf1 <-  paste0(myoutd,"/immunotherapy_all_score_scatterpie_",kk,".pdf")
  myinf <- config %>% filter(CI_geneset == kk)
  
  # load data ----
  resdata <- data.frame()
  for (i in myinf$order) {
    load(myinf$location[i])
    tmp <- TimiFS(enrich_res,
                  condition = myinf$CI_condition[i],
                  cutoff = myinf$CI_cutoff[i])
    tmp$Dataset <- myinf$oldname[i]
    resdata <- rbind(resdata,tmp)
    rm(enrich_res,tmp)
  } 
  
  # plot features ----
  resdata$total <- resdata$Favorable.Score+resdata$Unfavorable.Score
  resdata$radius <- sqrt((resdata$total/pi/50))
  
  # dataset display order -----
  ds_order <- myinf$oldname[nrow(myinf):1] 
  
  # rank cell order:favorble-->unfavorable----------
  
  colnames(resdata)
  fa <- resdata %>% 
    mutate(Percentage = resdata$Favorable.Score/resdata$total*100) %>%
    filter(Percentage>50) %>%
    select(Cell.Type) %>%
    table() %>%
    data.frame() %>%
    arrange(Freq) 
  colnames(fa) <- c("cell","fa")
  unf <-  resdata %>% 
    mutate(Percentage = resdata$Unfavorable.Score/resdata$total*100) %>%
    filter(Percentage<50) %>%
    select(Cell.Type) %>%
    table() %>%
    data.frame(stringsAsFactors = F) %>%
    arrange(Freq) 
  colnames(unf) <- c("cell","unf")
  
  
  cel.order <-  merge(fa,unf,by="cell") %>%
    arrange(-fa,unf) 
  
  
  cel.order <-cel.order$cell
  
  
  # scatterpie plot
  
  resdata$x <- resdata$Cell.Type %>% 
    factor( levels=cel.order)%>% 
    as.numeric()
  resdata$y <- resdata$Dataset %>% 
    factor(levels = ds_order ) %>% 
    as.numeric()
  
  colnames(resdata)
  
  
  p1<- ggplot() + 
    geom_scatterpie(aes(x=x, y=y, r=radius), data=resdata,
                    cols=c("Favorable.Score","Unfavorable.Score"), 
                    legend_name = "Cell",
                    color=NA)  + 
    scale_fill_manual(values=c("#f46d43","#3288bd")) +
    coord_equal() +
    scale_x_discrete(limits = factor(seq(1,length(cel.order))),
                     labels = cel.order) +
    scale_y_discrete(limits = factor(seq(1,length(ds_order)+2)),
                     breaks = factor(seq(1,length(ds_order))),
                     labels = c(ds_order)) +
    theme_bw(base_size = 15, base_family = "serif") + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) +
    labs(x="Cell Type",y="Dataset",
         title=kk,
         base_size = 30, base_family = "serif",face="bold") +
    theme(
      legend.position = c(0.4,0.97),
      legend.box = "horizontal",
      legend.direction= "horizontal",
      panel.grid.major = element_line(color = "#d9d9d9",
                                      linewidth = 0.5,
                                      linetype = 2),
      panel.grid=element_blank(),
      legend.key.width = unit(0.5,"cm"),
      legend.title = element_text(face="bold", color="black",family = "serif", size=30),
      legend.text= element_text(face="bold", color="black",family = "serif", size=30),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(face="bold", color="black", size=25, 
                                 angle = 40,  vjust = 1, hjust=1),
      axis.text.y = element_text(face="bold", color="black", size=25),
      axis.title.x = element_text(face="bold", color="black", size=30),
      axis.title.y = element_text(face="bold",color="black", size=30))+
    geom_scatterpie_legend(resdata$radius, x=length(cel.order)-3, y=length(ds_order)+1.5) 
  
  # export plot
  pdf(myoutf1,
      width = 25, height = 23)
  print(p1)
  
  dev.off()
  
}



