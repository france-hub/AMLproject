#This script plots gene signature enrichments from Zheng et al.
#Fig. 1A

rm(list = ls())

library(rstudioapi)
library(tidyverse)
library(reshape2)
library(magrittr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# read data 
data<-readxl::read_xlsx("zhenggsea.xlsx")
colnames(data) <- gsub("_", " ", colnames(data)) #remove underscores
colnames(data) <- gsub("(.*) ", "\\1_", colnames(data)) #put underscore at the last space

data <- reshape2::melt(data) %>% separate(variable, c("group_id", "variable"), sep = "_") #melt and separate FDR and NES
data <- reshape2::dcast(data, geneset+group_id~variable) %>% 
  set_colnames(c("geneset", "group_id", "NES", "p_adj")) #cast and set colnames
data <- data %>% dplyr::mutate(FDR = case_when(p_adj < 0.05 ~ 10, TRUE ~ 1)) #add FDR column with < 0.1 = 10 otherwise 1 (to scale the size of dots)
data <- data %>% mutate(condition = case_when(grepl("AML", group_id)~ "AML bas vs HD",
                grepl("PostCRvsNR", group_id)~ "Res vs NonRes post",
                grepl("PreCRvsNR", group_id)~"Res vs NonRes bas"))
data$geneset <- gsub("Zheng_Sci2021_CD8_|Zheng_Sci2021_CD8.c03_", "", data$geneset)
#define palette
pal_exp <- colorRampPalette(c("blue", "white", "red"))(256)

data$FDR <- factor(data$FDR) #factorize to have only 2 FDR dots
levels(data$FDR) <- c(">= 0.05", "< 0.05") #change the levels to label the dots

tiff("./DotPlot.tiff", width = 5*380, height = 5*270, res = 300, pointsize = 5)     
ggplot(data, aes(x= geneset, y= condition, color=NES, group=condition)) +  geom_point(aes(size=FDR))+
  labs(x = NULL, y = NULL) + theme_classic() +  
  guides(size = guide_legend(nrow = 2, byrow = T))+
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.4, "cm")) +
  scale_color_gradientn(colours = pal_exp) 
dev.off()
