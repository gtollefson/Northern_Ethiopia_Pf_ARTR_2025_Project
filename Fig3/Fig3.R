################################################################################
#
# Script for Figure 3: Upset plots for prevalence of the Pfk13 
# 622I mutation and other drug-resistance mutation combinations 
# in monoclonal P. falciparum isolates related to lumefantrine 
# susceptibility
#
# Author: Abebe A Fola
# June 2, 2025
#
################################################################################

################################################################################
# Fig3A - MDR1_622I mutations prevalence
################################################################################
# Load required libraries
library (ggplot2)
library(ggpubr)
library(UpSetR)

getwd() # To check the current directory
setwd("C:/Users/PC/Desktop/Upsetplot") # change this to your current directory

Fig3A <- read.csv("Fig3Ainputfile.csv", header=TRUE, sep="," ) # Upload your file 
head(Fig3A)

#plot upsetplot
bar_cols1 <- c("black", "#999999", "#56B4E9", "#009E73", "#E69F00") # you can modify colors 
upset(Fig3A, nsets = 5, nintersects = 30, mb.ratio = c(0.5, 0.5),
      text.scale = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5),
      order.by = c("freq", "degree"), decreasing = c(FALSE, FALSE),
      keep.order = T, main.bar.color = bar_cols1,
      point.size = 3.5, line.size = 1)

# Save plot
ggsave("Fig3A.svg", dpi=600, width=7.5, height=6)
ggsave("Fig3A.pdf", dpi=600, width=7.5, height=6)

#################
# Fig3B - CRT_622I mutations prevalence
#################
# Load required libraries
library (ggplot2)
library(ggpubr)
library(UpSetR)

getwd() # To check the current directory
setwd("C:/Users/PC/Desktop/Upsetplot") #change this to your current directory

Fig3B <- read.csv("Fig3Binputfile.csv", header=TRUE, sep="," ) # Upload your file 
head(Fig3B)

# plot upsetplot
bar_cols1 <- c("black", "#999999", "#56B4E9", "#009E73", "#E69F00", "#F0E442", "#D55E00") # you can modify colors as you want
upset(Fig3B, nsets = 5, nintersects = 30, mb.ratio = c(0.5, 0.5),
      text.scale = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5),
      order.by = c("freq", "degree"), decreasing = c(FALSE, FALSE),
      keep.order = T, main.bar.color = bar_cols1,
      point.size = 3.5, line.size = 1)

# Save plot
ggsave("Fig3B.svg", dpi=600, width=7.5, height=6)
ggsave("Fig3B.pdf", dpi=600, width=7.5, height=6)
