######################################## Script for Creating ADMIXTURE Plots  ########################################################

#using pophelper library as described in Francis 2016 - Molecular Ecology Resources
#manual found at: royfrancis.github.io/pophelper/#1_introduction
#ADMIXTURE plots with no missing data & LD-pruned

#CVs: 1000 x 1000

#################################################################################################################################################

######## Set-up ########

getwd() #check working directory
remove(list = ls())

#load libraries
#library(devtools) #v.2.4.5
#devtools::install_github('royfrancis/pophelper')
library(tidyverse) #v.2.0.0
library(pophelper) #v.2.3.1

###############################################################################################################################################################
  
######## Gmi ADMIXTURE plots ########

#read in data
preHWE_Gmi_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
  preHWE_Gmi_alist <- readQ(files = preHWE_Gmi_afiles)
Gmi_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD/", full.names = TRUE)
  Gmi_alist <- readQ(files = Gmi_afiles)
Gmi_Ham_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_Ham/", full.names = TRUE)
  Gmi_Ham_alist <- readQ(files = Gmi_Ham_afiles)
Gmi_A_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popA/", full.names = TRUE)
  Gmi_A_alist <- readQ(files = Gmi_A_afiles)
Gmi_A_Ham_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popA_Ham", full.names = TRUE)
  Gmi_A_Ham_alist <- readQ(files = Gmi_A_Ham_afiles)
Gmi_A_nohighhet_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popA_nohighhet", full.names = TRUE)
  Gmi_A_nohighhet_alist <- readQ(files = Gmi_A_nohighhet_afiles)
Gmi_B_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popB/", full.names = TRUE)
  Gmi_B_alist <- readQ(files = Gmi_B_afiles)
Gmi_A_Ham_nohighhet_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popA_Ham_nohighhet", full.names = TRUE)
  Gmi_A_Ham_nohighhet_alist <- readQ(files = Gmi_A_Ham_nohighhet_afiles)
Gmi_A_Bas_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_popA_Bas", full.names = TRUE)
  Gmi_A_Bas_alist <- readQ(files = Gmi_A_Bas_afiles)
  
#### preHWE ####
  
grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 30),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##
#found in *log.out files
#lowest CV error = best structure (most likely # of demes)
#K1, K2, K3, K4, K5

CVs <- c(0.48439, 0.18692, 0.19556, 0.20465, 0.20720)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for preHWE, only going up to K = 2 (to catch cryptic structure)

K2 <- plotQ(preHWE_Gmi_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

#### postHWE ####

grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 30),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.39575, 0.20567, 0.21540, 0.21424, 0.22204)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K3 <- plotQ(Gmi_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K4 <- plotQ(Gmi_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K5 <- plotQ(Gmi_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

#### postHWE Ham ####

grplab <- c(rep("Hist - Hamilo", 30), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.45593, 0.21492, 0.23158, 0.23004, 0.24711)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_Ham_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, Ham")

#### population A ####

grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 19),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 64)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.22062, 0.22156, 0.23430, 0.24612, 0.25743)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_A_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A")

K3 <- plotQ(Gmi_A_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A")

K4 <- plotQ(Gmi_A_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A")

K5 <- plotQ(Gmi_A_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A")

#### population A Ham ####

grplab <- c(rep("Hist - Hamilo", 19), rep("Contemp - Hamilo", 64)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.24576, 0.25925, 0.28455, 0.31205, 0.34407)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_A_Ham_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, species A, Ham")


#### population A nohighhet ####

grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 18),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 61)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.21607, 0.22429, 0.23638, 0.25181, 0.26277)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_A_nohighhet_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, nohighhet")

K3 <- plotQ(Gmi_A_nohighhet_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, nohighhet")

K4 <- plotQ(Gmi_A_nohighhet_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, nohighhet")

K5 <- plotQ(Gmi_A_nohighhet_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, nohighhet")

#### population B ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 30)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.25708, 0.22549, 0.21229, 0.23234, 0.26028)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_B_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species B")

K3 <- plotQ(Gmi_B_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species B")

K4 <- plotQ(Gmi_B_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species B")

K5 <- plotQ(Gmi_B_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species B")

#### population A Ham nohighhet ####

grplab <- c(rep("Hist - Hamilo", 18), rep("Contemp - Hamilo", 61)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.23858, 0.26493, 0.29628, 0.32688, 0.35243)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_A_Ham_nohighhet_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, species A, Ham nohighhet")

K3 <- plotQ(Gmi_A_Ham_nohighhet_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Ham nohighhet")

K4 <- plotQ(Gmi_A_Ham_nohighhet_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Ham nohighhet")

K5 <- plotQ(Gmi_A_Ham_nohighhet_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Ham nohighhet")

#### population A Bas ####

grplab <- c(rep("Hist - Basud", 12), rep("Contemp - Basud", 78)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.22387, 0.24424, 0.26895, 0.29168, 0.31514)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Gmi_A_Bas_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Bas")

K3 <- plotQ(Gmi_A_Bas_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Bas")

K4 <- plotQ(Gmi_A_Bas_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Bas")

K5 <- plotQ(Gmi_A_Bas_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, Bas")


###############################################################################################################################################################

######## Ela ADMIXTURE plots ########

#read in data
preHWE_Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
  preHWE_Ela_alist <- readQ(files = preHWE_Ela_afiles)
Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD/", full.names = TRUE)
  Ela_alist <- readQ(files = Ela_afiles)
Ela_Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Ela/", full.names = TRUE)
  Ela_Ela_alist <- readQ(files = Ela_Ela_afiles)
Ela_Ela_nohighhet_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Ela_nohighhet/", full.names = TRUE)
  Ela_Ela_nohighhet_alist <- readQ(files = Ela_Ela_nohighhet_afiles)
Ela_Lle_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Lle/", full.names = TRUE)
  Ela_Lle_alist <- readQ(files = Ela_Lle_afiles)

#### preHWE ####

grplab <- c(rep("Hist - Hamilo", 30), rep("Contemp - Hamilo", 95)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.34568, 0.26414, 0.28227, 0.30269, 0.30993)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for preHWE, only going up to K = 2 (to catch cryptic structure)

K2 <- plotQ(preHWE_Ela_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

#### postHWE ####

grplab <- c(rep("Hist - Hamilo", 30), rep("Contemp - Hamilo", 95)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.31801, 0.24937, 0.26883, 0.28747, 0.29566)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Ela_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#16537e", "#8aa9be"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs")

K3 <- plotQ(Ela_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K4 <- plotQ(Ela_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K5 <- plotQ(Ela_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

#### Ela species ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.26160, 0.28282, 0.30614, 0.33166, 0.35005)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Ela_Ela_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#16537e", "#8aa9be"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, Ela species")

K3 <- plotQ(Ela_Ela_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species")

K4 <- plotQ(Ela_Ela_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species")

K5 <- plotQ(Ela_Ela_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species")

#### Ela nohighhet species ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 92)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.26188, 0.28546, 0.30789, 0.33534, 0.35679)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Ela_Ela_nohighhet_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#16537e", "#8aa9be"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, Ela species, nohighhet")

K3 <- plotQ(Ela_Ela_nohighhet_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species, nohighhet")

K4 <- plotQ(Ela_Ela_nohighhet_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species, nohighhet")

K5 <- plotQ(Ela_Ela_nohighhet_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Ela species, nohighhet")

#### Lle species ####

grplab <- c(rep("Hist - Hamilo", 19), rep("Contemp - Hamilo", 1)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.21137, 0.26199, 0.31604, 0.36322, 0.42057)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Ela_Lle_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Lle species")

K3 <- plotQ(Ela_Lle_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Lle species")

K4 <- plotQ(Ela_Lle_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Lle species")

K5 <- plotQ(Ela_Lle_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Lle species")

###############################################################################################################################################################

######## Aeb ADMIXTURE plots ########

#read in data
preHWE_Aen_afiles <- list.files(path = "Data/Aen_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
  preHWE_Aen_alist <- readQ(files = preHWE_Aen_afiles)
Aen_afiles <- list.files(path = "Data/Aen_Ham/ADMIXTURE/noLD/", full.names = TRUE)
  Aen_alist <- readQ(files = Aen_afiles)
Aen_A_afiles <- list.files(path = "Data/Aen_Ham/ADMIXTURE/noLD_popA/", full.names = TRUE)
  Aen_A_alist <- readQ(files = Aen_A_afiles)
Aen_B_afiles <- list.files(path = "Data/Aen_Ham/ADMIXTURE/noLD_popB/", full.names = TRUE)
  Aen_B_alist <- readQ(files = Aen_B_afiles)

#### preHWE ####

grplab <- c(rep("Hist - Hamilo", 57), rep("Contemp - Hamilo", 95)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.22404, 0.18235, 0.17828, 0.19661, 0.14966)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#going up to 5 here since the CV is weird

K2 <- plotQ(preHWE_Aen_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K3 <- plotQ(preHWE_Aen_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K4 <- plotQ(preHWE_Aen_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K5 <- plotQ(preHWE_Aen_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

#### postHWE ####

grplab <- c(rep("Hist - Hamilo", 57), rep("Contemp - Hamilo", 95)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.22983, 0.18280, 0.18622, 0.19139, 0.20288)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Aen_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#a25505", "#e3ccb4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs")

K3 <- plotQ(Aen_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K4 <- plotQ(Aen_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

K5 <- plotQ(Aen_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

#### pop A ####

grplab <- c(rep("Hist - Hamilo", 55), rep("Contemp - Hamilo", 92)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.12015, 0.11453, 0.12527, 0.12715, 0.13940)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line(linewidth = 2) + 
  geom_point(size = 8) + 
  theme_bw() + 
  labs(y = "Cross-Validation Error", x = "K") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Aen_A_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#a25505", "#e3ccb4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, pop A")

K3 <- plotQ(Aen_A_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop A")

K4 <- plotQ(Aen_A_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop A")

K5 <- plotQ(Aen_A_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop A")

#### pop B ####
#this is weird, bc only 5 individuals

grplab <- c(rep("Hist - Hamilo", 2), rep("Contemp - Hamilo", 3)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(1.93109, 1.73175, 1.59797, 0.17989, 0.04782)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

#CV plot
CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", linewidth = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c( .5, .5, .5, .5), "cm"))
CV_plot

## ADMIXTURE plots ##
#for postHWE, go up to K = 5

K2 <- plotQ(Aen_B_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop B")

K3 <- plotQ(Aen_B_alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#9999FF", "#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop B")

K4 <- plotQ(Aen_B_alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop B")

K5 <- plotQ(Aen_B_alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, pop B")
