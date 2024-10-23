######################################## Script for Creating ADMIXTURE Plots  ########################################################

#Using pophelper library as described in Francis 2016 - Molecular Ecology Resources
#Manual found at: royfrancis.github.io/pophelper/#1_introduction
#ADMIXTURE plots with no missing data & LD-pruned

#CVs: 1000 x 1000

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
#library(devtools) #v.2.4.5
#devtools::install_github('royfrancis/pophelper')
library(tidyverse) #v.2.0.0
library(pophelper) #v.2.3.1
library(here) #v.1.0.1

setwd(here()) #set working directory
getwd() #check

###############################################################################################################################################################
  
######## Gmi ADMIXTURE plots ########

#read in data
preHWE_Gmi_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
  preHWE_Gmi_alist <- readQ(files = preHWE_Gmi_afiles)
Gmi_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD/", full.names = TRUE)
  Gmi_alist <- readQ(files = Gmi_afiles)
Gmi_Ham_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_Ham", full.names = TRUE)
  Gmi_Ham_alist <- readQ(files = Gmi_Ham_afiles)
Gmi_A_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_A/", full.names = TRUE)
  Gmi_A_alist <- readQ(files = Gmi_A_afiles)
Gmi_B_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_B/", full.names = TRUE)
  Gmi_B_alist <- readQ(files = Gmi_B_afiles)
Gmi_A_nohighhet_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_A_nohighhet", full.names = TRUE)
  Gmi_A_nohighhet_alist <- readQ(files = Gmi_A_nohighhet_afiles)
Gmi_A_Ham_nohighhet_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_A_nohighhet_Ham", full.names = TRUE)
  Gmi_A_Ham_nohighhet_alist <- readQ(files = Gmi_A_Ham_nohighhet_afiles)
Gmi_A_Ham_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_A_Ham", full.names = TRUE)
  Gmi_A_Ham_alist <- readQ(files = Gmi_A_Ham_afiles)
Gmi_A_Ham_nohighhet_unrelated_afiles <- list.files(path = "Data/Gmi_Ham/ADMIXTURE/noLD_A_nohighhet_Ham_unrelated", full.names = TRUE)
  Gmi_A_Ham_nohighhet_unrelated_alist <- readQ(files = Gmi_A_Ham_nohighhet_unrelated_afiles)
  
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

CVs <- c(0.48640, 0.18673, 0.19500, 0.20421, 0.20758)
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
            clustercol = c("#1c3b0e", "#afc8a4"),
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

CVs <- c(0.39831, 0.20585, 0.21609, 0.22534, 0.22204)
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

K2 <- plotQ(Gmi_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs")

#### Ham only ####

grplab <- c(rep("Hist - Hamilo", 30), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.45877, 0.21651, 0.23161, 0.23225, 0.24689)
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

K2 <- plotQ(Gmi_Ham_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, Ham")

#### population A ####

grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 19),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 66)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.22463, 0.22042, 0.23281, 0.24408, 0.25445)
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

K2 <- plotQ(Gmi_A_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A")

#### population B ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 28)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.21830, 0.20489, 0.20301, 0.22874, 0.25915)
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
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species B")

#### population A nohighhet ####

grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 18),
            rep("Contemp - Basud", 78), rep("Contemp - Hamilo", 63)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.21825, 0.22445, 0.23477, 0.24802, 0.26028)
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
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, species A, nohighhet")

#### population A Ham nohighhet ####

grplab <- c(rep("Hist - Hamilo", 18), rep("Contemp - Hamilo", 63)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.24366, 0.26473, 0.29437, 0.32086, 0.34861)
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
#for final filtered dataset, go up to K = 5

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

#### population A Ham ####

grplab <- c(rep("Hist - Hamilo", 19), rep("Contemp - Hamilo", 66)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.25189, 0.25137, 0.27735, 0.30620, 0.33059)
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

K2 <- plotQ(Gmi_A_Ham_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#1c3b0e", "#afc8a4"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, species A, Ham")

###############################################################################################################################################################

######## Ela ADMIXTURE plots ########

#read in data
preHWE_Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
  preHWE_Ela_alist <- readQ(files = preHWE_Ela_afiles)
Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD/", full.names = TRUE)
  Ela_alist <- readQ(files = Ela_afiles)
Ela_Ela_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Ela/", full.names = TRUE)
  Ela_Ela_alist <- readQ(files = Ela_Ela_afiles)
Ela_Lle_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Lle/", full.names = TRUE)
  Ela_Lle_alist <- readQ(files = Ela_Lle_afiles)
Ela_Ela_nohighhet_afiles <- list.files(path = "Data/Ela_Ham/ADMIXTURE/noLD_Ela_nohighhet/", full.names = TRUE)
  Ela_Ela_nohighhet_alist <- readQ(files = Ela_Ela_nohighhet_afiles)

#### preHWE ####

grplab <- c(rep("Hist - Hamilo", 30), rep("Contemp - Hamilo", 95)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.33754, 0.25972, 0.28040, 0.29909, 0.31035)
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
            clustercol = c("#8aa9be", "#16537e"),
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

CVs <- c(0.31401, 0.24902, 0.26915, 0.28819, 0.29737)
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

K2 <- plotQ(Ela_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#8aa9be", "#16537e"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs")

#### Ela species ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 94)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.26139, 0.28555, 0.30825, 0.33280, 0.34896)
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

K2 <- plotQ(Ela_Ela_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#16537e", "#8aa9be"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 1, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1, 
            showtitle = FALSE, titlelab = "ADMIXTURE plot", showsubtitle = FALSE, subtitlelab = "LD-pruned SNPs, Ela species")

#### Lle species ####

grplab <- c(rep("Hist - Hamilo", 19), rep("Contemp - Hamilo", 1)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.20447, 0.25217, 0.31195, 0.35439, 0.42894)
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

K2 <- plotQ(Ela_Lle_alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath = "Plots/PCAs_ADMIXTURE",
            clustercol = c("#8aa9be", "#16537e"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "LD-pruned SNPs, Lle species")

#### Ela nohighhet species ####

grplab <- c(rep("Hist - Hamilo", 11), rep("Contemp - Hamilo", 92)) #creates group label (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
  meta.data$loc <- as.character(meta.data$loc)

## ADMIXTURE cross-validation scores ##

CVs <- c(0.26176, 0.28598, 0.31130, 0.33626, 0.35778)
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
#for final filtered dataset, go up to K = 5

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