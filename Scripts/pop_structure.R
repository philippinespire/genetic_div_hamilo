#################################################### Script for Creating PCA & ADMIXTURE Plots  ########################################################

#adjust paths as needed

#################################################################################################################################################

######## Set-up ########

##### set working directory ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

remove(list = ls()) #clear working env

##### load libraries ####
#library(devtools)
#devtools::install_github('royfrancis/pophelper') #to install pophelper the first time
library(tidyverse)
library(gridExtra)
library(pophelper)

#### import PCA data ####
#read in data
eigenval <- read.csv("PIRE.Aen.Ham.HWE.eigenval", header = FALSE, sep = " ") #eigenvalues
data_PC <- read.csv("PIRE.Aen.Ham.HWE.eigenvec", header = FALSE, sep = " ") #eigenvectors
head(data_PC) #double check data is read in correctly 

#add column names
colnames(data_PC) <-c ('Population', 'Individual', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 
                       'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20')
data_PC$Population <- as.character(data_PC$Population) #make sure values in population column are characters

#add columns for Location & Era
#fill out data for each column based on specific conditions, then rearrange the data to have Location & Era on the first 2 columns
data_PC <- data_PC %>%
  mutate(Location =
           case_when(endsWith(Population, "Bas") ~ "Basud",
                     endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% #add lines as needed, depending on dataset
  mutate(Era =
           case_when(endsWith(Population, "ABas") ~ "Albatross",
                     endsWith(Population, "AHam") ~ "Albatross",
                     endsWith(Population,"CBat") ~ "Contemporary",
                     endsWith(Population, "CBas") ~ "Contemporary" #add lines as needed, depending on dataset
           )) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

View(data_PC) #check if your data frame was created correctly, with correct values & without NAs

#### import ADMIXTURE data ####
#Make sure your .Q files are in one folder
afiles <- list.files(path = "../PIRE_Gmi_Ham/ADMIXTURE/preHWE/", full.names = TRUE)
alist <- readQ(files = afiles) 
lapply(alist, attributes) #check to make sure read in properly

################################################################################################################################################

######## PCA ########

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
#then divide 1st eigenvalue by that sum to get % variance explained by PC 1, etc.
#values for varPC1 and varPC2 will be for your plot below!
varPC1 <- (eigenval[1,1] / sum(eigenval$V1))*100
varPC2 <- (eigenval[2,1] / sum(eigenval$V1))*100

#with one location
PCA_12 <- ggplot(data = data_PC, aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 16, stroke = 5) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 70.78% of total variance)", y = "PC2 (explains 12.70% of total variance)") #this is where varPC1 and varPC2 go
PCA_12_annotated <- PCA_12 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 12)), fill=guide_legend(override.aes=list(shape=21)))
PCA_12_annotated

#with multiple locations
PCA_12 <- ggplot(data = data_PC, aes(x = PC1, y = PC2, color = Era, shape = Location)) + 
  geom_point(size = 16, stroke = 5) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 70.78% of total variance)", y = "PC2 (explains 12.70% of total variance)") #this is where varPC1 and varPC2 go
PCA_12_annotated <- PCA_12 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 12)), fill=guide_legend(override.aes=list(shape=21)))
PCA_12_annotated

################################################################################################################################################

######## ADMIXTURE ########

#find out how many individuals there are for each Location & Era
table(data_PC$Location, data_PC$Era) #if you only have one site to compare, use table(data_PC$Era) only

#create group labels based on the results of your table query above
grplab <- c(rep("Hist - Basud", 12), rep("Hist - Hamilo", 29), 
            rep("Contemp - Basud", 20), rep("Contemp - Hamilo", 33)) #create group labels (# individuals in each population)

#put group labels into meta.data for plots
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)

#### ADMIXTURE cross-validation scores #####
#found in *log.out files
#lowest CV error = best structure (most likely # of demes)

#K1: 0.58307
#K2: 0.22519
#K3: 0.25069
#K4: 0.27119
#K5: 0.26734

CVs <- c(0.58307, 0.22519, 0.25069, 0.27119, 0.26734)
Ks <- c(1, 2, 3, 4, 5)

CV_df <- as.data.frame(cbind(CVs, Ks))

CV_plot <- ggplot(data = CV_df, aes(x=Ks, y = CVs)) + 
  geom_line() + 
  geom_point()
CV_plot_annotated <- CV_plot + theme_bw() + 
  labs(title = "Cross-validation error plot", y = "Cross-validation error", x = "K") + 
  theme(axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        legend.text = element_text(size = 30), legend.title = element_text(size = 30))
CV_plot_annotated

#### ADMIXTURE plots ####
#Often only K2 is needed (to catch cryptic structure) --> should at least create plots up to the K with the lowest CV score
#plots will be exported to your working directory

K2 <- plotQ(alist[2], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath=getwd(),
            clustercol = c("#FF9329", "#2121D9"),
            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K3 <- plotQ(alist[3], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath=getwd(),
            clustercol = c("#9999FF", "#2121D9", "#FF9329"),
            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K4 <- plotQ(alist[4], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath=getwd(),
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")

K5 <- plotQ(alist[5], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, exportpath=getwd(),
            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6, 
            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
            showtitle = TRUE, titlelab = "ADMIXTURE plot", showsubtitle = TRUE, subtitlelab = "preHWE SNPs")