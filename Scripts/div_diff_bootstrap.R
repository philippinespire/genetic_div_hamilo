############################## Script for Diversity Diff Permutations  #######################################################

#randomly assigns individuals to two time points (10,000X)
#calculates difference in Ho & He
#written for Aen

#################################################################################################################################################

######## Set-up ########

#set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#getwd()

remove(list = ls())

#load libraries
library(adegenet) #v.2.1.10
library(hierfstat) #v.0.5.11
library(boot) #v.1.3.28
library(vcfR) #v.1.14.0

#read in data
vcf <- read.vcfR("Data/Aen_Ham/Aen.renamed.noLD.A.nohighhet.vcf")
  genind <- vcfR2genind(vcf) #convert to genind object for analyses

################################################################################################################################################

######## Ho & He permutation estimates ########

#add population level data
#will resample this to create simulated datasets
pop <- c(rep(1, times = 30), rep(2, times = 92)) #1 = Albatross, 2 = Contemporary

#create vectors to populate
permutation_Ho <- c()
permutation_He <- c()
permutation_Fis <- c()

for (i in 1:10000) { 
  list <- sample(pop, replace = FALSE)
  pop(genind) <- list
  sum_stats <- basic.stats(genind)
  
  #calculate Ho & He diff
  permutation_Ho[i] <- diff(colMeans(sum_stats$Ho))
  permutation_He[i] <- diff(colMeans(sum_stats$Hs))
  
  #pull out Ho & He
  Hist_Ho <- sum_stats$Ho[,1]
  Contemp_Ho <- sum_stats$Ho[,2]
  Hist_He <- sum_stats$Hs[,1]
  Contemp_He <- sum_stats$Hs[,2]
  
  #bind into a dataframe
  Hist_div <- as.data.frame(cbind(Hist_Ho, Hist_He))
    colnames(Hist_div) <- c("Ho", "He")
  Contemp_div <- as.data.frame(cbind(Contemp_Ho, Contemp_He))
    colnames(Contemp_div) <- c("Ho", "He")
  
  #calculate Fis
  Hist_div$Fis <- 1 - (Hist_div$Ho/Hist_div$He)
    Hist_div$Fis[is.nan(Hist_div$Fis)] <- 0
  Contemp_div$Fis <- 1 - (Contemp_div$Ho/Contemp_div$He)
    Contemp_div$Fis[is.nan(Contemp_div$Fis)] <- 0
  
  Hist_mean_Fis <- mean(Hist_div$Fis) #0.006
  Contemp_mean_Fis <- mean(Contemp_div$Fis) #0.006
  
  #calculate Fis diff
  permutation_Fis[i] <- Contemp_mean_Fis - Hist_mean_Fis
}

permutation_diff_div <- as.data.frame(permutation_Ho)
  permutation_diff_div$He <- permutation_He
  permutation_diff_div$Fis <- permutation_Fis
  permutation_diff_div$iteration <- 1:10000
  colnames(permutation_diff_div) <- c("Ho", "He", "Fis", "iteration")

write.csv(permutation_diff_div, "Data/Aen_Ham/Aen_permutation_div.csv")