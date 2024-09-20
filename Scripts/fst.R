####################################### Script for Fst ##################################################

#Calculates pairwise Fst using the hierfstat package
#Identifies high Fst windows from pixy output

#pixy script
#from https://pixy.readthedocs.io/en/latest/plotting.html

##################################################################################################

######## Set-up ########

getwd()
remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(adegenet) #v.2.1.10
library(hierfstat) #v.0.5.11
library(vcfR) #v.1.14.0
library(boot) #v.1.3.28
library(here) #v.1.0.1

#read in VCFs
Gmi_vcf <- read.vcfR("Data/Gmi_Ham/Gmi.renamed.noLD.A.nohighhet.Ham.vcf")
  Gmi_genind <- vcfR2genind(Gmi_vcf)
Gmi_cryptic_vcf <- read.vcfR("Data/Gmi_Ham/Gmi.renamed.noLD.Ham.vcf")
  Gmi_cryptic_genind <- vcfR2genind(Gmi_cryptic_vcf)
Ela_vcf <- read.vcfR("Data/Ela_Ham/Ela.renamed.noLD.Ela.nohighhet.vcf")
  Ela_genind <- vcfR2genind(Ela_vcf)
Ela_cryptic_vcf <- read.vcfR("Data/Ela_Ham/Ela.renamed.noLD.vcf")
  Ela_cryptic_genind <- vcfR2genind(Ela_cryptic_vcf)
Aen_vcf <- read.vcfR("Data/Aen_Ham/Aen.renamed.noLD.A.nohighhet.vcf")
  Aen_genind <- vcfR2genind(Aen_vcf) 
Aen_unrelated_vcf <- read.vcfR("Data/Aen_Ham/Aen.renamed.noLD.A.nohighhet.unrelated.vcf")
  Aen_unrelated_genind <- vcfR2genind(Aen_unrelated_vcf)

#read in pixy fst windows 
Gmi_fst <- read.table(here("Data/Gmi_Ham/pixy", "pixy_fst.txt"), header = TRUE)
Ela_fst <- read.table(here("Data/Ela_Ham/pixy", "pixy_fst.txt"), header = TRUE)
Aen_fst <- read.table(here("Data/Aen_Ham/pixy", "pixy_fst.txt"), header = TRUE)
Aen_unrelated_fst <- read.table(here("Data/Aen_Ham/pixy", "pixy_unrelated_fst.txt"), 
                                header = TRUE)

#read in BayeScan
Gmi_BayeScan <- read.csv(here("Data/Gmi_Ham", "Gmi_A_Ham_nohighhet_SNPs_pr10_fst.txt"), header = TRUE, sep = " ")
  Gmi_BayeScan <- Gmi_BayeScan[, -2] #cleaning up bc didn't read columns in correctly
  Gmi_BayeScan <- Gmi_BayeScan[, 1:6]
  colnames(Gmi_BayeScan) <- c("index", "prob", "log10(PO)", "qval", "alpha", "fst")
Ela_BayeScan <- read.csv(here("Data/Ela_Ham", "Ela_nohighhet_pr10_fst.txt"), header = TRUE, sep = " ")
  Ela_BayeScan <- Ela_BayeScan[, -2]
  Ela_BayeScan <- Ela_BayeScan[, 1:6]
  colnames(Ela_BayeScan) <- c("index", "prob", "log10(PO)", "qval", "alpha", "fst")
Aen_BayeScan <- read.csv(here("Data/Aen_Ham", "Aen_A_SNPs_pr10_fst.txt"), header = TRUE, sep = " ")
  Aen_BayeScan <- Aen_BayeScan[, -2]
  Aen_BayeScan <- Aen_BayeScan[, 1:6]
  colnames(Aen_BayeScan) <- c("index", "prob", "log10(PO)", "qval", "alpha", "fst")

#########################################################################################################################

######## Gmi Fst #######

#### BayeScan output ####
  
## identify any outliers ##
#qvalue < 0.05 --> FDR 5% (SNPs with this q-value have 5% prob of being a false positive)
Gmi_pot_outliers <- subset(Gmi_BayeScan, qval <= 0.05) #3 outliers

#plot
Gmi_BS_plot <- ggplot(data = Gmi_BayeScan, aes(y = qval, x = fst)) + 
  geom_point(color = "#1c3b0e", alpha = 0.5, size = 16) + 
  geom_hline(aes(yintercept = 0.05), color = "black", linewidth = 4, linetype = "dashed") + 
  xlab(bquote(~F[sT])) + ylab("q-value") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Gmi_BS_plot

#### Calculate pairwise-Fst ####

## Calculate pairwise-Fst for filtered dataset ####
#add population level data
pop <- c(rep(1, times = 18), rep(2, times = 63))
  pop(Gmi_genind) <- pop #1 = Albatross, 2 = Contemporary
  Gmi_genind

Gmi_hierf <- genind2hierfstat(Gmi_genind) #convert to hierfstat db for pairwise analyses
Gmi_pairwise_fst <- genet.dist(Gmi_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#0.00467

## bootstrap pairwise_fst for 95% CI ##
#need to convert pop character to numeric for bootstrap to work
Gmi_hierf$pop <- as.numeric(Gmi_hierf$pop)
  class(Gmi_hierf$pop) #check to make sure numeric

#bootstrap pairwise estimates
Gmi_pairwise_boot <- boot.ppfst(dat = Gmi_hierf, nboot = 1000, 
                                quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
Gmi_ci_upper <- Gmi_pairwise_boot$ul #0.00617
Gmi_ci_lower <- Gmi_pairwise_boot$ll #0.00339

## Calculate pairwise-Fst for cryptic species dataset ####
#add population level data
pop <- c(rep(1, times = 7), 2, rep(1, times = 7), 2, 1, 2, 1, 1, 1, 2, 1, rep(2, times = 7),
         3, 3, 3, 4, 3, 3, 4, rep(3, times = 8), 4, 3, 4, rep(3, times = 6), 
         4, 3, 3, 4, 3, 4, 4, 3, 4, rep(3, times = 6), 4, 3, 3, 4, 3, 3, 4, 3, 4, 4, 3, 3, 3, 4, 3, 4, 
         rep(3, times = 6), 4, 4, 3, 3, 3, 4, 3, 4, 3, 4, rep(3, times = 5), 4, 3, 3, 3, 3, 4, 3, 3, 4, 
         4, 3, 4, 4, 3, 3, 3, 3, 4)
pop(Gmi_cryptic_genind) <- pop #1 = Albatross - A, 2 = Albatross - B, 3 = Contemporary - A, 4 = Contemporary - B
Gmi_cryptic_genind

Gmi_cryptic_hierf <- genind2hierfstat(Gmi_cryptic_genind) #convert to hierfstat db for pairwise analyses
Gmi_cryptic_pairwise_fst <- genet.dist(Gmi_cryptic_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#1-2: 0.66023, 1-3: 0.00348, 1-4: 0.72117, 2-3: 0.66341, 2-4: 0.05885, 3-4: 0.69433 

#### pixy fst windows ####

#remove Na
Gmi_fst <- Gmi_fst[!is.na(Gmi_fst$avg_wc_fst), ] #lost 0 windows

#pull out windows with fst > 0.15
Gmi_Ham_A_pixy_pot_outliers <- subset(Gmi_fst, avg_wc_fst >= 0.15) #10
  Gmi_Ham_A_pixy_pot_outliers <- Gmi_Ham_A_pixy_pot_outliers[order(Gmi_Ham_A_pixy_pot_outliers$avg_wc_fst, 
                                                                   decreasing = TRUE),]

#plot windows
Gmi_fst$NUM <- c(1:2484) #just assigning each window a number so can plot

Gmi_pixy_fst_window_plot <- ggplot(data = Gmi_fst, aes(x = NUM, y = avg_wc_fst)) + 
  geom_point(color = "#134f5c", size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0.15, color = "black", linewidth = 2, linetype = "dashed") + 
  ggtitle("Gmi pixy window outliers") + labs(y = "avg fst", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
Gmi_pixy_fst_window_plot

#############################################################################################################################

######## Ela Fst #######

#### BayeScan output ####

## identify any outliers ##
#qvalue < 0.05 --> FDR 5% (SNPs with this q-value have 5% prob of being a false positive)
Ela_pot_outliers <- subset(Ela_BayeScan, qval <= 0.05) #1 outlier

#plot
Ela_BS_plot <- ggplot(data = Ela_BayeScan, aes(y = qval, x = fst)) + 
  geom_point(color = "#16537e", alpha = 0.5, size = 16) + 
  geom_hline(aes(yintercept = 0.05), color = "black", linewidth = 4, linetype = "dashed") + 
  xlab(bquote(~F[sT])) + ylab("q-value") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ela_BS_plot

#### Calculate pairwise-Fst ####

## Calculate pairwise-Fst for filtered dataset ####
#add population level data
pop <- c(rep(1, times = 11), rep(2, times = 92))
  pop(Ela_genind) <- pop #1 = Albatross, 2 = Contemporary
  Ela_genind

Ela_hierf <- genind2hierfstat(Ela_genind) #convert to hierfstat db for pairwise analyses
Ela_pairwise_fst <- genet.dist(Ela_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#0.00653

## bootstrap pairwise_fst for 95% CI ##
#need to convert pop character to numeric for bootstrap to work
Ela_hierf$pop <- as.numeric(Ela_hierf$pop)

#bootstrap pairwise estimates
Ela_pairwise_boot <- boot.ppfst(dat = Ela_hierf, nboot = 1000, 
                                quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
Ela_ci_upper <- Ela_pairwise_boot$ul #0.00721
Ela_ci_lower <- Ela_pairwise_boot$ll #0.00583

## Calculate pairwise-Fst for cryptic species dataset ####
#add population level data
pop <- c(rep(2, times = 4), rep(1, times = 5), rep(2, times = 7), rep(1, times = 4), rep(2, times = 8), 1, 1, 
         rep(3, times = 42), 4, rep(3, times = 52))
pop(Ela_cryptic_genind) <- pop #1 = Albatross - Ela, 2 = Albatross - Lle, 3 = Contemporary - Ela, 4 = Contemporary - Ela
Ela_cryptic_genind

Ela_cryptic_hierf <- genind2hierfstat(Ela_cryptic_genind) #convert to hierfstat db for pairwise analyses
Ela_cryptic_pairwise_fst <- genet.dist(Ela_cryptic_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#1-2: 0.49223, 1-3: 0.00654, 1-4: 0.28935, 2-3: 0.43207, 2-4: 0.05352, 3-4: 0.31759


#### pixy fst windows ####

#remove Na
Ela_fst <- Ela_fst[!is.na(Ela_fst$avg_wc_fst), ] #lost 1 window

#pull out windows with fst > 0.15
Ela_pixy_pot_outliers <- subset(Ela_fst, avg_wc_fst >= 0.15) #241
Ela_pixy_pot_outliers <- Ela_pixy_pot_outliers[order(Ela_pixy_pot_outliers$avg_wc_fst, 
                                                                 decreasing = TRUE),]

#plot windows
Ela_fst$NUM <- c(1:7972) #just assigning each window a number so can plot

Ela_pixy_fst_window_plot <- ggplot(data = Ela_fst, aes(x = NUM, y = avg_wc_fst)) + 
  geom_point(color = "#134f5c", size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0.15, color = "black", linewidth = 2, linetype = "dashed") + 
  ggtitle("Ela pixy window outliers") + labs(y = "avg fst", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
Ela_pixy_fst_window_plot

#############################################################################################################################

######## Aen Fst #######

#### BayeScan output ####

## identify any outliers ##
#qvalue < 0.05 --> FDR 5% (SNPs with this q-value have 5% prob of being a false positive)
Aen_pot_outliers <- subset(Aen_BayeScan, qval <= 0.05) #0 outliers

#plot
Aen_BS_plot <- ggplot(data = Aen_BayeScan, aes(y = qval, x = fst)) + 
  geom_point(color = "#a25505", alpha = 0.5, size = 16) + 
  geom_hline(aes(yintercept = 0.05), color = "black", linewidth = 4, linetype = "dashed") + 
  xlab(bquote(~F[sT])) + ylab("q-value") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_BS_plot

#### Calculate pairwise-Fst ####

#add population level data
pop <- c(rep(1, times = 30), rep(2, times = 92))
  pop(Aen_genind) <- pop #1 = Albatross, 2 = Contemporary
  Aen_genind

Aen_hierf <- genind2hierfstat(Aen_genind) #convert to hierfstat db for pairwise analyses
Aen_pairwise_fst <- genet.dist(Aen_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#0.05743

## bootstrap pairwise_fst for 95% CI ##

#need to convert pop character to numeric for bootstrap to work
Aen_hierf$pop <- as.numeric(Aen_hierf$pop)

#bootstrap pairwise estimates
Aen_pairwise_boot <- boot.ppfst(dat = Aen_hierf, nboot = 1000, 
                                quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
Aen_ci_upper <- Aen_pairwise_boot$ul #0.06508
Aen_ci_lower <- Aen_pairwise_boot$ll #0.05020

#### pixy fst windows ####

#remove Na
Aen_fst <- Aen_fst[!is.na(Aen_fst$avg_wc_fst), ] #lost 6 windows

#pull out windows with fst > 0.15
Aen_pixy_pot_outliers <- subset(Aen_fst, avg_wc_fst >= 0.15) #260
  Aen_pixy_pot_outliers <- Aen_pixy_pot_outliers[order(Aen_pixy_pot_outliers$avg_wc_fst, 
                                                     decreasing = TRUE),]

#plot windows
Aen_fst$NUM <- c(1:1065) #just assigning each window a number so can plot

Aen_pixy_fst_window_plot <- ggplot(data = Aen_fst, aes(x = NUM, y = avg_wc_fst)) + 
  geom_point(color = "#a25505", size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0.15, color = "black", linewidth = 2, linetype = "dashed") + 
  ggtitle("Aen pixy window outliers") + labs(y = "avg fst", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
Aen_pixy_fst_window_plot

#############################################################################################################################

######## Aen unrelated Fst #######

#### BayeScan output ####

## identify any outliers ##
#qvalue < 0.05 --> FDR 5% (SNPs with this q-value have 5% prob of being a false positive)
Aen_pot_outliers <- subset(Aen_BayeScan, qval <= 0.05) #0 outliers

#plot
Aen_BS_plot <- ggplot(data = Aen_BayeScan, aes(y = qval, x = fst)) + 
  geom_point(color = "#a25505", alpha = 0.5, size = 16) + 
  geom_hline(aes(yintercept = 0.05), color = "black", linewidth = 4, linetype = "dashed") + 
  xlab(bquote(~F[sT])) + ylab("q-value") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_BS_plot

#### Calculate pairwise-Fst ####

#add population level data
pop <- c(rep(1, times = 30), rep(2, times = 84))
  pop(Aen_unrelated_genind) <- pop #1 = Albatross, 2 = Contemporary
  Aen_unrelated_genind

Aen_unrelated_hierf <- genind2hierfstat(Aen_unrelated_genind) #convert to hierfstat db for pairwise analyses
Aen_unrelated_pairwise_fst <- genet.dist(Aen_unrelated_hierf, method = "WC84") #calculates Weir & Cockerham's Fst
#0.05669

## bootstrap pairwise_fst for 95% CI ##

#need to convert pop character to numeric for bootstrap to work
Aen_unrelated_hierf$pop <- as.numeric(Aen_unrelated_hierf$pop)

#bootstrap pairwise estimates
Aen_unrelated_pairwise_boot <- boot.ppfst(dat = Aen_unrelated_hierf, nboot = 1000, 
                                          quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
Aen_unrelated_ci_upper <- Aen_unrelated_pairwise_boot$ul #0.06532
Aen_unrelated_ci_lower <- Aen_unrelated_pairwise_boot$ll #0.04896

#### pixy fst windows ####

#remove Na
Aen_unrelated_fst <- Aen_unrelated_fst[!is.na(Aen_unrelated_fst$avg_wc_fst), ] #lost 6 windows

#pull out windows with fst > 0.15
Aen_unrelated_pixy_pot_outliers <- subset(Aen_unrelated_fst, avg_wc_fst >= 0.15) #260
Aen_unrelated_pixy_pot_outliers <- Aen_unrelated_pixy_pot_outliers[order(Aen_unrelated_pixy_pot_outliers$avg_wc_fst, 
                                                                         decreasing = TRUE),]

#plot windows
Aen_unrelated_fst$NUM <- c(1:1062) #just assigning each window a number so can plot

Aen_pixy_fst_window_plot <- ggplot(data = Aen_unrelated_fst, aes(x = NUM, y = avg_wc_fst)) + 
  geom_point(color = "#a25505", size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0.15, color = "black", linewidth = 2, linetype = "dashed") + 
  ggtitle("Aen unrelated pixy window outliers") + labs(y = "avg fst", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
Aen_pixy_fst_window_plot
