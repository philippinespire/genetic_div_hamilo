################################## Script for Genetic Diversity Estimates  #######################################################

#Calculates Ho, He & Fis
#Bootstrapping across sites
#Also reads in permuted data to compare point estimates with

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(adegenet) #v.2.1.10
library(hierfstat) #v.0.5.11
library(tidyverse) #v.2.0.0
library(boot) #v.1.3.28
library(vcfR) #v.1.14.0
library(here) #v.1.0.1

#write function to calculate mean
samp_mean <- function(x, i) {
  mean(x[i])
}

#read in vcfs 
Gmi_vcf <- read.vcfR("Data/Gmi_Ham/Gmi.renamed.noLD.A.nohighhet.Ham.vcf")
  Gmi_genind <- vcfR2genind(Gmi_vcf)
Ela_vcf <- read.vcfR("Data/Ela_Ham/Ela.renamed.noLD.Ela.nohighhet.vcf")
  Ela_genind <- vcfR2genind(Ela_vcf) #convert to genind object for analyses
  
#read in permutation data
Gmi_Ham_permutation <- read.csv(here("Data/Gmi_Ham", "GmiHam_permutation_div.csv"))
Ela_permutation <- read.csv(here("Data/Ela_Ham", "Ela_permutation_div.csv")) #using w highhet for now, since w/o still running

################################################################################################################################################
  
######## Gmi diversity estimates ########
  
#add population level data
pop <- c(rep(1, times = 18), rep(2, times = 63))
  pop(Gmi_genind) <- pop #1 = Albatross - Hamilo, 2 = Contemp = Hamilo
  Gmi_genind #check to make sure pops read right (should range from 18-63 individuals)
  
#calculate diversity metrics w/in pops
sum_stats_Gmi <- basic.stats(Gmi_genind)
  
#### Calculate mean Ho, He ####
  
#mean Ho & He for each population
Gmi_Ho_means <- colMeans(sum_stats_Gmi$Ho) #1: 0.10796, 2: 0.10165
Gmi_Hs_means <- colMeans(sum_stats_Gmi$Hs) #1: 0.10218, 2: 0.09606

#Ho % change: 5.84% decrease
#He % change: 5.99% decrease
  
Gmi_Ho_means_diff <- mean(sum_stats_Gmi$Ho[,2]) - mean(sum_stats_Gmi$Ho[,1])  #-0.00631
Gmi_Hs_means_diff <- mean(sum_stats_Gmi$Hs[,2]) - mean(sum_stats_Gmi$Hs[,1])  #-0.00612
  
##### Pull out per-marker estimates ####
  
#pull out Ho for each population
Hist_Gmi_Ho <- sum_stats_Gmi$Ho[,1]
Contemp_Gmi_Ho <- sum_stats_Gmi$Ho[,2]
  
#pull out He for each population
Hist_Gmi_He <- sum_stats_Gmi$Hs[,1]
Contemp_Gmi_He <- sum_stats_Gmi$Hs[,2]
  
#bind Ho & He into dataframe for each pop
Hist_Gmi_div <- as.data.frame(cbind(Hist_Gmi_Ho, Hist_Gmi_He))
  colnames(Hist_Gmi_div) <- c("Ho", "He")

Contemp_Gmi_div <- as.data.frame(cbind(Contemp_Gmi_Ho, Contemp_Gmi_He))
  colnames(Contemp_Gmi_div) <- c("Ho", "He")
  
##### Calculate mean Fis ####
  
#calculate Fis for each population
Hist_Gmi_div$Fis <- 1 - (Hist_Gmi_div$Ho/Hist_Gmi_div$He)
  Hist_Gmi_div$Fis[is.nan(Hist_Gmi_div$Fis)] <- 0

Contemp_Gmi_div$Fis <- 1 - (Contemp_Gmi_div$Ho/Contemp_Gmi_div$He)
  Contemp_Gmi_div$Fis[is.nan(Contemp_Gmi_div$Fis)] <- 0
  
#calculate mean Fis
Hist_Gmi_mean_Fis <- mean(Hist_Gmi_div$Fis) #-0.01125
Contemp_Gmi_mean_Fis <- mean(Contemp_Gmi_div$Fis) #-0.01127
  
Contemp_Gmi_mean_Fis - Hist_Gmi_mean_Fis #-0.00002
  
#write out
#write.csv(Hist_Gmi_div, "Data/Gmi_Ham/Hist_Gmi_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Contemp_Gmi_div, "Data/Gmi_Ham/Contemp_Gmi_div.csv", quote = FALSE, row.names = TRUE)
  
#### Ho & He bootstrapping ####
  
#read in if running separately
#Hist_Gmi_div <- read.csv("Data/Gmi_Ham/Hist_Gmi_div.csv")
#  colnames(Hist_Gmi_div) <- c("locus", "Ho", "He", "Fis")
#Contemp_Gmi_div <- read.csv("Data/Gmi_Ham/Contemp_Gmi_div.csv")
#  colnames(Contemp_Gmi_div) <- c("locus", "Ho", "He", "Fis")
  
## Ho bootstrapping ##
#historical Ho
boot_Hist_Gmi_Ho <- boot(data = Hist_Gmi_div$Ho, statistic = samp_mean, R = 1000)
  Hist_Gmi_Ho_95ci <- boot.ci(boot_Hist_Gmi_Ho, conf = 0.95, type = "norm")
  Hist_Gmi_Ho_95ci_normal <- Hist_Gmi_Ho_95ci$normal

#contemporary Ho
boot_contemp_Gmi_Ho <- boot(data = Contemp_Gmi_div$Ho, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Ho_95ci <- boot.ci(boot_contemp_Gmi_Ho, conf = 0.95, type = "norm")
  Contemp_Gmi_Ho_95ci_normal <- Contemp_Gmi_Ho_95ci$normal 
 
## He bootstrapping ## 
#historical He
boot_Hist_Gmi_He <- boot(data = Hist_Gmi_div$He, statistic = samp_mean, R = 1000)
  Hist_Gmi_He_95ci <- boot.ci(boot_Hist_Gmi_He, conf = 0.95, type = "norm")
  Hist_Gmi_He_95ci_normal <- Hist_Gmi_He_95ci$normal
  
#contemporary He
boot_contemp_Gmi_He <- boot(data = Contemp_Gmi_div$He, statistic = samp_mean, R = 1000)
  Contemp_Gmi_He_95ci <- boot.ci(boot_contemp_Gmi_He, conf = 0.95, type = "norm")
  Contemp_Gmi_He_95ci_normal <- Contemp_Gmi_He_95ci$normal 
  
## Fis bootstrapping ##
#historical Fis
boot_Hist_Gmi_Fis <- boot(data = Hist_Gmi_div$Fis, statistic = samp_mean, R = 1000)
  Hist_Gmi_Fis_95ci <- boot.ci(boot_Hist_Gmi_Fis, conf = 0.95, type = "norm")
  Hist_Gmi_Fis_95ci_normal <- Hist_Gmi_Fis_95ci$normal
  
#contemporary Fis
boot_contemp_Gmi_Fis <- boot(data = Contemp_Gmi_div$Fis, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Fis_95ci <- boot.ci(boot_contemp_Gmi_Fis, conf = 0.95, type = "norm")
  Contemp_Gmi_Fis_95ci_normal <- Contemp_Gmi_Fis_95ci$normal 
  
## write out Ho, He, Fis summary tables ##
Ho_mean <- as.data.frame(c(0.10796, 0.10165))
He_mean <- as.data.frame(c(0.10218, 0.09606))
Fis_mean <- as.data.frame(c(-0.01125, -0.01127))
  
#Ho summary table
Ho_ci <- rbind(Hist_Gmi_Ho_95ci_normal, Contemp_Gmi_Ho_95ci_normal) #combine df w/ci info for each pop into one dataframe
Ho_sum <- cbind(Ho_mean, Ho_ci)
  colnames(Ho_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Ho_sum$Pop <- c("Hist", "Contemp")
  
Ho_sum$diff_lower <- Ho_sum$mean - Ho_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
Ho_sum$diff_upper <- Ho_sum$`97.5_per` - Ho_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization
  
#He summary table
He_ci <- rbind(Hist_Gmi_He_95ci_normal, Contemp_Gmi_He_95ci_normal)
He_sum <- cbind(He_mean, He_ci)
  colnames(He_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  He_sum$Pop <- c("Hist", "Contemp")
  
He_sum$diff_lower <- He_sum$mean - He_sum$`2.5_per`
He_sum$diff_upper <- He_sum$`97.5_per` - He_sum$mean
  
#Fis summary table
Fis_ci <- rbind(Hist_Gmi_Fis_95ci_normal, Contemp_Gmi_Fis_95ci_normal)
Fis_sum <- cbind(Fis_mean, Fis_ci)
  colnames(Fis_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Fis_sum$Pop <- c("Hist", "Contemp")
  
Fis_sum$diff_lower <- Fis_sum$mean - Fis_sum$`2.5_per`
Fis_sum$diff_upper <- Fis_sum$`97.5_per` - Fis_sum$mean
  
#combine dataframes
div_sum_all <- as.data.frame(rbind(Ho_sum, He_sum, Fis_sum))
div_sum_all$metric <- c("Ho", "Ho", "He", "He", "Fis", "Fis")
  
#write out data
#write.csv(div_sum_all, "Data/Gmi_Ham/Gmi_div_cis.csv")

## Calculate 95% for diff through time ##
#Ho diff
Gmi_Ho_boot_diff <- as.data.frame(boot_Hist_Gmi_Ho$t)
  Gmi_Ho_boot_diff$Contemp_Ho <- boot_contemp_Gmi_Ho$t
  colnames(Gmi_Ho_boot_diff) <- c("Hist_Ho", "Contemp_Ho")

Gmi_Ho_boot_diff$Ho_diff <- Gmi_Ho_boot_diff$Contemp_Ho - Gmi_Ho_boot_diff$Hist_Ho
  Gmi_Ho_boot_diff_CIs <- quantile(Gmi_Ho_boot_diff$Ho_diff, c(0.025, 0.975)) #-0.01118, -0.00122

#He diff
Gmi_He_boot_diff <- as.data.frame(boot_Hist_Gmi_He$t)
  Gmi_He_boot_diff$Contemp_He <- boot_contemp_Gmi_He$t
  colnames(Gmi_He_boot_diff) <- c("Hist_He", "Contemp_He")
  
Gmi_He_boot_diff$He_diff <- Gmi_He_boot_diff$Contemp_He - Gmi_He_boot_diff$Hist_He
  Gmi_He_boot_diff_CIs <- quantile(Gmi_He_boot_diff$He_diff, c(0.025, 0.975)) #-0.01033, -0.00184

#Fis diff
Gmi_Fis_boot_diff <- as.data.frame(boot_Hist_Gmi_Fis$t)
  Gmi_Fis_boot_diff$Contemp_Fis <- boot_contemp_Gmi_Fis$t
  colnames(Gmi_Fis_boot_diff) <- c("Hist_Fis", "Contemp_Fis")
  
Gmi_Fis_boot_diff$Fis_diff <- Gmi_Fis_boot_diff$Contemp_Fis - Gmi_Fis_boot_diff$Hist_Fis
  Gmi_Fis_boot_diff_CIs <- quantile(Gmi_Fis_boot_diff$Fis_diff, c(0.025, 0.975)) #-0.00505, 0.00496

## Calculate empirical p-value ##
#Ho
Ho_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, Ho < -0.00631))) #number of permutations with value greater than observed (really less than)
  (Ho_perm_greater + 1)/10001 #0.0914

#He
He_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, He < -0.00612)))
  (He_perm_greater + 1)/10001 #0.0455

#Fis
Fis_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, Fis < -0.00002)))
  (Fis_perm_greater + 1)/10001 #0.5796
  
#### Gmi Ho, He & Fis visualization ####
  
#read in data (if running separately)
#div_sum_all <- read.csv("Data/Gmi_Ham/Gmi_div_cis.csv", row.names = 1)
  
#add factor level
div_sum_all$Pop <- factor(div_sum_all$Pop, levels = c("Hist", "Contemp"))
  
#het only dataframe
het_sum_all <- subset(div_sum_all, div_sum_all$metric != "Fis")
  
#fis dataframe
Fis_sum <- subset(div_sum_all, div_sum_all$metric == "Fis")
  
#plot of mean Ho & He w/95% CI error bars
Het_plot <- ggplot(data = het_sum_all, aes(x = metric, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#1c3b0e","#afc8a4")) +
  scale_shape_manual(values = c(19, 15)) +
  ylim(0.05, 0.15) + ggtitle("Gmi Ho & He 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
Het_plot
  
#plot of mean Fis w/95% CI error bars
Fis_plot <- ggplot(data = Fis_sum, aes(x = metric, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 2, color = "black") + 
  scale_color_manual(values = c("#1c3b0e","#afc8a4")) +
  scale_shape_manual(values = c(19, 15)) +
  ggtitle("Gmi Fis 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
Fis_plot

#### Null distribution plots ####

#plot of null Ho diff distribution
Ho_null_plot <- ggplot(Gmi_Ham_permutation, aes(x = Ho)) + 
  geom_density(color = "#bac4b6", fill = "#bac4b6") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00887), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00973), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00631), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~H[o])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
Ho_null_plot

#plot of null He diff distribution
He_null_plot <- ggplot(Gmi_Ham_permutation, aes(x = He)) + 
  geom_density(color = "#bac4b6", fill = "#bac4b6") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00704), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00761), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00612), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~H[e])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
He_null_plot

#plot of null Fis diff distribution
Fis_null_plot <- ggplot(Gmi_Ham_permutation, aes(x = Fis)) + 
  geom_density(color = "#bac4b6", fill = "#bac4b6") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.01296), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00927), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00002), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~F[IS])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
Fis_null_plot
  
################################################################################################################################################
  
######## Ela diversity estimates ########

#add population level data
pop <- c(rep(1, times = 11), rep(2, times = 92))
  pop(Ela_genind) <- pop #1 = Albatross, 2 = Contemporary
  Ela_genind #check to make sure pops read right (should range from 11-92 individuals)

#calculate diversity metrics w/in pops
sum_stats_Ela <- basic.stats(Ela_genind)
 
#### Calculate mean Ho, He ####
  
#mean Ho & He for each population
Ela_Ho_means <- colMeans(sum_stats_Ela$Ho) #1: 0.11900, 2: 0.11577
Ela_Hs_means <- colMeans(sum_stats_Ela$Hs) #1: 0.12109, 2: 0.11408

#Ho % change: 2.74% decrease
#He % change: 5.79% decrease

Ela_Ho_means_diff <- diff(colMeans(sum_stats_Ela$Ho)) #-0.00324
Ela_Hs_means_diff <- diff(colMeans(sum_stats_Ela$Hs)) #-0.00701

##### Pull out per-marker estimates ####

#pull out Ho for each population
Hist_Ela_Ho <- sum_stats_Ela$Ho[,1]
Contemp_Ela_Ho <- sum_stats_Ela$Ho[,2]

#pull out He for each population
Hist_Ela_He <- sum_stats_Ela$Hs[,1]
Contemp_Ela_He <- sum_stats_Ela$Hs[,2]

#bind Ho & He into dataframe for each pop
Hist_Ela_div <- as.data.frame(cbind(Hist_Ela_Ho, Hist_Ela_He))
  head(Hist_Ela_div) #check content is right
  colnames(Hist_Ela_div) <- c("Ho", "He")

Contemp_Ela_div <- as.data.frame(cbind(Contemp_Ela_Ho, Contemp_Ela_He))
  colnames(Contemp_Ela_div) <- c("Ho", "He")

##### Calculate mean Fis ####
  
#calculate Fis for each population
Hist_Ela_div$Fis <- 1 - (Hist_Ela_div$Ho/Hist_Ela_div$He)
  Hist_Ela_div$Fis[is.nan(Hist_Ela_div$Fis)] <- 0
  
Contemp_Ela_div$Fis <- 1 - (Contemp_Ela_div$Ho/Contemp_Ela_div$He)
  Contemp_Ela_div$Fis[is.nan(Contemp_Ela_div$Fis)] <- 0

#calculate mean Fis
Hist_Ela_mean_Fis <- mean(Hist_Ela_div$Fis) #0.00846
Contemp_Ela_mean_Fis <- mean(Contemp_Ela_div$Fis) #0.00972

Contemp_Ela_mean_Fis - Hist_Ela_mean_Fis #0.00126

#write out
#write.csv(Hist_Ela_div, "Data/Ela_Ham/Hist_Ela_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Contemp_Ela_div, "Data/Ela_Ham/Contemp_Ela_div.csv", quote = FALSE, row.names = TRUE)

#### Ho & He bootstrapping ####

#read in if running separately
#Hist_Ela_div <- read.csv("Data/Ela_Ham/Hist_Ela_div.csv")
#  colnames(Hist_Ela_div) <- c("locus", "Ho", "He", "Fis")
#Contemp_Ela_div <- read.csv("Data/Ela_Ham/Contemp_Ela_div.csv")
#  colnames(Contemp_Ela_div) <- c("locus", "Ho", "He", "Fis")
  
## Ho bootstrapping ##
#historical Ho
boot_Hist_Ela_Ho <- boot(data = Hist_Ela_div$Ho, statistic = samp_mean, R = 1000) #1000 permutations of Ho
  Hist_Ela_Ho_95ci <- boot.ci(boot_Hist_Ela_Ho, conf = 0.95, type = "norm") #get 95% CI for Ho
  Hist_Ela_Ho_95ci_normal <- Hist_Ela_Ho_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for Ho

#contemporary Ho
boot_Contemp_Ela_Ho <- boot(data = Contemp_Ela_div$Ho, statistic = samp_mean, R = 1000)
  Contemp_Ela_Ho_95ci <- boot.ci(boot_Contemp_Ela_Ho, conf = 0.95, type = "norm")
  Contemp_Ela_Ho_95ci_normal <- Contemp_Ela_Ho_95ci$normal
  
## He bootstrapping ##
#historical He
boot_Hist_Ela_He <- boot(data = Hist_Ela_div$He, statistic = samp_mean, R = 1000)
  Hist_Ela_He_95ci <- boot.ci(boot_Hist_Ela_He, conf = 0.95, type = "norm")
  Hist_Ela_He_95ci_normal <- Hist_Ela_He_95ci$normal
  
#contemporary He
boot_Contemp_Ela_He <- boot(data = Contemp_Ela_div$He, statistic = samp_mean, R = 1000)
  Contemp_Ela_He_95ci <- boot.ci(boot_Contemp_Ela_He, conf = 0.95, type = "norm")
  Contemp_Ela_He_95ci_normal <- Contemp_Ela_He_95ci$normal

## Fis bootstrapping ##
#historical Fis
boot_Hist_Ela_Fis <- boot(data = Hist_Ela_div$Fis, statistic = samp_mean, R = 1000)
  Hist_Ela_Fis_95ci <- boot.ci(boot_Hist_Ela_Fis, conf = 0.95, type = "norm")
  Hist_Ela_Fis_95ci_normal <- Hist_Ela_Fis_95ci$normal
  
#contemporary Fis
boot_Contemp_Ela_Fis <- boot(data = Contemp_Ela_div$Fis, statistic = samp_mean, R = 1000)
  Contemp_Ela_Fis_95ci <- boot.ci(boot_Contemp_Ela_Fis, conf = 0.95, type = "norm")
  Contemp_Ela_Fis_95ci_normal <- Contemp_Ela_Fis_95ci$normal

## write out Ho, He, Fis summary tables ##
Ho_mean <- as.data.frame(c(0.11900, 0.11577))
He_mean <- as.data.frame(c(0.12109, 0.11408))
Fis_mean <- as.data.frame(c(0.00846, 0.00972))

#Ho summary table
Ho_ci <- rbind(Hist_Ela_Ho_95ci_normal, Contemp_Ela_Ho_95ci_normal) #combine df w/ci info for each pop into one dataframe
Ho_sum <- cbind(Ho_mean, Ho_ci)
  colnames(Ho_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Ho_sum$Pop <- c("Hist", "Contemp")

Ho_sum$diff_lower <- Ho_sum$mean - Ho_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
Ho_sum$diff_upper <- Ho_sum$`97.5_per` - Ho_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#He summary table
He_ci <- rbind(Hist_Ela_He_95ci_normal, Contemp_Ela_He_95ci_normal)
He_sum <- cbind(He_mean, He_ci)
  colnames(He_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  He_sum$Pop <- c("Hist", "Contemp")

He_sum$diff_lower <- He_sum$mean - He_sum$`2.5_per`
He_sum$diff_upper <- He_sum$`97.5_per` - He_sum$mean

#Fis summary table
Fis_ci <- rbind(Hist_Ela_Fis_95ci_normal, Contemp_Ela_Fis_95ci_normal)
Fis_sum <- cbind(Fis_mean, Fis_ci)
  colnames(Fis_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Fis_sum$Pop <- c("Hist", "Contemp")

Fis_sum$diff_lower <- Fis_sum$mean - Fis_sum$`2.5_per`
Fis_sum$diff_upper <- Fis_sum$`97.5_per` - Fis_sum$mean

#combine dataframes
div_sum_all <- as.data.frame(rbind(Ho_sum, He_sum, Fis_sum))
  div_sum_all$metric <- c("Ho", "Ho", "He", "He", "Fis", "Fis")

#write out data
#write.csv(div_sum_all, "Data/Ela_Ham/Ela_div_cis.csv")

## Calculate 95% for diff through time ##
#Ho diff
Ela_Ho_boot_diff <- as.data.frame(boot_Hist_Ela_Ho$t)
  Ela_Ho_boot_diff$Contemp_Ho <- boot_Contemp_Ela_Ho$t
  colnames(Ela_Ho_boot_diff) <- c("Hist_Ho", "Contemp_Ho")
  
Ela_Ho_boot_diff$Ho_diff <- Ela_Ho_boot_diff$Contemp_Ho - Ela_Ho_boot_diff$Hist_Ho
  Ela_Ho_boot_diff_CIs <- quantile(Ela_Ho_boot_diff$Ho_diff, c(0.025, 0.975)) #-0.00529, -0.00119
  
#He diff
Ela_He_boot_diff <- as.data.frame(boot_Hist_Ela_He$t)
  Ela_He_boot_diff$Contemp_He <- boot_Contemp_Ela_He$t
  colnames(Ela_He_boot_diff) <- c("Hist_He", "Contemp_He")
  
Ela_He_boot_diff$He_diff <- Ela_He_boot_diff$Contemp_He - Ela_He_boot_diff$Hist_He
  Ela_He_boot_diff_CIs <- quantile(Ela_He_boot_diff$He_diff, c(0.025, 0.975)) #-0.00886, -0.00517
  
#Fis diff
Ela_Fis_boot_diff <- as.data.frame(boot_Hist_Ela_Fis$t)
  Ela_Fis_boot_diff$Contemp_Fis <- boot_Contemp_Ela_Fis$t
  colnames(Ela_Fis_boot_diff) <- c("Hist_Fis", "Contemp_Fis")
  
Ela_Fis_boot_diff$Fis_diff <- Ela_Fis_boot_diff$Contemp_Fis - Ela_Fis_boot_diff$Hist_Fis
  Ela_Fis_boot_diff_CIs <- quantile(Ela_Fis_boot_diff$Fis_diff, c(0.025, 0.975)) #-0.00090, 0.00346
  
## Calculate empirical p-value ##
#Ho
Ho_perm_greater <- as.numeric(nrow(subset(Ela_permutation, Ho < -0.00324))) #number of permutations with value greater than observed (really less than)
  (Ho_perm_greater + 1)/10001 #0.0024
  
#He
He_perm_greater <- as.numeric(nrow(subset(Ela_permutation, He < -0.00701)))
  (He_perm_greater + 1)/10001 #0.0001
  
#Fis
Fis_perm_greater <- as.numeric(nrow(subset(Ela_permutation, Fis > 0.00126)))
  (Fis_perm_greater + 1)/10001 #0.1134
  
#### Ela Ho, He & Fis visualization ####

#read in data (if running separately)
#div_sum_all <- read.csv("Data/Ela_Ham/Ela_div_cis.csv", row.names = 1)

#add factor level
div_sum_all$Pop <- factor(div_sum_all$Pop, levels = c("Hist", "Contemp"))
  
#het only dataframe
het_sum_all <- subset(div_sum_all, div_sum_all$metric != "Fis")

#fis dataframe
Fis_sum <- subset(div_sum_all, div_sum_all$metric == "Fis")

#plot of mean Ho & He w/95% CI error bars
Het_plot <- ggplot(data = het_sum_all, aes(x = metric, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  scale_shape_manual(values = c(19, 15)) +
  ylim(0.05, 0.15) + ggtitle("Ela Ho & He 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
Het_plot

#plot of mean Fis w/95% CI error bars
Fis_plot <- ggplot(data = Fis_sum, aes(x = metric, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  scale_shape_manual(values = c(19, 15)) +
  ylim(0, 0.015) + ggtitle("Ela Fis 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
Fis_plot

#### Null distribution plots ####
Ela_permutation_nona <- na.omit(Ela_permutation) #not sure why this happened?

#plot of null Ho diff distribution
Ho_null_plot <- ggplot(Ela_permutation, aes(x = Ho)) + 
  geom_density(color = "#b9cbd8", fill = "#b9cbd8") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00229), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00229), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00324), linewidth = 6, color = "#16537e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~H[o])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
Ho_null_plot

#plot of null He diff distribution
He_null_plot <- ggplot(Ela_permutation, aes(x = He)) + 
  geom_density(color = "#b9cbd8", fill = "#b9cbd8") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00186), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00195), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00701), linewidth = 6, color = "#16537e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~H[e])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
He_null_plot

#plot of null Fis diff distribution
Fis_null_plot <- ggplot(Ela_permutation, aes(x = Fis)) + 
  geom_density(color = "#b9cbd8", fill = "#b9cbd8") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.01814), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.01494), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = 0.00126), linewidth = 6, color = "#16537e") + #"real" mean diff
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab(bquote("Temporal Change in"~F[IS])) + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
Fis_null_plot