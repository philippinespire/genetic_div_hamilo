################################## Script for Genetic Diversity Estimates  #######################################################

#Calculates Ho, He & Fis
#Bootstrapping across sites
#Also reads in permuted across time points data to compare point estimates with

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
Gmi_vcf <- read.vcfR("Data/Gmi_Ham/Gmi.rename.noLD.Apop.nohighhet.vcf")
  Gmi_genind <- vcfR2genind(Gmi_vcf)
Ela_vcf <- read.vcfR("Data/Ela_Ham/Ela.rename.sorted.noLD.Ela.nohighhet.vcf")
  Ela_genind <- vcfR2genind(Ela_vcf) #convert to genind object for analyses
Aen_vcf <- read.vcfR("Data/Aen_Ham/Aen.rename.noLD.A.vcf")
  Aen_genind <- vcfR2genind(Aen_vcf) #convert to genind object for analyses
  
#read in permutation data
Gmi_Ham_permutation <- read.csv(here("Data/Gmi_Ham", "GmiHam_permutation_div.csv"))
Ela_permutation <- read.csv(here("Data/Ela_Ham", "Ela_permutation_div.csv")) #using w highhet for now, since w/o still running
Aen_permutation <- read.csv(here("Data/Aen_Ham", "Aen_permutation_div.csv"))

################################################################################################################################################
  
######## Gmi diversity estimates ########
  
#add population level data
pop <- c(rep(1, times = 12), rep(2, times = 18), 
         rep(3, times = 78), rep(4, times = 61))
  pop(Gmi_genind) <- pop #1 = Albatross - Basud, 2 = Albatross - Hamilo, 3 = Contemp - Basud, 4 = Contemp = Hamilo
  Gmi_genind #check to make sure pops read right (should range from 12-78 individuals)
  
#calculate diversity metrics w/in pops
sum_stats_Gmi <- basic.stats(Gmi_genind)
  
#### Calculate mean Ho, He ####
  
#mean Ho & He for each population
Gmi_Ho_means <- colMeans(sum_stats_Gmi$Ho) #1: 0.10647, 2: 0.10784, 3: 0.08964, 4: 0.09901 
Gmi_Hs_means <- colMeans(sum_stats_Gmi$Hs) #1: 0.10545, 2: 0.10195, 3: 0.08800, 4: 0.09350
  
Gmi_Bas_Ho_means_diff <- mean(sum_stats_Gmi$Ho[,3]) - mean(sum_stats_Gmi$Ho[,1])  #-0.01683
Gmi_Ham_Ho_means_diff <- mean(sum_stats_Gmi$Ho[,4]) - mean(sum_stats_Gmi$Ho[,2])  #-0.00883

Gmi_Bas_Hs_means_diff <- mean(sum_stats_Gmi$Hs[,3]) - mean(sum_stats_Gmi$Hs[,1])  #-0.01745
Gmi_Ham_Hs_means_diff <- mean(sum_stats_Gmi$Hs[,4]) - mean(sum_stats_Gmi$Hs[,2])  #-0.00826
  
##### Pull out per-marker estimates ####
  
#pull out Ho for each population
Hist_Gmi_Bas_Ho <- sum_stats_Gmi$Ho[,1]
Hist_Gmi_Ham_Ho <- sum_stats_Gmi$Ho[,2]
Contemp_Gmi_Bas_Ho <- sum_stats_Gmi$Ho[,3]
Contemp_Gmi_Ham_Ho <- sum_stats_Gmi$Ho[,4]
  
#pull out He for each population
Hist_Gmi_Bas_He <- sum_stats_Gmi$Hs[,1]
Hist_Gmi_Ham_He <- sum_stats_Gmi$Hs[,2]
Contemp_Gmi_Bas_He <- sum_stats_Gmi$Hs[,3]
Contemp_Gmi_Ham_He <- sum_stats_Gmi$Hs[,4]
  
#bind Ho & He into dataframe for each pop
Hist_Gmi_Bas_div <- as.data.frame(cbind(Hist_Gmi_Bas_Ho, Hist_Gmi_Bas_He))
  head(Hist_Gmi_Bas_div) #check content is right
  colnames(Hist_Gmi_Bas_div) <- c("Ho", "He")
  
Hist_Gmi_Ham_div <- as.data.frame(cbind(Hist_Gmi_Ham_Ho, Hist_Gmi_Ham_He))
  colnames(Hist_Gmi_Ham_div) <- c("Ho", "He")
  
Contemp_Gmi_Bas_div <- as.data.frame(cbind(Contemp_Gmi_Bas_Ho, Contemp_Gmi_Bas_He))
  colnames(Contemp_Gmi_Bas_div) <- c("Ho", "He")
  
Contemp_Gmi_Ham_div <- as.data.frame(cbind(Contemp_Gmi_Ham_Ho, Contemp_Gmi_Ham_He))
  colnames(Contemp_Gmi_Ham_div) <- c("Ho", "He")
  
##### Calculate mean Fis ####
  
#calculate Fis for each population
Hist_Gmi_Bas_div$Fis <- 1 - (Hist_Gmi_Bas_div$Ho/Hist_Gmi_Bas_div$He)
  Hist_Gmi_Bas_div$Fis[is.nan(Hist_Gmi_Bas_div$Fis)] <- 0
  
Hist_Gmi_Ham_div$Fis <- 1 - (Hist_Gmi_Ham_div$Ho/Hist_Gmi_Ham_div$He)
  Hist_Gmi_Ham_div$Fis[is.nan(Hist_Gmi_Ham_div$Fis)] <- 0
  
Contemp_Gmi_Bas_div$Fis <- 1 - (Contemp_Gmi_Bas_div$Ho/Contemp_Gmi_Bas_div$He)
  Contemp_Gmi_Bas_div$Fis[is.nan(Contemp_Gmi_Bas_div$Fis)] <- 0
  
Contemp_Gmi_Ham_div$Fis <- 1 - (Contemp_Gmi_Ham_div$Ho/Contemp_Gmi_Ham_div$He)
  Contemp_Gmi_Ham_div$Fis[is.nan(Contemp_Gmi_Ham_div$Fis)] <- 0
  
#calculate mean Fis
Hist_Gmi_Bas_mean_Fis <- mean(Hist_Gmi_Bas_div$Fis) #0.00195
Hist_Gmi_Ham_mean_Fis <- mean(Hist_Gmi_Ham_div$Fis) #-0.01153
Contemp_Gmi_Bas_mean_Fis <- mean(Contemp_Gmi_Bas_div$Fis) #0.01332
Contemp_Gmi_Ham_mean_Fis <- mean(Contemp_Gmi_Ham_div$Fis) #-0.00946
  
Contemp_Gmi_Bas_mean_Fis - Hist_Gmi_Bas_mean_Fis #0.01137
Contemp_Gmi_Ham_mean_Fis - Hist_Gmi_Ham_mean_Fis #0.00207
  
#write out
#write.csv(Hist_Gmi_Bas_div, "Data/Gmi_Ham/Hist_Gmi_Bas_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Hist_Gmi_Ham_div, "Data/Gmi_Ham/Hist_Gmi_Ham_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Contemp_Gmi_Bas_div, "Data/Gmi_Ham/Contemp_Gmi_Bas_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Contemp_Gmi_Ham_div, "Data/Gmi_Ham/Contemp_Gmi_Ham_div.csv", quote = FALSE, row.names = TRUE)
  
#### Ho & He bootstrapping ####
  
#read in if running separately
#Hist_Gmi_Bas_div <- read.csv("Data/Gmi_Ham/Hist_Gmi_Bas_div.csv")
#  colnames(Hist_Gmi_div) <- c("locus", "Ho", "He", "Fis")
#Hist_Gmi_Ham_div <- read.csv("Data/Gmi_Ham/Hist_Gmi_Ham_div.csv")
#  colnames(Hist_Gmi_Ham_div) <- c("locus", "Ho", "He", "Fis")
#Contemp_Gmi_Bas_div <- read.csv("Data/Gmi_Ham/Contemp_Gmi_Bas_div.csv")
#  colnames(Contemp_Gmi_Bas_div) <- c("locus", "Ho", "He", "Fis")
#Contemp_Gmi_Ham_div <- read.csv("Data/Gmi_Ham/Contemp_Gmi_Ham_div.csv")
#  colnames(Contemp_Gmi_Ham_div) <- c("locus", "Ho", "He", "Fis")
  
## Ho bootstrapping ##
#historical Bas Ho
boot_Hist_Gmi_Bas_Ho <- boot(data = Hist_Gmi_Bas_div$Ho, statistic = samp_mean, R = 1000) #1000 permutations of Ho
  Hist_Gmi_Bas_Ho_95ci <- boot.ci(boot_Hist_Gmi_Bas_Ho, conf = 0.95, type = "norm") #get 95% CI for Ho
  Hist_Gmi_Bas_Ho_95ci_normal <- Hist_Gmi_Bas_Ho_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for Ho
 
#historical Ham Ho
boot_Hist_Gmi_Ham_Ho <- boot(data = Hist_Gmi_Ham_div$Ho, statistic = samp_mean, R = 1000)
  Hist_Gmi_Ham_Ho_95ci <- boot.ci(boot_Hist_Gmi_Ham_Ho, conf = 0.95, type = "norm")
  Hist_Gmi_Ham_Ho_95ci_normal <- Hist_Gmi_Ham_Ho_95ci$normal
  
#contemporary Bas Ho
boot_contemp_Gmi_Bas_Ho <- boot(data = Contemp_Gmi_Bas_div$Ho, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Bas_Ho_95ci <- boot.ci(boot_contemp_Gmi_Bas_Ho, conf = 0.95, type = "norm")
  Contemp_Gmi_Bas_Ho_95ci_normal <- Contemp_Gmi_Bas_Ho_95ci$normal   
  
#contemporary Ham Ho
boot_contemp_Gmi_Ham_Ho <- boot(data = Contemp_Gmi_Ham_div$Ho, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Ham_Ho_95ci <- boot.ci(boot_contemp_Gmi_Ham_Ho, conf = 0.95, type = "norm")
  Contemp_Gmi_Ham_Ho_95ci_normal <- Contemp_Gmi_Ham_Ho_95ci$normal 
 
## He bootstrapping ## 
#historical Bas He
boot_Hist_Gmi_Bas_He <- boot(data = Hist_Gmi_Bas_div$He, statistic = samp_mean, R = 1000)
  Hist_Gmi_Bas_He_95ci <- boot.ci(boot_Hist_Gmi_Bas_He, conf = 0.95, type = "norm")
  Hist_Gmi_Bas_He_95ci_normal <- Hist_Gmi_Bas_He_95ci$normal
  
#historical Ham He
boot_Hist_Gmi_Ham_He <- boot(data = Hist_Gmi_Ham_div$He, statistic = samp_mean, R = 1000)
  Hist_Gmi_Ham_He_95ci <- boot.ci(boot_Hist_Gmi_Ham_He, conf = 0.95, type = "norm")
  Hist_Gmi_Ham_He_95ci_normal <- Hist_Gmi_Ham_He_95ci$normal
  
#contemporary Bas He
boot_contemp_Gmi_Bas_He <- boot(data = Contemp_Gmi_Bas_div$He, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Bas_He_95ci <- boot.ci(boot_contemp_Gmi_Bas_He, conf = 0.95, type = "norm")
  Contemp_Gmi_Bas_He_95ci_normal <- Contemp_Gmi_Bas_He_95ci$normal   
  
#contemporary Ham He
boot_contemp_Gmi_Ham_He <- boot(data = Contemp_Gmi_Ham_div$He, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Ham_He_95ci <- boot.ci(boot_contemp_Gmi_Ham_He, conf = 0.95, type = "norm")
  Contemp_Gmi_Ham_He_95ci_normal <- Contemp_Gmi_Ham_He_95ci$normal 
  
## Fis bootstrapping ##
#historical Bas Fis
boot_Hist_Gmi_Bas_Fis <- boot(data = Hist_Gmi_Bas_div$Fis, statistic = samp_mean, R = 1000)
  Hist_Gmi_Bas_Fis_95ci <- boot.ci(boot_Hist_Gmi_Bas_Fis, conf = 0.95, type = "norm")
  Hist_Gmi_Bas_Fis_95ci_normal <- Hist_Gmi_Bas_Fis_95ci$normal
  
#historical Ham Fis
boot_Hist_Gmi_Ham_Fis <- boot(data = Hist_Gmi_Ham_div$Fis, statistic = samp_mean, R = 1000)
  Hist_Gmi_Ham_Fis_95ci <- boot.ci(boot_Hist_Gmi_Ham_Fis, conf = 0.95, type = "norm")
  Hist_Gmi_Ham_Fis_95ci_normal <- Hist_Gmi_Ham_Fis_95ci$normal
  
#contemporary Bas Fis
boot_contemp_Gmi_Bas_Fis <- boot(data = Contemp_Gmi_Bas_div$Fis, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Bas_Fis_95ci <- boot.ci(boot_contemp_Gmi_Bas_Fis, conf = 0.95, type = "norm")
  Contemp_Gmi_Bas_Fis_95ci_normal <- Contemp_Gmi_Bas_Fis_95ci$normal   
  
#contemporary Ham Fis
boot_contemp_Gmi_Ham_Fis <- boot(data = Contemp_Gmi_Ham_div$Fis, statistic = samp_mean, R = 1000)
  Contemp_Gmi_Ham_Fis_95ci <- boot.ci(boot_contemp_Gmi_Ham_Fis, conf = 0.95, type = "norm")
  Contemp_Gmi_Ham_Fis_95ci_normal <- Contemp_Gmi_Ham_Fis_95ci$normal 
  
## write out Ho, He, Fis summary tables ##
Ho_mean <- as.data.frame(c(0.10647, 0.10784, 0.08964, 0.09901))
He_mean <- as.data.frame(c(0.10545, 0.10195, 0.08800, 0.09350))
Fis_mean <- as.data.frame(c(0.00195, -0.01153, 0.01332, -0.00946))
  
#Ho summary table
Ho_ci <- rbind(Hist_Gmi_Bas_Ho_95ci_normal, Hist_Gmi_Ham_Ho_95ci_normal, 
               Contemp_Gmi_Bas_Ho_95ci_normal, Contemp_Gmi_Ham_Ho_95ci_normal) #combine df w/ci info for each pop into one dataframe
Ho_sum <- cbind(Ho_mean, Ho_ci)
  colnames(Ho_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Ho_sum$Pop <- c("Hist - Bas", "Hist - Ham", "Contemp - Bas", "Contemp - Ham")
  
Ho_sum$diff_lower <- Ho_sum$mean - Ho_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
Ho_sum$diff_upper <- Ho_sum$`97.5_per` - Ho_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization
  
#He summary table
He_ci <- rbind(Hist_Gmi_Bas_He_95ci_normal, Hist_Gmi_Ham_He_95ci_normal, 
               Contemp_Gmi_Bas_He_95ci_normal, Contemp_Gmi_Ham_He_95ci_normal)
He_sum <- cbind(He_mean, He_ci)
  colnames(He_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  He_sum$Pop <- c("Hist - Bas", "Hist - Ham", "Contemp - Bas", "Contemp - Ham")
  
He_sum$diff_lower <- He_sum$mean - He_sum$`2.5_per`
He_sum$diff_upper <- He_sum$`97.5_per` - He_sum$mean
  
#Fis summary table
Fis_ci <- rbind(Hist_Gmi_Bas_Fis_95ci_normal, Hist_Gmi_Ham_Fis_95ci_normal, 
                Contemp_Gmi_Bas_Fis_95ci_normal, Contemp_Gmi_Ham_Fis_95ci_normal)
Fis_sum <- cbind(Fis_mean, Fis_ci)
  colnames(Fis_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Fis_sum$Pop <- c("Hist - Bas", "Hist - Ham", "Contemp - Bas", "Contemp - Ham")
  
Fis_sum$diff_lower <- Fis_sum$mean - Fis_sum$`2.5_per`
Fis_sum$diff_upper <- Fis_sum$`97.5_per` - Fis_sum$mean
  
#combine dataframes
div_sum_all <- as.data.frame(rbind(Ho_sum, He_sum, Fis_sum))
div_sum_all$metric <- c("Ho", "Ho", "Ho", "Ho", 
                        "He", "He", "He", "He", 
                        "Fis", "Fis", "Fis", "Fis")
  
#write out data
#write.csv(div_sum_all, "Data/Gmi_Ham/Gmi_A_div_cis.csv")

## Calculate 95% for diff through time ##
#Ho diff
Gmi_Ham_Ho_boot_diff <- as.data.frame(boot_Hist_Gmi_Ham_Ho$t)
  Gmi_Ham_Ho_boot_diff$Contemp_Ho <- boot_contemp_Gmi_Ham_Ho$t
  colnames(Gmi_Ham_Ho_boot_diff) <- c("Hist_Ho", "Contemp_Ho")

Gmi_Ham_Ho_boot_diff$Ho_diff <- Gmi_Ham_Ho_boot_diff$Contemp_Ho - Gmi_Ham_Ho_boot_diff$Hist_Ho
  Gmi_Ham_Ho_boot_diff_CIs <- quantile(Gmi_Ham_Ho_boot_diff$Ho_diff, c(0.025, 0.975)) #-0.01373, -0.00375

#He diff
Gmi_Ham_He_boot_diff <- as.data.frame(boot_Hist_Gmi_Ham_He$t)
  Gmi_Ham_He_boot_diff$Contemp_He <- boot_contemp_Gmi_Ham_He$t
  colnames(Gmi_Ham_He_boot_diff) <- c("Hist_He", "Contemp_He")
  
Gmi_Ham_He_boot_diff$He_diff <- Gmi_Ham_He_boot_diff$Contemp_He - Gmi_Ham_He_boot_diff$Hist_He
  Gmi_Ham_He_boot_diff_CIs <- quantile(Gmi_Ham_He_boot_diff$He_diff, c(0.025, 0.975)) #-0.01258, -0.00387

#Fis diff
Gmi_Ham_Fis_boot_diff <- as.data.frame(boot_Hist_Gmi_Ham_Fis$t)
  Gmi_Ham_Fis_boot_diff$Contemp_Fis <- boot_contemp_Gmi_Ham_Fis$t
  colnames(Gmi_Ham_Fis_boot_diff) <- c("Hist_Fis", "Contemp_Fis")
  
Gmi_Ham_Fis_boot_diff$Fis_diff <- Gmi_Ham_Fis_boot_diff$Contemp_Fis - Gmi_Ham_Fis_boot_diff$Hist_Fis
  Gmi_Ham_Fis_boot_diff_CIs <- quantile(Gmi_Ham_Fis_boot_diff$Fis_diff, c(0.025, 0.975)) #-0.00275, 0.00674

## Calculate empirical p-value ##
#Ho
Ho_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, Ho < -0.00883))) #number of permutations with value greater than observed (really less than)
  (Ho_perm_greater + 1)/10001 #0.0055

#He
He_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, He < -0.00826)))
  (He_perm_greater + 1)/10001 #0.0006

#Fis
Fis_perm_greater <- as.numeric(nrow(subset(Gmi_Ham_permutation, Fis > 0.00207)))
  (Fis_perm_greater + 1)/10001 #0.2556
  
#### Gmi Ho, He & Fis visualization ####
  
#read in data (if running separately)
div_sum_all <- read.csv("Data/Gmi_Ham/Gmi_A_div_cis.csv", row.names = 1)
  
#add factor level
div_sum_all$Pop <- factor(div_sum_all$Pop, levels = c("Hist - Bas", "Hist - Ham", 
                                                      "Contemp - Bas", "Contemp - Ham"))
  
#het only dataframe
het_sum_all <- subset(div_sum_all, div_sum_all$metric != "Fis")
  
#fis dataframe
Fis_sum <- subset(div_sum_all, div_sum_all$metric == "Fis")
  
#plot of mean Ho & He w/95% CI error bars
Het_plot <- ggplot(data = het_sum_all, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#a25505", "#16537e", "#e1bf9b", "#8aa9be")) +
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
Fis_plot <- ggplot(data = Fis_sum, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 2, color = "black") + 
scale_color_manual(values = c("#a25505", "#16537e", "#e1bf9b", "#8aa9be")) +
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
  geom_vline(aes(xintercept = -0.01373), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.00375), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00883), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  annotate("text", x = -0.013, y = 97, label = "E", size = 40) + 
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
  geom_vline(aes(xintercept = -0.01258), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.00387), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00826), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  annotate("text", x = -0.012, y = 140, label = "B", size = 40) + 
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
  geom_vline(aes(xintercept = -0.00275), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00674), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = 0.00207), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  annotate("text", x = -0.021, y = 68, label = "K", size = 40) + 
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
Ela_Ho_means <- colMeans(sum_stats_Ela$Ho) #1: 0.12544, 2: 0.12091
Ela_Hs_means <- colMeans(sum_stats_Ela$Hs) #1: 0.12589, 2: 0.11705

Ela_Ho_means_diff <- diff(colMeans(sum_stats_Ela$Ho)) #-0.00453
Ela_Hs_means_diff <- diff(colMeans(sum_stats_Ela$Hs)) #-0.00884

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
Hist_Ela_mean_Fis <- mean(Hist_Ela_div$Fis) #0.00557
Contemp_Ela_mean_Fis <- mean(Contemp_Ela_div$Fis) #0.00815

Contemp_Ela_mean_Fis - Hist_Ela_mean_Fis #0.000258

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
Ho_mean <- as.data.frame(c(0.12544, 0.12091))
He_mean <- as.data.frame(c(0.12589, 0.11705))
Fis_mean <- as.data.frame(c(0.00556, 0.00815))

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
#write.csv(div_sum_all, "PIRE_Ela_Ham/Ela_div_cis.csv")

## Calculate 95% for diff through time ##
#Ho diff
Ela_Ho_boot_diff <- as.data.frame(boot_Hist_Ela_Ho$t)
  Ela_Ho_boot_diff$Contemp_Ho <- boot_Contemp_Ela_Ho$t
  colnames(Ela_Ho_boot_diff) <- c("Hist_Ho", "Contemp_Ho")
  
Ela_Ho_boot_diff$Ho_diff <- Ela_Ho_boot_diff$Contemp_Ho - Ela_Ho_boot_diff$Hist_Ho
  Ela_Ho_boot_diff_CIs <- quantile(Ela_Ho_boot_diff$Ho_diff, c(0.025, 0.975)) #-0.00764, -0.00161
  
#He diff
Ela_He_boot_diff <- as.data.frame(boot_Hist_Ela_He$t)
  Ela_He_boot_diff$Contemp_He <- boot_Contemp_Ela_He$t
  colnames(Ela_He_boot_diff) <- c("Hist_He", "Contemp_He")
  
Ela_He_boot_diff$He_diff <- Ela_He_boot_diff$Contemp_He - Ela_He_boot_diff$Hist_He
  Ela_He_boot_diff_CIs <- quantile(Ela_He_boot_diff$He_diff, c(0.025, 0.975)) #-0.01160, -0.00618
  
#Fis diff
Ela_Fis_boot_diff <- as.data.frame(boot_Hist_Ela_Fis$t)
  Ela_Fis_boot_diff$Contemp_Fis <- boot_Contemp_Ela_Fis$t
  colnames(Ela_Fis_boot_diff) <- c("Hist_Fis", "Contemp_Fis")
  
Ela_Fis_boot_diff$Fis_diff <- Ela_Fis_boot_diff$Contemp_Fis - Ela_Fis_boot_diff$Hist_Fis
  Ela_Fis_boot_diff_CIs <- quantile(Ela_Fis_boot_diff$Fis_diff, c(0.025, 0.975)) #-0.00078, 0.00588
  
## Calculate empirical p-value ##
#Ho
Ho_perm_greater <- as.numeric(nrow(subset(Ela_permutation, Ho < -0.00453))) #number of permutations with value greater than observed (really less than)
  (Ho_perm_greater + 1)/10001 #0.0137
  
#He
He_perm_greater <- as.numeric(nrow(subset(Ela_permutation, He < -0.00884)))
  (He_perm_greater + 1)/10001 #0.0001
  
#Fis
Fis_perm_greater <- as.numeric(nrow(subset(Ela_permutation, Fis > 0.00259)))
  (Fis_perm_greater + 1)/10001 #0.1009
  
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
Het_plot <- ggplot(data = het_sum_all, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
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
Fis_plot <- ggplot(data = Fis_sum, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
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
  geom_vline(aes(xintercept = -0.00764), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.00161), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00453), linewidth = 6, color = "#16537e") + #"real" mean diff
  annotate("text", x = -0.0115, y = 195, label = "F", size = 40) + 
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
  geom_vline(aes(xintercept = -0.0116), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.00618), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00884), linewidth = 6, color = "#16537e") + #"real" mean diff
  annotate("text", x = -0.0113, y = 270, label = "C", size = 40) + 
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
  geom_vline(aes(xintercept = -0.00078), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00588), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = 0.00259), linewidth = 6, color = "#16537e") + #"real" mean diff
  annotate("text", x = -0.036, y = 75, label = "L", size = 40) + 
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

######## Aen Ho, He & Fis estimates ########

#add population level data
pop <- c(rep(1, times = 55), rep(2, times = 92))
  pop(Aen_genind) <- pop #1 = Albatross, 2 = Contemporary
  Aen_genind #check to make sure pops read right (should range from 55-92 individuals)

#calculate diversity metrics w/in pops
sum_stats_Aen <- basic.stats(Aen_genind)

#### Calculate mean Ho, He ####

#mean Ho & He for each population
Aen_Ho_means <- colMeans(sum_stats_Aen$Ho) #1: 0.04544, 2: 0.04021
Aen_Hs_means <- colMeans(sum_stats_Aen$Hs) #1: 0.05798, 2: 0.03557

Aen_Ho_means_diff <- diff(colMeans(sum_stats_Aen$Ho)) #-0.00534
Aen_Hs_means_diff <- diff(colMeans(sum_stats_Aen$Hs)) #-0.02241

##### Pull out per-marker estimates ####

#pull out Ho for each population
Hist_Aen_Ho <- sum_stats_Aen$Ho[,1]
Contemp_Aen_Ho <- sum_stats_Aen$Ho[,2]

#pull out He for each population
Hist_Aen_He <- sum_stats_Aen$Hs[,1]
Contemp_Aen_He <- sum_stats_Aen$Hs[,2]

#bind Ho & He into dataframe for each pop
Hist_Aen_div <- as.data.frame(cbind(Hist_Aen_Ho, Hist_Aen_He))
  head(Hist_Aen_div) #check content is right
  colnames(Hist_Aen_div) <- c("Ho", "He")

Contemp_Aen_div <- as.data.frame(cbind(Contemp_Aen_Ho, Contemp_Aen_He))
  colnames(Contemp_Aen_div) <- c("Ho", "He")

##### Calculate mean Fis ####

#calculate Fis for each population
Hist_Aen_div$Fis <- 1 - (Hist_Aen_div$Ho/Hist_Aen_div$He)
  Hist_Aen_div$Fis[is.nan(Hist_Aen_div$Fis)] <- 0

Contemp_Aen_div$Fis <- 1 - (Contemp_Aen_div$Ho/Contemp_Aen_div$He)
  Contemp_Aen_div$Fis[is.nan(Contemp_Aen_div$Fis)] <- 0

#calculate mean Fis
Hist_Aen_mean_Fis <- mean(Hist_Aen_div$Fis) #0.09523
Contemp_Aen_mean_Fis <- mean(Contemp_Aen_div$Fis) #-0.10618

Contemp_Aen_mean_Fis - Hist_Aen_mean_Fis #-0.10618

#write out
#write.csv(Hist_Aen_div, "Data/Aen_Ham/Hist_Aen_div.csv", quote = FALSE, row.names = TRUE)
#write.csv(Contemp_Aen_div, "Data/Aen_Ham/Contemp_Aen_div.csv", quote = FALSE, row.names = TRUE)

#### Ho & He bootstrapping ####

#read in if running separately
#Hist_Aen_div <- read.csv("Data/Aen_Ham/Hist_Aen_div.csv")
  #colnames(Hist_Aen_div) <- c("locus", "Ho", "He", "Fis")
#Contemp_Aen_div <- read.csv("Data/Aen_Ham/Contemp_Aen_div.csv")
  #colnames(Contemp_Aen_div) <- c("locus", "Ho", "He", "Fis")

## Ho bootstrapping ##
#historical Ho
boot_Hist_Aen_Ho <- boot(data = Hist_Aen_div$Ho, statistic = samp_mean, R = 1000) #1000 permutations of Ho
  Hist_Aen_Ho_95ci <- boot.ci(boot_Hist_Aen_Ho, conf = 0.95, type = "norm") #get 95% CI for Ho
  Hist_Aen_Ho_95ci_normal <- Hist_Aen_Ho_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for Ho

#contemporary Ho
boot_Contemp_Aen_Ho <- boot(data = Contemp_Aen_div$Ho, statistic = samp_mean, R = 1000)
  Contemp_Aen_Ho_95ci <- boot.ci(boot_Contemp_Aen_Ho, conf = 0.95, type = "norm")
  Contemp_Aen_Ho_95ci_normal <- Contemp_Aen_Ho_95ci$normal

## He bootstrapping ##
#historical He
boot_Hist_Aen_He <- boot(data = Hist_Aen_div$He, statistic = samp_mean, R = 1000)
  Hist_Aen_He_95ci <- boot.ci(boot_Hist_Aen_He, conf = 0.95, type = "norm")
  Hist_Aen_He_95ci_normal <- Hist_Aen_He_95ci$normal

#contemporary He
boot_Contemp_Aen_He <- boot(data = Contemp_Aen_div$He, statistic = samp_mean, R = 1000)
  Contemp_Aen_He_95ci <- boot.ci(boot_Contemp_Aen_He, conf = 0.95, type = "norm")
  Contemp_Aen_He_95ci_normal <- Contemp_Aen_He_95ci$normal

## Fis bootstrapping ##
#historical Fis
boot_Hist_Aen_Fis <- boot(data = Hist_Aen_div$Fis, statistic = samp_mean, R = 1000)
  Hist_Aen_Fis_95ci <- boot.ci(boot_Hist_Aen_Fis, conf = 0.95, type = "norm")
  Hist_Aen_Fis_95ci_normal <- Hist_Aen_Fis_95ci$normal

#contemporary Fis
boot_Contemp_Aen_Fis <- boot(data = Contemp_Aen_div$Fis, statistic = samp_mean, R = 1000)
  Contemp_Aen_Fis_95ci <- boot.ci(boot_Contemp_Aen_Fis, conf = 0.95, type = "norm")
  Contemp_Aen_Fis_95ci_normal <- Contemp_Aen_Fis_95ci$normal

## write out Ho, He, Fis summary tables ##
Ho_mean <- as.data.frame(c(0.0454, 0.0402))
He_mean <- as.data.frame(c(0.0580, 0.0356))
Fis_mean <- as.data.frame(c(0.0952, -0.0109))

#Ho summary table
Ho_ci <- rbind(Hist_Aen_Ho_95ci_normal, Contemp_Aen_Ho_95ci_normal) #combine df w/ci info for each pop into one dataframe
Ho_sum <- cbind(Ho_mean, Ho_ci)
  colnames(Ho_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Ho_sum$Pop <- c("Hist", "Contemp")

Ho_sum$diff_lower <- Ho_sum$mean - Ho_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
Ho_sum$diff_upper <- Ho_sum$`97.5_per` - Ho_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#He summary table
He_ci <- rbind(Hist_Aen_He_95ci_normal, Contemp_Aen_He_95ci_normal)
He_sum <- cbind(He_mean, He_ci)
  colnames(He_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  He_sum$Pop <- c("Hist", "Contemp")

He_sum$diff_lower <- He_sum$mean - He_sum$`2.5_per`
He_sum$diff_upper <- He_sum$`97.5_per` - He_sum$mean

#Fis summary table
Fis_ci <- rbind(Hist_Aen_Fis_95ci_normal, Contemp_Aen_Fis_95ci_normal)
Fis_sum <- cbind(Fis_mean, Fis_ci)
  colnames(Fis_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  Fis_sum$Pop <- c("Hist", "Contemp")

Fis_sum$diff_lower <- Fis_sum$mean - Fis_sum$`2.5_per`
Fis_sum$diff_upper <- Fis_sum$`97.5_per` - Fis_sum$mean

#combine dataframes
div_sum_all <- as.data.frame(rbind(Ho_sum, He_sum, Fis_sum))
  div_sum_all$metric <- c("Ho", "Ho", "He", "He", "Fis", "Fis")
    
#write out data
#write.csv(div_sum_all, "Data/Aen_Ham/Aen_div_cis.csv")

## Calculate 95% for diff through time ##
#Ho diff
Aen_Ho_boot_diff <- as.data.frame(boot_Hist_Aen_Ho$t)
  Aen_Ho_boot_diff$Contemp_Ho <- boot_Contemp_Aen_Ho$t
  colnames(Aen_Ho_boot_diff) <- c("Hist_Ho", "Contemp_Ho")
  
Aen_Ho_boot_diff$Ho_diff <- Aen_Ho_boot_diff$Contemp_Ho - Aen_Ho_boot_diff$Hist_Ho
  Aen_Ho_boot_diff_CIs <- quantile(Aen_Ho_boot_diff$Ho_diff, c(0.025, 0.975)) #-0.00764, -0.00161
  
#He diff
Aen_He_boot_diff <- as.data.frame(boot_Hist_Aen_He$t)
  Aen_He_boot_diff$Contemp_He <- boot_Contemp_Aen_He$t
  colnames(Aen_He_boot_diff) <- c("Hist_He", "Contemp_He")
  
Aen_He_boot_diff$He_diff <- Aen_He_boot_diff$Contemp_He - Aen_He_boot_diff$Hist_He
  Aen_He_boot_diff_CIs <- quantile(Aen_He_boot_diff$He_diff, c(0.025, 0.975)) #-0.01160, -0.00618
  
#Fis diff
Aen_Fis_boot_diff <- as.data.frame(boot_Hist_Aen_Fis$t)
  Aen_Fis_boot_diff$Contemp_Fis <- boot_Contemp_Aen_Fis$t
  colnames(Aen_Fis_boot_diff) <- c("Hist_Fis", "Contemp_Fis")
  
Aen_Fis_boot_diff$Fis_diff <- Aen_Fis_boot_diff$Contemp_Fis - Aen_Fis_boot_diff$Hist_Fis
  Aen_Fis_boot_diff_CIs <- quantile(Aen_Fis_boot_diff$Fis_diff, c(0.025, 0.975)) #-0.00078, 0.00588

## Calculate empirical p-value ##
#Ho
Ho_perm_greater <- as.numeric(nrow(subset(Aen_permutation, Ho < -0.00534))) #number of permutations with value greater than observed (really less than)
  (Ho_perm_greater + 1)/10001 #0.0002
  
#He
He_perm_greater <- as.numeric(nrow(subset(Aen_permutation, He < -0.02241)))
  (He_perm_greater + 1)/10001 #0.0001
  
#Fis
Fis_perm_greater <- as.numeric(nrow(subset(Aen_permutation, Fis < -0.10618)))
  (Fis_perm_greater + 1)/10001 #0.0001
  
#### Aen Ho, He & Fis visualization ####

#read in data (if running separately)
#div_sum_all <- read.csv("PIRE_Aen_Ham/Aen_div_cis.csv", row.names = 1)

#add factor level
div_sum_all$Pop <- factor(div_sum_all$Pop, levels = c("Hist", "Contemp"))

#het only dataframe
het_sum_all <- subset(div_sum_all, div_sum_all$metric != "Fis")

#fis dataframe
Fis_sum <- subset(div_sum_all, div_sum_all$metric == "Fis")

#plot of mean Ho & He w/95% CI error bars
Het_plot <- ggplot(data = het_sum_all, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ylim(0, 0.15) + ggtitle("Aen Ho & He 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
Het_plot

#plot of mean Fis w/95% CI error bars
Fis_plot <- ggplot(data = Fis_sum, aes(x = metric, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 2, color = "black") + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ggtitle("Aen Fis 95% CIs") + theme_bw() + 
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
Ho_null_plot <- ggplot(Aen_permutation, aes(x = Ho)) + 
  geom_density(color = "#e3ccb4", fill = "#e3ccb4") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.01363), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00304), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00534), linewidth = 6, color = "#a25505") + #"real" mean diff
  annotate("text", x = -0.013, y = 242, label = "D", size = 40) + 
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
He_null_plot <- ggplot(Aen_permutation, aes(x = He)) + 
  geom_density(color = "#e3ccb4", fill = "#e3ccb4") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.03060), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.01465), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.02241), linewidth = 6, color = "#a25505") + #"real" mean diff
  annotate("text", x = -0.0295, y = 185, label = "A", size = 40) + 
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
Fis_null_plot <- ggplot(Aen_permutation, aes(x = Fis)) + 
  geom_density(color = "#e3ccb4", fill = "#e3ccb4") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.12066), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.25 CI lower
  geom_vline(aes(xintercept = -0.09312), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.10618), linewidth = 6, color = "#a25505") + #"real" mean diff
  annotate("text", x = -0.1155, y = 13.5, label = "J", size = 40) + 
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

###########################################################################################################################

#to convert genepop file to relatedness input file
#remove header rows first in unix
library(tidyverse)
Lle_genepop <- read.table("GENEPOP3.gen", sep = " ")
Lle_genepop2 <- Lle_genepop[,-2] #do twice
dim(Lle_genepop2)
test_2[test_2=="0"] <- "000000"
test_2 <- unite(test, all, sep = "", 2:40, remove = TRUE)
colnames <- c("ind", "genotypes")
colnames(test_2) <- colnames
n1 <- nchar(as.character(test_2$genotypes[2])) - 3
s1 <- seq(3, n1, by = 3)
nm1 <- paste0("x", seq_len(length(s1) +1))
test_3 <- test_2 %>% separate(genotypes, into = nm1, sep = s1)
write.table(Lle_genepop5, "relatednessformat.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#relatedness stuff
library(related)
rel_info <- readgenotypedata("relatednessformat.txt")

#calculate point estimates of relatedness
rel_output <- coancestry("relatednessformat.txt", wang = 2)

relatedness <- rel_output$relatedness
rel_Alb<- subset(relatedness, grepl("^Lle-AHam", ind1.id) & grepl("^Lle-AHam", ind2.id))
rel_Con<- subset(relatedness, grepl("^Lle-CNas", ind1.id) & grepl("^Lle-CNas", ind2.id))

rel_inpop_df <- rbind(rel_Alb, rel_Con)

################################################################################################################################################
