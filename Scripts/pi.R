############################################### Script to calculate pi ################################################

#from https://pixy.readthedocs.io/en/latest/plotting.html
#using all sites VCF (so know which are monomorphic and which are missing bc of quality control)
#calculated in 1000 bp windows
#bootstrapping across sites

#probe targets are ~120 bp in length, so filter to windows with at least 100 genotyped sites = most of that probe region is called

#################################################################################################################################

remove(list = ls())

#load libraries
library(here) #v.1.0.1
library(tidyverse) #v.2.0.0
library(boot) #v.1.3.28
library(scales) #v.1.2.1

#write function to calculate mean
samp_mean <- function(x, i) {
  mean(x[i])
}

#read in pixy data
Aen_pi <- read.table(here("Data/Aen_Ham/pixy", "pixy_pi.txt"), header = TRUE)
Gmi_pi <- read.table(here("Data/Gmi_Ham/pixy/", "pixy_pi.txt"), header = TRUE)
Ela_pi <- read.table(here("Data/Ela_Ham/pixy/", "pixy_pi.txt"), header = TRUE)

#read in permutation data
Aen_pi_permutation <- read.csv(here("Data/Aen_Ham", "Aen_permutation_pi.csv"))
Gmi_pi_permutation <- read.csv(here("Data/Gmi_Ham", "Gmi_Ham_permutation_pi.csv"))
Ela_pi_permutation <- read.csv(here("Data/Ela_Ham", "Ela_permutation_pi.csv"))

##########################################################################################################################

######## Gmi pi estimates ########

#subsetting dataframe to windows with >100 genotyped sites
Gmi_pi_100bp <- subset(Gmi_pi, no_sites >= 100 & 
                         pop != "Gmi-AHam-B" & pop != "Gmi-CBat-B") #left with 10324 windows (2581/time point)

#create dataframe for each pop
Gmi_ABas_pi_100bp <- subset(Gmi_pi_100bp, pop == "Gmi-ABas-A")
Gmi_AHam_pi_100bp <- subset(Gmi_pi_100bp, pop == "Gmi-AHam-A")
Gmi_CBas_pi_100bp <- subset(Gmi_pi_100bp, pop == "Gmi-CBas-A")
Gmi_CHam_pi_100bp <- subset(Gmi_pi_100bp, pop == "Gmi-CBat-A")

#### calculate mean & median pi/pop ####

#calculate mean pi
Gmi_ABas_pi_100bp_mean <- mean(Gmi_ABas_pi_100bp$avg_pi) #0.00317
Gmi_AHam_pi_100bp_mean <- mean(Gmi_AHam_pi_100bp$avg_pi) #0.00290
Gmi_CBas_pi_100bp_mean <- mean(Gmi_CBas_pi_100bp$avg_pi) #0.00242
Gmi_CHam_pi_100bp_mean <- mean(Gmi_CHam_pi_100bp$avg_pi) #0.00279

#calculate median pi
Gmi_ABas_pi_100bp_median <- median(Gmi_ABas_pi_100bp$avg_pi) #0.00186
Gmi_AHam_pi_100bp_median <- median(Gmi_AHam_pi_100bp$avg_pi) #0.00167
Gmi_CBas_pi_100bp_median <- median(Gmi_CBas_pi_100bp$avg_pi) #0.00122
Gmi_CHam_pi_100bp_median <- median(Gmi_CHam_pi_100bp$avg_pi) #0.00153

#### pi bootstrapping ####

#bootstrap for Gmi Bas Albatross
boot_Gmi_ABas_pi_100bp <- boot(data = Gmi_ABas_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_ABas_pi_100bp_95ci <- boot.ci(boot_Gmi_ABas_pi_100bp, conf = 0.95, type = "norm")
  Gmi_ABas_pi_100bp_95ci_normal <- Gmi_ABas_pi_100bp_95ci$normal

#bootstrap for Gmi Ham Albatross
boot_Gmi_AHam_pi_100bp <- boot(data = Gmi_AHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_AHam_pi_100bp_95ci <- boot.ci(boot_Gmi_AHam_pi_100bp, conf = 0.95, type = "norm")
  Gmi_AHam_pi_100bp_95ci_normal <- Gmi_AHam_pi_100bp_95ci$normal

#bootstrap for Gmi Bas Contemporary
boot_Gmi_CBas_pi_100bp <- boot(data = Gmi_CBas_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_CBas_pi_100bp_95ci <- boot.ci(boot_Gmi_CBas_pi_100bp, conf = 0.95, type = "norm")
  Gmi_CBas_pi_100bp_95ci_normal <- Gmi_CBas_pi_100bp_95ci$normal

#bootstrap for Gmi Ham Contemporary
boot_Gmi_CHam_pi_100bp <- boot(data = Gmi_CHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_CHam_pi_100bp_95ci <- boot.ci(boot_Gmi_CHam_pi_100bp, conf = 0.95, type = "norm")
  Gmi_CHam_pi_100bp_95ci_normal <- Gmi_CHam_pi_100bp_95ci$normal

## write out pi summary table ##
pi_mean <- as.data.frame(c(0.00317, 0.00290, 0.00242, 0.00279))

#pi summary table
pi_ci <- rbind(Gmi_ABas_pi_100bp_95ci_normal, Gmi_AHam_pi_100bp_95ci_normal, 
               Gmi_CBas_pi_100bp_95ci_normal, Gmi_CHam_pi_100bp_95ci_normal) #combine df w/ci info for each pop into one dataframe
pi_sum <- cbind(pi_mean, pi_ci)
  colnames(pi_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  pi_sum$Pop <- c("Hist - Bas", "Hist - Ham", "Contemp - Bas", "Contemp - Ham")

pi_sum$diff_lower <- pi_sum$mean - pi_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_sum$diff_upper <- pi_sum$`97.5_per` - pi_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
#write.csv(pi_sum, "Data/Gmi_Ham/Gmi_A_pi_cis.csv")

## Calculate 95% for diff through time ##
Gmi_pi_boot_diff <- as.data.frame(boot_Gmi_AHam_pi_100bp$t)
  Gmi_pi_boot_diff$Contemp_pi <- boot_Gmi_CHam_pi_100bp$t
  colnames(Gmi_pi_boot_diff) <- c("Hist_pi", "Contemp_pi")

Gmi_pi_boot_diff$pi_diff <- Gmi_pi_boot_diff$Contemp_pi - Gmi_pi_boot_diff$Hist_pi
  Gmi_pi_boot_diff_CIs <- quantile(Gmi_pi_boot_diff$pi_diff, c(0.025, 0.975)) #-0.00031, 0.00010

## Create null distribution of change through time ##
#remove Bas from dataframe to resample
Gmi_Ham_pi_100bp <- subset(Gmi_pi_100bp, pop != "Gmi-ABas-A" & pop != "Gmi-CBas-A") #left with 5162 windows
  
#create vector to populate
permutation_pi <- c()
  
#for loop to randomly sample df 10 000x
for (i in 1:10000) { 
  permutation_AHam <- Gmi_Ham_pi_100bp[sample(nrow(Gmi_Ham_pi_100bp), 2581), ] #random subset = # Alb rows
  permutation_CHam <- Gmi_Ham_pi_100bp[sample(nrow(Gmi_Ham_pi_100bp), 2581), ] #random subset = # Contemp rows
    
  #calculate means
  permutation_AHam_mean <- mean(permutation_AHam$avg_pi)
  permutation_CHam_mean <- mean(permutation_CHam$avg_pi)
    
  #calculate diff
  permutation_pi[i] <- permutation_CHam_mean - permutation_AHam_mean
}
  
#write out
permutation_diff_pi <- as.data.frame(permutation_pi)
  permutation_diff_pi$iteration <- 1:10000
  
#write.csv(permutation_diff_pi, "Data/Gmi_Ham/Gmi_Ham_permutation_pi.csv")
  
## Calculate empirical p-value ##
#number of permutations with value greater than observed (really less than)
pi_perm_greater <- as.numeric(nrow(subset(Gmi_pi_permutation, permutation_pi < -0.00011))) 
  (pi_perm_greater + 1)/10001 #0.0642
  
#### Visualize results ####

#ordering x-axis
pi_sum$Pop <- factor(pi_sum$Pop, levels = c("Hist - Bas", "Contemp - Bas", "Hist - Ham", "Contemp - Ham"))

#plot of mean pi w/95% CI error bars
pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#a25505", "#e1bf9b", "#16537e", "#8aa9be")) +
  ggtitle("Gmi A pi 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 1), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
pi_plot

#### Null distribution plots ####

quantile(Gmi_pi_permutation$permutation_pi, c(0.025, 0.975))

#plot of null pi diff distribution
pi_null_plot <- ggplot(Gmi_pi_permutation, aes(x = permutation_pi)) + 
  geom_density(color = "#bac4b6", fill = "#bac4b6") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00031), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.0001), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00011), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
  annotate("text", x = -0.0003, y = 5200, label = "H", size = 40) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab("Temporal Change in π") + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
pi_null_plot

#### Calculate change through time ####
#selective sweeps --> regions with declining pi could be swept/linked to sweep

## merge Alb & Contemp ##
pi_Alb <- as.data.frame(Gmi_AHam_pi_100bp$chromosome)
  pi_Alb$pi <- Gmi_AHam_pi_100bp$avg_pi
  colnames(pi_Alb) <- c("chromosome", "avg_pi")
pi_Contemp <- as.data.frame(Gmi_CHam_pi_100bp$chromosome)
  pi_Contemp$pi <- Gmi_CHam_pi_100bp$avg_pi
  colnames(pi_Contemp) <- c("chromosome", "avg_pi")

pi_merge <- merge(pi_Alb, pi_Contemp, by = "chromosome")
  colnames(pi_merge) <- c("chromosome", "Alb_avg_pi", "Contemp_avg_pi")
  pi_merge$delta <- pi_merge$Contemp_avg_pi - pi_merge$Alb_avg_pi

#calculate mean and median delta  
mean(na.omit(pi_merge$delta)) #-0.00011 --> losing diversity through time
median(na.omit(pi_merge$delta)) #0 
  
#identify ones that became negative through time
#poss selective sweeps
pot_outliers <- subset(pi_merge, delta < 0)

#diff scatter plot
pi_diff_scatter_plot <- ggplot(data = pi_merge, aes(x = Alb_avg_pi, y = Contemp_avg_pi)) + 
  geom_point(size = 10, color = "#7e81a7", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  ggtitle("Gmi delta pi") + labs(y = "Contemporary", x = "Albatross") + 
theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
pi_diff_scatter_plot

##################################################################################################################################

######## Ela pi estimates ########

#subsetting dataframe to windows with >100 genotyped sites
Ela_pi_100bp <- subset(Ela_pi, no_sites >= 100 & pop != "Ela-AHam") #left with 3646 windows (1823/time point)

#create dataframe for each pop
Ela_AHam_pi_100bp <- subset(Ela_pi_100bp, pop == "Ela-AHam")
Ela_CHam_pi_100bp <- subset(Ela_pi_100bp, pop == "Ela-CNas")

#### calculate mean & median pi/pop ####

#calculate mean pi
Ela_AHam_pi_100bp_mean <- mean(Ela_AHam_pi_100bp$avg_pi) #0.00395
Ela_CHam_pi_100bp_mean <- mean(Ela_CHam_pi_100bp$avg_pi) #0.00374

#calculate median pi
Ela_AHam_pi_100bp_median <- median(Ela_AHam_pi_100bp$avg_pi) #0.00275
Ela_CHam_pi_100bp_median <- median(Ela_CHam_pi_100bp$avg_pi) #0.00263

#### pi bootstrapping ####

#bootstrap for Ela Albatross
boot_Ela_AHam_pi_100bp <- boot(data = Ela_AHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Ela_AHam_pi_100bp_95ci <- boot.ci(boot_Ela_AHam_pi_100bp, conf = 0.95, type = "norm")
  Ela_AHam_pi_100bp_95ci_normal <- Ela_AHam_pi_100bp_95ci$normal

#bootstrap for Ela Contemporary
boot_Ela_CHam_pi_100bp <- boot(data = Ela_CHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Ela_CHam_pi_100bp_95ci <- boot.ci(boot_Ela_CHam_pi_100bp, conf = 0.95, type = "norm")
  Ela_CHam_pi_100bp_95ci_normal <- Ela_CHam_pi_100bp_95ci$normal

## write out pi summary table ##
pi_mean <- as.data.frame(c(0.00395, 0.00374))

#pi summary table
pi_ci <- rbind(Ela_AHam_pi_100bp_95ci_normal, Ela_CHam_pi_100bp_95ci_normal) #combine df w/ci info for each pop into one dataframe
pi_sum <- cbind(pi_mean, pi_ci)
  colnames(pi_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  pi_sum$Pop <- c("Hist", "Contemp")
  
pi_sum$diff_lower <- pi_sum$mean - pi_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_sum$diff_upper <- pi_sum$`97.5_per` - pi_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
#write.csv(pi_sum, "Data/Ela_Ham/Ela_pi_cis.csv")

## Calculate 95% for diff through time ##
Ela_pi_boot_diff <- as.data.frame(boot_Ela_AHam_pi_100bp$t)
  Ela_pi_boot_diff$Contemp_pi <- boot_Ela_CHam_pi_100bp$t
  colnames(Ela_pi_boot_diff) <- c("Hist_pi", "Contemp_pi")

Ela_pi_boot_diff$pi_diff <- Ela_pi_boot_diff$Contemp_pi - Ela_pi_boot_diff$Hist_pi
  Ela_pi_boot_diff_CIs <- quantile(Ela_pi_boot_diff$pi_diff, c(0.025, 0.975)) #-0.00048, 0.00006

## Create null distribution of change through time ##
#create vector to populate
permutation_pi <- c()
  
#for loop to randomly sample df 10 000x
for (i in 1:10000) { 
  permutation_AHam <- Ela_pi_100bp[sample(nrow(Ela_pi_100bp), 1823), ] #random subset = # Alb rows
  permutation_CHam <- Ela_pi_100bp[sample(nrow(Ela_pi_100bp), 1823), ] #random subset = # Contemp rows
  
  #calculate means
  permutation_AHam_mean <- mean(permutation_AHam$avg_pi)
  permutation_CHam_mean <- mean(permutation_CHam$avg_pi)
  
  #calculate diff
  permutation_pi[i] <- permutation_CHam_mean - permutation_AHam_mean
}
  
#write out
permutation_diff_pi <- as.data.frame(permutation_pi)
  permutation_diff_pi$iteration <- 1:10000
  
#write.csv(permutation_diff_pi, "Data/Ela_Ham/Ela_permutation_pi.csv")

## Calculate empirical p-value ##
#number of permutations with value greater than observed (really less than)
pi_perm_greater <- as.numeric(nrow(subset(Ela_pi_permutation, permutation_pi < -0.00021))) 
  (pi_perm_greater + 1)/10001 #0.0162
  
#### Visualize results ####

#ordering x-axis
pi_sum$Pop <- factor(pi_sum$Pop, levels = c("Hist", "Contemp"))

#plot of mean pi w/95% CI error bars
pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ggtitle("Ela pi 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
pi_plot

#### Null distribution plots ####

quantile(Ela_pi_permutation$permutation_pi, c(0.025, 0.975))

#plot of null pi diff distribution
pi_null_plot <- ggplot(Ela_pi_permutation, aes(x = permutation_pi)) + 
  geom_density(color = "#b9cbd8", fill = "#b9cbd8") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00048), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00006), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00021), linewidth = 6, color = "#16537e") + #"real" mean diff
  annotate("text", x = -0.0005, y = 3950, label = "I", size = 40) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab("Temporal Change in π") + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
pi_null_plot

#### Calculate change through time ####
#selective sweeps --> regions with declining pi could be swept/linked to sweep

## merge Alb & Contemp ##
pi_Alb <- as.data.frame(Ela_AHam_pi_100bp$chromosome)
  pi_Alb$window_pos_1 <- Ela_AHam_pi_100bp$window_pos_1
  pi_Alb$window_pos_2 <- Ela_AHam_pi_100bp$window_pos_2
  pi_Alb$pi <- Ela_AHam_pi_100bp$avg_pi
  colnames(pi_Alb) <- c("chromosome", "window_pos_1", "window_pos_2", "avg_pi")
pi_Contemp <- as.data.frame(Ela_CHam_pi_100bp$chromosome)
  pi_Contemp$window_pos_1 <- Ela_CHam_pi_100bp$window_pos_1
  pi_Contemp$window_pos_2 <- Ela_CHam_pi_100bp$window_pos_2
  pi_Contemp$pi <- Ela_CHam_pi_100bp$avg_pi
  colnames(pi_Contemp) <- c("chromosome", "window_pos_1", "window_pos_2", "avg_pi")

pi_merge <- merge(pi_Alb, pi_Contemp, by = c("chromosome", "window_pos_1", "window_pos_2"))
  colnames(pi_merge) <- c("chromosome", "window_pos_1", "window_pos_2", "Alb_avg_pi", "Contemp_avg_pi")
  pi_merge$delta <- pi_merge$Contemp_avg_pi - pi_merge$Alb_avg_pi

#calculate mean and median delta  
mean(na.omit(pi_merge$delta)) #-0.00021 --> losing diversity through time
median(na.omit(pi_merge$delta)) #0

#identify ones that became negative through time
#poss selective sweeps
pot_outliers <- subset(pi_merge, delta < 0)

#diff scatter plot
pi_diff_scatter_plot <- ggplot(data = pi_merge, aes(x = Alb_avg_pi, y = Contemp_avg_pi)) + 
  geom_point(size = 10, color = "#7e81a7", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  ggtitle("Ela delta pi") + labs(y = "Contemporary", x = "Albatross") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
pi_diff_scatter_plot

###########################################################################################################################

######## Aen pi estimates ########

#subsetting dataframe to windows with >100 genotyped sites
Aen_pi_100bp <- subset(Aen_pi, no_sites >= 100 & 
                         pop != "Aen-AHam-B" & pop != "Aen-Cbat-B") #left with 76 windows (38/time point)

#create dataframe for each pop
Aen_AHam_pi_100bp <- subset(Aen_pi_100bp, pop == "Aen-AHam-A")
Aen_CHam_pi_100bp <- subset(Aen_pi_100bp, pop == "Aen-Cbat-A")

#### calculate mean & median pi/pop ####

#calculate mean pi
Aen_AHam_pi_100bp_mean <- mean(Aen_AHam_pi_100bp$avg_pi) #0.00061
Aen_CHam_pi_100bp_mean <- mean(Aen_CHam_pi_100bp$avg_pi) #0.00037

#calculate median pi
Aen_AHam_pi_100bp_median <- median(Aen_AHam_pi_100bp$avg_pi) #0.00028
Aen_CHam_pi_100bp_median <- median(Aen_CHam_pi_100bp$avg_pi) #0.00000

#### pi bootstrapping ####

#bootstrap for Aen Albatross
boot_Aen_AHam_pi_100bp <- boot(data = Aen_AHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Aen_AHam_pi_100bp_95ci <- boot.ci(boot_Aen_AHam_pi_100bp, conf = 0.95, type = "norm")
  Aen_AHam_pi_100bp_95ci_normal <- Aen_AHam_pi_100bp_95ci$normal #0.00016, 0.00104

#bootstrap for Aen Contemporary
boot_Aen_CHam_pi_100bp <- boot(data = Aen_CHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Aen_CHam_pi_100bp_95ci <- boot.ci(boot_Aen_CHam_pi_100bp, conf = 0.95, type = "norm")
  Aen_CHam_pi_100bp_95ci_normal <- Aen_CHam_pi_100bp_95ci$normal #0.00005, 0.00069

## write out pi summary table ##
pi_mean <- as.data.frame(c(0.00061, 0.00037))

#pi summary table
pi_ci <- rbind(Aen_AHam_pi_100bp_95ci_normal, Aen_CHam_pi_100bp_95ci_normal) #combine df w/ci info for each pop into one dataframe
pi_sum <- cbind(pi_mean, pi_ci)
  colnames(pi_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  pi_sum$Pop <- c("Hist", "Contemp")

pi_sum$diff_lower <- pi_sum$mean - pi_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_sum$diff_upper <- pi_sum$`97.5_per` - pi_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
#write.csv(pi_sum, "Data/Aen_Ham/Aen_pi_cis.csv")

## Calculate 95% for diff through time ##
Aen_pi_boot_diff <- as.data.frame(boot_Aen_AHam_pi_100bp$t)
  Aen_pi_boot_diff$Contemp_pi <- boot_Aen_CHam_pi_100bp$t
  colnames(Aen_pi_boot_diff) <- c("Hist_pi", "Contemp_pi")

Aen_pi_boot_diff$pi_diff <- Aen_pi_boot_diff$Contemp_pi - Aen_pi_boot_diff$Hist_pi
  Aen_pi_boot_diff_CIs <- quantile(Aen_pi_boot_diff$pi_diff, c(0.025, 0.975)) #-0.00078, 0.00025
  
## Create null distribution of change through time ##
#create vector to populate
permutation_pi <- c()
  
#for loop to randomly sample df 10 000x
for (i in 1:10000) { 
  permutation_AHam <- Aen_pi_100bp[sample(nrow(Aen_pi_100bp), 38), ] #random subset = # Alb rows
  permutation_CHam <- Aen_pi_100bp[sample(nrow(Aen_pi_100bp), 38), ] #random subset = # Contemp rows
    
  #calculate means
  permutation_AHam_mean <- mean(permutation_AHam$avg_pi)
  permutation_CHam_mean <- mean(permutation_CHam$avg_pi)
    
  #calculate diff
  permutation_pi[i] <- permutation_CHam_mean - permutation_AHam_mean
}

#write out
permutation_diff_pi <- as.data.frame(permutation_pi)
  permutation_diff_pi$iteration <- 1:10000
  
#write.csv(permutation_diff_pi, "Data/Aen_Ham/Aen_permutation_pi.csv")

## Calculate empirical p-value ##
#number of permutations with value greater than observed (really less than)
pi_perm_greater <- as.numeric(nrow(subset(Aen_pi_permutation, permutation_pi < -0.00024))) 
  (pi_perm_greater + 1)/10001 #0.1118
  
#### Visualize results ####

#ordering x-axis
pi_sum$Pop <- factor(pi_sum$Pop, levels = c("Hist", "Contemp"))

#plot of mean pi w/95% CI error bars
pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean, color = Pop)) + 
  geom_point(aes(), size = 10) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ggtitle("Aen pi 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
pi_plot

#### Null distribution plots ####

quantile(Aen_pi_permutation$permutation_pi, c(0.025, 0.975))

#plot of null pi diff distribution
pi_null_plot <- ggplot(Aen_pi_permutation, aes(x = permutation_pi)) + 
  geom_density(color = "#e3ccb4", fill = "#e3ccb4") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00078), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00025), linewidth = 6, linetype = "dashed", color = "#bd8850") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00024), linewidth = 6, color = "#a25505") + #"real" mean diff
  annotate("text", x = -0.00075, y = 1900, label = "G", size = 40) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(labels = scales::label_number()) + 
  xlab("Temporal Change in π") + ylab("Density") + 
  theme_classic() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(1,1,1,1), "cm"),)
pi_null_plot

#### Calculate change through time ####
#selective sweeps --> regions with declining pi could be swept/linked to sweep

## merge Alb & Contemp ##
pi_Alb <- as.data.frame(Aen_AHam_pi_100bp$chromosome)
  pi_Alb$pi <- Aen_AHam_pi_100bp$avg_pi
  colnames(pi_Alb) <- c("chromosome", "avg_pi")
pi_Contemp <- as.data.frame(Aen_CHam_pi_100bp$chromosome)
  pi_Contemp$pi <- Aen_CHam_pi_100bp$avg_pi
  colnames(pi_Contemp) <- c("chromosome", "avg_pi")

pi_merge <- merge(pi_Alb, pi_Contemp, by = "chromosome")
colnames(pi_merge) <- c("chromosome", "Alb_avg_pi", "Contemp_avg_pi")
pi_merge$delta <- pi_merge$Contemp_avg_pi - pi_merge$Alb_avg_pi

#calculate mean and median delta  
mean(na.omit(pi_merge$delta)) #-0.00024 --> losing diversity through time
median(na.omit(pi_merge$delta)) #-0.00013

#identify ones that became negative through time
#poss selective sweeps
pot_outliers <- subset(pi_merge, delta < 0)

#diff scatter plot
pi_diff_scatter_plot <- ggplot(data = pi_merge, aes(x = Alb_avg_pi, y = Contemp_avg_pi)) + 
  geom_point(size = 10, color = "#7e81a7", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  ggtitle("Aen delta pi") + labs(y = "Contemporary", x = "Albatross") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"))
pi_diff_scatter_plot

###########################################################################################################################

######## All species comparisons ########
#designed to be run separately

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(gridExtra) #v.2.3
library(here) #v.1.0.1

#read in data
all_data <- read.csv(here("Data", "allsp_div_cis.csv"), header = TRUE)

#plot He
He_plot <- ggplot(data = all_data[all_data$metric == "He",], 
                  aes(x = factor(species, levels = c("Aen", "Gmi", "Ela")), y = mean, color = era_species, shape = era)) + 
  geom_point(size = 16, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, 
                linewidth = 4, position = position_dodge(width = 0.2)) + 
  annotate("text", x = 0.6, y = 0.145, label = "A", size = 40) + 
  annotate("text", x = 2.2, y = 0.019, label = "A. end: p > 0.001", size = 15, hjust = 0, color = "#a25505") + 
  annotate("text", x = 2.2, y = 0.012, label = "G. min: p > 0.001", size = 15, hjust = 0, color = "#1c3b0e") + 
  annotate("text", x = 2.2, y = 0.005, label = "E. lat: p > 0.001", size = 15, hjust = 0, color = "#16537e") + 
  scale_color_manual(values = c("#e3ccb4", "#8aa9be", "#afc8a4", "#a25505", "#16537e", "#1c3b0e")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("A.\nendrachtensis", "G.\nminuta", "E.\nlaterofenestra")) + 
  ylim(0, 0.15) +
  ylab(bquote(H[e]~"mean")) + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_blank(),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
He_plot

#plot Ho
Ho_plot <- ggplot(data = all_data[all_data$metric == "Ho",], 
                  aes(x = factor(species, levels = c("Aen", "Gmi", "Ela")), y = mean, color = era_species, shape = era)) + 
  geom_point(size = 16, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, 
                linewidth = 4, position = position_dodge(width = 0.2)) + 
  annotate("text", x = 0.6, y = 0.145, label = "B", size = 40) + 
  annotate("text", x = 2.2, y = 0.019, label = "A. end: p > 0.001", size = 15, hjust = 0, color = "#a25505") + 
  annotate("text", x = 2.2, y = 0.012, label = "G. min: p > 0.01", size = 15, hjust = 0, color = "#1c3b0e") + 
  annotate("text", x = 2.2, y = 0.005, label = "E. lat: p > 0.05", size = 15, hjust = 0, color = "#16537e") + 
  scale_color_manual(values = c("#e3ccb4", "#8aa9be", "#afc8a4", "#a25505", "#16537e", "#1c3b0e")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("A.\nendrachtensis", "G.\nminuta", "E.\nlaterofenestra")) + 
  ylim(0, 0.15) +
  ylab(bquote(H[o]~"mean")) + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_blank(),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ho_plot

#plot pi
Pi_plot <- ggplot(data = all_data[all_data$metric == "pi",], 
                  aes(x = factor(species, levels = c("Aen", "Gmi", "Ela")), y = mean, color = era_species, shape = era)) + 
  geom_point(size = 16, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, 
                linewidth = 4, position = position_dodge(width = 0.2)) + 
  annotate("text", x = 0.6, y = 0.004, label = "C", size = 40) + 
  annotate("text", x = 2.2, y = 0.0006, label = "A. end: p = 0.1118", size = 15, hjust = 0, color = "#a25505") + 
  annotate("text", x = 2.2, y = 0.0004, label = "G. min: p = 0.0642", size = 15, hjust = 0, color = "#1c3b0e") + 
  annotate("text", x = 2.2, y = 0.0002, label = "E. lat: p > 0.05", size = 15, hjust = 0, color = "#16537e") + 
  scale_color_manual(values = c("#e3ccb4", "#8aa9be", "#afc8a4", "#a25505", "#16537e", "#1c3b0e")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("A.\nendrachtensis", "G.\nminuta", "E.\nlaterofenestra")) + 
  ylab("π mean") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_blank(),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Pi_plot

#plot fis
Fis_plot <- ggplot(data = all_data[all_data$metric == "Fis",], 
                  aes(x = factor(species, levels = c("Aen", "Gmi", "Ela")), y = mean, color = era_species, shape = era)) + 
  geom_point(size = 16, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, 
                linewidth = 4, position = position_dodge(width = 0.2)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 4, color = "black") + 
  annotate("text", x = 0.6, y = 0.105, label = "D", size = 40) + 
  annotate("text", x = 2.2, y = 0.108, label = "A. end: p > 0.0001", size = 15, hjust = 0, color = "#a25505") + 
  annotate("text", x = 2.2, y = 0.102, label = "G. min: p = 0.2556", size = 15, hjust = 0, color = "#1c3b0e") + 
  annotate("text", x = 2.2, y = 0.096, label = "E. lat: p = 0.1009", size = 15, hjust = 0, color = "#16537e") + 
  scale_color_manual(values = c("#e3ccb4", "#8aa9be", "#afc8a4", "#a25505", "#16537e", "#1c3b0e")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("A.\nendrachtensis", "G.\nminuta", "E.\nlaterofenestra")) + 
  ylab(bquote(F[IS]~"mean")) + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_blank(),
        legend.position = "none", #not showing legend when putting all plots together
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Fis_plot

#plot all together
div_all_plot <- grid.arrange(He_plot, Ho_plot, fis_plot, pi_plot, ncol = 2)
div_all_plot