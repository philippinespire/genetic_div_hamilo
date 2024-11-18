############################################### Script to calculate pi ################################################

#From https://pixy.readthedocs.io/en/latest/plotting.html
#Using all sites VCF (so know which are monomorphic and which are missing bc of quality control)
#Calculated in 1Kb windows
#Bootstrapping across sites

#Probe targets are ~120 bp in length, so filter to windows with at least 100 genotyped sites = most of that probe region is called

#Also code to plot all of the diversity estimates together for each species

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
Gmi_pi <- read.table(here("Data/Gmi_Ham/pixy/", "pixy_pi.txt"), header = TRUE)
Ela_pi <- read.table(here("Data/Ela_Ham/pixy/", "pixy_pi.txt"), header = TRUE)

#read in permutation data
Gmi_pi_permutation <- read.csv(here("Data/Gmi_Ham", "Gmi_permutation_pi.csv"))
Ela_pi_permutation <- read.csv(here("Data/Ela_Ham", "Ela_permutation_pi.csv"))

##########################################################################################################################

######## Gmi pi estimates ########

#subsetting dataframe to windows with >100 genotyped sites
Gmi_pi_100bp <- subset(Gmi_pi, no_sites >= 100) #left with 5244 windows (2622/time point)

#create dataframe for each pop
Gmi_AHam_pi_100bp <- subset(Gmi_pi_100bp, pop == "AHam") #2622 windows
Gmi_CHam_pi_100bp <- subset(Gmi_pi_100bp, pop == "CBat") #2622 windows

#### calculate mean & median pi/pop ####

#calculate mean pi
Gmi_AHam_pi_100bp_mean <- mean(Gmi_AHam_pi_100bp$avg_pi) #0.00178
Gmi_CHam_pi_100bp_mean <- mean(Gmi_CHam_pi_100bp$avg_pi) #0.00172

#calculate median pi
Gmi_AHam_pi_100bp_median <- median(Gmi_AHam_pi_100bp$avg_pi) #0.00063
Gmi_CHam_pi_100bp_median <- median(Gmi_CHam_pi_100bp$avg_pi) #0.00073

#### pi bootstrapping ####

#bootstrap for Gmi Ham Albatross
boot_Gmi_AHam_pi_100bp <- boot(data = Gmi_AHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_AHam_pi_100bp_95ci <- boot.ci(boot_Gmi_AHam_pi_100bp, conf = 0.95, type = "norm")
  Gmi_AHam_pi_100bp_95ci_normal <- Gmi_AHam_pi_100bp_95ci$normal

#bootstrap for Gmi Ham Contemporary
boot_Gmi_CHam_pi_100bp <- boot(data = Gmi_CHam_pi_100bp$avg_pi, statistic = samp_mean, R = 1000) #1000 permutations of pi
  Gmi_CHam_pi_100bp_95ci <- boot.ci(boot_Gmi_CHam_pi_100bp, conf = 0.95, type = "norm")
  Gmi_CHam_pi_100bp_95ci_normal <- Gmi_CHam_pi_100bp_95ci$normal

## write out pi summary table ##
pi_mean <- as.data.frame(c(0.00178, 0.00172))

#pi summary table
pi_ci <- rbind(Gmi_AHam_pi_100bp_95ci_normal, Gmi_CHam_pi_100bp_95ci_normal) #combine df w/ci info for each pop into one dataframe
pi_sum <- cbind(pi_mean, pi_ci)
  colnames(pi_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  pi_sum$Pop <- c("Hist", "Contemp")

pi_sum$diff_lower <- pi_sum$mean - pi_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_sum$diff_upper <- pi_sum$`97.5_per` - pi_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
#write.csv(pi_sum, "Data/Gmi_Ham/Gmi_pi_cis.csv")

## Calculate 95% for diff through time ##
Gmi_pi_boot_diff <- as.data.frame(boot_Gmi_AHam_pi_100bp$t)
  Gmi_pi_boot_diff$Contemp_pi <- boot_Gmi_CHam_pi_100bp$t
  colnames(Gmi_pi_boot_diff) <- c("Hist_pi", "Contemp_pi")

Gmi_pi_boot_diff$pi_diff <- Gmi_pi_boot_diff$Contemp_pi - Gmi_pi_boot_diff$Hist_pi
  Gmi_pi_boot_diff_CIs <- quantile(Gmi_pi_boot_diff$pi_diff, c(0.025, 0.975)) #-0.00020, 0.00008

## Create null distribution of change through time ##
#ONLY NEED TO RUN ONCE
#create vector to populate
permutation_pi <- c()
  
#for loop to randomly sample df 10 000x
for (i in 1:10000) { 
  permutation_AHam <- Gmi_pi_100bp[sample(nrow(Gmi_pi_100bp), 2622), ] #random subset = # Alb rows
  permutation_CHam <- Gmi_pi_100bp[sample(nrow(Gmi_pi_100bp), 2622), ] #random subset = # Contemp rows
    
  #calculate means
  permutation_AHam_mean <- mean(permutation_AHam$avg_pi)
  permutation_CHam_mean <- mean(permutation_CHam$avg_pi)
    
  #calculate diff
  permutation_pi[i] <- permutation_CHam_mean - permutation_AHam_mean
}
  
#write out
permutation_diff_pi <- as.data.frame(permutation_pi)
  permutation_diff_pi$iteration <- 1:10000
  
write.csv(permutation_diff_pi, "Data/Gmi_Ham/Gmi_permutation_pi.csv")
  
## Calculate empirical p-value ##
#number of permutations with value greater than observed (really less than)
pi_perm_greater <- as.numeric(nrow(subset(permutation_diff_pi, permutation_pi < -0.00006))) 
  (pi_perm_greater + 1)/10001 #0.1166
  
#### Visualize results ####

#ordering x-axis
pi_sum$Pop <- factor(pi_sum$Pop, levels = c("Hist", "Contemp"))

#plot of mean pi w/95% CI error bars
pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#1c3b0e","#afc8a4")) +
  scale_shape_manual(values = c(19, 15)) +
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

quantile(permutation_diff_pi$permutation_pi, c(0.025, 0.975))

#plot of null pi diff distribution
pi_null_plot <- ggplot(permutation_diff_pi, aes(x = permutation_pi)) + 
  geom_density(color = "#bac4b6", fill = "#bac4b6") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.0001), linewidth = 6, linetype = "dashed", color = "#607556") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.0001), linewidth = 6, linetype = "dashed", color = "#607556") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00006), linewidth = 6, color = "#1c3b0e") + #"real" mean diff
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
mean(na.omit(pi_merge$delta)) #-0.00006 --> losing diversity through time
median(na.omit(pi_merge$delta)) #0 
  
#identify ones that became negative through time (lost diversity)
#poss selective sweeps
pot_outliers <- subset(pi_merge, delta < 0) #937 windows

#diff scatter plot
pi_diff_scatter_plot <- ggplot(data = pi_merge, aes(x = Alb_avg_pi, y = Contemp_avg_pi)) + 
  geom_point(size = 14, color = "#1c3b0e", alpha = 0.5) + 
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
Ela_pi_100bp <- subset(Ela_pi, no_sites >= 100) #left with 3948 windows (1974/time point)

#create dataframe for each pop
Ela_AHam_pi_100bp <- subset(Ela_pi_100bp, pop == "AHam") #1974 windows
Ela_CHam_pi_100bp <- subset(Ela_pi_100bp, pop == "CNas") #1974 windows

#### calculate mean & median pi/pop ####

#calculate mean pi
Ela_AHam_pi_100bp_mean <- mean(Ela_AHam_pi_100bp$avg_pi) #0.00408
Ela_CHam_pi_100bp_mean <- mean(Ela_CHam_pi_100bp$avg_pi) #0.00387

#calculate median pi
Ela_AHam_pi_100bp_median <- median(Ela_AHam_pi_100bp$avg_pi) #0.00319
Ela_CHam_pi_100bp_median <- median(Ela_CHam_pi_100bp$avg_pi) #0.00299

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
pi_mean <- as.data.frame(c(0.00408, 0.00387))

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
  Ela_pi_boot_diff_CIs <- quantile(Ela_pi_boot_diff$pi_diff, c(0.025, 0.975)) #-0.00046, 0.00001

## Create null distribution of change through time ##
#ONLY NEED TO RUN ONCE
#create vector to populate
permutation_pi <- c()
  
#for loop to randomly sample df 10 000x
for (i in 1:10000) { 
  permutation_AHam <- Ela_pi_100bp[sample(nrow(Ela_pi_100bp), 1974), ] #random subset = # Alb rows
  permutation_CHam <- Ela_pi_100bp[sample(nrow(Ela_pi_100bp), 1974), ] #random subset = # Contemp rows
  
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
pi_perm_greater <- as.numeric(nrow(subset(permutation_diff_pi, permutation_pi < -0.00021))) 
  (pi_perm_greater + 1)/10001 #0.007
  
#### Visualize results ####

#ordering x-axis
pi_sum$Pop <- factor(pi_sum$Pop, levels = c("Hist", "Contemp"))

#plot of mean pi w/95% CI error bars
pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  scale_shape_manual(values = c(19, 15)) +
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

quantile(permutation_diff_pi$permutation_pi, c(0.025, 0.975))

#plot of null pi diff distribution
pi_null_plot <- ggplot(permutation_diff_pi, aes(x = permutation_pi)) + 
  geom_density(color = "#b9cbd8", fill = "#b9cbd8") + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black") +
  geom_vline(aes(xintercept = -0.00017), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.25 CI lower
  geom_vline(aes(xintercept = 0.00017), linewidth = 6, linetype = "dashed", color = "#5b86a4") + #0.95 CI upper
  geom_vline(aes(xintercept = -0.00021), linewidth = 6, color = "#16537e") + #"real" mean diff
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
median(na.omit(pi_merge$delta)) #-0.000002

#identify ones that became negative through time (lost diversity)
#poss selective sweeps
pot_outliers <- subset(pi_merge, delta < 0) #990 windows

#diff scatter plot
pi_diff_scatter_plot <- ggplot(data = pi_merge, aes(x = Alb_avg_pi, y = Contemp_avg_pi)) + 
  geom_point(size = 14, color = "#16537e", alpha = 0.5) + 
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

######## All species comparisons ########
#designed to be run separately

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(gridExtra) #v.2.3
library(here) #v.1.0.1

#read in data
all_data <- read.csv(here("Data", "2sp_div_cis.csv"), header = TRUE)

#plot heterozygosity
Het_plot <- ggplot(data = all_data[all_data$metric == "He" | all_data$metric == "Ho", ],
                   aes(x = factor(era, levels = c("Hist", "Contemp")), y = mean, 
                       group = species_metric, color = species, shape = metric, linetype = metric)) + 
  geom_line(linewidth = 4) + 
  geom_point(size = 16) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0, linewidth = 4) + 
  scale_color_manual(values = c("#16537e", "#afc8a4")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("Hist.", "Contemp.")) + 
  ylab("Mean Heterozygosity") + 
  ylim(0.092, 0.125) +
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
Het_plot

#plot pi
Pi_plot <- ggplot(data = all_data[all_data$metric == "pi", ],
                   aes(x = factor(era, levels = c("Hist", "Contemp")), y = mean, 
                       group = species_metric, color = species, shape = metric, linetype = metric)) + 
  geom_line(linewidth = 4) + 
  geom_point(size = 16) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0, linewidth = 4) + 
  scale_color_manual(values = c("#16537e", "#afc8a4")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("Hist.", "Contemp.")) + 
  ylab("Mean π") + 
  theme_bw() + 
  ylim(0.0010, 0.0045) +
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
Fis_plot <- ggplot(data = all_data[all_data$metric == "Fis", ],
                  aes(x = factor(era, levels = c("Hist", "Contemp")), y = mean, 
                      group = species_metric, color = species, shape = metric, linetype = metric)) + 
  geom_line(linewidth = 4) + 
  geom_point(size = 16) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0, linewidth = 4) + 
  scale_color_manual(values = c("#16537e", "#afc8a4")) + 
  scale_shape_manual(values = c(19, 15)) +
  scale_x_discrete(labels = c("Hist.", "Contemp.")) + 
  ylab(bquote("Mean"~F[IS])) + 
  theme_bw() + 
  #ylim(0.0010, 0.0045) +
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