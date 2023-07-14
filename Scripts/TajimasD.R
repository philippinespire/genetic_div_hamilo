###################################### Script for Tajima D  #######################################################

#Created for transcriptome project
#Uses Tajima's D estimates from VCFtools --> calculated in 10kb windows (no sliding steps)

#Calculates change through time
#plots distribution shift through time

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(plotrix) #v.3.8.2
library(here) #v.1.0.1
library(ggridges) #v.0.5.4
library(purrr) #v.1.0.1

#write function to calculate mean
samp_mean <- function(x, i) {
  mean(x[i])
}

#read in data
Gmi_Alb_TD <- read.table(here("Data/Gmi_Ham", "Gmi_Apop_Ham_Alb_nohighhet.Tajima.D"), header = TRUE)
  Gmi_Contemp_TD <- read.table(here("Data/Gmi_Ham", "Gmi_Apop_Ham_Contemp_nohighhet.Tajima.D"), header = TRUE)
  Gmi_All_TD <- read.table(here("Data/Gmi_Ham", "Gmi_Apop_Ham_All_nohighhet.Tajima.D"), header = TRUE)
Ela_Alb_TD <- read.table(here("Data/Ela_Ham", "Ela_Alb.Tajima.D"), header = TRUE) #no highhet Alb individs
  Ela_Contemp_TD <- read.table(here("Data/Ela_Ham", "Ela_Contemp_nohighhet.Tajima.D"), header = TRUE)
  Ela_All_TD <- read.table(here("Data/Ela_Ham", "Ela_All_nohighhet.Tajima.D"), header = TRUE)
Aen_Alb_TD <- read.table(here("Data/Aen_Ham", "Aen_A_Alb.Tajima.D"), header = TRUE)
  Aen_Contemp_TD <- read.table(here("Data/Aen_Ham", "Aen_A_Contemp.Tajima.D"), header = TRUE)
  Aen_All_TD <- read.table(here("Data/Aen_Ham", "Aen_A_All.Tajima.D"), header = TRUE)

#clean data
Gmi_Alb_TD <- Gmi_Alb_TD[which(Gmi_Alb_TD$TajimaD != 'NaN'),]
  Gmi_Contemp_TD <- Gmi_Contemp_TD[which(Gmi_Contemp_TD$TajimaD != "NaN"),]
  Gmi_All_TD <- Gmi_All_TD[which(Gmi_All_TD$TajimaD != 'NaN'),]
Ela_Alb_TD <- Ela_Alb_TD[which(Ela_Alb_TD$TajimaD != 'NaN'),]
  Ela_Contemp_TD <- Ela_Contemp_TD[which(Ela_Contemp_TD$TajimaD != 'NaN'),]
  Ela_All_TD <- Ela_All_TD[which(Ela_All_TD$TajimaD != 'NaN'),] #bc cryptic species .. so use this as baseline for # "real" sequences w/SNPs to expect
Aen_Alb_TD <- Aen_Alb_TD[which(Aen_Alb_TD$TajimaD != 'NaN'),]
  Aen_Contemp_TD <- Aen_Contemp_TD[which(Aen_Contemp_TD$TajimaD != 'NaN'),]
  Aen_All_TD <- Aen_All_TD[which(Aen_All_TD$TajimaD != 'NaN'),] #bc some windows with no SNPs at all

#coercing Tajima D to numeric
Gmi_Alb_TD$TajimaD <- as.numeric(Gmi_Alb_TD$TajimaD)
  Gmi_Contemp_TD$TajimaD <- as.numeric(Gmi_Contemp_TD$TajimaD)
  Gmi_All_TD$TajimaD <- as.numeric(Gmi_All_TD$TajimaD)
Ela_Alb_TD$TajimaD <- as.numeric(Ela_Alb_TD$TajimaD)
  Ela_Contemp_TD$TajimaD <- as.numeric(Ela_Contemp_TD$TajimaD)
  Ela_All_TD$TajimaD <- as.numeric(Ela_All_TD$TajimaD)
Aen_Alb_TD$TajimaD <- as.numeric(Aen_Alb_TD$TajimaD)
  Aen_Contemp_TD$TajimaD <- as.numeric(Aen_Contemp_TD$TajimaD)
  Aen_All_TD$TajimaD <- as.numeric(Aen_All_TD$TajimaD)

#################################################################################################################################################

######## Gmi Tajima's D estimates ########
  
## Adding monomorphic sequences back in ##
#since SNPs called with cryptic species -- species-level diff could be called as SNP, when w/in species is not
#means have some monomorphic sequences that aren't "real" -- have no SNP in this species, shouldn't be in pool
#also, if SNPs disappear across time points will be missing sequences through time
#creates diff num of sequences across pops, want these to be all the same
  
##### set-up dfs for comparison ####
  
##create unique seq for comparison ##
Gmi_Alb_TD$uniqseq <- paste(Gmi_Alb_TD$CHROM, "-", Gmi_Alb_TD$BIN_START)
  Gmi_Contemp_TD$uniqseq <- paste(Gmi_Contemp_TD$CHROM, "-", Gmi_Contemp_TD$BIN_START)
  Gmi_All_TD$uniqseq <- paste(Gmi_All_TD$CHROM, "-", Gmi_All_TD$BIN_START)
  
## Alb pop calculations ##
#number of monomorphic sequences not included
nrow(Gmi_All_TD) - nrow(Gmi_Alb_TD) #470 --> need to add this many instances of 0 into dataframe
  
#create dataframe of sequences to add back
contigs_needed_Alb <- as.data.frame(setdiff(Gmi_All_TD$uniqseq, Gmi_Alb_TD$uniqseq))
  colnames(contigs_needed_Alb) <- c("uniqseq")
  
contigs_needed_Alb <- separate(contigs_needed_Alb, uniqseq, 
                               c("CHROM", "BIN_START"), sep = "-", remove = FALSE) #separate uniqseq out to original columns
  contigs_needed_Alb$N_SNPS <- c(rep(0, times = 470)) #add N_SNPS column
  contigs_needed_Alb$TajimaD <- c(rep("NA", times = 470)) #add TajimaD column
  
#create full Tajima D df
full_AlbTD <- rbind(Gmi_Alb_TD, contigs_needed_Alb)
  full_AlbTD <- full_AlbTD[order(full_AlbTD$uniqseq), ] #reorder by contig
  nrow(full_AlbTD) #check equals 2548 (# total sequences)
full_AlbTD$Pop <- c(rep("Albatross", times = 2548))
full_AlbTD$NUM <- c(1:2548)
  
## Contemp pop calculations ##
#number of monomorphic sequences not included
nrow(Gmi_All_TD) - nrow(Gmi_Contemp_TD) #85 --> need to add this many instances of 0 into dataframe
  
#create dataframe of sequences to add back
contigs_needed_Contemp <- as.data.frame(setdiff(Gmi_All_TD$uniqseq, Gmi_Contemp_TD$uniqseq))
  colnames(contigs_needed_Contemp) <- c("uniqseq")
  
contigs_needed_Contemp <- separate(contigs_needed_Contemp, uniqseq, 
                                   c("CHROM", "BIN_START"), sep = "-", remove = FALSE)
  contigs_needed_Contemp$N_SNPS <- c(rep(0, times = 85)) #add N_SNPS column
  contigs_needed_Contemp$TajimaD <- c(rep("NA", times = 85)) #add TajimaD column
  
#create full Tajima D df
full_ContempTD <- rbind(Gmi_Contemp_TD, contigs_needed_Contemp)
  full_ContempTD <- full_ContempTD[order(full_ContempTD$uniqseq), ]
  nrow(full_ContempTD) #check equals 2548 (# total sequences)
full_ContempTD$Pop <- c(rep("Contemporary", times = 2548))
full_ContempTD$NUM <- c(1:2548)
  
## combine pop data ##
#set-up All DF
full_AllTD <- Gmi_All_TD[order(Gmi_All_TD$uniqseq), ]
  nrow(full_AllTD) #check equals 2548 (# total sequences)
full_AllTD$Pop <- c(rep("All", times = 2548))
full_AllTD$NUM <- c(1:2548)
  
#combine TD for all seqs
TD_only_all <- rbind(full_AlbTD, full_ContempTD, full_AllTD)
  TD_only_all$TajimaD <- as.numeric(TD_only_all$TajimaD)
  
#write out
#write.csv(TD_only_all, "Data/Gmi_Ham/Gmi_popA_Ham_TajimasD_combined_full.csv")
  
#### Calculate means by pop ####
  
#calculate means
Alb_mean_TD_all <- mean(Gmi_Alb_TD$TajimaD) #-0.53290
Contemp_mean_TD_all <- mean(Gmi_Contemp_TD$TajimaD) #-0.58844
All_mean_TD_all <- mean(Gmi_All_TD$TajimaD) #-0.56168

Contemp_mean_TD_all - Alb_mean_TD_all #-0.05354
  
#### TD bootstrapping ####
  
#historical TD
boot_Hist_Gmi_TD <- boot(data = Gmi_Alb_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Hist_Gmi_TD_95ci <- boot.ci(boot_Hist_Gmi_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Hist_Gmi_TD_95ci_normal <- Hist_Gmi_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
#contemporary TD
boot_Contemp_Gmi_TD <- boot(data = Gmi_Contemp_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Contemp_Gmi_TD_95ci <- boot.ci(boot_Contemp_Gmi_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Contemp_Gmi_TD_95ci_normal <- Contemp_Gmi_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
#all TD
boot_All_Gmi_TD <- boot(data = Gmi_All_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  All_Gmi_TD_95ci <- boot.ci(boot_All_Gmi_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  All_Gmi_TD_95ci_normal <- All_Gmi_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD

## write out TD summary tables ##
TD_mean <- as.data.frame(c(-0.53290, -0.58844, -0.56168))
  
#TD summary table
TD_ci <- rbind(Hist_Gmi_TD_95ci_normal, Contemp_Gmi_TD_95ci_normal, 
               All_Gmi_TD_95ci_normal) #combine df w/ci info for each into one dataframe
TD_sum <- cbind(TD_mean, TD_ci)
  colnames(TD_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  TD_sum$Pop <- c("Hist", "Contemp", "All")
  
TD_sum$diff_lower <- TD_sum$mean - TD_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
TD_sum$diff_upper <- TD_sum$`97.5_per` - TD_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization
  
#write out data
#write.csv(TD_sum, "Data/Gmi_Ham/Gmi_popA_Ham_TD_cis.csv")

## Calculate 95% for diff through time ##
Gmi_TD_boot_diff <- as.data.frame(boot_Hist_Gmi_TD$t)
  Gmi_TD_boot_diff$Contemp_TD <- boot_Contemp_Gmi_TD$t
  colnames(Gmi_TD_boot_diff) <- c("Hist_TD", "Contemp_TD")

Gmi_TD_boot_diff$TD_diff <- Gmi_TD_boot_diff$Contemp_TD - Gmi_TD_boot_diff$Hist_TD
  Gmi_TD_boot_diff_CIs <- quantile(Gmi_TD_boot_diff$TD_diff, c(0.025, 0.975)) #-0.09793, -0.00937

#### Mann-Whitney U test ####
  
#subset data
#need Alb & Contemp TD in two separate vectors to compare
Alb_df <- TD_only_all[which(TD_only_all$Pop == "Albatross"), ]
  Alb_df_noNA <- Alb_df[which(Alb_df$TajimaD!='NA'),]
  Alb_TD <- Alb_df_noNA$TajimaD
Contemp_df <- TD_only_all[which(TD_only_all$Pop == "Contemporary"), ]
  Contemp_df_noNA <- Contemp_df[which(Contemp_df$TajimaD!='NA'),]
  Contemp_TD <- Contemp_df_noNA$TajimaD
  
#run Mann-Whitney U-test
TD_MUtest <- wilcox.test(Alb_TD, Contemp_TD) #W = 2586434, p = 0.534 #distributions equal
  
#### Visualize data ####
  
#read in data
TD_only_all <- read.csv(here("Data/Gmi_Ham", "Gmi_popA_Ham_TajimasD_combined_full.csv"), header = TRUE)
  TD_only_noall <- subset(TD_only_all, Pop != "All")
  
## Scatter plot ##
TD_scatter_plot <- ggplot(data = TD_only_noall, aes(x = NUM, y = TajimaD, color = Pop)) + 
  geom_point(size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 2, linetype = "dashed") + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ggtitle("Gmi Tajima's D") + labs(y = "Tajima's D", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
TD_scatter_plot
  
## Density plots ##
#ordering x-axis
TD_only_noall$Pop[TD_only_noall$Pop == "Albatross"] <- "Historical"
  TD_only_noall$Pop <- factor(TD_only_noall$Pop, levels = c("Historical", "Contemporary")) #ordering X-axis
  
#first density plot to get coordinates
TD_density_plot_1 <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, aes(x = TajimaD, y = Pop), scale = 1.2) 
  
#pull out density plot information
ingredients <- ggplot_build(TD_density_plot_1) %>% purrr::pluck("data", 1)
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()
  
#add mean coordinates
density_lines$mean <- c(-0.53290, -0.58844) #to get lines to plot at mean instead of elsewhere
  density_lines$lower_ci <- c(-0.57049, -0.61906)
  density_lines$upper_ci <- c(-0.49971, -0.55753)
density_lines$pop <- c("Historical", "Contemporary") #for coloring purposes
  
#density
TD_density_plot <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, 
                       aes(x = TajimaD, y = Pop, color = Pop, fill = Pop), 
                       alpha = 0.5, scale = 1.2) + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black", linetype = "dotted") + 
  geom_segment(data = density_lines, 
               aes(x = lower_ci, y = ymin, xend = lower_ci, yend = ymin+density*scale*iscale), 
               color = "#607556", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = upper_ci, y = ymin, xend = upper_ci, yend = ymin+density*scale*iscale), 
               color = "#607556", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = mean, y = ymin, xend = mean, yend = ymin+density*scale*iscale), 
               color = "black", linewidth = 4) +
  annotate("text", x = 3, y = 3, label = "B", size = 40) + 
  scale_color_manual(values = c("#76896e", "#bac4b6")) +
  scale_fill_manual(values = c("#76896e", "#bac4b6")) + 
  scale_y_discrete(expand = expansion(mult = c(0.01, .7))) + #bring y axis down x
  labs(x = "Tajima's D") +  
  theme_ridges() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black"), 
        axis.text.x = element_text(size = 55, color = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 55, color = "black", vjust = -0.5, hjust = 0.5),
        legend.position = "none", 
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))
TD_density_plot

#### Identify change through time ####

#calculate diff through time
TD_diff <- as.data.frame(cbind(Alb_df$uniqseq, Alb_df$TajimaD, Contemp_df$TajimaD))
  colnames(TD_diff) <- c("uniqseq", "Alb_TD", "Contemp_TD")
  TD_diff$Alb_TD <- as.numeric(as.character(TD_diff$Alb_TD)) #bc was factor
  TD_diff$Contemp_TD <- as.numeric(as.character(TD_diff$Contemp_TD))
  TD_diff$delta <- TD_diff$Contemp_TD - TD_diff$Alb_TD

#calculate mean and median delta  
mean(na.omit(TD_diff$delta)) #0.00984 --> individual sites getting slightly more positive through time
median(na.omit(TD_diff$delta)) #0.04522

#identify ones that became negative through time
#poss selective sweeps
TD_pot_selection <- subset(TD_diff, delta < 0)

#diff scatter plot
TD_diff_scatter_plot <- ggplot(data = TD_diff, aes(x = Alb_TD, y = Contemp_TD)) + 
  geom_point(size = 16, color = "#1c3b0e", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  annotate("text", x = -1.8, y = 2.8, label = "B", size = 40) + 
  labs(y = "Contemporary", x = "Albatross") + 
  ylim(-2, 3) + xlim(-2, 3) +
  theme_bw() + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
TD_diff_scatter_plot
  
#################################################################################################################################################
   
######## Ela Tajima's D estimates ########
  
##### set-up dfs for comparison ####
  
##create unique seq for comparison ##
Ela_Alb_TD$uniqseq <- paste(Ela_Alb_TD$CHROM, "-", Ela_Alb_TD$BIN_START)
  Ela_Contemp_TD$uniqseq <- paste(Ela_Contemp_TD$CHROM, "-", Ela_Contemp_TD$BIN_START)
  Ela_All_TD$uniqseq <- paste(Ela_All_TD$CHROM, "-", Ela_All_TD$BIN_START)
  
## Alb pop calculations ##
#number of monomorphic sequences not included
nrow(Ela_All_TD) - nrow(Ela_Alb_TD) #257 --> need to add this many instances of 0 into dataframe
  
#create dataframe of sequences to add back
contigs_needed_Alb <- as.data.frame(setdiff(Ela_All_TD$uniqseq, Ela_Alb_TD$uniqseq))
  colnames(contigs_needed_Alb) <- c("uniqseq")
  
contigs_needed_Alb <- separate(contigs_needed_Alb, uniqseq, 
                               c("CHROM", "BIN_START"), sep = "-", remove = FALSE) #separate uniqseq out to original columns
  contigs_needed_Alb$N_SNPS <- c(rep(0, times = 257)) #add N_SNPS column
  contigs_needed_Alb$TajimaD <- c(rep("NA", times = 257)) #add TajimaD column
  
#create full Tajima D df
full_AlbTD <- rbind(Ela_Alb_TD, contigs_needed_Alb)
  full_AlbTD <- full_AlbTD[order(full_AlbTD$uniqseq), ] #reorder by contig
  nrow(full_AlbTD) #check equals 3921 (# total sequences)
full_AlbTD$Pop <- c(rep("Albatross", times = 3921))
full_AlbTD$NUM <- c(1:3921)
  
## Contemp pop calculations ##
#number of monomorphic sequences not included
nrow(Ela_All_TD) - nrow(Ela_Contemp_TD) #66 --> need to add this many instances of 0 into dataframe
  
#create dataframe of sequences to add back
contigs_needed_Contemp <- as.data.frame(setdiff(Ela_All_TD$uniqseq, Ela_Contemp_TD$uniqseq))
  colnames(contigs_needed_Contemp) <- c("uniqseq")
  
contigs_needed_Contemp <- separate(contigs_needed_Contemp, uniqseq, 
                                   c("CHROM", "BIN_START"), sep = "-", remove = FALSE)
  contigs_needed_Contemp$N_SNPS <- c(rep(0, times = 66)) #add N_SNPS column
  contigs_needed_Contemp$TajimaD <- c(rep("NA", times = 66)) #add TajimaD column
  
#create full Tajima D df
full_ContempTD <- rbind(Ela_Contemp_TD, contigs_needed_Contemp)
  full_ContempTD <- full_ContempTD[order(full_ContempTD$uniqseq), ]
  nrow(full_ContempTD) #check equals 3921 (# total sequences)
full_ContempTD$Pop <- c(rep("Contemporary", times = 3921))
full_ContempTD$NUM <- c(1:3921)
  
## combine pop data ##
#set-up All DF
full_AllTD <- Ela_All_TD[order(Ela_All_TD$uniqseq), ]
  nrow(full_AllTD) #check equals 3921 (# total sequences)
full_AllTD$Pop <- c(rep("All", times = 3921))
full_AllTD$NUM <- c(1:3921)
  
#combine TD for all seqs
TD_only_all <- rbind(full_AlbTD, full_ContempTD, full_AllTD)
  TD_only_all$TajimaD <- as.numeric(TD_only_all$TajimaD)
  
#write out
#write.csv(TD_only_all, "Data/Ela_Ham/Ela_TajimasD_combined_full.csv")
  
#### Calculate means by pop ####
  
#calculate means
Alb_mean_TD_all <- mean(Ela_Alb_TD$TajimaD) #-0.36573
Contemp_mean_TD_all <- mean(Ela_Contemp_TD$TajimaD) #-0.30136
All_mean_TD_all <- mean(Ela_All_TD$TajimaD) #-0.32220

Contemp_mean_TD_all - Alb_mean_TD_all #0.06436
  
#### TD bootstrapping ####
  
#historical TD
boot_Hist_Ela_TD <- boot(data = Ela_Alb_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Hist_Ela_TD_95ci <- boot.ci(boot_Hist_Ela_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Hist_Ela_TD_95ci_normal <- Hist_Ela_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
#contemporary TD
boot_Contemp_Ela_TD <- boot(data = Ela_Contemp_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Contemp_Ela_TD_95ci <- boot.ci(boot_Contemp_Ela_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Contemp_Ela_TD_95ci_normal <- Contemp_Ela_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
#all TD
boot_All_Ela_TD <- boot(data = Ela_All_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  All_Ela_TD_95ci <- boot.ci(boot_All_Ela_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  All_Ela_TD_95ci_normal <- All_Ela_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
## write out TD summary tables ##
TD_mean <- as.data.frame(c(-0.36573, -0.30136, -0.32220))
  
#TD summary table
TD_ci <- rbind(Hist_Ela_TD_95ci_normal, Contemp_Ela_TD_95ci_normal, 
               All_Ela_TD_95ci_normal) #combine df w/ci info for each into one dataframe
TD_sum <- cbind(TD_mean, TD_ci)
  colnames(TD_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  TD_sum$Pop <- c("Hist", "Contemp", "All")
  
TD_sum$diff_lower <- TD_sum$mean - TD_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
TD_sum$diff_upper <- TD_sum$`97.5_per` - TD_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization
  
#write out data
#write.csv(TD_sum, "Data/Ela_Ham/Ela_TD_cis.csv")
  
## Calculate 95% for diff through time ##
Ela_TD_boot_diff <- as.data.frame(boot_Hist_Ela_TD$t)
  Ela_TD_boot_diff$Contemp_TD <- boot_Contemp_Ela_TD$t
  colnames(Ela_TD_boot_diff) <- c("Hist_TD", "Contemp_TD")

Ela_TD_boot_diff$TD_diff <- Ela_TD_boot_diff$Contemp_TD - Ela_TD_boot_diff$Hist_TD
  Ela_TD_boot_diff_CIs <- quantile(Ela_TD_boot_diff$TD_diff, c(0.025, 0.975)) #0.02293, 0.10735

#### Mann-Whitney U test ####
  
#subset data
#need Alb & Contemp TD in two separate vectors to compare
Alb_df <- TD_only_all[which(TD_only_all$Pop == "Albatross"), ]
  Alb_df_noNA <- Alb_df[which(Alb_df$TajimaD!='NA'),]
  Alb_TD <- Alb_df_noNA$TajimaD
Contemp_df <- TD_only_all[which(TD_only_all$Pop == "Contemporary"), ]
  Contemp_df_noNA <- Contemp_df[which(Contemp_df$TajimaD!='NA'),]
  Contemp_TD <- Contemp_df_noNA$TajimaD
  
#run Mann-Whitney U-test
TD_MUtest <- wilcox.test(Alb_TD, Contemp_TD) #W = 6881087, p = 0.054 #distributions not equal
  
#### Visualize data ####
  
#read in data
TD_only_all <- read.csv(here("Data/Ela_Ham", "Ela_TajimasD_combined_full.csv"), header = TRUE)
  TD_only_noall <- subset(TD_only_all, Pop != "All")
  
## Scatter plot ##
TD_scatter_plot <- ggplot(data = TD_only_noall, aes(x = NUM, y = TajimaD, color = Pop)) + 
  geom_point(size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 2, linetype = "dashed") + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  ggtitle("Ela Tajima's D") + labs(y = "Tajima's D", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
TD_scatter_plot
  
## Density plots ##
#ordering x-axis
TD_only_noall$Pop[TD_only_noall$Pop == "Albatross"] <- "Historical"
  TD_only_noall$Pop <- factor(TD_only_noall$Pop, levels = c("Historical", "Contemporary")) #ordering X-axis

#first density plot to get coordinates
TD_density_plot_1 <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, aes(x = TajimaD, y = Pop), scale = 1.2) 

#pull out density plot information
ingredients <- ggplot_build(TD_density_plot_1) %>% purrr::pluck("data", 1)
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()

#add mean coordinates
density_lines$mean <- c(-0.36573, -0.30136) #to get lines to plot at mean instead of elsewhere
density_lines$lower_ci <- c(-0.3959, -0.33372)
density_lines$upper_ci <- c(-0.33566, -0.26905)
density_lines$pop <- c("Historical", "Contemporary") #for coloring purposes

#density
TD_density_plot <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, 
                       aes(x = TajimaD, y = Pop, color = Pop, fill = Pop), 
                       alpha = 0.5, scale = 1.2) + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black", linetype = "dotted") + 
  geom_segment(data = density_lines, 
               aes(x = lower_ci, y = ymin, xend = lower_ci, yend = ymin+density*scale*iscale), 
               color = "#607685", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = upper_ci, y = ymin, xend = upper_ci, yend = ymin+density*scale*iscale), 
               color = "#607685", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = mean, y = ymin, xend = mean, yend = ymin+density*scale*iscale), 
               color = "black", linewidth = 4) +
  annotate("text", x = 3, y = 3, label = "C", size = 40) + 
  scale_color_manual(values = c("#7c98ab", "#b8cbd8")) +
  scale_fill_manual(values = c("#7c98ab", "#b8cbd8")) + 
  scale_y_discrete(expand = expansion(mult = c(0.01, .7))) + #bring y axis down x
  labs(x = "Tajima's D") +  
  theme_ridges() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black"), 
        axis.text.x = element_text(size = 55, color = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 55, color = "black", vjust = -0.5, hjust = 0.5),
        legend.position = "none", 
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))
TD_density_plot

#### Identify change through time ####

#calculate diff through time
TD_diff <- as.data.frame(cbind(Alb_df$uniqseq, Alb_df$TajimaD, Contemp_df$TajimaD))
  colnames(TD_diff) <- c("uniqseq", "Alb_TD", "Contemp_TD")
  TD_diff$Alb_TD <- as.numeric(as.character(TD_diff$Alb_TD)) #bc was factor
  TD_diff$Contemp_TD <- as.numeric(as.character(TD_diff$Contemp_TD))
TD_diff$delta <- TD_diff$Contemp_TD - TD_diff$Alb_TD

#calculate mean and median delta  
mean(na.omit(TD_diff$delta)) #0.09308 --> individual sites getting slightly more positive through time
median(na.omit(TD_diff$delta)) #0.12961 

#identify ones that became negative through time
#poss selective sweeps
TD_pot_selection <- subset(TD_diff, delta < 0)

#diff scatter plot
TD_diff_scatter_plot <- ggplot(data = TD_diff, aes(x = Alb_TD, y = Contemp_TD)) + 
  geom_point(size = 16, color = "#16537e", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  annotate("text", x = -1.8, y = 2.8, label = "C", size = 40) + 
  labs(y = "Contemporary", x = "Albatross") + 
  ylim(-2, 3) + xlim(-2, 3) +
  theme_bw() + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
TD_diff_scatter_plot

#################################################################################################################################################
  
######## Aen Tajima's D estimates ########

##### set-up dfs for comparison ####

##create unique seq for comparison ##
Aen_Alb_TD$uniqseq <- paste(Aen_Alb_TD$CHROM, "-", Aen_Alb_TD$BIN_START)
  Aen_Contemp_TD$uniqseq <- paste(Aen_Contemp_TD$CHROM, "-", Aen_Contemp_TD$BIN_START)
  Aen_All_TD$uniqseq <- paste(Aen_All_TD$CHROM, "-", Aen_All_TD$BIN_START)

## Alb pop calculations ##
#number of monomorphic sequences not included
nrow(Aen_All_TD) - nrow(Aen_Alb_TD) #5 --> need to add this many instances of 0 into dataframe

#create dataframe of sequences to add back
contigs_needed_Alb <- as.data.frame(setdiff(Aen_All_TD$uniqseq, Aen_Alb_TD$uniqseq))
  colnames(contigs_needed_Alb) <- c("uniqseq")

contigs_needed_Alb <- separate(contigs_needed_Alb, uniqseq, 
                               c("CHROM", "BIN_START"), sep = "-", remove = FALSE) #separate uniqseq out to original columns
  contigs_needed_Alb$N_SNPS <- c(rep(0, times = 5)) #add N_SNPS column
  contigs_needed_Alb$TajimaD <- c(rep("NA", times = 5)) #add TajimaD column

#create full Tajima D df
full_AlbTD <- rbind(Aen_Alb_TD, contigs_needed_Alb)
  full_AlbTD <- full_AlbTD[order(full_AlbTD$uniqseq), ] #reorder by contig
  nrow(full_AlbTD) #check equals 391 (# total sequences)
full_AlbTD$Pop <- c(rep("Albatross", times = 391))
full_AlbTD$NUM <- c(1:391)

## Contemp pop calculations ##
#number of monomorphic sequences not included
nrow(Aen_All_TD) - nrow(Aen_Contemp_TD) #244 --> need to add this many instances of 0 into dataframe

#create dataframe of sequences to add back
contigs_needed_Contemp <- as.data.frame(setdiff(Aen_All_TD$uniqseq, Aen_Contemp_TD$uniqseq))
  colnames(contigs_needed_Contemp) <- c("uniqseq")

contigs_needed_Contemp <- separate(contigs_needed_Contemp, uniqseq, 
                                   c("CHROM", "BIN_START"), sep = "-", remove = FALSE)
  contigs_needed_Contemp$N_SNPS <- c(rep(0, times = 244)) #add N_SNPS column
  contigs_needed_Contemp$TajimaD <- c(rep("NA", times = 244)) #add TajimaD column

#create full Tajima D df
full_ContempTD <- rbind(Aen_Contemp_TD, contigs_needed_Contemp)
  full_ContempTD <- full_ContempTD[order(full_ContempTD$uniqseq), ]
  nrow(full_ContempTD) #check equals 391 (# total sequences)
full_ContempTD$Pop <- c(rep("Contemporary", times = 391))
full_ContempTD$NUM <- c(1:391)

## combine pop data ##
#set-up All DF
full_AllTD <- Aen_All_TD[order(Aen_All_TD$uniqseq), ]
  nrow(full_AllTD) #check equals 391 (# total sequences)
full_AllTD$Pop <- c(rep("All", times = 391))
full_AllTD$NUM <- c(1:391)

#combine TD for all seqs
TD_only_all <- rbind(full_AlbTD, full_ContempTD, full_AllTD)
  TD_only_all$TajimaD <- as.numeric(TD_only_all$TajimaD)

#write out
#write.csv(TD_only_all, "Data/Aen_Ham/Ham_A_TajimasD_combined_full.csv")

#### Calculate means by pop ####

#calculate means
Alb_mean_TD_all <- mean(Aen_Alb_TD$TajimaD) #-0.765
Contemp_mean_TD_all <- mean(Aen_Contemp_TD$TajimaD) #0.410
All_mean_TD_all <- mean(Aen_All_TD$TajimaD) #-0.739

Contemp_mean_TD_all - Alb_mean_TD_all #1.17434

#### TD bootstrapping ####

#historical TD
boot_Hist_Aen_TD <- boot(data = Aen_Alb_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Hist_Aen_TD_95ci <- boot.ci(boot_Hist_Aen_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Hist_Aen_TD_95ci_normal <- Hist_Aen_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD

#contemporary TD
boot_Contemp_Aen_TD <- boot(data = Aen_Contemp_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  Contemp_Aen_TD_95ci <- boot.ci(boot_Contemp_Aen_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  Contemp_Aen_TD_95ci_normal <- Contemp_Aen_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
 
#all TD
boot_All_Aen_TD <- boot(data = Aen_All_TD$TajimaD, statistic = samp_mean, R = 1000) #1000 permutations of TD
  All_Aen_TD_95ci <- boot.ci(boot_All_Aen_TD, conf = 0.95, type = "norm") #get 95% CI for TD
  All_Aen_TD_95ci_normal <- All_Aen_TD_95ci$normal #pull out normal distribution 2.5 & 97.5 percentiles for TD
  
## write out TD summary tables ##
TD_mean <- as.data.frame(c(-0.765, 0.410, -0.739))

#TD summary table
TD_ci <- rbind(Hist_Aen_TD_95ci_normal, Contemp_Aen_TD_95ci_normal, 
               All_Aen_TD_95ci_normal) #combine df w/ci info for each into one dataframe
TD_sum <- cbind(TD_mean, TD_ci)
  colnames(TD_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
  TD_sum$Pop <- c("Hist", "Contemp", "All")

TD_sum$diff_lower <- TD_sum$mean - TD_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
TD_sum$diff_upper <- TD_sum$`97.5_per` - TD_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
#write.csv(TD_sum, "Data/Aen_Ham/Aen_TD_cis.csv")

## Calculate 95% for diff through time ##
Aen_TD_boot_diff <- as.data.frame(boot_Hist_Aen_TD$t)
  Aen_TD_boot_diff$Contemp_TD <- boot_Contemp_Aen_TD$t
  colnames(Aen_TD_boot_diff) <- c("Hist_TD", "Contemp_TD")

Aen_TD_boot_diff$TD_diff <- Aen_TD_boot_diff$Contemp_TD - Aen_TD_boot_diff$Hist_TD
  Aen_TD_boot_diff_CIs <- quantile(Aen_TD_boot_diff$TD_diff, c(0.025, 0.975)) #0.98900, 1.36678

#### Mann-Whitney U test ####

#subset data
#need Alb & Contemp TD in two separate vectors to compare
Alb_df <- TD_only_all[which(TD_only_all$Pop == "Albatross"), ]
  Alb_df_noNA <- Alb_df[which(Alb_df$TajimaD!='NA'),]
  Alb_TD <- Alb_df_noNA$TajimaD
Contemp_df <- TD_only_all[which(TD_only_all$Pop == "Contemporary"), ]
  Contemp_df_noNA <- Contemp_df[which(Contemp_df$TajimaD!='NA'),]
  Contemp_TD <- Contemp_df_noNA$TajimaD

#run Mann-Whitney U-test
TD_MUtest <- wilcox.test(Alb_TD, Contemp_TD) #W = 8776, p < 2.2e-16 #distributions not equal

#### Visualize data ####

#read in data
TD_only_all <- read.csv(here("Data/Aen_Ham", "Ham_A_TajimasD_combined_full.csv"), header = TRUE)
  TD_only_noall <- subset(TD_only_all, Pop != "All")

## Scatter plot ##
TD_scatter_plot <- ggplot(data = TD_only_noall, aes(x = NUM, y = TajimaD, color = Pop)) + 
  geom_point(size = 10) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 2, linetype = "dashed") + 
  scale_color_manual(values = c("#a25505", "#e1bf9b")) +
  ggtitle("Aen Tajima's D") + labs(y = "Tajima's D", x = "position") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())
TD_scatter_plot

## Density plots ##
#ordering x-axis
TD_only_noall$Pop[TD_only_noall$Pop == "Albatross"] <- "Historical"
TD_only_noall$Pop <- factor(TD_only_noall$Pop, levels = c("Historical", "Contemporary")) #ordering X-axis

#first density plot to get coordinates
TD_density_plot_1 <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, aes(x = TajimaD, y = Pop), scale = 1.2) 

#pull out density plot information
ingredients <- ggplot_build(TD_density_plot_1) %>% purrr::pluck("data", 1)
density_lines <- ingredients %>%
  group_by(group) %>% filter(density == max(density)) %>% ungroup()

#add mean coordinates
density_lines$mean <- c(-0.765, 0.41) #to get lines to plot at mean instead of elsewhere
density_lines$lower_ci <- c(-0.83807, 0.2337)
density_lines$upper_ci <- c(-0.68706, 0.57152)
density_lines$pop <- c("Historical", "Contemporary") #for coloring purposes

#density
TD_density_plot <- ggplot() + 
  geom_density_ridges2(data = TD_only_noall, 
                       aes(x = TajimaD, y = Pop, color = Pop, fill = Pop), 
                       alpha = 0.5, scale = 1.2) + 
  geom_vline(aes(xintercept = 0), linewidth = 4, color = "black", linetype = "dotted") + 
  geom_segment(data = density_lines, 
               aes(x = lower_ci, y = ymin, xend = lower_ci, yend = ymin+density*scale*iscale), 
               color = "#bd8850", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = upper_ci, y = ymin, xend = upper_ci, yend = ymin+density*scale*iscale), 
               color = "#bd8850", linewidth = 4, linetype = "solid") +
  geom_segment(data = density_lines, 
               aes(x = mean, y = ymin, xend = mean, yend = ymin+density*scale*iscale), 
               color = "black", linewidth = 4) +
  annotate("text", x = 3, y = 2.5, label = "A", size = 40) + 
  scale_color_manual(values = c("#d0ab84", "#ebdbca")) +
  scale_fill_manual(values = c("#d0ab84", "#ebdbca")) + 
  scale_y_discrete(expand = expansion(mult = c(0.01, .7))) + #bring y axis down x
  labs(x = "Tajima's D") +  
  theme_ridges() + 
  theme(plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black"), 
        axis.text.x = element_text(size = 55, color = "black"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 55, color = "black", vjust = -0.5, hjust = 0.5),
        legend.position = "none", 
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))
TD_density_plot

#### Identify change through time ####

#calculate diff through time
TD_diff <- as.data.frame(cbind(Alb_df$uniqseq, Alb_df$TajimaD, Contemp_df$TajimaD))
  colnames(TD_diff) <- c("uniqseq", "Alb_TD", "Contemp_TD")
  TD_diff$Alb_TD <- as.numeric(as.character(TD_diff$Alb_TD)) #bc was factor
  TD_diff$Contemp_TD <- as.numeric(as.character(TD_diff$Contemp_TD))
TD_diff$delta <- TD_diff$Contemp_TD - TD_diff$Alb_TD

#calculate mean and median delta  
mean(na.omit(TD_diff$delta)) #0.67844 --> individual sites getting more positive through time
median(na.omit(TD_diff$delta)) #0.59955 

#identify ones that became negative through time
#poss selective sweeps
TD_pot_selection <- subset(TD_diff, delta < 0)

#diff scatter plot
TD_diff_scatter_plot <- ggplot(data = TD_diff, aes(x = Alb_TD, y = Contemp_TD)) + 
  geom_point(size = 16, color = "#a25505", alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 4, linetype = "dashed", color = "black") + 
  annotate("text", x = -1.8, y = 2.8, label = "A", size = 40) + 
  labs(y = "Contemporary", x = "Albatross") + 
  ylim(-2, 3) + xlim(-2, 3) +
  theme_bw() + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1),
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
TD_diff_scatter_plot

#################################################################################################################################################

### extra stuff ####

Alb_contemp <- as.data.frame(cbind(Alb_df$uniqseq, Alb_df$TajimaD, Contemp_df$TajimaD))
  colnames <- c("uniqseq", "Alb_TD", "Contemp_TD")
  colnames(Alb_contemp) <- colnames
  Alb_contemp$Alb_TD <- as.character(Alb_contemp$Alb_TD) #bc was factor and factor --> numeric doesn't work
  Alb_contemp$Alb_TD <- as.numeric(Alb_contemp$Alb_TD)
  Alb_contemp$Contemp_TD <- as.character(Alb_contemp$Contemp_TD)
  Alb_contemp$Contemp_TD <- as.numeric(Alb_contemp$Contemp_TD)

Alb_contemp$delta <- Alb_contemp$Contemp_TD - Alb_contemp$Alb_TD
class(Alb_contemp$Alb_TD)
delta <- na.omit(Alb_contemp$delta)
mean(delta)
median(delta)

#subset by starting tD
Alb_TD_lessneg1 <- subset(Alb_contemp, Alb_TD <= -1)
  Alb_TD_lessneg1$ID <- c(rep("less than -1", 159))
Alb_TD_0neg1 <- subset(Alb_contemp, Alb_TD <= 0 & Alb_TD > -1)
  Alb_TD_0neg1$ID <- c(rep("0 to -1", 175))
Alb_TD_0pos1 <- subset(Alb_contemp, Alb_TD > 0 & Alb_TD < 1)
  Alb_TD_0pos1$ID <- c(rep("0 to 1", 31))
Alb_TD_greatpos1 <- subset(Alb_contemp, Alb_TD >= 1)
  Alb_TD_greatpos1$ID <- c(rep("great than 1", 21))
  
Alb_contemp_cat <- as.data.frame(rbind(Alb_TD_lessneg1, Alb_TD_0neg1, Alb_TD_0pos1, Alb_TD_greatpos1))
  Alb_contemp_cat$ID2 <- factor(Alb_contemp_cat$ID, levels = c("less than -1", "0 to -1", "0 to 1", "great than 1")) #ordering X-axis

TD_delta_boxplot <- ggplot(data = Alb_contemp, aes(y = delta)) + 
  geom_boxplot(lwd = 1.5) + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#000000")
TD_delta_boxplot_annotated <- TD_delta_boxplot + theme_bw() + 
  labs(x = "Sampling Location", y = "Tajima's D") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 28, color = "black"),
        axis.title = element_text(size = 30), legend.position = "top",
        legend.text = element_text(size = 30), legend.title = element_blank())
TD_delta_boxplot_annotated

TD_delta_cat_boxplot <- ggplot(data = Alb_contemp_cat, aes(x = ID2, y = delta)) + 
  geom_boxplot(lwd = 1.5) + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#000000") +
  annotate("text", x = 3.75, y = 2, label = "Mean overall = 0.040", size = 10) + 
  annotate("text", x = 3.75, y = 1.5, label = "Median overall = 0.087", size = 10)
TD_delta_cat_boxplot_annotated <- TD_delta_cat_boxplot + theme_bw() + 
  labs(x = "Albatross TD", y = "Delta (Contemp TD - Alb TD)") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 28, color = "black"),
        axis.title = element_text(size = 30), legend.position = "top",
        legend.text = element_text(size = 30), legend.title = element_blank())
TD_delta_cat_boxplot_annotated

##########################################################################################################################################

