################################################### Script for Relatedness  #######################################################

#Generates relatedness estimates w/relatedness r package (calls coancestry fortran program)

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load packages
#related had to be installed with Rtools -- download tar file and install
library(tidyverse) #2.0.0
library(related) #1.0
library(boot) #v.1.3.28

#mean function for bootstrapping
samp_mean <- function(x, i) {
  mean(x[i])
}

########################################################################################################################################################

######## Gmi relatedness estimates ########

#### Format data ####
#read in genepop data to transform
#removed header rows first in unix (including pop rows)
Gmi_genepop <- read.table("Data/Gmi_Ham/Gmi.renamed.noLD.A.nohighhet.Ham.geno",
                          sep = " ") #removed header rows first in unix (including pop rows)
  Gmi_genepop <- Gmi_genepop[,-2]
  Gmi_genepop <- Gmi_genepop[,-2] #doing twice, to remove blank columns
  dim(Gmi_genepop) #check dim

Gmi_genepop[Gmi_genepop == "0"] <- "000000" #force zeros to read properly

#transform columns
Gmi_genepop2 <- unite(Gmi_genepop, all, sep = "", 
                      2:7135, remove = TRUE) #merge all genotype columns together
  colnames <- c("ind", "genotypes") #rename columns
  colnames(Gmi_genepop2) <- colnames

#create number of columns to split genotypes into
n1 <- nchar(as.character(Gmi_genepop2$genotypes[2])) - 3
s1 <- seq(3, n1, by = 3)
nm1 <- paste0("x", seq_len(length(s1) +1))

#separate genotypes out into columns (each SNP gets 2 columns, 1 per allele)
Gmi_related_format <- Gmi_genepop2 %>% 
  separate(genotypes, into = nm1, sep = s1)

#write out
write.table(Gmi_related_format, "Data/Gmi_Ham/Gmi.relatedness.format.txt", sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#### Calculate point estimates of relatedness ####
#use wang estimator bc those seem to be least biased w/small sample sizes (Wang 2017) --> although if unmodified, may slightly underestimate relatedness
rel_output <- coancestry("Data/Gmi_Ham/Gmi.relatedness.format.txt", 
                         wang = 2) #calculates Wang point estimates

#pull out data by population
relatedness <- rel_output$relatedness
  relatedness_wang <- relatedness[, c("pair.no", "ind1.id", "ind2.id", "group", "wang")]

rel_Alb<- subset(relatedness_wang, grepl("^Gmi-AHam", ind1.id) & grepl("^Gmi-AHam", ind2.id))
  rel_Alb$Era <- "Albatross"
rel_Con<- subset(relatedness_wang, grepl("^Gmi-CBat", ind1.id) & grepl("^Gmi-CBat", ind2.id))
  rel_Con$Era <- "Contemporary"

rel_inpop_df <- rbind(rel_Alb, rel_Con)

#plot by Era
rel_violin_plot <- ggplot(data = rel_inpop_df, 
                          aes(x = Era, y = wang, color = Era, fill = Era)) + 
  geom_violin() + 
  scale_color_manual(values = c("#1c3b0e","#afc8a4")) +
  scale_fill_manual(values = c("#1c3b0e","#afc8a4")) +
  labs(y = "Pair-wise relatedness (wang)", x = "Era") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
rel_violin_plot

#### Calculate mean pairwise relatedness w/in pop ####

#calculate mean pairwise relatedness w/in pops
Alb_rel_mean <- mean(rel_Alb$wang) #0.0060
Contemp_rel_mean <- mean(rel_Con$wang) #0.0612

#calculate median pairwise relatedness w/in pops
Alb_rel_median <- median(rel_Alb$wang) #0.0111
Contemp_rel_median <- median(rel_Con$wang) #0.0863

rel_mean_df <- data.frame(Alb_rel_mean, Contemp_rel_mean)

#### Bootstrap for 95% CIs ####

#bootstrap Alb population
boot_Alb_wang <- boot(data = rel_Alb$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
  Alb_95ci_wang <- boot.ci(boot_Alb_wang, 
                           conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#bootstrap Contemp population
boot_Contemp_wang <- boot(data = rel_Con$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
  Contemp_95ci_wang <- boot.ci(boot_Contemp_wang, 
                               conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#summary table
t_rel_mean_df <- data.frame(t(rel_mean_df)) #transpose mean df
wang_mean_rel <- t_rel_mean_df[, 1] #pull out wang mean relatedness
estimator_vector <- c("wang", "wang") #create column with estimator name

wang_mean_rel <- data.frame(wang_mean_rel, estimator_vector) #combine means & estimator vector
  colnames(wang_mean_rel) <- c("mean", "estimator")
  rownames(wang_mean_rel) <- c("Albatross", "Contemporary")

Alb_95ci_wang_normal <- Alb_95ci_wang$normal #pull out normal distribution  2.5 & 97.5 percentiles for wang ci
  Contemp_95ci_wang_normal <- Contemp_95ci_wang$normal

wang_norm_ci <- rbind(Alb_95ci_wang_normal, Contemp_95ci_wang_normal) #combine df w/ci info for each pop into one dataframe
  colnames(wang_norm_ci) <- c("ci", "2.5_per", "97.5_per")
  rownames(wang_norm_ci) <- c("Albatross", "Contemporary")

#merge dataframes
mean_rel <- cbind(wang_mean_rel, wang_norm_ci) #combine wang ci info and sample mean info into one dataframe
  pop_vector <- c("Albatross", "Contemporary")
  
mean_rel <- cbind(mean_rel, pop_vector)
  colnames(mean_rel) <- c("mean", "estimator", "CI", "2.5_per", "97.5_per", "Pop")

mean_rel$diff_lower <- mean_rel$mean - mean_rel$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
mean_rel$diff_upper <- mean_rel$`97.5_per` - mean_rel$mean # calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out
#write.csv(mean_rel, "Data/Gmi_Ham/Gmi_wang_relatedness_cis.csv")
#write.csv(relatedness, "Data/Gmi_Ham/Gmi_relatedness_raw.csv")

#plot means
mean_rel_plot <- ggplot(data = mean_rel, aes(x = Pop, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#1c3b0e","#afc8a4")) +
  scale_shape_manual(values = c(19, 15)) +
  ggtitle("Gmi Wang pair-wise relatedness 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
mean_rel_plot

########################################################################################################################################################

######## Ela relatedness estimates ########

#### Format data ####
#read in genepop data to transform
#removed header rows first in unix (including pop rows)
Ela_genepop <- read.table("Data/Ela_Ham/Lle.renamed.noLD.Ela.nohighhet.geno",
                          sep = " ") #removed header rows first in unix (including pop rows)
  Ela_genepop <- Ela_genepop[,-2]
  Ela_genepop <- Ela_genepop[,-2] #doing twice, to remove blank columns
  dim(Ela_genepop) #check dim

Ela_genepop[Ela_genepop == "0"] <- "000000" #force zeros to read properly

#transform columns
Ela_genepop2 <- unite(Ela_genepop, all, sep = "", 
                      2:47064, remove = TRUE) #merge all genotype columns together
  colnames <- c("ind", "genotypes") #rename columns
  colnames(Ela_genepop2) <- colnames

#create number of columns to split genotypes into
n1 <- nchar(as.character(Ela_genepop2$genotypes[2])) - 3
s1 <- seq(3, n1, by = 3)
nm1 <- paste0("x", seq_len(length(s1) +1))

#separate genotypes out into columns (each SNP gets 2 columns, 1 per allele)
Ela_related_format <- Ela_genepop2 %>% 
  separate(genotypes, into = nm1, sep = s1)

#write out
write.table(Ela_related_format, "Data/Ela_Ham/Ela.relatedness.format.txt", sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#### Calculate point estimates of relatedness ####
#use wang estimator bc those seem to be least biased w/small sample sizes (Wang 2017) --> although if unmodified, may slightly underestimate relatedness
rel_output <- coancestry("Data/Ela_Ham/Ela.relatedness.format.txt", 
                         wang = 2) #calculates Wang point estimates

#pull out data by population
relatedness <- rel_output$relatedness
  relatedness_wang <- relatedness[, c("pair.no", "ind1.id", "ind2.id", "group", "wang")]

rel_Alb<- subset(relatedness_wang, grepl("^Lle-AHam", ind1.id) & grepl("^Lle-AHam", ind2.id))
  rel_Alb$Era <- "Albatross"
rel_Con<- subset(relatedness_wang, grepl("^Lle-CNas", ind1.id) & grepl("^Lle-CNas", ind2.id))
  rel_Con$Era <- "Contemporary"

rel_inpop_df <- rbind(rel_Alb, rel_Con)

#plot by Era
rel_violin_plot <- ggplot(data = rel_inpop_df, 
                          aes(x = Era, y = wang, color = Era, fill = Era)) + 
  geom_violin() + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  scale_fill_manual(values = c("#16537e", "#8aa9be")) +
  labs(y = "Pair-wise relatedness (wang)", x = "Era") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
rel_violin_plot

#### Calculate mean pairwise relatedness w/in pop ####

#calculate mean pairwise relatedness w/in pops
Alb_rel_mean <- mean(rel_Alb$wang) #-0.0508
Contemp_rel_mean <- mean(rel_Con$wang) #0.0223

#calculate median pairwise relatedness w/in pops
Alb_rel_median <- median(rel_Alb$wang) #-0.0504
Contemp_rel_median <- median(rel_Con$wang) #0.0228

rel_mean_df <- data.frame(Alb_rel_mean, Contemp_rel_mean)

#### Bootstrap for 95% CIs ####

#bootstrap Alb population
boot_Alb_wang <- boot(data = rel_Alb$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
  Alb_95ci_wang <- boot.ci(boot_Alb_wang, 
                           conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#bootstrap Contemp population
boot_Contemp_wang <- boot(data = rel_Con$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
  Contemp_95ci_wang <- boot.ci(boot_Contemp_wang, 
                               conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#summary table
t_rel_mean_df <- data.frame(t(rel_mean_df)) #transpose mean df
  wang_mean_rel <- t_rel_mean_df[, 1] #pull out wang mean relatedness
estimator_vector <- c("wang", "wang") #create column with estimator name

wang_mean_rel <- data.frame(wang_mean_rel, estimator_vector) #combine means & estimator vector
  colnames(wang_mean_rel) <- c("mean", "estimator")
  rownames(wang_mean_rel) <- c("Albatross", "Contemporary")

Alb_95ci_wang_normal <- Alb_95ci_wang$normal #pull out normal distribution  2.5 & 97.5 percentiles for wang ci
Contemp_95ci_wang_normal <- Contemp_95ci_wang$normal

wang_norm_ci <- rbind(Alb_95ci_wang_normal, Contemp_95ci_wang_normal) #combine df w/ci info for each pop into one dataframe
  colnames(wang_norm_ci) <- c("ci", "2.5_per", "97.5_per")
  rownames(wang_norm_ci) <- c("Albatross", "Contemporary")

#merge dataframes
mean_rel <- cbind(wang_mean_rel, wang_norm_ci) #combine wang ci info and sample mean info into one dataframe
  pop_vector <- c("Albatross", "Contemporary")
  
mean_rel <- cbind(mean_rel, pop_vector)
  colnames(mean_rel) <- c("mean", "estimator", "CI", "2.5_per", "97.5_per", "Pop")

mean_rel$diff_lower <- mean_rel$mean - mean_rel$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
mean_rel$diff_upper <- mean_rel$`97.5_per` - mean_rel$mean # calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out
#write.csv(mean_rel, "Data/Ela_Ham/Ela_wang_relatedness_cis.csv")
#write.csv(relatedness, "Data/Ela_Ham/Ela_relatedness_raw.csv")

#plot means
mean_rel_plot <- ggplot(data = mean_rel, aes(x = Pop, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#16537e", "#8aa9be")) +
  scale_shape_manual(values = c(19, 15)) +
  ggtitle("Ela Wang pair-wise relatedness 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
mean_rel_plot

########################################################################################################################################################

######## Aen relatedness estimates ########

#### Format data ####
#read in genepop data to transform
#removed header rows first in unix (including pop rows)
Aen_genepop <- read.table("Data/Aen_Ham/Aen.renamed.noLD.A.nohighhet.geno",
                          sep = " ") #removed header rows first in unix (including pop rows)
  Aen_genepop <- Aen_genepop[,-2]
  Aen_genepop <- Aen_genepop[,-2] #doing twice, to remove blank columns
  dim(Aen_genepop) #check dim

Aen_genepop[Aen_genepop == "0"] <- "000000" #force zeros to read properly
  
#transform columns
Aen_genepop2 <- unite(Aen_genepop, all, sep = "", 
                      2:3234, remove = TRUE) #merge all genotype columns together
  colnames <- c("ind", "genotypes") #rename columns
  colnames(Aen_genepop2) <- colnames

#create number of columns to split genotypes into
n1 <- nchar(as.character(Aen_genepop4$genotypes[2])) - 3
s1 <- seq(3, n1, by = 3)
nm1 <- paste0("x", seq_len(length(s1) +1))

#separate genotypes out into columns (each SNP gets 2 columns, 1 per allele)
Aen_related_format <- Aen_genepop2 %>% 
  separate(genotypes, into = nm1, sep = s1)

#write out
write.table(Aen_related_format, "Data/Aen_Ham/Aen.relatedness.format.txt", sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#### Calculate point estimates of relatedness ####
#use wang estimator bc those seem to be least biased w/small sample sizes (Wang 2017) --> although if unmodified, may slightly underestimate relatedness
rel_output <- coancestry("Data/Aen_Ham/Aen.relatedness.format.txt", 
                         wang = 2) #calculates Wang point estimates

#pull out data by population
relatedness <- rel_output$relatedness
  relatedness_wang <- relatedness[, c("pair.no", "ind1.id", "ind2.id", "group", "wang")]

rel_Alb<- subset(relatedness_wang, grepl("^Aen-AHam", ind1.id) & grepl("^Aen-AHam", ind2.id))
  rel_Alb$Era <- "Albatross"
rel_Con<- subset(relatedness_wang, grepl("^Aen-Cbat", ind1.id) & grepl("^Aen-Cbat", ind2.id))
  rel_Con$Era <- "Contemporary"

rel_inpop_df <- rbind(rel_Alb, rel_Con)

#plot by Era
rel_violin_plot <- ggplot(data = rel_inpop_df, 
                          aes(x = Era, y = wang, color = Era, fill = Era)) + 
  geom_violin() + 
  scale_color_manual(values = c("#a25505", "#e3ccb4")) +
  scale_fill_manual(values = c("#a25505", "#e3ccb4")) +
  labs(y = "Pair-wise relatedness (wang)", x = "Era") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
rel_violin_plot

#### Calculate mean pairwise relatedness w/in pop ####

#calculate mean pairwise relatedness w/in pops
Alb_rel_mean <- mean(rel_Alb$wang) #-0.5419
Contemp_rel_mean <- mean(rel_Con$wang) #0.0482

#calculate median pairwise relatedness w/in pops
Alb_rel_median <- median(rel_Alb$wang) #-0.5326
Contemp_rel_median <- median(rel_Con$wang) #0.0465

rel_mean_df <- data.frame(Alb_rel_mean, Contemp_rel_mean)

#### Bootstrap for 95% CIs ####

#bootstrap Alb population
boot_Alb_wang <- boot(data = rel_Alb$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
Alb_95ci_wang <- boot.ci(boot_Alb_wang, 
                       conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#bootstrap Contemp population
boot_Contemp_wang <- boot(data = rel_Con$wang, statistic = samp_mean, R = 1000) #1000 permutations of wang mean pairwise relatedness
Contemp_95ci_wang <- boot.ci(boot_Contemp_wang, 
                             conf = 0.95, type = "norm") #get 95% CI for wang pairwise relatedness

#summary table
t_rel_mean_df <- data.frame(t(rel_mean_df)) #transpose mean df
  wang_mean_rel <- t_rel_mean_df[, 1] #pull out wang mean relatedness
estimator_vector <- c("wang", "wang") #create column with estimator name

wang_mean_rel <- data.frame(wang_mean_rel, estimator_vector) #combine means & estimator vector
  colnames(wang_mean_rel) <- c("mean", "estimator")
  rownames(wang_mean_rel) <- c("Albatross", "Contemporary")

Alb_95ci_wang_normal <- Alb_95ci_wang$normal #pull out normal distribution  2.5 & 97.5 percentiles for wang ci
Contemp_95ci_wang_normal <- Contemp_95ci_wang$normal

wang_norm_ci <- rbind(Alb_95ci_wang_normal, Contemp_95ci_wang_normal) #combine df w/ci info for each pop into one dataframe
  colnames(wang_norm_ci) <- c("ci", "2.5_per", "97.5_per")
  rownames(wang_norm_ci) <- c("Albatross", "Contemporary")

#merge dataframes
mean_rel <- cbind(wang_mean_rel, wang_norm_ci) #combine wang ci info and sample mean info into one dataframe
  pop_vector <- c("Albatross", "Contemporary")
  
mean_rel <- cbind(mean_rel, pop_vector)
  colnames(mean_rel) <- c("mean", "estimator", "CI", "2.5_per", "97.5_per", "Pop")

mean_rel$diff_lower <- mean_rel$mean - mean_rel$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
mean_rel$diff_upper <- mean_rel$`97.5_per` - mean_rel$mean # calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out
#write.csv(mean_rel, "Data/Aen_Ham/Aen_wang_relatedness_cis.csv")
#write.csv(relatedness, "Data/Aen_Ham/Aen_relatedness_raw.csv")

#plot means
mean_rel_plot <- ggplot(data = mean_rel, aes(x = Pop, y = mean, color = Pop, shape = Pop)) + 
  geom_point(aes(), size = 10, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper), width = 0, linewidth = 2, 
                position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#a25505", "#e3ccb4")) +
  scale_shape_manual(values = c(19, 15)) +
  ggtitle("Aen Wang pair-wise relatedness 95% CIs") + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), 
        plot.title = element_text(size = 30, color = "black", face = "bold"),
        axis.title = element_text(size = 28, color = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 28, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.position = "top",
        legend.title = element_blank())
mean_rel_plot
