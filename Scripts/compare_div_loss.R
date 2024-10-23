############################################# Script for Diversity Loss Comparisons  #############################################################

#Compares loss in genetic diversity in this study to other observed losses
#Uses supp material from Leigh et al. (2019) - https://doi.org/10.1111/eva.12810

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(here) #v.1.0.1

#read in data
Leighetal_df <- read.csv(here("Data", "Leigh_et.al.data_used_in_synthesis.csv"))
  Leighetal_df <- Leighetal_df[1:99, 1:30] #loaded a bunch of empty columns and rows, so trimming to actual data

#################################################################################################################################################

######## Calculate average diversity loss in Leigh et al ########

#keep only classes that have >5 data points
Leighetal_great5_df <- subset(Leighetal_df, Leighetal_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == "Actinopterygii" |
                                Leighetal_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == "Aves" | 
                                Leighetal_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == "Insecta" | 
                                Leighetal_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == "Mammalia")
  
#rename to common taxa
Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life.[Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == 
                                                                       "Actinopterygii"] <- "Fish"
Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life.[Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == 
                                                                       "Aves"] <- "Birds"
Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life.[Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == 
                                                                       "Insecta"] <- "Insects"
Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life.[Leighetal_great5_df$Class..In.Species.2000...ITIS.Catalogue.of.Life. == 
                                                                       "Mammalia"] <- "Mammals"
  
## calculate percent changes in mean Ho & He over time ##
#positive values indicate a LOSS over time
Leighetal_great5_df$PercChange_Ho <- ((Leighetal_great5_df$Modern.Observed.Heterozygosity.Mean - Leighetal_great5_df$Historial.Observed.Heterozygosity.Mean)/
                                        Leighetal_great5_df$Historial.Observed.Heterozygosity.Mean)*100
Leighetal_great5_df$PercChange_He <- ((Leighetal_great5_df$Modern.Expected..Heterozygosity.Mean - Leighetal_great5_df$Historial.Expected..Heterozygosity.Mean)/
                                        Leighetal_great5_df$Historial.Expected..Heterozygosity.Mean)*100

#pull out relevant columns from Leighetal
Leighetal_HoHe_great5_df <- Leighetal_great5_df[ , grep("Class..In.Species.2000...ITIS.Catalogue.of.Life.|Historial.Observed.Heterozygosity.Mean|Historial.Expected..Heterozygosity.Mean|Modern.Observed.Heterozygosity.Mean|Modern.Expected..Heterozygosity.Mean|PercChange_Ho|PercChange_He", 
                                                        colnames(Leighetal_great5_df))]
  Leighetal_HoHe_great5_df$Dataset <- "Leighetal"
  Leighetal_HoHe_great5_df$n <- 1
  colnames(Leighetal_HoHe_great5_df) <- c("Species", "Ho_Hist", "He_Hist", "Ho_Contemp", "He_Contemp", 
                                          "PercChange_Ho","PercChange_He", "Dataset", "n")

#calculate averages across classes
avg_Ho_He_loss_by_class <-
  Leighetal_great5_df %>%
  group_by(Class..In.Species.2000...ITIS.Catalogue.of.Life.) %>%
  summarize(avg_Holoss <- mean(PercChange_Ho, na.rm = TRUE), #calculate average loss in Ho over time per tax group
            avg_Heloss <- mean(PercChange_He, na.rm = TRUE), #calculate average loss in He over time per tax group
            n()) #number of observations in each tax group

avg_Ho_He_loss_by_class$Dataset <- "Leighetal"
  colnames(avg_Ho_He_loss_by_class) <- c("Species", "PercChange_Ho", "PercChange_He", "n", "Dataset")

##################################################################################################################################################

######## Merge Leigh et al & current data ########
#Merging means and raw data separately --> Ela & Gmi data is the same regardless

#make dataframe from Ela & Gmi div data
Ho_Hist <- c(0.11900, 0.10796) #Ela, Gmi
He_Hist <- c(0.12109, 0.10218) #Ela, Gmi

Ho_Contemp <- c(0.11577, 0.10165) #Ela, Gmi
He_Contemp <- c(0.11408, 0.09606) #Ela, Gmi

combined_df <- data.frame(Ho_Hist, He_Hist, Ho_Contemp, He_Contemp)
  combined_df$Species <- c("E. lat", "G. min")
  combined_df$Dataset <- "original"
  combined_df$PercChange_Ho <- ((combined_df$Ho_Contemp - combined_df$Ho_Hist)/combined_df$Ho_Hist)*100
  combined_df$PercChange_He <- ((combined_df$He_Contemp - combined_df$He_Hist)/combined_df$He_Hist)*100
  combined_df$n <- 1

#merge raw df together
combined_all_raw_df <- rbind(combined_df, Leighetal_HoHe_great5_df)

#merge mean df together
combined_nopoint_df <- combined_df[ , -grep("Ho_|He_", colnames(combined_df))] #remove point estimate columns for merging
  
combined_all_mean_df <- rbind(combined_nopoint_df, avg_Ho_He_loss_by_class)

##################################################################################################################################################
  
######## Plot data ########

#violin plot of Ho loss by taxa
Ho <- ggplot() + 
  geom_point(data = combined_all_mean_df, 
             aes(x = factor(Species, levels = c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                 y = PercChange_Ho, color = Species), size = 0) + #setting this first to get right order on X-axis (but not really plotting anything)
  geom_violin(data = combined_all_raw_df[combined_all_raw_df$Dataset == "Leighetal", ], 
              aes(x = factor(Species, levels = c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                  y = PercChange_Ho, fill = Species), color = NA, lwd = 2) + 
  geom_point(data = combined_all_mean_df, 
             aes(x = factor(Species, c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                 y = PercChange_Ho, color = Species, shape = Dataset, size = Dataset)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 2, color = "black") +
  scale_color_manual(values = c("#757677", "#16537e", "#757677", "#afc8a4", "#757677", "#757677")) + 
  scale_fill_manual(values = c("#ada5d0", "#DDCC77", "#CC6677", "#B7E0F4")) +
  scale_shape_manual(values = c(18, 17)) +
  scale_size_manual(values = c(12, 14)) + 
  ylab(bquote("% Change in"~H[o])) + 
  xlab("Taxa") + 
  theme_bw() + 
  coord_flip() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 10)), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 10)), 
        axis.title.x = element_text(size = 65, color = "black", vjust = -1),
        axis.title.y = element_blank(),
        legend.position = "none", 
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ho

#violin plot of He loss by taxa
He <- ggplot() + 
  geom_point(data = combined_all_mean_df, 
             aes(x = factor(Species, levels = c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                 y = PercChange_He, color = Species), size = 0) + #setting this first to get right order on X-axis (but not really plotting anything)
  geom_violin(data = combined_all_raw_df[combined_all_raw_df$Dataset == "Leighetal", ], 
              aes(x = factor(Species, levels = c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                  y = PercChange_He, fill = Species), color = NA, lwd = 2) + 
  geom_point(data = combined_all_mean_df, 
             aes(x = factor(Species, c("G. min", "E. lat", "Fish", "Mammals", "Birds", "Insects")), 
                 y = PercChange_He, color = Species, shape = Dataset, size = Dataset)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 2, color = "black") +
  scale_color_manual(values = c("#757677", "#16537e", "#757677", "#afc8a4", "#757677", "#757677")) + 
  scale_fill_manual(values = c("#ada5d0", "#DDCC77", "#CC6677", "#B7E0F4")) +
  scale_shape_manual(values = c(18, 17)) +
  scale_size_manual(values = c(12, 14)) + 
  ylab(bquote("% Change in"~H[e])) + 
  xlab("Taxa") + 
  theme_bw() + 
  coord_flip() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 10)), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 10)), 
        axis.title.x = element_text(size = 65, color = "black", vjust = -1),
        axis.title.y = element_blank(),
        legend.position = "none", 
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
He
