######################################### Script for Creating PCA Plots  ########################################################

#PCA from eigenval & eigenvec files made with Plink
#plot size = 1500 x 1500

#################################################################################################################################################

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(here) #v.1.0.1

eigenvec_names <- c("Population", "Individual", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                    "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

###############################################################################################################################################################

######## Gmi PCAs ########
  
#read in data
preHWE_Gmi_Ham_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.rescaled.preHWE.eigenval"), header = FALSE, sep = " ")
  preHWE_Gmi_Ham_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.rescaled.preHWE.eigenvec"), header = FALSE, sep = " ")
  colnames(preHWE_Gmi_Ham_data) <- eigenvec_names
Gmi_Ham_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.eigenval"), header = FALSE, sep = " ") #noLD
  Gmi_Ham_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_data) <- eigenvec_names
Gmi_Ham_Ham_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.Ham.eigenval"), header = FALSE, sep = " ") #noLD
  Gmi_Ham_Ham_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.Ham.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_Ham_data) <- eigenvec_names
Gmi_Ham_A_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.eigenval"), header = FALSE, sep = " ")
  Gmi_Ham_A_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_A_data) <- eigenvec_names
Gmi_Ham_B_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.B.eigenval"), header = FALSE, sep = " ")
  Gmi_Ham_B_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.B.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_B_data) <- eigenvec_names
Gmi_Ham_A_nohighhet_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.nohighhet.eigenval"), header = FALSE, sep = " ")
  Gmi_Ham_A_nohighhet_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.nohighhet.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_A_nohighhet_data) <- eigenvec_names
Gmi_Ham_A_Ham_nohighhet_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.nohighhet.Ham.eigenval"), header = FALSE, sep = " ")
  Gmi_Ham_A_Ham_nohighhet_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.nohighhet.Ham.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_A_Ham_nohighhet_data) <- eigenvec_names
Gmi_Ham_A_Ham_eigenval <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.Ham.eigenval"), header = FALSE, sep = " ")
  Gmi_Ham_A_Ham_data <- read.csv(here("Data/Gmi_Ham/PCAs", "PIRE.Gmi.Ham.A.Ham.eigenvec"), header = FALSE, sep = " ")
  colnames(Gmi_Ham_A_Ham_data) <- eigenvec_names
  
#### preHWE PCA ####    
    
#add columns for Location & Era
preHWE_Gmi_Ham_data <- preHWE_Gmi_Ham_data %>%
  mutate(Location =
               case_when(endsWith(Population, "Bas") ~ "Basud",
                     endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "ABas") ~ "Historical",
                     endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary",
                     endsWith(Population, "CBas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (preHWE_Gmi_Ham_eigenval[1,1] / sum(preHWE_Gmi_Ham_eigenval$V1))*100
varPC2 <- (preHWE_Gmi_Ham_eigenval[2,1] / sum(preHWE_Gmi_Ham_eigenval$V1))*100
        
#preHWE PCA
preHWE_PCA_12 <- ggplot(data = preHWE_Gmi_Ham_data, aes(x = PC1, y = PC2, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
    scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 79.49% of total variance)", y = "PC2 (explains 9.24% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
preHWE_PCA_12

#### post HWE PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_data <- Gmi_Ham_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Bas") ~ "Basud",
                     endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "ABas") ~ "Historical",
                     endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary",
                     endsWith(Population, "CBas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_eigenval[1,1] / sum(Gmi_Ham_eigenval$V1))*100
varPC2 <- (Gmi_Ham_eigenval[2,1] / sum(Gmi_Ham_eigenval$V1))*100
varPC3 <- (Gmi_Ham_eigenval[3,1] / sum(Gmi_Ham_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_data, aes(x = PC1, y = PC2, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 66.69% of total variance)", y = "PC2 (explains 13.41% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_data, aes(x = PC1, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 66.69% of total variance)", y = "PC3 (explains 1.94% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_data, aes(x = PC2, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC2 (explains 13.41% of total variance)", y = "PC3 (explains 1.94% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Ham PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_Ham_data <- Gmi_Ham_Ham_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_Ham_eigenval[1,1] / sum(Gmi_Ham_Ham_eigenval$V1))*100
varPC2 <- (Gmi_Ham_Ham_eigenval[2,1] / sum(Gmi_Ham_Ham_eigenval$V1))*100
varPC3 <- (Gmi_Ham_Ham_eigenval[3,1] / sum(Gmi_Ham_Ham_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_Ham_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 61.17% of total variance)", y = "PC2 (explains 13.14% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_Ham_data, 
                 aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 61.17% of total variance)", y = "PC3 (explains 2.47% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_A_Ham_nohighhet_unrelated_data, 
                 aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 13.14% of total variance)", y = "PC3 (explains 2.47% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Species A PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_A_data <- Gmi_Ham_A_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Bas") ~ "Basud",
                     endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "ABas") ~ "Historical",
                     endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary",
                     endsWith(Population, "CBas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_A_eigenval[1,1] / sum(Gmi_Ham_A_eigenval$V1))*100
varPC2 <- (Gmi_Ham_A_eigenval[2,1] / sum(Gmi_Ham_A_eigenval$V1))*100
varPC3 <- (Gmi_Ham_A_eigenval[3,1] / sum(Gmi_Ham_A_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_A_data, aes(x = PC1, y = PC2, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 24.43% of total variance)", y = "PC2 (explains 5.81% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_A_data, aes(x = PC1, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 24.43% of total variance)", y = "PC3 (explains 5.40% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_A_data, aes(x = PC2, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC2 (explains 5.81% of total variance)", y = "PC3 (explains 5.40% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Species B PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_B_data <- Gmi_Ham_B_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_B_eigenval[1,1] / sum(Gmi_Ham_B_eigenval$V1))*100
varPC2 <- (Gmi_Ham_B_eigenval[2,1] / sum(Gmi_Ham_B_eigenval$V1))*100
varPC3 <- (Gmi_Ham_B_eigenval[3,1] / sum(Gmi_Ham_B_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_B_data, aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 28.96% of total variance)", y = "PC2 (explains 16.37% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_B_data, aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 28.96% of total variance)", y = "PC3 (explains 7.74% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_B_data, aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 16.37% of total variance)", y = "PC3 (explains 7.74% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Species A nohighhet PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_A_nohighhet_data <- Gmi_Ham_A_nohighhet_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Bas") ~ "Basud",
                     endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "ABas") ~ "Historical",
                     endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary",
                     endsWith(Population, "CBas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_A_nohighhet_eigenval[1,1] / sum(Gmi_Ham_A_nohighhet_eigenval$V1))*100
varPC2 <- (Gmi_Ham_A_nohighhet_eigenval[2,1] / sum(Gmi_Ham_A_nohighhet_eigenval$V1))*100
varPC3 <- (Gmi_Ham_A_nohighhet_eigenval[3,1] / sum(Gmi_Ham_A_nohighhet_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_A_nohighhet_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 12.69% of total variance)", y = "PC2 (explains 7.17% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_A_nohighhet_data, 
                 aes(x = PC1, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC1 (explains 12.69% of total variance)", y = "PC3 (explains 5.92% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_A_nohighhet_data, 
                 aes(x = PC2, y = PC3, color = Era, shape = Location)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Basud River", "Hamilo Cove")) + 
  labs(x = "PC2 (explains 7.17% of total variance)", y = "PC3 (explains 5.92% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Species A Ham nohighhet PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_A_Ham_nohighhet_data <- Gmi_Ham_A_Ham_nohighhet_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_A_Ham_nohighhet_eigenval[1,1] / sum(Gmi_Ham_A_Ham_nohighhet_eigenval$V1))*100
varPC2 <- (Gmi_Ham_A_Ham_nohighhet_eigenval[2,1] / sum(Gmi_Ham_A_Ham_nohighhet_eigenval$V1))*100
varPC3 <- (Gmi_Ham_A_Ham_nohighhet_eigenval[3,1] / sum(Gmi_Ham_A_Ham_nohighhet_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_A_Ham_nohighhet_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 11.16% of total variance)", y = "PC2 (explains 7.17% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_A_Ham_nohighhet_data, 
                 aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 11.16% of total variance)", y = "PC3 (explains 5.67% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_A_Ham_nohighhet_data, 
                 aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 7.17% of total variance)", y = "PC3 (explains 5.67% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Species A Ham PCA ####    
#noLD

#add columns for Location & Era
Gmi_Ham_A_Ham_data <- Gmi_Ham_A_Ham_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Bat") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population,"CBat") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Gmi_Ham_A_Ham_eigenval[1,1] / sum(Gmi_Ham_A_Ham_eigenval$V1))*100
varPC2 <- (Gmi_Ham_A_Ham_eigenval[2,1] / sum(Gmi_Ham_A_Ham_eigenval$V1))*100
varPC3 <- (Gmi_Ham_A_Ham_eigenval[3,1] / sum(Gmi_Ham_A_Ham_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Gmi_Ham_A_Ham_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 20.40% of total variance)", y = "PC2 (explains 5.76% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Gmi_Ham_A_Ham_data, 
                 aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 20.40% of total variance)", y = "PC3 (explains 5.49% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Gmi_Ham_A_Ham_nohighhet_data, 
                 aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 5.76% of total variance)", y = "PC3 (explains 5.49% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

###############################################################################################################################################################

######## Ela PCAs ########

#read in data
preHWE_Ela_Ham_eigenval <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.rescaled.preHWE.eigenval"), header = FALSE, sep = " ")
  preHWE_Ela_Ham_data <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.rescaled.preHWE.eigenvec"), header = FALSE, sep = " ")
  colnames(preHWE_Ela_Ham_data) <- eigenvec_names
Ela_Ham_eigenval <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.eigenval"), header = FALSE, sep = " ") #noLD
  Ela_Ham_data <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.eigenvec"), header = FALSE, sep = " ")
  colnames(Ela_Ham_data) <- eigenvec_names
Ela_Ham_Ela_eigenval <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Ela.eigenval"), header = FALSE, sep = " ")
  Ela_Ham_Ela_data <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Ela.eigenvec"), header = FALSE, sep = " ")
  colnames(Ela_Ham_Ela_data) <- eigenvec_names
Ela_Ham_Lle_eigenval <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Lle.eigenval"), header = FALSE, sep = " ")
  Ela_Ham_Lle_data <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Lle.eigenvec"), header = FALSE, sep = " ")
  colnames(Ela_Ham_Lle_data) <- eigenvec_names
Ela_Ham_Ela_nohighhet_eigenval <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Ela.nohighhet.eigenval"), header = FALSE, sep = " ")
  Ela_Ham_Ela_nohighhet_data <- read.csv(here("Data/Ela_Ham/PCAs", "PIRE.Ela.Ham.Ela.nohighhet.eigenvec"), header = FALSE, sep = " ")
  colnames(Ela_Ham_Ela_nohighhet_data) <- eigenvec_names

#### preHWE PCA ####    

#add columns for Location & Era
preHWE_Ela_Ham_data <- preHWE_Ela_Ham_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Nas") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population, "CNas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (preHWE_Ela_Ham_eigenval[1,1] / sum(preHWE_Ela_Ham_eigenval$V1))*100
varPC2 <- (preHWE_Ela_Ham_eigenval[2,1] / sum(preHWE_Ela_Ham_eigenval$V1))*100

#preHWE PCA
preHWE_PCA_12 <- ggplot(data = preHWE_Ela_Ham_data, aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 55.83% of total variance)", y = "PC2 (explains 3.58% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
preHWE_PCA_12

#### postHWE PCA ####    
#noLD

#add columns for Location & Era
Ela_Ham_data <- Ela_Ham_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Nas") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population, "CNas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Ela_Ham_eigenval[1,1] / sum(Ela_Ham_eigenval$V1))*100
varPC2 <- (Ela_Ham_eigenval[2,1] / sum(Ela_Ham_eigenval$V1))*100
varPC3 <- (Ela_Ham_eigenval[3,1] / sum(Ela_Ham_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Ela_Ham_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#8aa9be", "#16537e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 45.55% of total variance)", y = "PC2 (explains 4.44% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Ela_Ham_data, aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 45.55% of total variance)", y = "PC3 (explains 3.76% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Ela_Ham_data, aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 4.44% of total variance)", y = "PC3 (explains 3.76% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Ela PCA ####    
#noLD

#add columns for Location & Era
Ela_Ham_Ela_data <- Ela_Ham_Ela_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Nas") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population, "CNas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Ela_Ham_Ela_eigenval[1,1] / sum(Ela_Ham_Ela_eigenval$V1))*100
varPC2 <- (Ela_Ham_Ela_eigenval[2,1] / sum(Ela_Ham_Ela_eigenval$V1))*100
varPC3 <- (Ela_Ham_Ela_eigenval[3,1] / sum(Ela_Ham_Ela_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Ela_Ham_Ela_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#8aa9be", "#16537e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 7.19% of total variance)", y = "PC2 (explains 6.11% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Ela_Ham_Ela_data, aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 7.19% of total variance)", y = "PC3 (explains 5.94% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Ela_Ham_Ela_data, aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 6.11% of total variance)", y = "PC3 (explains 5.94% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Lle PCA ####    
#noLD

#add columns for Location & Era
Ela_Ham_Lle_data <- Ela_Ham_Lle_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Nas") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population, "CNas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Ela_Ham_Lle_eigenval[1,1] / sum(Ela_Ham_Lle_eigenval$V1))*100
varPC2 <- (Ela_Ham_Lle_eigenval[2,1] / sum(Ela_Ham_Lle_eigenval$V1))*100
varPC3 <- (Ela_Ham_Lle_eigenval[3,1] / sum(Ela_Ham_Lle_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Ela_Ham_Lle_data, aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC2") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 9.46% of total variance)", y = "PC2 (explains 8.48% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_12

PCA_13 <- ggplot(data = Ela_Ham_Lle_data, aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 9.46% of total variance)", y = "PC3 (explains 6.27% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Ela_Ham_Lle_data, aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 8.48% of total variance)", y = "PC3 (explains 6.27% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23

#### Ela nohighhet PCA ####    
#noLD

#add columns for Location & Era
Ela_Ham_Ela_nohighhet_data <- Ela_Ham_Ela_nohighhet_data %>%
  mutate(Location =
           case_when(endsWith(Population, "Ham") ~ "Hamilo",
                     endsWith(Population, "Nas") ~ "Hamilo")) %>% 
  mutate(Era =
           case_when(endsWith(Population, "AHam") ~ "Historical",
                     endsWith(Population, "CNas") ~ "Contemporary")) %>%
  relocate(Location, .before = Population) %>%
  relocate(Era, .before = Location)

#calculate % variance each PC explains by adding up all eigenvalues in *.eigenval file
varPC1 <- (Ela_Ham_Ela_nohighhet_eigenval[1,1] / sum(Ela_Ham_Ela_nohighhet_eigenval$V1))*100
varPC2 <- (Ela_Ham_Ela_nohighhet_eigenval[2,1] / sum(Ela_Ham_Ela_nohighhet_eigenval$V1))*100
varPC3 <- (Ela_Ham_Ela_nohighhet_eigenval[3,1] / sum(Ela_Ham_Ela_nohighhet_eigenval$V1))*100

#PCA
PCA_12 <- ggplot(data = Ela_Ham_Ela_nohighhet_data, 
                 aes(x = PC1, y = PC2, color = Era, shape = Era)) + 
  geom_point(size = 18, alpha = 0.5) +
  scale_color_manual(values = c("#8aa9be", "#16537e"), labels = c("Contemporary", "Historical")) + 
  scale_shape_manual(values = c(19, 15), labels = c("Contemporary", "Historical")) + 
  labs(x = "PC1 (explains 6.25% of total variance)", y = "PC2 (explains 6.11% of total variance)") + 
  scale_size(guide = "none") + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(linewidth = 4), 
        plot.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 4), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", vjust = 3),
        axis.title.x = element_text(size = 55, color = "black", vjust = -1), 
        legend.position = "none", 
        plot.margin = unit(c(0.5,1.5,1,1), "cm"))
PCA_12

PCA_13 <- ggplot(data = Ela_Ham_Ela_nohighhet_data, 
                 aes(x = PC1, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC1 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC1 (explains 6.25% of total variance)", y = "PC3 (explains 5.67% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_13

PCA_23 <- ggplot(data = Ela_Ham_Ela_nohighhet_data, 
                 aes(x = PC2, y = PC3, color = Era, shape = Era)) + 
  geom_point(size = 16) + ggtitle("PC2 v. PC3") + 
  scale_color_manual(values = c("#000000", "#999999"), labels = c("Historical", "Contemporary")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Historical", "Contemporary")) + 
  labs(x = "PC2 (explains 6.11% of total variance)", y = "PC3 (explains 5.67% of total variance)") + 
  scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = "center", axis.line = element_line(linewidth = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + 
  guides(color = guide_legend(override.aes = list(size = 16)))
PCA_23