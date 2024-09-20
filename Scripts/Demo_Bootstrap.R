######################## Script for calculating 95% CIs for momi models ########################################################

#reads in parameter estimates from momi SFS-bootstrap runs
#calculates 95% CIs
#creates plots of Ne change through time

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(here) #v.1.0.1
library(tidyverse) #v.2.0.0
library(scales)

#read in contemporary bootstrap data
Gmi_demo_contemp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                  "Gmi_contemp2changeexpg_bootstraps.csv"), 
                             header = FALSE)
  Gmi_demo_contemp <- as.data.frame(t(Gmi_demo_contemp))
    colnames(Gmi_demo_contemp) <- c("N_Hist", "N_Alb", "N_Cont", "T_Exp", "T_Bot")
    Gmi_demo_contemp <- Gmi_demo_contemp[-1,]
Ela_demo_contemp <- read.csv(here("Data/Ela_Ham/momi2", 
                                  "Ela_contemp2changeexpg_bootstraps.csv"), 
                             header = FALSE)
  Ela_demo_contemp <- as.data.frame(t(Ela_demo_temp))
    colnames(Ela_demo_contemp) <- c("N_Hist", "N_Alb", "N_Cont", "T_Exp", "T_Bot")
    Ela_demo_contemp <- Ela_demo_contemp[-1,]

#read in historical bootstrap data
Gmi_demo_temp <- read.csv(here("Data/Gmi_Ham/momi2", 
                               "Gmi_temponly2histchange_bootstraps.csv"), 
                          header = FALSE)
  Gmi_demo_temp <- as.data.frame(t(Gmi_demo_temp))
    colnames(Gmi_demo_temp) <- c("N_Hist", "N_Rec", "N_Cont", "T_Exp", "T_Rec")
    Gmi_demo_temp <- Gmi_demo_temp[-1,]
Ela_demo_temp <- read.csv(here("Data/Ela_Ham/momi2", 
                               "Ela_temponly2histchange_bootstraps.csv"), 
                          header = FALSE)
  Ela_demo_temp <- as.data.frame(t(Ela_demo_temp))
   colnames(Ela_demo_temp) <- c("N_Hist", "N_Rec", "N_Cont", "T_Exp", "T_Rec")
   Ela_demo_temp <- Ela_demo_temp[-1,]
  
#read in historical and contemporary bootstrap data
Gmi_demo_contemptemp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                      "Gmi_temp3change_bootstraps.csv"), 
                                 header = FALSE)
Gmi_demo_contemptemp <- as.data.frame(t(Gmi_demo_contemptemp))
  colnames(Gmi_demo_contemptemp) <- c("N_Hist", "N_Rec", "N_Alb", "N_Cont", "T_Exp", "T_Rec", "T_Bot")
  Gmi_demo_contemptemp <- Gmi_demo_contemptemp[-1,]
Ela_demo_contemptemp <- read.csv(here("Data/Ela_Ham/momi2", 
                                      "Ela_temp2histchange_bootstraps.csv"), 
                                 header = FALSE)
  Ela_demo_contemptemp <- as.data.frame(t(Ela_demo_contemptemp))
    colnames(Ela_demo_contemptemp) <- c("N_Hist", "N_Rec", "N_Cont", "T_Exp", "T_Rec")
    Ela_demo_contemptemp <- Ela_demo_contemptemp[-1,]

#####################################################################################################
  
######## Get 95% CIs ########
  
#### 95% CIs from contemporary data ####
#Gmi
Gmi_N_Hist <- quantile(as.numeric(Gmi_demo_contemp$N_Hist), c(0.025, 0.975))
Gmi_N_Alb <- quantile(as.numeric(Gmi_demo_contemp$N_Alb), c(0.025, 0.975))
Gmi_N_Cont <- quantile(as.numeric(Gmi_demo_contemp$N_Cont), c(0.025, 0.975))
Gmi_T_Exp <- quantile(as.numeric(Gmi_demo_contemp$T_Exp), c(0.025, 0.975))
Gmi_T_Bot <- quantile(as.numeric(Gmi_demo_contemp$T_Bot), c(0.025, 0.975))

#Ela
Ela_N_Hist <- quantile(as.numeric(Ela_demo_contemp$N_Hist), c(0.025, 0.975))
Ela_N_Alb <- quantile(as.numeric(Ela_demo_contemp$N_Alb), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_contemp$N_Cont), c(0.025, 0.975))
Ela_T_Exp <- quantile(as.numeric(Ela_demo_contemp$T_Exp), c(0.025, 0.975))
Ela_T_Bot <- quantile(as.numeric(Ela_demo_contemp$T_Bot), c(0.025, 0.975))

#### 95% CIs from historical data ####
#Gmi
Gmi_N_Hist <- quantile(as.numeric(Gmi_demo_temp$N_Hist), c(0.025, 0.975))
Gmi_N_Rec <- quantile(as.numeric(Gmi_demo_temp$N_Rec), c(0.025, 0.975))
Gmi_N_Cont <- quantile(as.numeric(Gmi_demo_temp$N_Cont), c(0.025, 0.975))
Gmi_T_Exp <- quantile(as.numeric(Gmi_demo_temp$T_Exp), c(0.025, 0.975))
Gmi_T_Rec <- quantile(as.numeric(Gmi_demo_temp$T_Rec), c(0.025, 0.975))

#Ela
Ela_N_Hist <- quantile(as.numeric(Ela_demo_temp$N_Hist), c(0.025, 0.975))
Ela_N_Rec <- quantile(as.numeric(Ela_demo_temp$N_Rec), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_temp$N_Cont), c(0.025, 0.975))
Ela_T_Exp <- quantile(as.numeric(Ela_demo_temp$T_Exp), c(0.025, 0.975))
Ela_T_Rec <- quantile(as.numeric(Ela_demo_temp$T_Rec), c(0.025, 0.975))

#### 95% CIs from contemporary and historical data ####
#Gmi
Gmi_N_Hist <- quantile(as.numeric(Gmi_demo_contemptemp$N_Hist), c(0.025, 0.975))
Gmi_N_Rec <- quantile(as.numeric(Gmi_demo_contemptemp$N_Rec), c(0.025, 0.975))
Gmi_N_Alb <- quantile(as.numeric(Gmi_demo_contemptemp$N_Alb), c(0.025, 0.975))
Gmi_N_Cont <- quantile(as.numeric(Gmi_demo_contemptemp$N_Cont), c(0.025, 0.975))
Gmi_T_Exp <- quantile(as.numeric(Gmi_demo_contemptemp$T_Exp), c(0.025, 0.975))
Gmi_T_Rec <- quantile(as.numeric(Gmi_demo_contemptemp$T_Rec), c(0.025, 0.975))
Gmi_T_Bot <- quantile(as.numeric(Gmi_demo_contemptemp$T_Bot), c(0.025, 0.975))

#Ela
Ela_N_Hist <- quantile(as.numeric(Ela_demo_contemptemp$N_Hist), c(0.025, 0.975))
Ela_N_Rec <- quantile(as.numeric(Ela_demo_contemptemp$N_Rec), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_contemptemp$N_Cont), c(0.025, 0.975))
Ela_T_Exp <- quantile(as.numeric(Ela_demo_contemptemp$T_Exp), c(0.025, 0.975))
Ela_T_Rec <- quantile(as.numeric(Ela_demo_contemptemp$T_Rec), c(0.025, 0.975))

################################################################################################

######## Make demographic bootstrap plots #########

#### Gmi plots #### 
## Contemporary plot ##
Gmi_Ne_boot_contemp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                     "gmi_momi2_contemp_boot_output_formatted.csv"), 
                                header = TRUE)

Gmi_Ne_contemp_boot_plot <- ggplot(data = Gmi_Ne_boot_contemp,
                           aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#1c3b0e")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Gmi_Ne_contemp_boot_plot

## Historical plot ##
Gmi_Ne_boot_temp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                  "gmi_momi2_temponly_boot_output_formatted.csv"), 
                             header = TRUE)

Gmi_Ne_temp_boot_plot <- ggplot(data = Gmi_Ne_boot_temp,
                           aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#1c3b0e")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Gmi_Ne_temp_boot_plot

## Contemporary & historical plot ##
Gmi_Ne_boot_contemp_temp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                          "gmi_momi2_temp_boot_output_formatted.csv"), 
                                     header = TRUE)

Gmi_Ne_contemp_temp_boot_plot <- ggplot(data = Gmi_Ne_boot_contemp_temp,
                                        aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#1c3b0e")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Gmi_Ne_contemp_temp_boot_plot

#### Ela plots ####
## Contemporary plot ##
Ela_Ne_boot_contemp <- read.csv(here("Data/Ela_Ham", 
                                     "ela_momi2_contemp_boot_output_formatted.csv"), 
                                header = TRUE)

Ela_Ne_contemp_boot_plot <- ggplot(data = Ela_Ne_boot_contemp,
                           aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#16537e")) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ela_Ne_contemp_boot_plot

## Historical plot ##
Ela_Ne_boot_temp <- read.csv(here("Data/Ela_Ham/momi2", 
                                  "ela_momi2_temponly_boot_output_formatted.csv"), 
                             header = TRUE)

Ela_Ne_temp_boot_plot <- ggplot(data = Ela_Ne_boot_temp,
                                aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#16537e")) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ela_Ne_temp_boot_plot

## Contemporary & historical plot ##
Ela_Ne_boot_contemp_temp <- read.csv(here("Data/Ela_Ham/momi2", 
                                          "ela_momi2_temp_boot_output_formatted.csv"), 
                                     header = TRUE)

Ela_Ne_contemp_temp_boot_plot <- ggplot(data = Ela_Ne_boot_contemp_temp,
                                        aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#BFC9CA", "#16537e")) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Ela_Ne_contemp_temp_boot_plot

#### Aen plots #### 

## Contemporary plot ##
Aen_Ne_boot_contemp <- read.csv(here("Data/Aen_Ham/momi2", 
                                     "aen_momi2_contemp_boot_output_formatted.csv"), 
                                header = TRUE)

Aen_Ne_contemp_boot_plot <- ggplot(data = Aen_Ne_boot_contemp,
                                   aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_Ne_contemp_boot_plot

## Temporal plot ##
Aen_Ne_boot_temp <- read.csv(here("Data/Aen_Ham/momi2", 
                                  "aen_momi2_temponly_boot_output_formatted.csv"), 
                                       header = TRUE)

Aen_Ne_temp_boot_plot <- ggplot(data = Aen_Ne_boot_temp,
                                aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_Ne_temp_boot_plot

## Contemporary & temporal plot ##
Aen_Ne_boot_contemp_temp <- read.csv(here("Data/Aen_Ham/momi2", 
                                          "aen_momi2_temp_boot_output_formatted.csv"), 
                                     header = TRUE)

Aen_Ne_contemp_temp_boot_plot <- ggplot(data = Aen_Ne_boot_contemp_temp,
                                                  aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_Ne_contemp_temp_boot_plot

#### Aen unrelated plots #### 

## Contemporary plot ##
Aen_unrelated_Ne_boot_contemp <- read.csv(here("Data/Aen_Ham/momi2", 
                                               "aen_unrelated_momi2_contemp_boot_output_formatted.csv"), 
                                          header = TRUE)

Aen_unrelated_Ne_contemp_boot_plot <- ggplot(data = Aen_unrelated_Ne_boot_contemp,
                                             aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_unrelated_Ne_contemp_boot_plot

## Temporal plot ##
Aen_unrelated_Ne_boot_temp <- read.csv(here("Data/Aen_Ham/momi2", 
                                            "aen_unrelated_momi2_temponly_boot_output_formatted.csv"), 
                                       header = TRUE)

Aen_unrelated_Ne_temp_boot_plot <- ggplot(data = Aen_unrelated_Ne_boot_temp,
                                          aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_unrelated_Ne_temp_boot_plot

## Contemporary & temporal plot ##
Aen_unrelated_Ne_boot_contemp_temp <- read.csv(here("Data/Aen_Ham/momi2", 
                                                    "aen_unrelated_momi2_temp_boot_output_formatted.csv"), 
                                               header = TRUE)

Aen_unrelated_Ne_contemp_temp_boot_plot <- ggplot(data = Aen_unrelated_Ne_boot_contemp_temp,
                                                  aes(x=year, y=Ne, group=run, color = type, size = type)) +
  geom_step() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_size_manual(values = c(1.5, 5)) +
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4), 
        plot.title = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 4), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 65, color = "black", vjust = 3),
        axis.title.x = element_text(size = 65, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,1,1), "cm"),)
Aen_unrelated_Ne_contemp_temp_boot_plot
