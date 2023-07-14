########################## Script for calculating 95% CIs for momi models ########################################################

#reads in parameter estimates from momi SFS-bootstrap runs
#calculates 95% CIs

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(here) #v.1.0.1
library(tidyverse) #v.2.0.0

#read in bootstrap data
Aen_demo_temp <- read.csv(here("Data/Aen_Ham", "Aen_temporal2change_bootstraps.csv"), 
                          header = FALSE)
  Aen_demo_temp <- as.data.frame(t(Aen_demo_temp))
    colnames(Aen_demo_temp) <- c("N_Hist", "N_Alb", "N_Cont", "T_exp", "T_bot")
    Aen_demo_temp <- Aen_demo_temp[-1,]
Gmi_demo_temp <- read.csv(here("Data/Gmi_Ham", "Gmi_temporal2change_bootstraps.csv"), 
                          header = FALSE)
  Gmi_demo_temp <- as.data.frame(t(Gmi_demo_temp))
  colnames(Gmi_demo_temp) <- c("N_Hist", "N_Alb", "N_Cont", "T_exp", "T_bot")
  Gmi_demo_temp <- Gmi_demo_temp[-1,]
Ela_demo_temp <- read.csv(here("Data/Ela_Ham", "Ela_temporal2change_bootstraps.csv"), 
                          header = FALSE)
  Ela_demo_temp <- as.data.frame(t(Ela_demo_temp))
  colnames(Ela_demo_temp) <- c("N_Hist", "N_Alb", "N_Cont", "T_exp", "T_bot")
  Ela_demo_temp <- Ela_demo_temp[-1,]

#####################################################################################################
  
#### Get 95% CIs ####
  
#Aen
Aen_N_Hist <- quantile(as.numeric(Aen_demo_temp$N_Hist), c(0.025, 0.975))
Aen_N_Alb <- quantile(as.numeric(Aen_demo_temp$N_Alb), c(0.025, 0.975))
Aen_N_Cont <- quantile(as.numeric(Aen_demo_temp$N_Cont), c(0.025, 0.975))
Aen_T_exp <- quantile(as.numeric(Aen_demo_temp$T_exp), c(0.025, 0.975))
Aen_T_bot <- quantile(as.numeric(Aen_demo_temp$T_bot), c(0.025, 0.975))

#Gmi
Gmi_N_Hist <- quantile(as.numeric(Gmi_demo_temp$N_Hist), c(0.025, 0.975))
Gmi_N_Alb <- quantile(as.numeric(Gmi_demo_temp$N_Alb), c(0.025, 0.975))
Gmi_N_Cont <- quantile(as.numeric(Gmi_demo_temp$N_Cont), c(0.025, 0.975))
Gmi_T_exp <- quantile(as.numeric(Gmi_demo_temp$T_exp), c(0.025, 0.975))
Gmi_T_bot <- quantile(as.numeric(Gmi_demo_temp$T_bot), c(0.025, 0.975))

#Ela
Ela_N_Hist <- quantile(as.numeric(Ela_demo_temp$N_Hist), c(0.025, 0.975))
Ela_N_Alb <- quantile(as.numeric(Ela_demo_temp$N_Alb), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_temp$N_Cont), c(0.025, 0.975))
Ela_T_exp <- quantile(as.numeric(Ela_demo_temp$T_exp), c(0.025, 0.975))
Ela_T_bot <- quantile(as.numeric(Ela_demo_temp$T_bot), c(0.025, 0.975))
