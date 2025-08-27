######################## Script for calculating 95% CIs for momi models ########################################################

#Reads in parameter estimates from momi SFS-bootstrap runs
#Calculates 95% CIs
#Creates plots of Ne change through time (size: 2500 x 1500)

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(here) #v.1.0.1
library(tidyverse) #v.2.0.0
library(scales) #v.1.2.1

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
  Ela_demo_contemp <- as.data.frame(t(Ela_demo_contemp))
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
                               "Ela_temponlyhistchange_bootstraps.csv"), 
                          header = FALSE)
  Ela_demo_temp <- as.data.frame(t(Ela_demo_temp))
   colnames(Ela_demo_temp) <- c("N_Hist", "N_Cont", "T_Exp")
   Ela_demo_temp <- Ela_demo_temp[-1,]
  
#read in historical and contemporary bootstrap data
Gmi_demo_contemptemp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                      "Gmi_temp3change_bootstraps.csv"), 
                                 header = FALSE)
Gmi_demo_contemptemp <- as.data.frame(t(Gmi_demo_contemptemp))
  colnames(Gmi_demo_contemptemp) <- c("N_Hist", "N_Rec", "N_Alb", "N_Cont", "T_Exp", "T_Rec", "T_Bot")
  Gmi_demo_contemptemp <- Gmi_demo_contemptemp[-1,]

Ela_demo_contemptemp <- read.csv(here("Data/Ela_Ham/momi2", 
                                      "Ela_temp3change_bootstraps.csv"), 
                                 header = FALSE)
  Ela_demo_contemptemp <- as.data.frame(t(Ela_demo_contemptemp))
    colnames(Ela_demo_contemptemp) <- c("N_Hist", "N_Rec", "N_Alb", "N_Cont", "T_Exp", "T_Rec", "T_Bot")
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
Gmi_T_Exp <- quantile(as.numeric(Gmi_demo_temp$T_Exp), c(0.025, 0.975)) + 110 #to get in actual YBP, these models end at 1908
Gmi_T_Rec <- quantile(as.numeric(Gmi_demo_temp$T_Rec), c(0.025, 0.975)) + 110

#Ela
Ela_N_Hist <- quantile(as.numeric(Ela_demo_temp$N_Hist), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_temp$N_Cont), c(0.025, 0.975))
Ela_T_Exp <- quantile(as.numeric(Ela_demo_temp$T_Exp), c(0.025, 0.975)) + 110 #to get in actual YBP, these models end at 1908

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
Ela_N_Alb <- quantile(as.numeric(Ela_demo_contemptemp$N_Alb), c(0.025, 0.975))
Ela_N_Cont <- quantile(as.numeric(Ela_demo_contemptemp$N_Cont), c(0.025, 0.975))
Ela_T_Exp <- quantile(as.numeric(Ela_demo_contemptemp$T_Exp), c(0.025, 0.975))
Ela_T_Rec <- quantile(as.numeric(Ela_demo_contemptemp$T_Rec), c(0.025, 0.975))
Ela_T_Bot <- quantile(as.numeric(Ela_demo_contemptemp$T_Bot), c(0.025, 0.975))

################################################################################################

######## Make demographic bootstrap plots #########
#for gradual (exponential) change models, need to predict values and plot
#for instantaneous change models, will use "formatted" files to plot
#for gradual change models, will also calculate Ne at 1908 (and 95% CIs)

######## Gmi plots ######## 

#### Contemporary plot ####
#exponential change model is best fit

## add empirical estimates to df for bootstrapping ##
Gmi_empirical_df <- data.frame(N_Hist = 63954.62473358225,
                               N_Alb = 20046711.84161918,
                               N_Cont = 60.349407868520636, 
                               T_Exp = 191147.55090317182, 
                               T_Bot = 74.85833552215723)
Gmi_demo_contemp_boot <- rbind(Gmi_demo_contemp, Gmi_empirical_df)

## predict values based on models ##
#using model structure to predict Ne back 1 MYA

Gmi_contemp_predict_df <- data.frame() #empty df to fill with predicted values
Gmi_contemp_Ne_1908_df <- data.frame() #empty df to fill with predicted 1908 values

#for loop over bootstrap runs
for (i in 1:nrow(Gmi_demo_contemp_boot)) {
  
  #track progress
  pct <- round(i/101 * 100)
  cat("Progress:", pct, "%\n")
  
  #assign variables
  run <- Gmi_demo_contemp_boot[i, ]
  
  n_contemp <- as.numeric(run$N_Cont)
  n_bot <- as.numeric(run$N_Alb)
  n_ancient <- as.numeric(run$N_Hist)
  t_onset_bot <- as.numeric(run$T_Bot)
  t_onset_ancient <- as.numeric(run$T_Exp)
  
  #calculate variables
  t_ancient <- t_onset_ancient - t_onset_bot
  r_bot <- log(n_contemp/n_bot)/t_onset_bot
  r_anc <- log(n_bot/n_ancient)/t_ancient
  lambda_bot <- exp(r_bot)
  lambda_anc <- exp(r_anc)
  
  #create vector to fill with estimated Nes
  pop_size <- c()
  pop_size[1] <- n_ancient #seed with n_ancient
  
  #estimate Ne from onset of ancient pop growth to start of recent bottleneck
  for(x in 2:(t_ancient)) {
    pop_size[x] <- pop_size[x-1]*lambda_anc
  }
  
  #force Ne at start of recent bottleneck to equal n_bot (should be, but make sure)
  pop_size[t_ancient] <- n_bot #n_bot
  
  #estimate Ne from onset of recent bottleneck to contemp
  for(x in (t_ancient + 1):(t_ancient + t_onset_bot)){
    pop_size[x] <- pop_size[x-1]*lambda_bot
  }
  
  #force Ne at contemp to equal n_contemp (should be, but make sure)
  pop_size[t_ancient + t_onset_bot] <- n_contemp
  
  #add time
  pop_size_df <- as.data.frame(pop_size)
    pop_size_df$time <- nrow(pop_size_df):1 #going back in time
  
  #grab Ne @ 1908 (t = 110)
  Ne_1908_run <- subset(pop_size_df, pop_size_df$time == 110)
    Ne_1908_run$run <- i
    Gmi_contemp_Ne_1908_df <- rbind(Gmi_contemp_Ne_1908_df, Ne_1908_run)
    
  #subset to make plotting easier
  pop_size_df_sub <- pop_size_df[seq(1, (nrow(pop_size_df)-100), by = 5), ] #subset by 5
    pop_size_df_contemp <- pop_size_df[(nrow(pop_size_df)-100):nrow(pop_size_df), ] #bc will plot on log scale, get every point from 1:100
    pop_size_df <- rbind(pop_size_df_sub, pop_size_df_contemp)
    
  #add in starting Ne back to 1 MYA
  pop_size_start_df <- data.frame(pop_size = n_ancient, time = 1000000) #just need 1 point to get straight line
    pop_size_all_df <- rbind(pop_size_start_df, pop_size_df) #combine together
    pop_size_all_df$run <- i #add in row
  
  Gmi_contemp_predict_df <- rbind(Gmi_contemp_predict_df, pop_size_all_df) #append to bootstrap dataset
}

## clean up bootstrap_df ##
colnames(Gmi_contemp_predict_df) <- c("pop_size", "time", "run")
Gmi_contemp_predict_df$pop_size[Gmi_contemp_predict_df$pop_size <= 1] <- 0 #ensure no negative Ne values (shouldn't be)

#specify which is from "real SFS" and which is from bootstrapped SFS
Gmi_contemp_predict_df$type[Gmi_contemp_predict_df$run == 101] <- "real"
  Gmi_contemp_predict_df$type[Gmi_contemp_predict_df$run != 101] <- "bootstrap"

## calculate 1908 Ne ##
#Ne from observed SFS
Gmi_contemp_Ne_1908 <- Gmi_contemp_Ne_1908_df[101, 1]

#95% CI from bootstrapped runs
Gmi_contemp_Ne_1908_boot <- subset(Gmi_contemp_Ne_1908_df, Gmi_contemp_Ne_1908_df$run != 101)
  Gmi_contemp_Ne_1908_CI <- quantile(Gmi_contemp_Ne_1908_boot$pop_size, c(0.025, 0.975))
  
## plot ##
Gmi_Ne_contemp_boot_plot <- ggplot() + 
  geom_line(data = Gmi_contemp_predict_df[Gmi_contemp_predict_df$type == "bootstrap", ], 
            aes(x = time, y = pop_size, group = run), 
                color = "#afc8a4", linewidth = 3, alpha = 0.3) +
  geom_line(data = Gmi_contemp_predict_df[Gmi_contemp_predict_df$type == "real", ], 
            aes(x = time, y = pop_size, group = run), 
                color = "#1c3b0e", linewidth = 5) +
  geom_vline(xintercept = 110, linetype = "dotted", color = "black", linewidth = 3) + #for historical time point, in years before present
  scale_x_log10(limits = c(1, 1e6), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Gmi_Ne_contemp_boot_plot

#### Historical plot ####
#instantaneous change model is best fit

## read in data ##
Gmi_Ne_boot_temp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                  "gmi_momi2_temponly_boot_output_formatted.csv"), 
                             header = TRUE)
  colnames(Gmi_Ne_boot_temp) <- c("type", "run", "time", "pop_size", "year_meaning", "Ne_meaning")
  Gmi_Ne_boot_temp$adjust_time <- Gmi_Ne_boot_temp$time + 110 #adjusting to make same scale as with contemp data
    Gmi_Ne_boot_temp$adjust_time[Gmi_Ne_boot_temp$adjust_time > 1000000] <- 1000000 #fix 1 MYA time point

## plot ##
Gmi_Ne_temp_boot_plot <- ggplot() + 
  geom_line(data = Gmi_Ne_boot_temp[ Gmi_Ne_boot_temp$type == "bootstrap", ], 
            aes(x = adjust_time, y = pop_size, group = run), 
            color = "#afc8a4", linewidth = 3, alpha = 0.3) +
  geom_line(data = Gmi_Ne_boot_temp[Gmi_Ne_boot_temp$type == "real", ], 
            aes(x = adjust_time, y = pop_size, group = run), 
            color = "#1c3b0e", linewidth = 5) +
  scale_x_log10(limits = c(110, 1e6), #want to start where actually have data (110)
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Gmi_Ne_temp_boot_plot

#### Contemporary and Historical plot ####
#instantaneous change model is best fit

## read in data ##
Gmi_Ne_boot_contemp_temp <- read.csv(here("Data/Gmi_Ham/momi2", 
                                          "gmi_momi2_temp_boot_output_formatted.csv"), 
                                     header = TRUE)
  colnames(Gmi_Ne_boot_contemp_temp) <- c("type", "run", "time", "pop_size", "year_meaning", "Ne_meaning")

## plot ##
Gmi_Ne_contemp_temp_boot_plot <- ggplot() + 
  geom_step(data = Gmi_Ne_boot_contemp_temp[Gmi_Ne_boot_contemp_temp$type == "bootstrap", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#afc8a4", linewidth = 3, alpha = 0.3) +
  geom_step(data = Gmi_Ne_boot_contemp_temp[Gmi_Ne_boot_contemp_temp$type == "real", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#1c3b0e", linewidth = 5) +
  geom_vline(xintercept = 110, linetype = "dotted", color = "black", linewidth = 3) + #for historical time point, in years before present
  scale_x_log10(limits = c(1, 1e6), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Gmi_Ne_contemp_temp_boot_plot

######## Ela plots ######## 

#### Contemporary plot ####
#exponential change model is best fit

## add empirical estimates to df for bootstrapping ##
Ela_empirical_df <- data.frame(N_Hist = 678667.77,
                               N_Alb = 42360343.53783952,
                               N_Cont = 10.333042154447817, 
                               T_Exp = 804608.0920460667, 
                               T_Bot = 16.80318729771402)
Ela_demo_contemp_boot <- rbind(Ela_demo_contemp, Ela_empirical_df)

## predict values based on models ##
#using model structure to predict Ne back 1 MYA

Ela_contemp_predict_df <- data.frame() #empty df to fill with predicted values
Ela_contemp_Ne_1908_df <- data.frame() #empty df to fill with predicted 1908 values

#for loop over bootstrap runs
for (i in 1:nrow(Ela_demo_contemp_boot)) {
  
  #track progress
  pct <- round(i/101 * 100)
  cat("Progress:", pct, "%\n")
  
  #assign variables
  run <- Ela_demo_contemp_boot[i, ]
  
  n_contemp <- as.numeric(run$N_Cont)
  n_bot <- as.numeric(run$N_Alb)
  n_ancient <- as.numeric(run$N_Hist)
  t_onset_bot <- as.numeric(run$T_Bot)
  t_onset_ancient <- as.numeric(run$T_Exp)
  
  #calculate variables
  t_ancient <- t_onset_ancient - t_onset_bot
  r_bot <- log(n_contemp/n_bot)/t_onset_bot
  r_anc <- log(n_bot/n_ancient)/t_ancient
  lambda_bot <- exp(r_bot)
  lambda_anc <- exp(r_anc)
  
  #create vector to fill with estimated Nes
  pop_size <- c()
  pop_size[1] <- n_ancient #seed with n_ancient
  
  #estimate Ne from onset of ancient pop growth to start of recent bottleneck
  for(x in 2:(t_ancient)) {
    pop_size[x] <- pop_size[x-1]*lambda_anc
  }
  
  #force Ne at start of recent bottleneck to equal n_bot (should be, but make sure)
  pop_size[t_ancient] <- n_bot #n_bot
  
  #estimate Ne from onset of recent bottleneck to contemp
  for(x in (t_ancient + 1):(t_ancient + t_onset_bot)){
    pop_size[x] <- pop_size[x-1]*lambda_bot
  }
  
  #force Ne at contemp to equal n_contemp (should be, but make sure)
  pop_size[t_ancient + t_onset_bot] <- n_contemp
  
  #add time
  pop_size_df <- as.data.frame(pop_size)
  pop_size_df$time <- nrow(pop_size_df):1 #going back in time
  
  #grab Ne @ 1908 (t = 110)
  Ne_1908_run <- subset(pop_size_df, pop_size_df$time == 110)
  Ne_1908_run$run <- i
  Ela_contemp_Ne_1908_df <- rbind(Ela_contemp_Ne_1908_df, Ne_1908_run)
  
  #subset to make plotting easier
  pop_size_df_sub <- pop_size_df[seq(1, (nrow(pop_size_df)-100), by = 5), ] #subset by 5
  pop_size_df_contemp <- pop_size_df[(nrow(pop_size_df)-100):nrow(pop_size_df), ] #bc will plot on log scale, get every point from 1:100
  pop_size_df <- rbind(pop_size_df_sub, pop_size_df_contemp)
  
  #add in starting Ne back to 1 MYA
  pop_size_start_df <- data.frame(pop_size = n_ancient, time = 1000000) #just need 1 point to get straight line
  pop_size_all_df <- rbind(pop_size_start_df, pop_size_df) #combine together
  pop_size_all_df$run <- i #add in row
  
  Ela_contemp_predict_df <- rbind(Ela_contemp_predict_df, pop_size_all_df) #append to bootstrap dataset
}

## clean up bootstrap_df ##
colnames(Ela_contemp_predict_df) <- c("pop_size", "time", "run")
Ela_contemp_predict_df$pop_size[Ela_contemp_predict_df$pop_size <= 1] <- 0 #ensure no negative Ne values (shouldn't be)

#specify which is from "real SFS" and which is from bootstrapped SFS
Ela_contemp_predict_df$type[Ela_contemp_predict_df$run == 101] <- "real"
Ela_contemp_predict_df$type[Ela_contemp_predict_df$run != 101] <- "bootstrap"

## calculate 1908 Ne ##
#Ne from observed SFS
Ela_contemp_Ne_1908 <- Ela_contemp_Ne_1908_df[101, 1]

#95% CI from bootstrapped runs
Ela_contemp_Ne_1908_boot <- subset(Ela_contemp_Ne_1908_df, Ela_contemp_Ne_1908_df$run != 101)
Ela_contemp_Ne_1908_CI <- quantile(Ela_contemp_Ne_1908_boot$pop_size, c(0.025, 0.975))

## plot ##
Ela_Ne_contemp_boot_plot <- ggplot() + 
  geom_line(data = Ela_contemp_predict_df[Ela_contemp_predict_df$type == "bootstrap", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#8aa9be", linewidth = 3, alpha = 0.3) +
  geom_line(data = Ela_contemp_predict_df[Ela_contemp_predict_df$type == "real", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#16537e", linewidth = 5) +
  geom_vline(xintercept = 110, linetype = "dotted", color = "black", linewidth = 3) + #for historical time point, in years before present
  scale_x_log10(limits = c(1, 1e6), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Ela_Ne_contemp_boot_plot

#### Historical plot ####
#instantaneous change model is best fit

## read in data ##
Ela_Ne_boot_temp <- read.csv(here("Data/Ela_Ham/momi2", 
                                  "ela_momi2_temponly_boot_output_formatted.csv"), 
                             header = TRUE)
colnames(Ela_Ne_boot_temp) <- c("type", "run", "time", "pop_size", "year_meaning", "Ne_meaning")
Ela_Ne_boot_temp$adjust_time <- Ela_Ne_boot_temp$time + 110 #adjusting to make same scale as with contemp data
Ela_Ne_boot_temp$adjust_time[Ela_Ne_boot_temp$adjust_time > 1000000] <- 1000000 #fix 1 MYA time point

## plot ##
Ela_Ne_temp_boot_plot <- ggplot() + 
  geom_line(data = Ela_Ne_boot_temp[Ela_Ne_boot_temp$type == "bootstrap", ], 
            aes(x = adjust_time, y = pop_size, group = run), 
            color = "#8aa9be", linewidth = 3, alpha = 0.3) +
  geom_line(data = Ela_Ne_boot_temp[Ela_Ne_boot_temp$type == "real", ], 
            aes(x = adjust_time, y = pop_size, group = run), 
            color = "#16537e", linewidth = 5) +
  scale_x_log10(limits = c(110, 1e6), #want to start where actually have data (110)
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Ela_Ne_temp_boot_plot

#### Contemporary and Historical plot ####
#instantaneous change model is best fit

## read in data ##
Ela_Ne_boot_contemp_temp <- read.csv(here("Data/Ela_Ham/momi2", 
                                          "ela_momi2_temp_boot_output_formatted.csv"), 
                                     header = TRUE)
colnames(Ela_Ne_boot_contemp_temp) <- c("type", "run", "time", "pop_size", "year_meaning", "Ne_meaning")

## plot ##
Ela_Ne_contemp_temp_boot_plot <- ggplot() + 
  geom_step(data = Ela_Ne_boot_contemp_temp[Ela_Ne_boot_contemp_temp$type == "bootstrap", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#8aa9be", linewidth = 3, alpha = 0.3) +
  geom_step(data = Ela_Ne_boot_contemp_temp[Ela_Ne_boot_contemp_temp$type == "real", ], 
            aes(x = time, y = pop_size, group = run), 
            color = "#16537e", linewidth = 5) +
  geom_vline(xintercept = 110, linetype = "dotted", color = "black", linewidth = 3) + #for historical time point, in years before present
  scale_x_log10(limits = c(1, 1e6), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 1e11), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0)) + 
  ylab(bquote(N[e])) + 
  xlab("Years before present") +
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 4),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 55, color = "black", margin = margin(r = 20)), 
        axis.text.x = element_text(size = 55, color = "black", margin = margin(t = 20)), 
        axis.title.y = element_text(size = 55, color = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 55, color = "black", margin = margin(t = 30)),
        legend.position = "none",
        plot.margin = unit(c(2,2,1,1), "cm"),)
Ela_Ne_contemp_temp_boot_plot
