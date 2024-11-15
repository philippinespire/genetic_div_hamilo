#################################################### Script for Ne Estimates  ###################################################################

#Calculates harmonic mean of Ne (over time) based on het loss equation

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

################################################################################################################################################

######## Ne estimates ########

## create heterozygosity dataframe ##
hist_He <- c(0.10218, 0.12109) #Gmi, Ela
cont_He <- c(0.09606, 0.11408)
hist_2.5_He <- c(0.09904, 0.11965)
cont_2.5_He <- c(0.09305, 0.11285)
hist_97.5_He <- c(0.10537, 0.12249)
cont_97.5_He <- c(0.09896, 0.11531)

het_df <- as.data.frame(cbind(hist_He, cont_He, hist_2.5_He, cont_2.5_He, hist_97.5_He, cont_97.5_He))

## calculate Ne ##
#het loss equation: Ht/Ho = (1 - 1/(2Ne))^t
#Ne = 1/(2*(1-((Ht/Ho)^(1/t))))
#t = # elapsed generations --> 1 generation = 1 year, t = 110

het_df$temp_Ne <- 1/(2*(1-((het_df$cont_He/het_df$hist_He)^(1/110))))

# calculate 95% confidence intervals #
#want largest and smallest differences
het_df$temp_2.5_Ne <- 1/(2*(1-((het_df$cont_2.5_He/het_df$hist_97.5_He)^(1/110)))) #"smallest" contemp and "biggest" historical -> less loss through time = higher Ne
het_df$temp_97.5_Ne <- 1/(2*(1-((het_df$cont_97.5_He/het_df$hist_2.5_He)^(1/110)))) #"biggest" contemp and "smallest" historical -> more loss through time = lower Ne
