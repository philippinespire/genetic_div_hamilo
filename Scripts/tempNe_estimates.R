#################################################### Script for Ne Estimates  ###################################################################

#Estimate generation time for species based on taxonomy
#Calculates harmonic mean of Ne (over time) based on het loss equation

#################################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
#remotes::install_github("ropensci/rfishbase")
#remotes::install_github("James-Thorson/FishLife", force = TRUE)
library(FishLife) #v.3.0.0

################################################################################################################################################

######## Estimate generation time ########
#use FishLife to estimate generation time

## Gazza minuta ##
Gmi_params <- Plot_taxa(Search_species(Genus = "Gazza", Species = "minuta")$match_taxonomy[1]) #G = 1.39376929

## Equulites laterofenestra ##
Ela_params <- Plot_taxa(Search_species(Genus = "Equulites", Species = "laterofenestra")$match_taxonomy[1]) #G = 1.39353794

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
#t = # elapsed generations --> 1 generation = 1.394 years (estimated from FishLife), t = 110 years/1.394 years/gen = 78.9 generations

het_df$temp_Ne <- 1/(2*(1-((het_df$cont_He/het_df$hist_He)^(1/78.9))))

# calculate 95% confidence intervals #
#want largest and smallest differences
het_df$temp_2.5_Ne <- 1/(2*(1-((het_df$cont_2.5_He/het_df$hist_97.5_He)^(1/78.9)))) #"smallest" contemp and "biggest" historical -> less loss through time = higher Ne
het_df$temp_97.5_Ne <- 1/(2*(1-((het_df$cont_97.5_He/het_df$hist_2.5_He)^(1/78.9)))) #"biggest" contemp and "smallest" historical -> more loss through time = lower Ne
