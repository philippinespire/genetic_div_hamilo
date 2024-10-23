####################################################### Script for Diploid Filtering ##############################################################

#Modified from Chris's version (https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/indvAlleleBalance.R)
#Calculates read depth v. heterozygosity stats to look for issues with reference bias and under-powered genotype calls

###################################################################################################################################################

######## Set-up ########

#set working directory
getwd()

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(magrittr) #v.2.0.3
library(janitor) #v.2.2.0
library(data.table) #v.1.14.8

#set user defined variables
alleledepthFILE = 'Data/Gmi_Ham/Gmi.A.rad.RAW-10-10-rescaled.Fltr17.2.recode.renamed.AD.tsv'
genoFILE = 'Data/Gmi_Ham/Gmi.A.rad.RAW-10-10-rescaled.Fltr17.2.recode.renamed.GT.tsv'
indvPATTERN = "gmi_"
modIdPATTERN = "_gm_.*"
het_cutoff_pop1 = 0.2125
het_cutoff_pop3 = 0.275

#calculate variables
outDIR=dirname(genoFILE)
outFilePREFIX <-
  basename(genoFILE) %>%
  str_remove("GT.tsv")

#define functions
getMode <- function(x) {
  keys <- 0:100 / 100
  keys[which.max(tabulate(match(round(x,
                                      2), 
                                keys)))]
}

#read in data
allele_depths <-
  read_tsv(alleledepthFILE, 
           col_types=cols(.default = "c")) %>%
  clean_names() %>%
  pivot_longer(cols=contains(indvPATTERN),
               names_to="id") %>%
  separate(col=value,
           into=(c("num_reads_ref",
                   "num_reads_alt")),
           convert=TRUE) %>%
  mutate(num_reads = num_reads_alt + num_reads_ref)

genotypes <-
  read_tsv(genoFILE,
           col_types=cols(.default = "c")) %>%
  clean_names() %>%
  pivot_longer(cols=contains(indvPATTERN),
               names_to="id") %>%
  separate(col=value,
           into=(c("state_1",
                   "state_2"))) %>%
 mutate(state_1 = na_if(state_1,
                         ""),
         state_2 = na_if(state_2,
                         ""),
         genotype = case_when(state_1 != state_2 ~ "hetero",
                              is.na(state_1) & is.na(state_2) ~ NA_character_,
                              state_1 == ref ~ "homo_ref",
                              state_1 == alt ~ "homo_alt",
                              TRUE ~ "error"))

#######################################################################################################################################################

######### Calculate heterozygosity by individual & locus ########

#combine data
num_loci <-
  genotypes %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  pull(chrom_pos) %>%
  unique() %>%
  length()

all_data <-
  allele_depths %>%
  left_join(genotypes) %>%
  left_join(genotypes %>%
              group_by(id,
                       genotype) %>%
              summarize(num_pos = n()) %>%
              pivot_wider(names_from = "genotype",
                          values_from = "num_pos") %>%
              mutate(across(hetero:homo_alt,
                            ~replace_na(.,
                                        0)),
                     heterozygosity_obs_ind = hetero/(num_loci - `NA`),
                     pop = case_when(heterozygosity_obs_ind >= het_cutoff_pop3 ~ 3,
                                     heterozygosity_obs_ind >= het_cutoff_pop1 ~ 1,
                                     TRUE ~ 2)) %>%
              select(id,
                     heterozygosity_obs_ind,
                     pop)) %>%
  separate(id,
           into=(c(NA,
                   "era",
                   NA,
                   NA,
                   NA,
                   NA)),
           remove = FALSE) %>%
  left_join(genotypes %>%
              mutate(chrom_pos = str_c(chrom,
                                       pos,
                                       sep="_")) %>%
              separate(id,
                       into=(c(NA,
                               "era",
                               NA,
                               NA,
                               NA,
                               NA)),
                       remove = FALSE) %>%
              group_by(chrom,
                       pos,
                       chrom_pos,
                       genotype,
                       era) %>%
              summarize(num_pos = n()) %>%
              pivot_wider(names_from = "genotype",
                          values_from = "num_pos") %>% 
              mutate(across(hetero:homo_alt,
                            ~replace_na(.,
                                        0)),
                     heterozygosity_obs_locus = hetero/(hetero + homo_ref + homo_alt)) %>% 
              select(chrom,
                     pos,
                     chrom_pos,
                     era,
                     heterozygosity_obs_locus)) %>%
  mutate(qual = as.double(qual),
         across(contains(c("num_reads",
                           "like",
                           "pos")),
                ~as.numeric(.)),
         num_reads = num_reads_alt + num_reads_ref,
         allele_balance = case_when(genotype == "hetero" ~ num_reads_alt / num_reads,
                                    TRUE ~ NA_real_),
         id = str_remove(id,
                         modIdPATTERN),
         chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) 

#save as RDS
saveRDS(all_data, 
       file=str_c(outDIR,
                  "/",
                  "all_data_diploid.rds",
                  sep=""))

####################################################################################################################

######## Investigate missing data patterns ########

#read in if running separately
all_data <- readRDS("Data/Gmi_Ham/all_data_diploid.rds")

#make column that indicates whether genotype is called or not
all_data <- all_data %>% 
  mutate(genotype_call = ifelse(is.na(all_data$genotype), 0, 1)) #if missing, code as 0, otherwise code as 1

#remove basud river individuals (only for Gmi)
all_data <- all_data[!grepl("^gmi_a_bas_*", all_data$id) & 
                       !grepl("^gmi_c_bas_*", all_data$id), ]

all_data_order <- all_data[order(all_data$id), ] #order by individual id

#plot missing data as heatmap
test <- ggplot(all_data_order, 
               aes(x = id, y = chrom_pos, 
                   fill = era, alpha = genotype_call)) + 
  geom_tile() + 
  scale_alpha_identity(guide = "none") +
  scale_fill_manual(values = c("#1c3b0e","#afc8a4")) + 
  xlab("Individual") + ylab("Locus") + 
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_text(size = 55, color = "black"),
        legend.position = "none")
test

#####################################################################################################################

######## Investigate read depth - het relationship ########
#two-fold reason for this:
#1. want to make sure lower read depth in historical data doesn't artificially reduce heterozygosity
#2. look for issues with reference/mapping bias (do homozygous calls have lower read depth than het -- does this vary by era)

#read in if running separately
all_data <- readRDS("Data/Gmi_Ham/all_data_diploid.rds")
  all_data$era[all_data$era == "c"] <- "Contemp"
  all_data$era[all_data$era == "a"] <- "Hist"

#### Calculate average AB per locus by locus, genotype, era ####
#calculate average AB by locus, era, genotype
avg_AB_by_locus_era_genotype <-
  all_data %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(mean_AB <- mean(allele_balance, na.rm = TRUE))

avg_AB_by_locus_era_genotype <- avg_AB_by_locus_era_genotype[!is.na(avg_AB_by_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads
avg_AB_by_locus_era_genotype <- subset(avg_AB_by_locus_era_genotype, avg_AB_by_locus_era_genotype$het_homo == "hetero")

avg_AB_by_locus_era_hetero <- avg_AB_by_locus_era_genotype[, -2]  

colnames(avg_AB_by_locus_era_hetero) <- c("chrom_pos", "era", "mean_AB")

#### Calculate average depth per locus by locus, genotype, era ####
#calculate average read depth by locus, era, genotype
avg_read_depth_by_locus_era_genotype <-
  all_data %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(mean_depth <- mean(num_reads, na.rm = TRUE))

avg_read_depth_by_locus_era_genotype <- avg_read_depth_by_locus_era_genotype[!is.na(avg_read_depth_by_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads

colnames(avg_read_depth_by_locus_era_genotype) <- c("chrom_pos", "het_homo", "era", "mean_depth")

#pull out by genotype
avg_read_depth_by_locus_era_genotype_het <- subset(avg_read_depth_by_locus_era_genotype, 
                                                   avg_read_depth_by_locus_era_genotype$het_homo == "hetero")
avg_read_depth_by_locus_era_genotype_homo <- subset(avg_read_depth_by_locus_era_genotype, 
                                                    avg_read_depth_by_locus_era_genotype$het_homo == "homo") 

#merge dataframes
avg_read_depth_combined <- merge(avg_read_depth_by_locus_era_genotype_het, avg_read_depth_by_locus_era_genotype_homo, 
                                 by = c("chrom_pos", "era")) #1175 loci that only have homozygous calls, no heterozygous calls -> not quite right, bc double counting across eras
  colnames(avg_read_depth_combined) <- c("chrom_pos", "era", "het", "mean_depth_het", "homo", "mean_depth_homo")
  avg_read_depth_combined <- avg_read_depth_combined[, -3] #removing het & homo columns (not needed)
  avg_read_depth_combined <- avg_read_depth_combined[, -4]
  
#calculate read ratio
avg_read_depth_combined$read_ratio <- avg_read_depth_combined$mean_depth_homo/avg_read_depth_combined$mean_depth_het

#### Sum number of reads by locus, era ####
#sum reads by locus, era
tot_num_reads_by_locus_era <-
  all_data %>%
  group_by(chrom_pos,
           era) %>%
  summarize(tot_reads <- sum(num_reads, na.rm = TRUE))

colnames(tot_num_reads_by_locus_era) <- c("chrom_pos", "era", "tot_depth")

#count by locus, era, genotype
count_locus_era_genotype <-
  all_data %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(n())

count_locus_era_genotype <- count_locus_era_genotype[!is.na(count_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads

colnames(count_locus_era_genotype) <- c("chrom_pos", "het_homo", "era", "num_genotypes")

#count by locus, era
count_locus_era <-
  count_locus_era_genotype %>%
  group_by(chrom_pos,
           era) %>%
  summarize(tot_genotypes <- sum(num_genotypes, na.rm = TRUE))

colnames(count_locus_era) <- c("chrom_pos", "era", "tot_genotypes")

## merge data frames ##
full_count <- merge(tot_num_reads_by_locus_era, count_locus_era_genotype, by = c("chrom_pos", "era"))
full_count_wtot <- merge(full_count, count_locus_era, by = c("chrom_pos", "era"))

#calculate % het loci per locus
full_count_wtot$perc_hethomo <- (full_count_wtot$num_genotypes/full_count_wtot$tot_genotypes) * 100

#### Calculate number of loci by depth, era, genotype ####
#count loci by depth, era, genotype
num_reads_by_depth_era <-
  all_data %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(num_reads,
           het_homo,
           era) %>%
  summarize(n())

num_reads_by_depth_era <- num_reads_by_depth_era[!is.na(num_reads_by_depth_era$het_homo), ] #just removing NAs -- most of these don't have reads

## Albatross reads ##
#identify only albatross het reads
a_sum_het <- subset(num_reads_by_depth_era, 
                    num_reads_by_depth_era$era == "Hist" & num_reads_by_depth_era$het_homo == "hetero")
  colnames(a_sum_het) <- c("num_reads", "het_homo", "era", "het_genotypes") #this count is now # het genotypes at X read depth across all alb individuals
  a_sum_het <- a_sum_het[, -2] #remove extra column

#identify only albatross homo reads
a_sum_homo <- subset(num_reads_by_depth_era, 
                     num_reads_by_depth_era$era == "Hist" & num_reads_by_depth_era$het_homo == "homo")
  colnames(a_sum_homo) <- c("num_reads", "het_homo", "era", "homo_genotypes") #this count is now # homo genotypes at X read depth across all alb individuals
  a_sum_homo <- a_sum_homo[, -2] #remove extra column

#create missing rows for merging
homo_nohet <- setdiff(a_sum_homo$num_reads, a_sum_het$num_reads) #in homo, not het
  het_missing <- data.frame(num_reads = homo_nohet, era = "Hist", het_genotypes = 0)
het_nohomo <- setdiff(a_sum_het$num_reads, a_sum_homo$num_reads) #in het, not homo
  homo_missing <- data.frame(num_reads = het_nohomo, era = "Hist", homo_genotypes = 0)

#merge albatross dataframes
a_sum_het_full <- rbind(a_sum_het, het_missing)
a_sum_homo_full <- rbind(a_sum_homo, homo_missing)

a_genotypes <- merge(a_sum_het_full, a_sum_homo_full, by = c("num_reads", "era"))

## Contemporary reads ##
#identify only contemp het reads
c_sum_het <- subset(num_reads_by_depth_era, 
                    num_reads_by_depth_era$era == "Contemp" & num_reads_by_depth_era$het_homo == "hetero")
  colnames(c_sum_het) <- c("num_reads", "het_homo", "era", "het_genotypes") #this count is now # het genotypes at X read depth across all contemp individuals
  c_sum_het <- c_sum_het[, -2] #remove extra column

c_sum_homo <- subset(num_reads_by_depth_era, 
                     num_reads_by_depth_era$era == "Contemp" & num_reads_by_depth_era$het_homo == "homo")
  colnames(c_sum_homo) <- c("num_reads", "het_homo", "era", "homo_genotypes") #this count is now # homo genotypes at X read depth across all contemp individuals
  c_sum_homo <- c_sum_homo[, -2] #remove extra column

#create missing rows for merging
homo_nohet <- setdiff(c_sum_homo$num_reads, c_sum_het$num_reads) #in homo, not het
  het_missing <- data.frame(num_reads = homo_nohet, era = "Contemp", het_genotypes = 0)
het_nohomo <- setdiff(c_sum_het$num_reads, c_sum_homo$num_reads) #in het, not homo
  homo_missing <- data.frame(num_reads = het_nohomo, era = "Contemp", homo_genotypes = 0)

#merge contemporary dataframes
c_sum_het_full <- rbind(c_sum_het, het_missing)
c_sum_homo_full <- rbind(c_sum_homo, homo_missing)

c_genotypes <- merge(c_sum_het, c_sum_homo, by = c("num_reads", "era"))

## All reads together ##
#check for missing rows before merging
alb_nocontemp <- setdiff(a_genotypes$num_reads, c_genotypes$num_reads) #in alb, not contemp
  contemp_missing <- data.frame(num_reads = alb_nocontemp, era = "Contemp", 
                                het_genotypes = 0, homo_genotypes = 0)
contemp_noalb <- setdiff(c_genotypes$num_reads, a_genotypes$num_reads) #in contemp, not alb
  alb_missing <- data.frame(num_reads = contemp_noalb, era = "Hist", 
                           het_genotypes = 0, homo_genotypes = 0)

#merge all dataframes
alb_genotypes_full <- rbind(a_genotypes, alb_missing)
contemp_genotypes_full <- rbind(c_genotypes, contemp_missing)

all_genotypes <- rbind(alb_genotypes_full, contemp_genotypes_full)
  all_genotypes$tot_genotypes <- all_genotypes$het_genotypes + all_genotypes$homo_genotypes #calculate total genotypes at depth (for an era)
  all_genotypes$perc_het <- all_genotypes$het_genotypes/all_genotypes$tot_genotypes #calculate % genotypes that are het at depth (for an era)

################################################################################################
  
######## Visualize results ########

#percent het genotypes by read depth --> at any given depth, what % of calls are het?
all_genotypes %>% 
  ggplot(aes(x=num_reads, y = perc_het, color = era)) + 
  geom_point(size = 4) + 
  labs(y = "% heterozygous genotypes", 
       x = "read depth") + 
  scale_x_continuous(limits = c(1, 100)) + 
  scale_y_continuous(limits = c(0, 0.5)) + 
  scale_color_manual(values = c("#afc8a4", "#1c3b0e")) +
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())

#tot num genotypes by read depth --> at any given depth, what is tot # genotypes?
all_genotypes %>%
  ggplot(aes(x=num_reads,
             y=tot_genotypes,
             color=era, 
             fill=era)) +
  geom_bar(stat = "identity") + 
  labs(y = "total # genotypes", 
       x = "read depth") + 
  scale_x_continuous(limits = c(1, 100)) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e")) + 
  scale_fill_manual(values = c("#afc8a4", "#1c3b0e")) + 
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_text(size = 24, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        strip.text = element_text(size = 24, color = "black")) +
  facet_grid(era ~ .,
             scales = "free")

#histogram of average read depth --> does per locus mean read depth differ between het & homo genotypes?
avg_read_depth_by_locus_era_genotype[avg_read_depth_by_locus_era_genotype$era == "Contemp", ] %>%
  ggplot(aes(x=mean_depth)) +
  geom_histogram(color = "#afc8a4", 
                 fill = "#afc8a4") + 
  labs(y = "count", 
       x = "mean read depth") + 
  scale_x_continuous(limits = c(1, 100)) +
  theme_bw() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 24, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        strip.text = element_text(size = 24, color = "black")) +
  facet_grid(het_homo ~ .,
             scales = "free")

#percent het genotypes by tot reads at locus --> relationship between total # reads & percent heterozygote
full_count_wtot[full_count_wtot$het_homo == "hetero" & full_count_wtot$era == "Contemp", ] %>% 
  ggplot(aes(x=tot_depth, y = perc_hethomo, color = era)) + 
  geom_point(size = 4) + 
  labs(y = "% heterozygous genotypes", 
       x = "total read count") + 
  scale_color_manual(values = c("#afc8a4")) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank())

#read depth ratio across eras --> avg read depth homo / avg read depth het --> do homozygous calls have ~half the depth of het?
avg_read_depth_combined %>% 
  ggplot(aes(x=read_ratio, color = era, fill = era)) + 
  geom_density() + 
  labs(y = "density", 
       x = "read ratio (homozygous/heterozygous)") + 
  scale_x_continuous(limits = c(0, 2)) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e")) + 
  scale_fill_manual(values = c("#afc8a4", "#1c3b0e")) +  
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank()) + 
  facet_grid(era ~ .,
             scales = "free")

#allele balance
avg_AB_by_locus_era_hetero %>% 
  ggplot(aes(x=mean_AB, color = era, fill = era)) + 
  geom_density() + 
  geom_vline(aes(xintercept = 0.375), linetype = "dashed", linewidth = 1.5, color = "#666666") +
  geom_vline(aes(xintercept = 0.625), linetype = "dashed", linewidth = 1.5, color = "#666666") +
  labs(y = "density", 
       x = "mean allele balance") + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("#afc8a4", "#1c3b0e")) + 
  scale_fill_manual(values = c("#afc8a4", "#1c3b0e")) +  
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank()) + 
  facet_grid(era ~ .,
             scales = "free")