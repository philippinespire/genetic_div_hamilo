####################################################### Script for Diploid Filtering ##############################################################

#Modified from Chris's version (https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/indvAlleleBalance.R)
#Creates a "greenlist" of loci that match diploidy assumptions for downstream analyses
#Filtering follow HDPlot guidelines written by McKinney et al. (2017) doi:10.1111/1755-0998.12613
#Filters for missing data by population
#Calculates read depth v. heterozygosity stats to look for issues with reference bias and under-powered genotype calls
#Filters based on reference bias issues

###################################################################################################################################################

######## Set-up ########

#set working directory
getwd()

remove(list = ls())

#load libraries
library(tidyverse)
library(magrittr)
library(janitor)
library(data.table)

#set user defined variables
alleledepthFILE = 'Data/Aen_Ham/Aen.A.rad.RAW-6-6-rescaled.Fltr18.1.HWE.keeploci.renamed.AD.tsv'
genoFILE = 'Data/Aen_Ham/Aen.A.rad.RAW-6-6-rescaled.Fltr18.1.HWE.keeploci.renamed.GT.tsv'
indvPATTERN = "aen_"
modIdPATTERN = "_ae_.*"
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

#### modified for octoploidy (added states 3-8) ####
#genotypes <-
#  read_tsv(genoFILE,
#           col_types=cols(.default = "c")) %>%
#  clean_names() %>%
#  pivot_longer(cols=contains(indvPATTERN),
#               names_to="id") %>%
#  separate(col=value,
#           into=(c("state_1",
#                   "state_2", 
#                   "state_3", 
#                   "state_4",
#                   "state_5", 
#                   "state_6",
#                  "state_7",
#                   "state_8"))) %>%
#  mutate(state_1 = na_if(state_1,
#                         ""),
#         state_2 = na_if(state_2,
#                         ""),
#         state_3 = na_if(state_3,
#                         ""),
#         state_4 = na_if(state_4,
#                         ""),
#         state_5 = na_if(state_5,
#                         ""),
#         state_6 = na_if(state_6,
#                         ""),
#         state_7 = na_if(state_7,
#                         ""),
#         state_8 = na_if(state_8,
#                         ""),
#         genotype = case_when(state_1 != state_2 ~ "hetero", state_1 != state_3 ~ "hetero", state_1 != state_4 ~ "hetero", state_1 != state_5 ~ "hetero", state_1 != state_6 ~ "hetero", state_1 != state_7 ~ "hetero", state_1 != state_8 ~ "hetero",
#                              state_2 != state_3 ~ "hetero", state_2 != state_4 ~ "hetero", state_2 != state_5 ~ "hetero", state_2 != state_6 ~ "hetero", state_2 != state_7 ~ "hetero", state_2 != state_8 ~ "hetero",
#                              state_3 != state_4 ~ "hetero", state_3 != state_5 ~ "hetero", state_3 != state_6 ~ "hetero", state_3 != state_7 ~ "hetero", state_3 != state_8 ~ "hetero",
#                              state_4 != state_5 ~ "hetero", state_4 != state_6 ~ "hetero", state_4 != state_7 ~ "hetero", state_4 != state_8 ~ "hetero",
#                              state_5 != state_6 ~ "hetero", state_5 != state_7 ~ "hetero", state_5 != state_8 ~ "hetero",
#                              state_6 != state_7 ~ "hetero", state_6 != state_8 ~ "hetero",
#                              state_7 != state_8 ~ "hetero",
#                              is.na(state_1) & is.na(state_2) ~ NA_character_,
#                              state_1 == ref ~ "homo_ref",
#                              state_1 == alt ~ "homo_alt",
#                              TRUE ~ "error"))

#######################################################################################################################################################

######### Calculate heterozygosity by individual & locus ########
#should be done for both octoploid and diploid data

#combine data
num_loci <-
  genotypes %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  pull(chrom_pos) %>%
  unique() %>%
  length()

#calculate heterozygosity
genotypes %>%
  group_by(id,
           genotype) %>%
  summarize(num_pos = n()) %>%
  pivot_wider(names_from = "genotype",
              values_from = "num_pos") %>%
  mutate(heterozygosity_obs_ind = hetero/(hetero + homo_ref + homo_alt),
         pop = case_when(heterozygosity_obs_ind >= het_cutoff_pop3 ~ 3,
                         heterozygosity_obs_ind >= het_cutoff_pop1 ~ 1,
                         TRUE ~ 2)) %>%
  select(id,
         heterozygosity_obs_ind,
         pop) %>%
  ggplot(aes(x=id,
             y=heterozygosity_obs_ind,
             color=pop)) +
  geom_point()

all_data <-
  allele_depths %>%
  #left_join(geno_likes) %>%
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

#for diploid, some of the loci were missing heterozygosity_obs_locus (bc missing homo_ref column)
#calculate separately
missing_het_loc <- all_data[is.na(all_data$heterozygosity_obs_locus), ] #pull rows with missing het_locus data

missing_het_loc <- missing_het_loc %>%  #count number of each genotype per locus/era
  group_by(chrom_pos, genotype, era) %>%
  summarize(num_pos = n()) %>%
  pivot_wider(names_from = "genotype", 
              values_from = "num_pos")

colnames(missing_het_loc) <- c("chrom_pos", "era", "hetero", "homo_alt", "homo_ref") #rename columns, one was given an NA --> likely what caused earlier calculation to fail

missing_het_loc <- missing_het_loc %>% #calculate het_locus
  mutate(across(hetero:homo_ref, ~replace_na(.,0)), )

missing_het_loc$heterozygosity_obs_locus <- 
  missing_het_loc$hetero/(missing_het_loc$hetero + 
                            missing_het_loc$homo_alt + 
                            missing_het_loc$homo_ref)

#merge two het_locus dataframes
missing_het_loc$chrom_pos_era <- paste0(missing_het_loc$chrom_pos, "_", missing_het_loc$era) #calculate column to merge by
  all_data$chrom_pos_era <- paste0(all_data$chrom_pos, "_", all_data$era)
all_data$heterozygosity_obs_locus2 <- NA #create new column to fill with newly-calculated het_loc data (just incase something is wrong)

all_data <- transform(all_data,
                      heterozygosity_obs_locus2 = 
                        missing_het_loc[match(chrom_pos_era, 
                                              missing_het_loc$chrom_pos_era), ]$heterozygosity_obs_locus) #basically, fills column in all_data with values from column in missing_het_loc, matching by chrom_pos_era

all_data <- all_data %>% #combine two het_locus columns into one
  mutate(across(heterozygosity_obs_locus, 
                coalesce, 
                heterozygosity_obs_locus2))

all_data <- subset(all_data, select = -c(chrom_pos_era, heterozygosity_obs_locus2)) #remove extra columns

#save as RDS
saveRDS(all_data, 
       file=str_c(outDIR,
                  "/",
                  "all_data_diploid.rds", #change depending on ploidy
                  sep=""))
 
 #######################################################################################################################
 
######## Summarize data by chrom, pos, era & pop ########
#REALLY ONLY NEEDED FOR AEN --> special diploid filtering
 
#read in if running separately
#all_data <- readRDS("Data/Aen_Ham/diploid_filtering/all_data_diploid.rds")

#summarize by chrom, pos, era & pop
hetero_data_chrom_pos_era <-
  all_data %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  filter(genotype == "hetero") %>% 
  group_by(chrom, 
           pos,
           chrom_pos,
           era) %>%
  summarize(total_num_reads = sum(num_reads),
            total_num_altreads = sum(num_reads_alt),
            total_num_refreads = sum(num_reads_ref),
            mean_read_depth = mean(num_reads),
            mode_read_depth = getMode(num_reads),
            median_read_depth = median(num_reads),
            mode_allele_balance = getMode(allele_balance),
            median_allele_balance = median(allele_balance),
            mean_allele_balance = mean(allele_balance),
            mean_heterozygosity_obs_locus = mean(heterozygosity_obs_locus), #does this really do anything? just one value bc by def calculated across locus
            n = n())

hetero_data_chrom_pos_era$tot_allelebalance <- 
  hetero_data_chrom_pos_era$total_num_altreads/hetero_data_chrom_pos_era$total_num_reads

#write out
write.csv(hetero_data_chrom_pos_era, file = "Data/Aen_Ham/meanAB_data_era_pos.csv")

######## Visualize data ########

#histogram of median allele balance by era
hetero_data_chrom_pos_era %>% 
  ggplot(aes(x=median_allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8,
                            1/6,
                            2/8,
                            2/6,
                            3/8,
                            4/8,
                            5/8,
                            4/6,
                            6/8,
                            5/6,
                            7/8),
             color="grey",
             linetype="dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Medians by Position") +
  facet_grid(era ~ .,
             scales = "free")

ggsave(paste(outDIR, 
             "/",
             'HIST-AB-medPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

#histogram of mode allele balance by era
hetero_data_chrom_pos_era %>%
  ggplot(aes(x=mode_allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8,
                            1/6,
                            2/8,
                            2/6,
                            3/8,
                            4/8,
                            5/8,
                            4/6,
                            6/8,
                            5/6,
                            7/8),
             color="grey",
             linetype="dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Modes by Position") +
  facet_grid(era ~ .,
             scales = "free")

ggsave(paste(outDIR, 
             "/",
             'HIST-AB-modePOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

#box plot of allele balance x individual
all_data %>%
  ggplot(aes(x=id,
             y=allele_balance,
             fill=era)) +
  geom_boxplot() +
  geom_hline(yintercept = c(1/8,
                            1/6,
                            1/5,
                            1/4),
             color="black",
             linetype="solid") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  labs(title = "Boxplots of Allele Balance",
       subtitle = "Individuals x Position, Most Granular, HLines at 1/8, 1/6, 1/5, 1/4",
       x = "Indiviudal ID") 

ggsave(paste(outDIR, 
             "/",
             outFilePREFIX,
             'BOXPL-AB-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

#box plot of read depth x individual
all_data %>%
  ggplot(aes(x=id,
             y=log10(num_reads),
             fill=era)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  labs(title = "Boxplots of Read Depth",
       subtitle = "Individuals x Position, Most Granular",
       x = "Indiviudal ID") 

ggsave(paste(outDIR, 
             "/",
             'BOXP-DP-INDxPOS.png', 
             sep = ""), 
       height = 9, 
       width = 6.5)

#bar plot of genotypes x individual
all_data %>%
  ggplot(aes(x=id,
             fill = genotype)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) 


ggsave(paste(outDIR, 
             "/",
             'BARPL-GT-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

##############################################################################################################

######### Filtering to green list of loci ########
#loci that meet diploid expectations
#again, REALLY ONLY NEED FOR AEN

#read in data if running separately
hetero_data_chrom_pos_era <- read.csv(file = "Data/Aen_Ham/diploid_filtering/meanAB_data_era_pos.csv", header = TRUE)
all_data <- readRDS("Data/Aen_Ham/diploid_filtering/all_data_diploid.rds")
all_data_octoploid <- readRDS("Data/Aen_Ham/diploid_filtering/all_data_octoploid.rds")

##### Calculate D #### 
# D = "read-ratio deviation" (statistical test for deviation from expected read ratio)
#equivalent to z-score
#Based on McKinney et al. 2017 doi:10.1111/1755-0998.12613

#add SD column
hetero_data_chrom_pos_era$read_SD <- sqrt(0.5*(1-0.5)*hetero_data_chrom_pos_era$total_num_reads)

#add z-score (equivalent to D from McKinney et al. 2017)
#describes deviation between observed and expected allelic-specific counts from a binomial distribution (w/ p = q = 0.5)
#higher z-scores = greater deviation (essentially, less likely to get these read counts if truly heterozygote)
hetero_data_chrom_pos_era$z_score <- 
  ((hetero_data_chrom_pos_era$total_num_reads/2) - 
  hetero_data_chrom_pos_era$total_num_altreads)/hetero_data_chrom_pos_era$read_SD

#look at z-score distribution
#this distribution will determine what filtering thresholds should be set at
#want z-score & num het cut-offs that capture the inner portion of the scatterplot and exclude the outer ring
hetero_data_chrom_pos_era %>%
  ggplot(aes(x=z_score,
             y = mean_heterozygosity_obs_locus,
             color = era)) +
  geom_point() +
  labs(y = "percent heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

ggsave(paste(outDIR,
             "/", 
             'zscore_dist.png', 
             sep = ""), 
       height = 8, 
       width = 9)

######## Filter to list that pass thresholds in both eras (het in both) ########
#filter in eras separately, so create era-specific dataframes
hetero_data_chrom_pos_contemp <- subset(hetero_data_chrom_pos_era, 
                                        hetero_data_chrom_pos_era$era == "cbat")

hetero_data_chrom_pos_alb <- subset(hetero_data_chrom_pos_era, 
                                    hetero_data_chrom_pos_era$era == "a")

#filter to loci w/het <=0.6 (based on McKinney et al. 2017)
hetero_data_chrom_pos_contemp_het0.6 <- 
  subset(hetero_data_chrom_pos_contemp, 
         hetero_data_chrom_pos_contemp$mean_heterozygosity_obs_locus <= 0.60000000)

hetero_data_chrom_pos_alb_het0.6 <- 
  subset(hetero_data_chrom_pos_alb, 
         hetero_data_chrom_pos_alb$mean_heterozygosity_obs_locus <= 0.60000000)

#filter by z_score (btwn -2.5 & 2.5 based on zscore plot and low depth of coverage in Alb)
hetero_data_chrom_pos_contemp_het0.6_zscore2.5 <- 
  subset(hetero_data_chrom_pos_contemp_het0.6, 
         hetero_data_chrom_pos_contemp_het0.6$z_score <= 2.5 & 
           hetero_data_chrom_pos_contemp_het0.6$z_score >= -2.5)

hetero_data_chrom_pos_alb_het0.6_zscore2.5 <- 
  subset(hetero_data_chrom_pos_alb_het0.6, 
         hetero_data_chrom_pos_alb_het0.6$z_score <= 2.5 & 
           hetero_data_chrom_pos_alb_het0.6$z_score >= -2.5)

#check filtering worked with scatterplot
#AB should center around 0.5, should be no % het greater than 0.6
hetero_data_chrom_pos_contemp_het0.6_zscore2.5 %>% 
  ggplot(aes(x=mean_allele_balance,
             y = mean_heterozygosity_obs_locus)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "percent heterozygotes")

hetero_data_chrom_pos_alb_het0.6_zscore2.5 %>% 
  ggplot(aes(x=mean_allele_balance,
             y = mean_heterozygosity_obs_locus)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "percent heterozygotes")

#list of loci that pass filters in both eras
greenlist_hetboth <- as.data.frame(intersect(hetero_data_chrom_pos_contemp_het0.6_zscore2.5$chrom_pos, 
                                             hetero_data_chrom_pos_alb_het0.6_zscore2.5$chrom_pos))
  colnames(greenlist_hetboth) <- "chrom_pos"

######## Filter to list that pass thresholds in one era but homozygous in other era ########
#loci that pass contemp filters but are NOT heterozygous in albatross
greenlist_hetcontemp_notalb <- as.data.frame(setdiff(hetero_data_chrom_pos_contemp_het0.6_zscore2.5$chrom_pos, 
                                                     hetero_data_chrom_pos_alb$chrom_pos))
  colnames(greenlist_hetcontemp_notalb) <- "chrom_pos"

#loci that pass albatross filters but are NOT heterozygous in contemp
greenlist_hetalb_notcontemp <- as.data.frame(setdiff(hetero_data_chrom_pos_alb_het0.6_zscore2.5$chrom_pos, 
                                                     hetero_data_chrom_pos_contemp$chrom_pos))
  colnames(greenlist_hetalb_notcontemp) <- "chrom_pos"

######## Filter to list that are homozygous only ########
#e.g. homozygous ref in one era and alt in the other
#list of all loci
all_loci_list <- as.data.frame(sort(unique(as.character(all_data$chrom_pos))))
  colnames(all_loci_list) <- "chrom_pos"

#list of only loci that have at least one heterozygote
het_loci_list <- as.data.frame(sort(unique(as.character(hetero_data_chrom_pos_era$chrom_pos))))
  colnames(het_loci_list) <- "chrom_pos"

#loci that don't have any het at all (either era)
greenlist_homo_only <- as.data.frame(setdiff(all_loci_list$chrom_pos, het_loci_list$chrom_pos))
  colnames(greenlist_homo_only) <- "chrom_pos"

#### Create list of homozygote loci where genotype calls differ based on ploidy level ####
#e.g. loci where all indv are homozygote at diploid level, but when called with octoploidy at least one indv switches to het
#want to remove these loci from final "greenlist"

#subset octoploid loci list to those that appear in diploid dataset
all_data_octoploid_subset <- subset(all_data_octoploid, chrom_pos %in% all_data$chrom_pos)
  
#create column to merge by (chrom-pos-individual ID)
all_data_octoploid_subset$chrom_pos_id <- paste(all_data_octoploid_subset$chrom_pos, 
                                                "_", all_data_octoploid_subset$id)
all_data$chrom_pos_id <- paste(all_data$chrom_pos, "_", all_data$id)

#merge dataframes
all_data_merged <- merge(all_data_octoploid_subset, all_data, by = "chrom_pos_id")
all_data_merged_nona <- all_data_merged[(!is.na(all_data_merged$genotype.y)), ] #bc some genotype x (octoploid) called when diploid wasn't --> often bc changed to missing during filtering bc read depth too low (<10X)
  
#create logic vector indicating whether or not all genotypes match for given locus
match_factor<-c(1:5070039) #length of all_data_merged_nona --> create vector to populate
  
for (i in 1:1815867) {
  if (all_data_merged_nona$genotype.x[i] == all_data_merged_nona$genotype.y[i]) {
    match_factor[i] <- "TRUE"
    }
  }
  
#add matching factor to merged dataframe
matching <- c(match_factor) #bc sometimes original vector doesn't add properly
  all_data_merged_nona$matching <- matching
  
#identify rows where genotypes don't match up (match_factor != TRUE)
genotype_mismatches <- subset(all_data_merged_nona, matching != TRUE)
  
#identify chrom & pos of loci where genotypes don't match up
genotype_mismatches_chrompos <- unique(sort((genotype_mismatches[, "chrom_pos.x"])))
  genotype_mismatches_chrompos_df <- as.data.frame(genotype_mismatches_chrompos)
    colnames(genotype_mismatches_chrompos_df) <- "chrom_bp"

##### Check mismatch genotype list against greenlist of homozygous loci from diploid data ####
#want to make sure not changing calls to hetero if called assuming octoploidy --> any that change calls at even one indv are thrown out

#filter original greenlist of homozygous sites to ones with no genotype mismatches
greenlist_homo_only_match <- as.data.frame(setdiff(greenlist_homo_only$chrom_pos, 
                                                          genotype_mismatches_chrompos_df$chrom_bp))
  colnames(greenlist_homo_only_match) <- "chrom_pos"

######## Merge all greenlists into one ########
#merge greenlists
greenlist_full <- as.data.frame(rbind(greenlist_hetboth, greenlist_hetcontemp_notalb, 
                                      greenlist_hetalb_notcontemp, greenlist_homo_only_match))
greenlist_full_sorted <- as.data.frame(sort(unique(as.character(greenlist_full$chrom_pos)))) #sort by chrom_bp
  colnames(greenlist_full_sorted) <- "chrom_pos"

#split chrom & pos into different columns
greenlist_full_split <- greenlist_full_sorted %>% 
  separate(chrom_pos, into = c("dDocent", "Contig", "contig_num", "pos"),
           sep = "_", remove = FALSE)

#create contig column
greenlist_full_split$contig <- paste(greenlist_full_split$dDocent, greenlist_full_split$Contig, 
                                     greenlist_full_split$contig_num, sep = "_")

#keep only contig & pos columns
greenlist_toprint <- greenlist_full_split[, c("contig", "pos")]

#check # of contigs kept
greenlist_full_contigs <- sort(unique(as.character(greenlist_toprint$contig)))
  
#write out list without header
write.table(greenlist_toprint, file = "Data/Aen_Ham/diploid_filtering/greenlist_loci_full_HD_2.5.txt", 
            col.names = FALSE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

######## Check distribution of diploid loci across contigs ########
#do these "diploid" loci cluster on the same contigs?
#want only loci on contigs with mostly diploid-passing loci (<80%) --> these most likely to be truly diploid

#read in if running separately
hetero_data_chrom_pos_era <- read.csv("Data/Aen_Ham/diploid_filtering/meanAB_data_era_pos.csv", header = TRUE)

pass_diploid <- read.table("Data/Aen_Ham/diploid_filtering/greenlist_loci_full_HD_2.5.txt", header = FALSE) #9251
#pass_diploid <- greenlist_toprint
  colnames(pass_diploid) <- c("chrom", "pos")
  pass_diploid$chrom_pos <- paste0(pass_diploid$chrom, "_", pass_diploid$pos) #add unique chrom_pos column
    pass_diploid <- pass_diploid[order(pass_diploid$chrom_pos), ]
  pass_diploid$status <- "pass" #add status (these loci passed HD diploid filters)

#subset het dataset to just loci on contigs with at least one passing locus
chrom_pass <- hetero_data_chrom_pos_era[hetero_data_chrom_pos_era$chrom %in% 
                                          pass_diploid$chrom, ]

#get just unique chrom_pos (not duplicates)
chrom_pass_unique <- as.data.frame(unique(chrom_pass$chrom_pos))
  colnames(chrom_pass_unique) <- c("chrom_pos")
  chrom_pass_unique <- separate(chrom_pass_unique, chrom_pos, 
                                c("CHROM", "contig", "num", "pos"), 
                                sep = "_", remove = FALSE)
  chrom_pass_unique$chrom <- paste0(chrom_pass_unique$CHROM, "_", 
                                    chrom_pass_unique$contig, "_", chrom_pass_unique$num)
  chrom_pass_unique <- subset(chrom_pass_unique, select = c("chrom_pos", "chrom", "pos"))

#### get list of loci on "diploid" chrom that failed HD filters ####
fail_diploid <- as.data.frame(setdiff(chrom_pass_unique$chrom_pos, 
                                      pass_diploid$chrom_pos)) #89645
  colnames(fail_diploid) <- c("chrom_pos")
  fail_diploid <- separate(fail_diploid, chrom_pos, 
                           c("CHROM", "contig", "num", "pos"), 
                           sep = "_", remove = FALSE)
  fail_diploid$chrom <- paste0(fail_diploid$CHROM, "_", 
                               fail_diploid$contig, "_", fail_diploid$num)
  fail_diploid <- subset(fail_diploid, select = c("chrom_pos", "chrom", "pos"))
  fail_diploid$status <- "fail" #add status (these loci failed HD diploid filters)

#### build full dataset of loci on contigs with at least one passing locus ####
full_loci_list <- as.data.table(rbind(pass_diploid, fail_diploid))
  full_loci_list <- full_loci_list[order(full_loci_list$chrom),] #149129

##### count status by contig ####
status_by_chrom <- full_loci_list[, .N, by = .(chrom, status)]

#converting to wide to get % loci that pass/contig
status_by_chrom_wide <- as.data.table(pivot_wider(data = status_by_chrom, 
                                                  names_from = status,
                                                  values_from = N))
status_by_chrom_wide[is.na(status_by_chrom_wide)] <- 0 #change NAs to 0

#calculate % loci that pass HD filters/contig
status_by_chrom_wide$perc_pass <- status_by_chrom_wide$pass/(status_by_chrom_wide$pass + 
                                                               status_by_chrom_wide$fail)

#### bin contigs by # total loci ####
#want to see, if contig is 100% passing, how many loci are actually on it (1 or 20?)
#calculate total num loci found on each contig
status_by_chrom_wide$numloci <- status_by_chrom_wide$pass + status_by_chrom_wide$fail

#create column to populate by # of loci bins (e.g. 1 locus, 2-5, 6-10, etc.)
status_by_chrom_wide$numloci_range <- "first" #create empty column to populate 
status_by_chrom_wide$numloci_range[status_by_chrom_wide$numloci == 1] <- "1"
status_by_chrom_wide$numloci_range[status_by_chrom_wide$numloci  >=2 & 
                                     status_by_chrom_wide$numloci <=5] <- "2-5"
status_by_chrom_wide$numloci_range[status_by_chrom_wide$numloci >= 6 & 
                                     status_by_chrom_wide$numloci <= 10] <- "6-10"
status_by_chrom_wide$numloci_range[status_by_chrom_wide$numloci >= 11 & 
                                     status_by_chrom_wide$numloci <= 20] <- "11-20"
status_by_chrom_wide$numloci_range[status_by_chrom_wide$numloci >20] <- ">20"

#histogram plot by % pass
perc_pass_plot <- ggplot(data = status_by_chrom_wide, 
                         aes(x = perc_pass, fill = numloci_range)) + 
  geom_histogram(binwidth = 0.1) + 
  theme_minimal()
perc_pass_plot #somewhat bimodal, makes sense

#get list of passing contigs (those with at least 80% of SNPs that pass original HD filters)
status_by_chrom_wide_0.8 <- subset(status_by_chrom_wide, 
                                   status_by_chrom_wide$perc_pass > 0.79000000000) #9547 contigs
  sum(status_by_chrom_wide_0.8$numloci) #23876 loci
  
#print out list of diploid contigs
#printing out as bed file, will subset VCF to this list of contigs
greenlist_contigs <- status_by_chrom_wide_0.8[, "chrom"]
  greenlist_contigs$chromStart <- 1 #adding pseudo-chromStart column
  greenlist_contigs$chromEnd <- 1000000 #adding pseudo-chromEnd column (needs to be large enough to include all loci)
  
#format for printing
greenlist_contigs <- format(greenlist_contigs, scientific = FALSE)

#write out list with header
write.table(greenlist_contigs, file = "Data/Aen_Ham/diploid_filtering/diploid_contigs.bed", 
            col.names = TRUE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

######## Visualize results ########
#subset original het_by_era to only chrom w/ >80% passing loci
chrom_pass_0.8 <- chrom_pass[chrom_pass$chrom %in% 
                               greenlist_contigs$chrom, ] #smaller --> prob bc alot of these loci that pass filters are homozygous across eras

#add diploid filter status
chrom_pass_0.8 <- merge(chrom_pass_0.8, 
                        full_loci_list[, c('chrom_pos', 'status')], all.x = TRUE)

#add SD column
chrom_pass_0.8$read_SD <- sqrt(0.5*(1-0.5)*chrom_pass_0.8$total_num_reads)

#add z-score (equivalent to D from McKinney et al. 2017)
chrom_pass_0.8$z_score <- 
  ((chrom_pass_0.8$total_num_reads/2) - 
     chrom_pass_0.8$total_num_altreads)/chrom_pass_0.8$read_SD

#HD plot
chrom_pass_0.8 %>%
  ggplot(aes(x=z_score,
             y = mean_heterozygosity_obs_locus,
             color = status)) +
  geom_point() +
  scale_x_continuous(limits = c(-25, 25)) +
  labs(y = "percent heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

#read depth v heterozygosity
hetero_data_chrom_pos_era %>% 
  ggplot(aes(x=mean_read_depth, y = mean_heterozygosity_obs_locus)) + 
  geom_point(size = 1.5) + 
  labs(y = "heterozygosity per locus", 
       x = "mean read depth per locus (heterozygote only)") + 
  scale_x_continuous(limits = c(0, 100)) + 
  theme_bw() +
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        strip.text = element_text(size = 28, color = "black")) +
  facet_grid(era ~ ., scales="free_y")

ggsave(paste("Plots",
             "/", 
             'Aen_depth_het.png', 
             sep = ""), 
       height = 8, 
       width = 9)

####################################################################################################################

######## Filtering for % missing within populations ########
#REALLY ONLY FOR AEN (or other species that have one population/species with very low individuals)
#When final Filter17 overfilters data
#Only focusing on A populations, since B populations have 1-2 individuals (and thus can have inflated missing data - also, these individuals will not be included in any of the downstream analyses)

#### read in data ####
AHam_miss <- as.data.frame(read.table("Data/Aen_Ham/diploid_filtering/Aen.A.rad.RAW-6-6-rescaled.Fltr17.2.Aen-AHam.lmiss", 
                                      header = TRUE)) #12307 loci (polymorphic), 577865 loci (monomorphic)
  AHam_miss$chrom_pos <- paste0(AHam_miss$CHR, "_", 
                                AHam_miss$POS)
CBat_miss <- as.data.frame(read.table("Data/Aen_Ham/diploid_filtering/Aen.A.rad.RAW-6-6-rescaled.Fltr17.2.Aen-Cbat.lmiss", 
                                      header = TRUE)) #12307 loci (polymorphic), 577865 loci (monomorphic)
  CBat_miss$chrom_pos <- paste0(CBat_miss$CHR, "_", 
                                CBat_miss$POS)

#### filter to loci missing in <50% of indv ####
  
AHam_miss_keep <- subset(AHam_miss, AHam_miss$F_MISS <= 0.5) #6620 loci (polymorphic), 239990 loci (monomorphic)
CBat_miss_keep <- subset(CBat_miss, CBat_miss$F_MISS <= 0.5) #12307 loci (polymorphic), 577865 loci (monomorphic), none were removed either time

#since no loci are missing from >50% of CBat pop, can just use AHam list as those to keep

#### write out loci to keep ####

keep_loci <- AHam_miss_keep[, c("CHR", "POS")]
  keep_loci_ordered <- keep_loci[order(keep_loci[,"POS"]), ]

write.table(keep_loci, file = "Data/Aen_Ham/diploid_filtering/loci_nomiss.txt", #change as needed
            col.names = FALSE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

#####################################################################################################################

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
  scale_fill_manual(values = c("#16537e", "#8aa9be")) + 
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
#2. look for issues with reference/mapping bias (do homozygous calls have lower read depth than het -- does this vary by era) ESP FOR AEN

#read in if running separately
all_data <- readRDS("Data/Ela_Ham/all_data_diploid.rds")
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

#### Visualize results ####
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
  scale_color_manual(values = c("#8aa9be", "#16537e")) + 
  scale_fill_manual(values = c("#8aa9be", "#16537e")) + 
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
  geom_histogram(color = "#8aa9be", 
                 fill = "#8aa9be") + 
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
  scale_color_manual(values = c("#8aa9be")) + 
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
  scale_color_manual(values = c("#8aa9be", "#16537e")) + 
  scale_fill_manual(values = c("#8aa9be", "#16537e")) +  
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
  scale_color_manual(values = c("#8aa9be", "#16537e")) + 
  scale_fill_manual(values = c("#8aa9be", "#16537e")) +  
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank()) + 
  facet_grid(era ~ .,
             scales = "free")

#######################################################################################################################################################

######## Check for Mapping/Reference Bias ########

#### Filtering for AB ####
#binomial test close to 0.5
#also keeping ones that have AB between 0.4 & 0.6 (bc binomial is too stringent when read depth is high)

## Sum number of reads by locus, era, genotype ##
#sum reads by locus, era
tot_num_reads_by_locus_era_genotype <-
  all_data %>%
  group_by(chrom_pos,
           genotype,
           era) %>%
  summarize(tot_reads <- sum(num_reads, na.rm = TRUE))

tot_num_reads_by_locus_era_genotype$era[tot_num_reads_by_locus_era_genotype$era == "a"] <- "Hist"
tot_num_reads_by_locus_era_genotype$era[tot_num_reads_by_locus_era_genotype$era == "cbat"] <- "Contemp"

tot_num_reads_by_locus_era_genotype_hetero <- subset(tot_num_reads_by_locus_era_genotype, tot_num_reads_by_locus_era_genotype$genotype == "hetero")
  tot_num_reads_by_locus_era_genotype_hetero <- tot_num_reads_by_locus_era_genotype_hetero[, -2]

colnames(tot_num_reads_by_locus_era_genotype_hetero) <- c("chrom_pos", "era", "tot_depth")

## Sum number of reads by locus, era, genotype ##
#sum reads by locus, era
tot_ref_reads_by_locus_era_genotype <-
  all_data %>%
  group_by(chrom_pos,
           genotype,
           era) %>%
  summarize(tot_ref_reads <- sum(num_reads_ref, na.rm = TRUE))

tot_ref_reads_by_locus_era_genotype$era[tot_ref_reads_by_locus_era_genotype$era == "a"] <- "Hist"
tot_ref_reads_by_locus_era_genotype$era[tot_ref_reads_by_locus_era_genotype$era == "cbat"] <- "Contemp"

tot_ref_reads_by_locus_era_genotype_hetero <- subset(tot_ref_reads_by_locus_era_genotype, tot_ref_reads_by_locus_era_genotype$genotype == "hetero")
  tot_ref_reads_by_locus_era_genotype_hetero <- tot_ref_reads_by_locus_era_genotype_hetero[, -2]

colnames(tot_ref_reads_by_locus_era_genotype_hetero) <- c("chrom_pos", "era", "tot_ref_depth")

#merge dataframes
reads_hetero <- merge(tot_ref_reads_by_locus_era_genotype_hetero, 
                      tot_num_reads_by_locus_era_genotype_hetero, by = c("chrom_pos", "era"))

## Run binomial test ##
reads_hetero$pval <- dbinom(reads_hetero$tot_ref_depth, reads_hetero$tot_depth, 0.5)
reads_hetero$AB <- reads_hetero$tot_ref_depth/reads_hetero$tot_depth

reads_hetero_hist <- subset(reads_hetero, reads_hetero$era == "Hist")

## Get passing list of loci ##
#only care about historical -- as this is where see reference bias
#essentially loci that either pass binomial test AND/OR have AB between 0.4 & 0.6
pass_AB_loci <- subset(reads_hetero_hist, reads_hetero_hist$pval >= 0.05) #1544 loci
tentpass_AB_loci <- subset(reads_hetero_hist, reads_hetero_hist$pval < 0.05 &
                          reads_hetero_hist$AB >= 0.4 & 
                          reads_hetero_hist$AB <= 0.6)

#merge list
green_AB_loci <- rbind(pass_AB_loci, tentpass_AB_loci)

#loci that failed
fail_AB_loci <- as.data.frame(setdiff(reads_hetero_hist$chrom_pos, green_AB_loci$chrom_pos))
  colnames(fail_AB_loci) <- c("chrom_pos")

#then need to add in homozygote only loci and ones that only show up in contemp
#any loci that show up in original dataset that don't show up in reads_hetero_hist
AB_loci_to_add <- as.data.frame(setdiff(tot_num_reads_by_locus_era_genotype$chrom_pos, reads_hetero_hist$chrom_pos))
  colnames(AB_loci_to_add) <- "chrom_pos"
  
combined_AB_loci_list <- c(green_AB_loci$chrom_pos, AB_loci_to_add$chrom_pos)

#### Filtering for Read Depth Ratio ####
#binomial test close to 1 (avg read depth in het = avg read depth in homo)
#also just hard filter to those above 0.8 RDR

## Calculate average depth per locus by locus, genotype, era ##
#calculate average read depth by locus, era, genotype
avg_read_depth_by_locus_era_genotype <-
  all_data %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(mean_depth <- mean(num_reads, na.rm = TRUE))

avg_read_depth_by_locus_era_genotype <- avg_read_depth_by_locus_era_genotype[!is.na(avg_read_depth_by_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads

avg_read_depth_by_locus_era_genotype$era[avg_read_depth_by_locus_era_genotype$era == "a"] <- "Hist"
avg_read_depth_by_locus_era_genotype$era[avg_read_depth_by_locus_era_genotype$era == "cbat"] <- "Contemp"

colnames(avg_read_depth_by_locus_era_genotype) <- c("chrom_pos", "het_homo", "era", "mean_depth")

#pull out by genotype
avg_read_depth_by_locus_era_genotype_het <- subset(avg_read_depth_by_locus_era_genotype, 
                                                   avg_read_depth_by_locus_era_genotype$het_homo == "hetero")
avg_read_depth_by_locus_era_genotype_homo <- subset(avg_read_depth_by_locus_era_genotype, 
                                                    avg_read_depth_by_locus_era_genotype$het_homo == "homo") 

#merge dataframes
avg_read_depth_combined <- merge(avg_read_depth_by_locus_era_genotype_het, avg_read_depth_by_locus_era_genotype_homo, 
                                 by = c("chrom_pos", "era"))
  colnames(avg_read_depth_combined) <- c("chrom_pos", "era", "het", "mean_depth_het", "homo", "mean_depth_homo")
  avg_read_depth_combined <- avg_read_depth_combined[, -3] #removing het & homo columns (not needed)
  avg_read_depth_combined <- avg_read_depth_combined[, -4]

#calculate read ratio
avg_read_depth_combined$read_ratio <- avg_read_depth_combined$mean_depth_homo/avg_read_depth_combined$mean_depth_het

## Run binomial test ##
#round avg depths off to make integer for binomial test
avg_read_depth_combined$mean_depth_het <- round(avg_read_depth_combined$mean_depth_het, 0)
avg_read_depth_combined$mean_depth_homo <- round(avg_read_depth_combined$mean_depth_homo, 0)

#subset to only ones with read ratio less than 0.8
#not worried about if higher (eg, depth at homo is ~ equal to or greater than het)
avg_read_depth_combined_small <- subset(avg_read_depth_combined, read_ratio < 0.8)

#binomial test
avg_read_depth_combined_small$pval <- dbinom(avg_read_depth_combined_small$mean_depth_homo, 
                                       avg_read_depth_combined_small$mean_depth_het, 1.0)

avg_read_depth_combined_small_hist <- subset(avg_read_depth_combined_small,
                                        avg_read_depth_combined_small$era == "Hist")

## Get failing list of loci ##
#only care about historical -- as this is where see reference bias
#essentially loci that either pass binomial test AND/OR have RDR >= 0.8
fail_RDR_loci <- subset(avg_read_depth_combined_small_hist, 
                        avg_read_depth_combined_small_hist$pval < 0.05) #1173 loci (only 1 passed)

#then need to add in homozygote only loci, ones that only show up in contemp, and ones that passed filter (high RDR or passed test)
#any loci that show up in original dataset that don't show up in fail list are fine
RDR_loci_pass <- as.data.frame(setdiff(avg_read_depth_by_locus_era_genotype$chrom_pos, 
                                       fail_RDR_loci$chrom_pos))
  colnames(RDR_loci_pass) <- "chrom_pos"

#### Merge two lists together ####
all_pass_loci <- subset(RDR_loci_pass, !(chrom_pos %in% fail_AB_loci$chrom_pos)) #7111 loci

## Format to print out and filter ##
all_pass_loci_sep <- 
  all_pass_loci %>% 
  separate("chrom_pos", into = c("dDocent", "Contig", "Contig_Number", "pos"), "_") %>%
  unite("contig", 1:2)

#order numerically
all_pass_loci_sep$Contig_Number <- as.numeric(all_pass_loci_sep$Contig_Number) #make numeric for ordering
all_pass_loci_sep_ordered <- all_pass_loci_sep[order(all_pass_loci_sep[, "Contig_Number"]), ]

all_pass_loci_sep_ordered <- unite(all_pass_loci_sep_ordered, "chrom", 1:2)

write.table(all_pass_loci_sep_ordered, 
            file = "Data/Aen_Ham/diploid_filtering/AB_RDR_loci.txt", #change as needed
            col.names = FALSE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

#### Re-plot read depth ratio and mean AB stats ####
#to see how well filtering "fixed" the problem

#subset all_loci df to just passing loci
all_data_pass <- subset(all_data, chrom_pos %in% all_pass_loci$chrom_pos) #931541

## Calculate average AB per locus by locus, genotype, era ##
#calculate average AB by locus, era, genotype
avg_AB_by_locus_era_genotype <-
  all_data_pass %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(mean_AB <- mean(allele_balance, na.rm = TRUE))

avg_AB_by_locus_era_genotype <- avg_AB_by_locus_era_genotype[!is.na(avg_AB_by_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads
avg_AB_by_locus_era_genotype <- subset(avg_AB_by_locus_era_genotype, avg_AB_by_locus_era_genotype$het_homo == "hetero")

avg_AB_by_locus_era_hetero <- avg_AB_by_locus_era_genotype[, -2]  

avg_AB_by_locus_era_hetero$era[avg_AB_by_locus_era_hetero$era == "a"] <- "Hist"
avg_AB_by_locus_era_hetero$era[avg_AB_by_locus_era_hetero$era == "cbat"] <- "Contemp"

colnames(avg_AB_by_locus_era_hetero) <- c("chrom_pos", "era", "mean_AB")

## Calculate average depth per locus by locus, genotype, era ##
#calculate average read depth by locus, era, genotype
avg_read_depth_by_locus_era_genotype <-
  all_data_pass %>%
  mutate(het_homo = ifelse(genotype == "hetero", "hetero", "homo")) %>% #creates "homo" column for homo_alt & homo_ref to group/sum by
  group_by(chrom_pos,
           het_homo,
           era) %>%
  summarize(mean_depth <- mean(num_reads, na.rm = TRUE))

avg_read_depth_by_locus_era_genotype <- avg_read_depth_by_locus_era_genotype[!is.na(avg_read_depth_by_locus_era_genotype$het_homo), ] #just removing NAs -- most of these don't have reads

avg_read_depth_by_locus_era_genotype$era[avg_read_depth_by_locus_era_genotype$era == "a"] <- "Hist"
avg_read_depth_by_locus_era_genotype$era[avg_read_depth_by_locus_era_genotype$era == "cbat"] <- "Contemp"

colnames(avg_read_depth_by_locus_era_genotype) <- c("chrom_pos", "het_homo", "era", "mean_depth")

#pull out by genotype
avg_read_depth_by_locus_era_genotype_het <- subset(avg_read_depth_by_locus_era_genotype, 
                                                   avg_read_depth_by_locus_era_genotype$het_homo == "hetero")
avg_read_depth_by_locus_era_genotype_homo <- subset(avg_read_depth_by_locus_era_genotype, 
                                                    avg_read_depth_by_locus_era_genotype$het_homo == "homo") 

#merge dataframes
avg_read_depth_combined <- merge(avg_read_depth_by_locus_era_genotype_het, avg_read_depth_by_locus_era_genotype_homo, 
                                 by = c("chrom_pos", "era"))
  colnames(avg_read_depth_combined) <- c("chrom_pos", "era", "het", "mean_depth_het", "homo", "mean_depth_homo")
  avg_read_depth_combined <- avg_read_depth_combined[, -3] #removing het & homo columns (not needed)
  avg_read_depth_combined <- avg_read_depth_combined[, -4]

#calculate read ratio
avg_read_depth_combined$read_ratio <- avg_read_depth_combined$mean_depth_homo/avg_read_depth_combined$mean_depth_het

## Plots ##
#read depth ratio across eras --> avg read depth homo / avg read depth het --> do homozygous calls have ~half the depth of het?
avg_read_depth_combined %>% 
  ggplot(aes(x=read_ratio, color = era, fill = era)) + 
  geom_density() + 
  labs(y = "density", 
       x = "read ratio (homozygous/heterozygous)") + 
  scale_x_continuous(limits = c(0, 2)) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_fill_manual(values = c("#e3ccb4", "#a25505")) +  
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
  labs(y = "density", 
       x = "mean allele balance") + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("#e3ccb4", "#a25505")) + 
  scale_fill_manual(values = c("#e3ccb4", "#a25505")) +  
  theme_bw() + 
  theme(axis.title = element_text(size = 28, color = "black"), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 28, color = "black"), 
        legend.title = element_blank()) + 
  facet_grid(era ~ .,
             scales = "free")