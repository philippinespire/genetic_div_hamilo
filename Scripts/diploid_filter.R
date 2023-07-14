#Modified from Chris's version (https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/indvAlleleBalance.R)
#Creates a "greenlist" of loci that match diploidy assumptions for downstream analyses
#Filtering follow HDPlot guidelines written by McKinney et al. (2017) doi:10.1111/1755-0998.12613

######## Set-up ########

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

remove(list = ls())

#load libraries
library(tidyverse)
library(magrittr)
library(janitor)
library(data.table)

#set user defined variables
alleledepthFILE = 'diploid_filtering/Aen.AB.rad.RAW-6-6.Fltr15.9.recode.AD.tsv'
genoFILE = 'diploid_filtering/Aen.AB.rad.RAW-6-6.Fltr15.9.recode.GT.tsv'
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
#                   "state_7",
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

######################################################################################################################

######### Calculate heterozygosity by individual ########
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

#heterozygosity by individual
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
  # left_join(geno_likes) %>%
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
                  outFilePREFIX,
                  "diploid_filtering/alldata_diploid.rds", #change depending on ploidy
                  sep=""))
 
 #######################################################################################################################
 
 ######## Summarize data by chrom, pos, era & pop ########
 
#read in if running separately
#all_data <- readRDS("diploid_filtering/alldata_diploid.rds")

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
write.csv(hetero_data_chrom_pos_era, file = "diploid_filtering/meanAB_data_era_pos.csv")

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
             outFilePREFIX,
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
             outFilePREFIX,
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
             outFilePREFIX,
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
             outFilePREFIX,
             'BARPL-GT-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

##############################################################################################################

######### Filtering to green list of loci ########
#loci that meet diploid expectations

#read in data if running separately
#hetero_data_chrom_pos_era <- read.csv(file = "diploid_filtering/meanAB_data_era_pos.csv", header = TRUE)
#all_data <- readRDS("diploid_filtering/alldata_diploid.rds")
#all_data_octoploid <- readRDS("diploid_filtering/alldata_octoploid.rds")

##### Calculate D #### 
# D = "read-ratio deviation" (statistical test for deviation from expected read ratio)
#equivalent to z-score
#Based on McKinney et al. 2017 doi:10.1111/1755-0998.12613

#add SD column
hetero_data_chrom_pos_era$read_SD <- sqrt(0.5*(1-0.5)*hetero_data_chrom_pos_era$total_num_reads)

#add z-score (equivalent to D from McKinney et al. 2017)
#describes deviation between observed and expected allelic-specific counts from a binomial distribtion (w/ p = q = 0.5)
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
  scale_x_continuous(limits = c(-25, 25)) +
  labs(y = "percent heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

ggsave(paste('zscore_dist.png', 
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
match_factor<-c(1:1815867) #length of all_data_merged_nona --> create vector to populate
  
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
greenlist_full_contigs <- sort(unique(as.character(greenlist_toprint2$contig)))
  
#write out list without header
write.table(greenlist_toprint, file = "diploid_filtering/greenlist_loci_full_HD_2.5.txt", 
            col.names = FALSE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

###############################################################################################

######## Check distribution of diploid loci across contigs ########
#do these "diploid" loci cluster on the same contigs?
#want only loci on contigs with mostly diploid-passing loci (<80%) --> these most likely to be truly diploid

#read in if running separately
hetero_data_chrom_pos_era <- read.csv("diploid_filtering/meanAB_data_era_pos.csv", header = TRUE)
pass_diploid <- read.table("diploid_filtering/greenlist_loci_full_HD_2.5.txt", header = FALSE) #9251
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
                                      pass_diploid$chrom_pos)) #23552
  colnames(fail_diploid) <- c("chrom_pos")
  fail_diploid <- separate(fail_diploid, chrom_pos, 
                           c("CHROM", "contig", "num", "pos"), 
                           sep = "_", remove = FALSE)
  fail_diploid$chrom <- paste0(fail_diploid$CHROM, "_", 
                               fail_diploid$contig, "_", fail_diploid$num)
  fail_diploid <- subset(fail_diploid, select = c("chrom_pos", "chrom", "pos"))
  fail_diploid$status <- "fail" #add status (these loci failed HD diploid filters)

#### build full dataset of loci on contigs with at leaast one passing locus ####
full_loci_list <- as.data.table(rbind(pass_diploid, fail_diploid))
  full_loci_list <- full_loci_list[order(full_loci_list$chrom),] #32803

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
                                   status_by_chrom_wide$perc_pass > 0.79000000000) #554 contigs
  sum(status_by_chrom_wide_0.8$numloci) #2264 loci
  
#print out list of diploid contigs
#printing out as bed file, will subset VCF to this list of contigs
greenlist_contigs <- status_by_chrom_wide_0.8[, "chrom"]
  greenlist_contigs$chromStart <- 1 #adding pseudo-chromStart column
  greenlist_contigs$chromEnd <- 1000000 #adding pseudo-chromEnd column (needs to be large enough to include all loci)

#write out list with header
write.table(greenlist_contigs, file = "diploid_filtering/diploid_contigs.bed", 
            col.names = TRUE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

####################################################################################################################

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