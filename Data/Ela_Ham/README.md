Contains files either read into R scripts for *Equulites laterofenestra*.

Directories:
 * **ADMIXTURE/:** Directories with ADMIXTURE output files generated with different sets of individuals (all sites & individuals prior to filtering for HWE, post HWE-filtering and LD-pruned sites: all individuals, all *Equulites laterofenestra* individuals, all *Leiognathus leuciscus* individuals, all *Equulites laterofenestra* individuals except those with high heterozygosity). Created by running `ADMIXTURE`. Read into `admixture.R`.
 * **PCAs/:** Contains eigenvalue and eigenvector files generated with different sets of individuals (all sites & individuals prior to filtering for HWE, post HWE-filtering and LD-pruned sites: all individuals, all *Equulites laterofenestra* individuals, all *Leiognathus leuciscus* individuals, all *Equulites laterofenestra* individuals except those with high heterozygosity, all *Equulites laterofenestra individuals except those with high heterozygosity and originally considered an 'outlier' from the first no-cryptids, no-highhet PCA). Created by running `plink`. Read into `PCAs.R`.
 * **momi2/:** Contains csv files with output of bootstrapped momi2 runs, formatted two different ways. Created by running `momi2`. Read into `Demo_Bootstrap.R`.
 * **pixy/:** Contains pi & fst output files (estimates of pi & fst in sliding windows). Created by running `pixy`. Read into `pi.R` & `fst.R`.

All other files:
 * **Ela.nohighhet.Alb.Tajima.D:** Tajima's D estimates (10000 bp windows) for the *Equulites laterofenestra* historical population. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Ela.nohighhet.All.Tajima.D:** Tajima's D estimates (10000 bp windows) for all *Equulites laterofenestra* individuals pooled together. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Ela.nohighhet.Contemp.Tajima.D:** Tajima's D estimates (10000 bp windows) for the *Equulites laterofenestra* contemporary population. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Ela_Ela_nohighhet_SNPs_pr100_fst.txt:** BayeScan output for *Equulites laterofenestra*. Created by running `BayeScan`. Read into `fst.R`.
 * **Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.2.recode.renamed.AD.tsv:** Allele depth information for all SNPs/individuals for *Equulites laterofenestra*, pulled from raw VCF. Read into `sequencing_stats.R`.
   **Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.2.recode.renamed.GT.tsv:** Genotype information for all SNPs/individuals for *Equulites laterofenestra*, pulled from raw VCF. Read into `sequencing_stats.R`.
