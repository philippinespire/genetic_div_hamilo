Contains files either read into R scripts for *Gazza minuta*.

Directories:
 * **ADMIXTURE/:** Directories with ADMIXTURE output files generated with different sets of individuals (all sites & individuals prior to filtering for HWE, post HWE-filtering and LD-pruned sites: all individuals, all Hamilo Cove individuals, all species A individuals, all species B individuals, all species A individuals only from Hamilo Cove, all species A individuals except those with high heterozygosity, all species A individuals from Hamilo Cove except those with high heterozygosity, all species A individuals from Hamilo Cove except those with high heterozygosity and originally considered an 'outlier' from the first no-cryptids, no-highhet, Hamilo-only PCA). Created by running `ADMIXTURE`. Read into `admixture.R`.
 * **PCAs/:** Contains eigenvalue and eigenvector files generated with different sets of individuals (all sites & individuals prior to filtering for HWE, post HWE-filtering and LD-pruned sites: all individuals, all Hamilo Cove individuals, all species A individuals, all species B individuals, all species A individuals only from Hamilo Cove, all species A individuals except those with high heterozygosity, all species A individuals from Hamilo Cove except those with high heterozygosity, all species A individuals from Hamilo Cove except those with high heterozygosity and originally considered an 'outlier' from the first no-cryptids, no-highhet, Hamilo-only PCA). Created by running `plink`. Read into `PCAs.R`.
 * **momi2/:** Contains csv files with output of bootstrapped momi2 runs for *Gazza minuta*, formatted two different ways. Created by running `momi2`. Read into `Demo_Bootstrap.R`.
 * **pixy/:** Contains pi & fst output files (estimates of pi & fst in sliding windows) for *Gazza minuta*. Created by running `pixy`. Read into `pi.R` & `fst.R`.

All other files:
 * **Gmi.A.nohighhet.Ham.Alb.Tajima.D:** Tajima's D estimates (10000 bp windows) for the *Gazza minuta* historical population. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Gmi.A.nohighhet.Ham.All.Tajima.D:** Tajima's D estimates (10000 bp windows) for all *Gazza minuta* individuals pooled together. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Gmi.A.nohighhet.Ham.Contemp.Tajima.D:** Tajima's D estimates (10000 bp windows) for the *Gazza minuta* contemporary population. Created by running `VCFtools`. Read into `TajimasD.R`.
 * **Ela_Ela_nohighhet_SNPs_pr100_fst.txt:** BayeScan output for *Gazza minuta*. Created by running `BayeScan`. Read into `fst.R`.
 * **Gmi.A.rad.RAW-10-10-rescaled.Fltr17.2.recode.renamed.AD.tsv:** Allele depth information for all SNPs/individuals for *Gazza minuta*, pulled from raw VCF. Read into `sequencing_stats.R`.
 * **Gmi.A.rad.RAW-10-10-rescaled.Fltr17.2.recode.renamed.GT.tsv:** Genotype information for all SNPs/individuals for *Gazza minuta*, pulled from raw VCF. Read into `sequencing_stats.R`.
 * **Gmi.A.rad.RAW-10-10-rescaled.Fltr17.2.recode.renamed.vcf:** Raw VCF for *Gazza minuta*. Read into `diversity.R` & `fst.R`.
 * **Gmi_A_nohighhet_HamNe.txt:** NeEstimator output for *Gazza minuta*. Data from this incorporated into `tempNe_estimates.R`.
 * **Gmi_A_nohighhet_Ham_SNPs_pr100_fst.txt:** BayeScan output for *Gazza minuta*. Created by running `BayeScan`. Read into `fst.R`.
