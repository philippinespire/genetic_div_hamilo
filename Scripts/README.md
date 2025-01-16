Scripts:
 * **admixture.R**: Reads in `ADMIXTURE` output files and creates output plots for visualization of `ADMIXTURE` results. Also creates CV plots to identify the "best" value of K.
 * **compare_div_loss.R**: Compares observed loss in genetic diversity in these species to other temporal losses (using supplemental material from [Leigh et al. (2019)](https://doi.org/10.1111/eva.12810).
 * **Demo_Bootstrap.R**: Reads in output from `momi2` and calculates 95% CIs for Ne estimates and timing of Ne size changes; creates plots of Ne change through time.
 * **div_diff_bootstrap.R**: Creates null distributions of genetic diversity estimates.
 * **diversity.R**: Calculates Ho, He, and Fis from allele frequencies along with bootstrapped 95% CIs. Compares point estimates to null distributions.
 * **fst.R**: Calculates pairwise Fst. Reads in `pixy` output to identify high Fst windows. Reads in `BayeScan` output to identify candidate SNPs under selection.
 * **PCAs.R**: Reads in eigenvector information from `plink` and creates PCA plots.
 * **pi.R**: Reads in `pixy` output and calculates mean pi per timepoint along with bootstrapped 95% CI. Also creates plots of diversity changes (Ho, He, Fis, pi) through time for both species together.
 * **relatedness.R**: Calculates pairwise relatedness (point estimates for all possible pairs, mean within- and across-time point relatedness, and bootstrapped 95% CIs).
 * **sequencing_stats.R**: Calculates read depth v. heterozygosity stats to look for issues with reference bias and under-powered genotype calls. Also generates various plots to analyze sequencing quality (amount of missing data, etc.).
 * **TajimasD.R**: Reads in .csv files containing Tajima's D output from `VCFtools` and calculates the mean (+ SE) Tajima's D in each timepoint and with both timepoints pooled. Also analyzes change in the distribution of Tajima's D through time.
 * **tempNe_estimates.R**: Calculates harmonic mean of Ne (over time) based on heterozygosity loss equation.

*If input files for a script were created by running code and/or calling programs on a remote workstation, details on the code used can be found in /Scripts/Upstream_Analyses.*
