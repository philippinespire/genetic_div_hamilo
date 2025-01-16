Contains code for upstream analyses that create input files for some of the R scripts:

Demographic modeling scripts:
 * **Ela_momi2_code.md:** Contains code (and model results) for running momi2 for *Equulites laterofenestra*. Creates input data for `Demo_Bootstrap.R`.
 * **Ela_momi_bootstrap_contemponly.py:** python script for running momi2 on bootstrapped contemporary SFS for *Equulites laterofenestra*. Creates input data for `Demo_Bootstrap.R`.
 * **Ela_momi_bootstrap_temp.py:** python script for running momi2 on bootstrapped contemporary and historical SFS for *Equulites laterofenestra*. Creates input data for `Demo_Bootstrap.R`.
 * **Ela_momi_bootstrap_temponly.py:** python script for running momi2 on bootstrapped historical SFS for *Equulites laterofenestra*. Creates input data for `Demo_Bootstrap.R`.
 * **Gmi_momi2_code.md:** Contains code (and model results) for running momi2 for *Gazza minuta*. Creates input data for `Demo_Bootstrap.R`.
 * **Gmi_momi_bootstrap_contemponly.py:** python script for running momi2 on bootstrapped contemporary SFS for *Gazza minuta*. Creates input data for `Demo_Bootstrap.R`.
 * **Gmi_momi_bootstrap_temp.py:** python script for running momi2 on bootstrapped contemporary and historical SFS for *Gazza minuta*. Creates input data for `Demo_Bootstrap.R`.
 * **Gmi_momi_bootstrap_temponly.py:** python script for running momi2 on bootstrapped historical SFS for *Gazza minuta*. Creates input data for `Demo_Bootstrap.R`.
 * **momi_bootstrap.sb:** sbatch script to run momi2 for bootstrapping (reads in python scripts). Creates input data for `Demo_Bootstrap.R`.

All other scripts:
 * **ADMIXTURE_code.md:** Contains code for running ADMIXTURE. Creates input files for `admixture.R`.
 * **BayeScan.sbatch:** sbatch script to run BayeScan. Creates input files for `fst.R`.
 * **GT_AD.md:** Contains code for pulling genotype and depth information from VCF files. Creates input files for `sequencing_stats.R`.
 * **PCA.md:** Contains code for running PCAs. Creates input files for `PCAs.R`.
 * **SFS_code.md:** Contains code for visualizing 2D-SFS with dadi.
 * **TajimasD_code.md:** Contains code for calculating Tajima's D with VCFtools. Creates input files for `TajimasD.R`.
 * **pixy.sbatch:** sbatch script to run pixy to generate sliding window estimates of pi and fst. Creates input files for `pi.R` & `fst.R`.
