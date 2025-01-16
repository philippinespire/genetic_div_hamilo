# Tajima's D code

Code for calculating Tajima's D using `VCFtools`. Run on an HPC cluster (Wahab at ODU) before transferring output files to local computer for downstream analyses and visualization in R.

`VCFtools` v.0.1.16 was already installed as a module on Wahab. To load in working environment, use following code:

``` bash
module load container_env bcftools
```

For Tajima's D, used the "all sites" VCF files that include both monomorphic and polymorphic sites. For each species, subsetted VCF file to either contemporary or historical individuals.

Calculated Tajima's D, both within each population and with the full dataset (all individuals included in VCF):

``` bash
#NOTE: don't grab interactive node for VCFtools jobs bc run in seconds, BUT if is more data-intensive, should get off log-in node to do this.

vcftools --vcf <VCF FILE> --TajimaD 10000 --out <OUT PREFIX>
```

Copied `*Tajima.D` files to local computer and read into R for downstream analyses (`Scripts/TajimasD.R`).
