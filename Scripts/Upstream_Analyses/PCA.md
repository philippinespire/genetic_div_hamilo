# PCA code

Code for calculating eigenvalues & eigenvectors with the program `plink`. Run on an HPC cluster (Wahab at ODU) before transferring output files to local computer for downstream analyses and visualization in R.

## Install plink

Install plink (v.1.9):

```sh
module load container_env python3

crun.python3 -c -s -p ~/.conda/envs/popgen #to create popgen conda env --> only need to do this once

#install plink
crun.python3 -p ~/.conda/envs/popgen conda install -c bioconda -c conda-forge plink
```

## Calculate eigenvalues & eigenvectors

```bash
#NOTE: probably best to grab an interactive node for this (don't run on log-in node).

module load container_env python3

crun.python3 -p ~/.conda/envs/popgen plink --vcf <VCF_FILE> --allow-extra-chr --pca var-wts --out <PIRE.SPECIES.LOC>
```

Copy `*.eigenvec` & `*.eigenval` files to local computer and read into R for downstream analysis/visualization (`Scripts/PCAs.R`).
