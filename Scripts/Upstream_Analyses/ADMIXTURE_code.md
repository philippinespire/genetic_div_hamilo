# ADMIXTURE code

Code for running `ADMIXTURE` to assess population structure. Requires input (`bim` & `bed` files) created with PLINK. Instructions for installing PLINK & ADMIXTURE can be found in `/popgen_analyses/README.md`.

## Install plink and admixture

Install `plink` (v.1.9):

```sh
module load container_env python3

crun.python3 -c -s -p ~/.conda/envs/popgen #to create popgen conda env --> only need to do this once

#install plink
crun.python3 -p ~/.conda/envs/popgen conda install -c bioconda -c conda-forge plink
```

Install `ADMIXTURE` (v.1.3): 

```sh
module load container_env python3

#install admixture
#assumes popgen conda environment has already been created
crun.python3 -p ~/.conda/envs/popgen conda install -c bioconda admixture
```

## Create input files for ADMIXTURE

Do this with `plink`.

```sh
#NOTE: probably best to grab an interactive node for this (don't run on log-in node).

module load container_env python3

#create input files
crun.python3 -p ~/.conda/envs/popgen plink --vcf <VCF_FILE> --allow-extra-chr --make-bed --out <PIRE.SPECIES.LOC>

#transform bim file so that can read into ADMIXTURE
awk '{$1=0;print $0}' PIRE.SPECIES.LOC.bim > PIRE.SPECIES.LOC.bim.tmp
mv PIRE.SPECIES.LOC.bim.tmp PIRE.SPECIES.LOC.bim
```

## Run ADMIXTURE

```sh
crun.python3 -p ~/.conda/envs/popgen admixture PIRE.SPECIES.LOC.bed 1 --cv > PIRE.SPECIES.LOC.log1.out
#run from 1-5
```

Copy `*.Q` files to local computer. Read into R for visualization (`popgen_analyses/pop_structure.R`).
