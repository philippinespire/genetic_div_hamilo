#!/bin/bash -l

#SBATCH --job-name=pixy
#SBATCH -o pixy-%j.out
#SBATCH -p main
#SBATCH -c 4

#to run use the following command:
#sbatch pixy.sbatch

module load container_env python3

#install pixy v.1.2.7 (if necessary)
crun.python3 -c -s -p ~/.conda/envs/popgen #create popgen env if not already done
crun.python3 -p ~/.conda/envs/popgen conda install -c conda-forge pixy
crun.python3 -p ~/.conda/envs/popgen conda install -c bioconda htslib

#code to run pixy on Wahab
bash
crun.python3 -p ~/.conda/envs/popgen \
  pixy --stats pi fst \
       --vcf ../gmi.all.A.nohighhet.Ham.vcf.gz \
       --populations allsites_popfile_nohighhet_Ham_A.txt \
       --window_size 1000 \
       --n_cores 4
