#!/bin/bash -l

#SBATCH --job-name=pixy
#SBATCH -o pixy-%j.out
#SBATCH -p main
#SBATCH -c 4

#to run use the following command:
#sbatch pixy.sbatch

module load container_env python3

#code to run pixy
bash
export SINGULARITY_BIND=/home/r3clark/.conda

crun.python3 -p ~/.conda/envs/popgen \
  pixy --stats pi fst \
       --vcf ../gmi.all.A.nohighhet.Ham.vcf.gz \
       --populations allsites_popfile_nohighhet_Ham_A.txt \
       --window_size 1000 \
       --n_cores 4
