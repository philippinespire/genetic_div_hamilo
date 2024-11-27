
#!/bin/bash -l

#SBATCH --job-name=bayescan
#SBATCH -o bayescan.%j.out
#SBATCH -p main
#SBATCH -c 4
##SBATCH --time=48:00:00

/home/rdclark/PIRE/BayeScan2.1/BayeScan2.1/binaries/BayeScan2.1_linux64bits /home/rdclark/PIRE/gazza_minuta/BayeScan/Gmi_SNPs.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 -o /home/rdclark/PIRE/gazza_minuta/BayeScan/Gmi_A_nohighhet_Ham_SNPs_pr100
