TajimasD\_upstream
================

Code for calculating Tajima's D using VCFtools. Run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses and visualization in R.

VCFtools v.1.16 was already installed as a module on Amarel. To load in working environment, use following code:

``` bash
module load VCFtools
```

For Tajima's D, used a VCF file unfiltered for minor allele count (originally 5,729 SNPs), and sub-setted to only synonymous SNPs (to avoid potential issues with purifying selection in coding & UTR regions). Synonymous SNPs were identified by mapping contigs to the *Amphiprion frenatus* genome and then performing a structural annotation on these mapped SNPs using SnpEff. The *A. frenatus* genome is available at <https://datadryad.org/stash/dataset/doi:10.5061/dryad.nv1sv>.

**To map contigs**

Performed BLASTn searches against the *A. frenatus* fasta. Each SNP was mapped to its corresponding bp and contig in *A. frenatus*, and these coordinates were used to create a new VCF file for input into SnpEff.

**To run SnpEff**

SnpEff was downloaded and unzipped in `~/Programs`. As `A.frenatus` is a non-model organism, had to create a new database from the `afrenatus.all_contigs.onlystandardFiltered.gff3` file (downloaded from Dryad).

Added database to `~/Programs/snpEff/snpEff.config`:

``` bash
#in Non-standard Databases section of config file, added following lines:

# Amphiprion frenatus
af1.genome : Amphiprion frenatus
```

Built database in `~/Programs/snpEff/data`:

``` bash
mkdir af1
cd ~/Programs/snpEff/data/af1

#copied afrenatus.all_contigs.onlystandardFiltered.gff3 here and renamed
mv afrenatus.all_contigs.onlystandardFiltered.gff3 genes.gff

#copied Amphiprion_frenatus_GenomeAssembly.fasta here and renamed
mv Amphiprion_frenatus_GenomeAssembly.fasta sequences.fa

#created database
cd ~/Programs/snpEff
java -jar snpEff.jar build -gff3 -v af1 #works if creates snpEffectPredictor.bin in data/af1
```

Created `snp.sbatch` file to run SnpEff:

``` bash
#!/bin/bash

#SBATCH --partition=p_mlp195 #p_mlp195 for pinsky lab partition, otherwise use main (or EOAS or E&E)
#SBATCH --job-name=snpeff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=2:00:00

#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=rene.clark@rutgers.edu

#script for running snpEff on Amarel

java -Xmx8g -jar snpEff.jar af1 data/af1/output.hicov2.snps.only_afrenatus_FULL.vcf > afrenatus_ann_FULL.vcf
```

Submitted job:

``` bash
sbatch snpEff.sbatch
```

Creates three files:

1.  `afrenatus_ann_FULL.vcf`: vcf file with the structural annotation included.
2.  `snpEff_genes_FULL.txt`: breakdown of \# and type of SNPs by contig.
3.  `snpEff_summary_FULL.html`: summary of \# and type of SNPs in entire VCF

Took `afrenatus_ann_FULL.vcf`, converted back to *A. clarkii* bp & contig coordinates, and subsetted to only synonymous SNPs.

**To calculate Tajima's D**

First, created sub-setted VCF file with only synonymous SNPs:

``` bash
#created synonymous_SNPs.txt (column w/contig & column w/bp of all synonymous SNPs)

#used VCFtools to filter output.hicov2.snps.only.mac1.vcf to only include synonymous SNPs
vcftools --vcf output.hicov2.snps.only.mac1.vcf --positions synonymous_SNPs.txt --recode --recode-INFO-all --out output.hicov2.snps.only.mac1.synonymous
```

Then, sub-setted `output.hicov2.snps.only.mac1.synonymous.vcf` for each sampling site (code for Japan below):

``` bash
#created J_individuals.txt (column of individuals in Japan)

#used VCFtools to filter output.hicov2.snps.only.mac1.synonymous.vcf to only include individuals from Japan
vcftools --vcf output.hicov2.snps.only.mac1.synonymous.vcf --keep J_individuals.txt --recode --recode-INFO-all --out Jmac1_SYN

#repeated for Philippines & Indonesia
```

Calculated Tajima's D, both within each population and with the full dataset (all individuals included in VCF):

``` bash
#NOTE: don't grab interactive node for VCFtools jobs bc run in seconds, BUT if is more data-intensive, should get off log-in node to do this.

vcftools --vcf Jmac1_SYN.vcf --TajimaD 10000 --out SYN_mac1_JTajimasD

#used bins of 10,000 bp so would calculate across each transcript individually (largest transcript <7000 bp)
#repeated for Philippines, Indonesia, and the original vcf file w/all individuals
```

Copied `*Tajima.D` files to local computer and opened in Excel to create .csv files that were then read into R for downstream analyses (`Scripts/TajimaD_script.R`).
