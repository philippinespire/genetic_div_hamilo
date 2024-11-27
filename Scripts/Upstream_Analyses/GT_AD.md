# Genotype & Depth Information

Code to pull genotype and depth information from VCF files.

```sh
module load container_env ddocent
bash

#set variables
VCFFILE=<VCF FILE>
JUNK_PATTERN=_.*-..-.*-.*_.*_L1_clmp_fp2_repr
NUM_CHR_ID=21
FILE_PREFIX=$(echo $VCFFILE | sed 's/vcf//')

#make header
paste <(echo -e 'chrom\tpos\tref\talt\tqual') <(crun vcf-query $VCFFILE -l | cut -c1-$NUM_CHR_ID | tr "\n" "\t" ) > individuals.tsv

#extract columns of info from VCF (test files first, limits records to 100 SNPs)
#format of resulting files is: CHROM, POS, REF, ALT, QUAL, IND1 ...
cat individuals.tsv <(crun vcf-query $VCFFILE -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%AD]\n' | head -n 100 ) | sed 's/\t$//' > ${FILE_PREFIX}AD.tsv #allele depth info
cat individuals.tsv <(crun vcf-query $VCFFILE -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' | head -n 100 ) | sed 's/\t$//' > ${FILE_PREFIX}GT.tsv #genotype info

#open one of the files up to make sure it looks okay

#extract columns of info from VCF (no limits)
#format of resulting files is: CHROM, POS, REF, ALT, QUAL, IND1 ...
cat individuals.tsv <(crun vcf-query $VCFFILE -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%AD]\n' | sed 's/\t$//' ) > ${FILE_PREFIX}AD.tsv #allele depth info
cat individuals.tsv <(crun vcf-query $VCFFILE -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' | sed 's/\t$//' ) > ${FILE_PREFIX}GT.tsv #genotype info
```

Copied `*tsv` files to local computer and read into R for downstream analyses (`Scripts/sequencing_stats.R`).
