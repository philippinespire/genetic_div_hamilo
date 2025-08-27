[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15085163.svg)](https://zenodo.org/doi/10.5281/zenodo.15085163)

This repository provides the data, scripts, and figures for the analyses in "Anthropocene genetic diversity loss in the marine tropics", which trackes genetic diversity loss in, and the demographic history of, populations of two tropical ponyfish from Hamilo Cove, Philippines (*Gazza minuta* and *Equulites laterofenestra*).

A complete list of all necessary software and packages (with version numbers) can be found at the bottom of this README.

This repository is organized into 3 folders:
1. **/Data:** Contains source data for R scripts, filtered VCF files, and metadata. Organized by species.
2. **/Plots:** Contains pngs of plots included in the manuscript.
3. **/Scripts:** Contains R scripts and code for analyses included in the manuscript.

*A record of the bioinformatics pipeline (from raw reads to filtered VCF files) for each species can be found in the relevant repository on the philippinespire GitHub organization page: either [**pire_gazza_minuta_cssl**](https://github.com/philippinespire/pire_gazza_minuta_cssl) or [**pire_leiognathus_leuciscus_cssl**](https://github.com/philippinespire/pire_leiognathus_leuciscus_cssl).* 

Please contact René Clark at rclark848[at]gmail.com with any questions.
_______________________________________________________

**Links to associated metadata:**
1. raw sequences (NCBI GenBank):
    * ***Gazza minuta***:
        * Historical samples: [BioProject PRJNA998057](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA998057/)
        * Contemporary samples: [BioProject PRJNA998814](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA998814/)
    * ***Equulites laterofenestra***:
        * Historical samples: [BioProject PRJNA998845](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA998845/)
        * Contemporary samples: [BioProject PRJNA999299](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA999299/)
2. metadata (GEOME):
   * ***Gazza minuta***:
        * Historical samples: https://n2t.net/ark:/21547/FMB2
        * Contemporary samples: https://n2t.net/ark:/21547/FMH2
   * ***Equulites laterofenestra***:
        * Historical samples: https://n2t.net/ark:/21547/FMR2
        * Contemporary samples: https://n2t.net/ark:/21547/FMX2
3. genome assembly used for mapping & genotype calls (Figshare):
    * ***Gazza minuta***: [DOI:10.6084/m9.figshare.29991607](https://doi.org/10.6084/m9.figshare.29991607.v1)
    * ***Equulites laterofenestra***: [DOI:10.6084/m9.figshare.29991661](https://doi.org/10.6084/m9.figshare.29991661.v1)
_______________________________________________________

**Necessary Programs/Software for Data Analysis**   
1. ADMIXTURE (v.1.3)
2. BayeScan (v.2.1)
3. dadi (v.2.0.5)
4. easySFS (v.1)
5. momi2 (v.2.1.19)
6. NeEstimator (v.2.0.1)
7. pixy (v.1.2.7)
8. PLINK (v.1.9)
9. R (v.4.2.2 or above)
10. RStudio (v.2023.06.1 or above)
11. VCFtools (v.0.1.14)

**Necessary R Packages**
1. adegenet (v.2.1.10)
2. boot (v.1.3.28)
3. data.table (v.1.14.8)
4. devtools (v.2.4.5)
5. FishLife (v.3.0.0)
6. ggridges (v.0.5.4)
7. here (v.1.0.1)
8. hierfstat (v.0.5.11)
9. janitor (v.2.2.0)
10. magrittr (v.2.0.3)
11. plotrix (v.3.8.2)
12. pophelper (v.2.3.1)
13. purrr (v.1.0.1)
14. related (v.1.0)
15. scales (v.1.2.1)
16. tidyverse (v.2.0.0)
17. vcfR (v.1.14.0)
