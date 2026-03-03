![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES?sort=semver)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES?style=social)](https://github.com/sebastian-gregoricchio/Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES/fork)


# Dauyey, Velds *et al.*<br>A spatial single-cell transcriptomic atlas of metastatic breast cancer progression

## Introduction
Here we deposited the codes used for the the analyses performed in [Dauyey, Velds *et al.*](https://doi.org/XYZ), "A spatial single-cell transcriptomic atlas of metastatic breast cancer progression" (Journal, 202X).

<br>

## Data access
Raw data are available at the BioImage Archive under accession number [S-BIAD2706](https://www.ebi.ac.uk/biostudies/studies/S-BIAD2706).

<br>


## System Requirements
All the codes are run on R (v4.4), therefore any computer can be used to run the analyses.
Installation will take about 30 min.

### R Dependencies
#### CRAN packages
- ggplot2 (>3.2)
- dplyr (1.1.4)
- data.table (1.18.0)
- tidyverse
- devtools


#### Bioconductor packages
- Seurat (5.3.0)
- tradeSeq (1.18.0)
- slingshot (2.12.0)
- clusterProfiler (4.12.6)
- msigdbr (25.1.1)

#### External packages
- spacexr

### Installation Guide:
For Bioconductor packages use:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("package_name")
```

For CRAN packages use:
```
install.packages("package_name")
```

External packages:
```
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

------------------------------------------

#### Contributors
![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES?size=50&padding=5&bots=true)
