This repository contains all data and code used to generate the analyses and figures presented in the manuscript:

*“Rising widespread Kelch13 R622I and Detection of C580Y of African origin in *Plasmodium falciparum* in northern Ethiopia further threaten artemisinin-based therapies”*

GT will add preprint author list here when uploaded and ultimately the citation to peer reviewed paper


## Repository Structure

```bash
.
├── Fig1/                 # Scripts and outputs for Figure 1 (Prevalence and geography)
├── Fig2/                 # Scripts and outputs for Figure 2 (Haplotype comparison)
├── Fig3/                 # Scripts and outputs for Figure 3 (Drug resistance combinations)
├── Fig4/                 # Scripts and outputs for Figure 4 (IBD and clonal network analyses)
├── Fig5/                 # Scripts and outputs for Figure 5 (Temporal trends and relatedness)
├── data/                 # Input data (genotype data, metadata, COI, IBD, rainfall, case counts, mapping geojson)
└── README.md             # Project overview and usage instructions

```

## Getting Started

1. Clone the repository
```bash
git clone https://github.com/gtollefson/Northern_Ethiopia_Pf_ARTR_2025_Project.git
cd Northern_Ethiopia_Pf_ARTR_2025_Project
```

2. Install dependencies

R (Figures 1, 3, 4, 5)
Ensure the following R packages are installed:
```r
install.packages(c(
  "dplyr", "ggplot2", "ggpubr", "ggsignif", "igraph", "ggraph",
  "tidygraph", "rstatix", "sf", "readr", "stringr", "tidyr", "viridis", "terra"
))
```

Python (Figure 2)
Python dependencies (e.g. pandas, matplotlib, scikit-learn) can be installed via:
```bash
pip install numpy==1.24.4 pandas==1.5.3 matplotlib==3.7.1 seaborn==0.12.2 scipy==1.10.1 tqdm==4.65.0

```

## Figures Overview

**Figure 1:**  
Geographic and statistical summary of K13 622I mutation prevalence and pfhrp2-based RDT negativity across districts, kebeles, and RDT result categories in Northwest Ethiopia.

**Figure 2:**  
Flanking haplotype comparison shows similarity between Ethiopian and Southeast Asian K13 C580Y mutants, highlighting potential shared ancestry or convergent evolution.

**Figure 3:**  
Upset plots display combinations of K13 622I with _mdr1_ and _crt_ resistance mutations in monoclonal P. falciparum isolates from Ethiopia.

**Figure 4:**  
Genetic relatedness of K13 mutant and wildtype parasites is illustrated through IBD sharing and network analysis stratified by RDT result, district, and kebele.

**Figure 5:**  
Temporal dynamics of K13 622I prevalence, RDT negativity, case counts, rainfall, and parasite relatedness across months and genotypes in Northwest Ethiopia.

## Data Sources
All necessary input files are stored in /data/ and include:

* Unique molecular identifying PCR-amplification bias corrected genotyping tables from MIPTools
* Sample Metadata
* Complexity of infection (COI) estimates
* Identity-by-descent (IBD) estimates
* Rainfall and case surveillance data
* Mapping file (geojson at country and woredas(district) level)

You will need to adjust file paths in the scripts to reflect your local installation.

## Citation
If using this repository or its outputs, please cite: Will add preprint citation here when uploaded and ultimately the citation to peer reviewed paper

## Contact
For questions about the codebase or contributions, contact the corresponding authors in the linked published mansucript.