# Moss qSIP Analysis Pipelines

This repository contains scripts and workflows for analyzing moss-associated nitrogen-fixing microorganisms using quantitative stable isotope probing (qSIP) and excess atom fraction (EAF) approaches. The pipelines support the identification of active diazotrophs, phylogenetic tree construction, and ranking of taxa by their isotopic enrichment signals.

## Contents

### 1. `moss_qsip_eaf_process`
Our data analysis includes the following main steps:
1. **Check and load required packages**  
   - Automatically check for missing R packages, install them if necessary, and load them (including `qsip`, `phyloseq`, `ggplot2`, `dplyr`, etc.).
2. **Load qSIP data tables and metadata**  
   - Import sample metadata (sample IDs, treatments, replicates, isotope labeling).  
   - Load qSIP OTU abundance–density data.  
   - Exclude non-target OTUs (Unassigned, Eukaryota, Archaea, mitochondria, chloroplast) to retain only bacterial OTUs.
3. **Perform data filtering and quality control**  
   - Retain OTUs present in at least **three qSIP fractions** of **two replicates** per treatment (control or warming, light and label).  
   - Generate a filtered dataset containing only reproducibly detected OTUs across treatments.  
   - Normalize abundances for downstream analysis.
4. **Calculate EAF (Excess Atom Fraction)**  
   - Compute OTU-specific ¹⁵N EAF values at the treatment and replicate levels.  
   - Perform calculations both **with and without 1,000 bootstrap resampling iterations** for robustness.  
   - Export results as `.csv` tables.
5. **Identify active OTUs with EAF > 0**  
   - Select OTUs with significantly positive EAF values in both bootstrap and non-bootstrap datasets.  
   - Apply significance filtering (e.g., median EAF > 0, p < 0.05).  
   - Export final active OTU table in both **long** and **wide** formats for downstream phylogenetic tree construction and ranking analyses.

### 2. Extract moss active OTU IDs for phylogenetic tree pipeline
Pipeline for extracting **active OTUs/ASVs** identified by qSIP-EAF analysis to be used in downstream phylogenetic analyses.
- Input: OTU/ASV table with EAF results
- Output: A filtered list of **active OTU IDs** for tree building
- Functions:
  - Threshold-based selection of active OTUs  
  - Export of FASTA sequences corresponding to active OTUs  

### 3. Phylogenetic tree construction
Pipeline for building phylogenetic trees of **active diazotrophs**.
- Input: FASTA sequences of active OTUs
- Output: Phylogenetic tree (Newick + annotated tree figure)
- Functions:
  - Multiple sequence alignment (MSA)  
  - Model selection and tree inference (e.g., IQ-TREE, RAxML)  
  - Visualization of active diazotroph placement in broader community context  

### 4. Ranked EAF OTUs with interval pipeline
Pipeline for ranking **active OTUs** based on their isotopic enrichment signal and confidence intervals.
- Input: qSIP-EAF results with bootstrapping
- Output: Ranked OTU list with confidence intervals (CSV + plots)
- Functions:
  - Calculation of mean EAF and 95% confidence intervals  
  - Ranking of OTUs by isotopic activity  
  - Visualization of ranked contribution to nitrogen fixation  

## Requirements
- R (≥4.0) with packages: `dplyr`, `tidyr`, `ggplot2`, `phyloseq`
- Python (≥3.8) with packages: `pandas`, `biopython`
- Tree building software: IQ-TREE / RAxML
- qSIP analysis tools  


## Usage
Each pipeline is organized in its own directory with scripts and example data.  
Refer to individual pipeline folders for detailed instructions.


## Citation
If you use this repository, please cite:  
*Zhou W., et al. (in prep). Moss-associated nitrogen fixation under global change revealed by quantitative stable isotope probing.

