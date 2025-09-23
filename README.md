# Moss qSIP Analysis Pipelines

This repository contains scripts and workflows for analysing moss-associated nitrogen-fixing microorganisms using quantitative stable isotope probing (qSIP) and excess atom fraction (EAF) approaches. The pipelines support the identification of active diazotrophs, phylogenetic tree construction, and ranking of taxa by their isotopic enrichment signals.  

All data analyses were conducted in **R** using the [`qSIP`package](https://github.com/bramstone/qsip). For detailed information on the underlying algorithms and functions, please refer to the `qSIP` documentation.

## Contents

### 1. `moss_qsip_eaf_process`
This data analysis process includes the following main steps:
1.1. **Check and load required packages**  
   - Automatically check for missing R packages, install them if necessary, and load them (including `qsip`, `phyloseq`, `ggplot2`, `dplyr`, etc.).
1.2. **Load qSIP data tables and metadata**  
   - Import sample metadata (ide_id, sampleIDs, fraction, timepoint, isotope, iso_trt, ecosystem, treatment, replicates).  
   - Load qSIP qPCR and OTU abundance–density data.  
   - Exclude non-target OTUs (Unassigned, Eukaryota, Archaea, mitochondria, chloroplast) to retain only bacterial OTUs.
1.3. **Perform data filtering and quality control**  
   - Retain OTUs present in at least **three qSIP fractions** of **two replicates** per treatment (control or warming, light and label).  
   - Generate a filtered dataset containing only reproducibly detected OTUs across treatments.  
   - Normalise abundances for downstream analysis.
1.4. **Calculate Excess Atom Fraction (EAF)**  
   - Compute OTU-specific ¹⁵N EAF values at the treatment and replicate levels.  
   - Perform calculations both **with and without 1,000 bootstrap resampling iterations** for robustness.  
   - Export results as `.csv` tables.
1.5. **Identify active OTUs with EAF > 0**  
   - Select OTUs with significantly positive EAF values in both bootstrap and non-bootstrap datasets.  
   - Export final active OTU table in both **long** and **wide** formats for downstream phylogenetic tree construction and ranking analyses.

### 2. Extract moss active OTU IDs for phylogenetic tree pipeline
This step prepares the list of moss-associated active OTUs (with significant EAF > 0) for downstream phylogenetic tree construction.
2.1. **Load filtered EAF dataset**  
   - Import the final filtered OTU table (`filtered_data_taxa_pro_signi_over_0_20250923.xlsx`).  
   - Convert `Treatment` to a factor (Control, Warming).  
   - Retain OTUs with positive EAF values from moss samples.
2.2. **Calculate treatment-specific mean EAF values**  
   - Group OTUs by `otu_id` and treatment.  
   - Compute the mean EAF per OTU for Control and Warming treatments.  
   - Reshape to wide format and replace missing values with zero.  
   - Calculate the difference in enrichment (`ΔEAF = Warming – Control`).
2.3. **Annotate with taxonomy**  
   - Merge the EAF table with taxonomic assignments (`taxa.csv`).  
   - Generate a final annotation table (`annotation_moss.csv`) containing OTU IDs, EAF values, ΔEAF, and taxonomy.
2.4. **Export OTU ID list**  
   - Extract active OTU IDs into a text file (`select_OTU_moss_509.txt`).  
   - This OTU list will be used to match sequences in `otus_taxa.fa` for phylogenetic tree construction.
   - 
### 3. Phylogenetic tree construction
This pipeline builds phylogenetic trees of **active diazotrophs** by extracting the **active** select_OTU_moss_509 sequences identified after SIP–EAF filtering, aligning them, and constructing a maximum likelihood tree. The resulting tree can be visualized and annotated in iTOL.
3.1. **Extract sequences** 
   - Selected OTUs from the full fasta file (`seqkit grep`).
3.2. **Align sequences**
   - Aligned with MUSCLE.
3.3. **Construct a phylogenetic tree**
   - Phylogenetic tree using IQ-TREE with 1000 bootstrap and ALRT replicates.
3.4. **Visualize and annotate the tree**
   -Visualised and annotated in [iTOL](http://itol.embl.de/) using taxonomy, abundance, and custom schemes generated with `table2itol.R`.

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

