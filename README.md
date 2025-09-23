# Moss qSIP Analysis Pipelines

This repository contains scripts and workflows for analysing moss-associated nitrogen-fixing microorganisms using quantitative stable isotope probing (qSIP) and excess atom fraction (EAF) approaches. The pipelines support the identification of active diazotrophs, phylogenetic tree construction, and ranking of taxa by their isotopic enrichment signals.  

All data analyses were conducted in **R** using the [`qSIP`package](https://github.com/bramstone/qsip). For detailed information on the underlying algorithms and functions, please refer to the `qSIP` documentation.

## Contents

### 1. `moss_qsip_eaf_process`
This data analysis process includes the following main steps:

1.1. **Check and load required packages**  
   - Automatically check for missing R packages, install them if necessary, and load (`qsip`, `phyloseq`, `ggplot2`, `dplyr`, etc.).

1.2. **Load qSIP data tables and metadata**  
   - Import metadata (`ide_id`, `sampleIDs`, `fraction`, `timepoint`, `isotope`, `iso_trt`, `ecosystem`, `treatment`, `replicates`).  
   - Load fraction density, qPCR, OTU abundance, and taxa data (`otu_id`, `Density.g.ml`, `avg_16S_g_soil`, `seq_abund`,  `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species` ).  
   - Exclude non-target OTUs at kingdom (Unassigned, Eukaryota, Archaea, mitochondria, chloroplast).

1.3. **Data filtering and quality control**  
   - Retain OTUs present in at least **three fractions** of **two replicates** per treatment (control or warming, light and label).  
   - Generate a filtered dataset with reproducibly detected OTUs.  
   - Normalize abundances for downstream analysis.

1.4. **Calculate excess atom fraction (EAF)**  
   - Compute OTU-specific ¹⁵N EAF values by treatment and replicate.  
   - Perform calculations **with and without 1,000 bootstrap resampling iterations**.  
   - Export results as `.csv` tables.

1.5. **Identify active OTUs (EAF > 0)**  
   - Select OTUs with significantly positive EAF values in both bootstrap and non-bootstrap datasets.  
   - Export the final active OTU table in **long** and **wide** formats for downstream tree construction and ranking analyses.

### 2. Extract moss active OTU IDs for phylogenetic tree pipeline
This step prepares the list of moss-associated active OTUs (EAF > 0) for phylogenetic tree construction.

2.1. **Load filtered dataset**  
   - Import the final filtered OTU table (`filtered_data_taxa_pro_signi_over_0_20250923.xlsx`).  
   - Convert `Treatment` to a factor (`Control`, `Warming`).  
   - Retain OTUs with positive EAF values in moss samples.

2.2. **Calculate treatment-specific mean EAF**  
   - Group OTUs by `otu_id` and treatment.  
   - Compute mean EAF per OTU for Control and Warming.  
   - Reshape to wide format, replace missing values with zero.  
   - Calculate enrichment difference (`ΔEAF = Warming – Control`).

2.3. **Annotate with taxonomy**  
   - Merge EAF values with taxonomic assignments (`taxa.csv`).  
   - Generate a final annotation table (`annotation_moss.csv`) including OTU IDs, EAF, ΔEAF, and taxonomy.

2.4. **Export OTU ID list**  
   - Extract OTU IDs to a text file (`select_OTU_moss_509.txt`).  
   - This OTU list will be used to extract sequences from `otus_taxa.fa` for tree construction.

### 3. Phylogenetic tree construction
This pipeline builds phylogenetic trees of **active diazotrophs** by extracting the `select_OTU_moss_509` sequences identified after SIP–EAF filtering, aligning them, and constructing a maximum likelihood tree. The resulting tree can be visualized and annotated in iTOL.

3.1. **Extract sequences**  
   - Select OTU sequences from the full fasta file using `seqkit grep`.

3.2. **Align sequences**  
   - Perform multiple sequence alignment with MUSCLE.

3.3. **Construct a phylogenetic tree**  
   - Build a maximum likelihood tree with IQ-TREE using 1000 bootstrap and ALRT replicates.

3.4. **Visualize and annotate**  
   - Upload the tree to [iTOL](http://itol.embl.de/) and annotate with taxonomy, abundance, and custom schemes generated using `table2itol.R`.

### 4. Ranked EAF OTUs with interval pipeline
Pipeline for ranking **active OTUs** based on their isotopic enrichment signal and confidence intervals.

This pipeline visualizes **ranked moss OTUs** based on their estimated ¹⁵N enrichment (EAF) values and confidence intervals, separately for control and warming treatments. The plots display OTUs ordered by rank within each phylum, with error bars representing uncertainty.

4.1. **Load and prepare data**  
   - Import taxonomic annotation (`taxa.csv`) and select relevant columns.  
   - Import EAF data calculated with 1,000 bootstrap resampling (`eaf_pro_1000.csv`).  
   - Merge EAF values with taxonomy to create a combined dataset.

4.2. **Subset treatments and rank OTUs**  
   - Filter OTUs for `Control` and `Warming` moss treatments separately.  
   - Order OTUs by **Phylum** and **median EAF (X50.)** within each group.  
   - Assign rank IDs to each OTU for plotting.

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
*Wei Zhou, Yuxuan Bai, Yuhong Xie, Bin Wei, Wolfgang Wanek, Kathrin Rousk, Genevieve Noyce, Dianye Zhang, Yunfeng Peng, and Yuanhe Yang. Key role of moss in satisfying the elevated plant nitrogen demand under warming in a permafrost ecosystem. 2025 

