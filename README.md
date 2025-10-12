# Moss qSIP Analysis Pipelines

This repository contains data and scripts for analysing moss-associated nitrogen-fixing microorganisms using quantitative stable isotope probing (qSIP) and excess atom fraction (EAF) approaches. The pipelines support the identification of active diazotrophs, phylogenetic tree construction, and ranking of taxa by their isotopic enrichment signals.  

All data analyses were conducted in **R** using the [`qSIP`package](https://github.com/bramstone/qsip). For detailed information on the underlying algorithms and functions, please refer to the `qSIP` documentation (https://github.com/bramstone/qsip).

## Structure
```

├── data
│   ├── qSIP_input
│   │   ├── qSIP_data_sheet_OTU_abundance_density.csv                  <-- Input dataset for qSIP calculations, including metadata, fraction density, qPCR, OTU abundance, and taxonomic information. Detailed description in Section 1.2: "Load qSIP data tables and metadata"
│   │   ├── group_sampleID_with_control_and_rep.xlsx                   <-- Metadata file with group, sample ID, and replicate information
│   │
│   ├── qSIP_output
│   │   ├── eaf_taxa_with_1000_bootstrap.csv                           <-- Taxa-specific 15N EAF values at the treatment level based on 1,000 bootstrap resampling iterations
│   │   ├── eaf_taxa_without_bootstrap.csv                             <-- Taxa-specific 15N EAF values at the replicate level without bootstrap resampling
│   │   ├── filtered_data_taxa_signi_over_0.xlsx                       <-- Taxa with EAF values consistently > 0 (with and without 1,000 bootstrap resampling)
│   │   └── wide_OTU_moss_result_after_eaf_over_0.csv                  <-- Matrix with each row representing an OTU and each column representing a sample
│   │
│   ├── Phylogenetic_tree
│   │   ├── annotation_moss.csv                                        <-- OTU × Sample EAF matrix for 509 active OTUs; rows represent OTU IDs, columns represent sample names
│   │   ├── taxa.csv                                                   <-- OTU-to-taxonomy reference table
│   │   ├── select_OTU_moss_509.txt                                    <-- List of 509 active OTU IDs
│   │   ├── all_taxa.fa                                                <-- Representative DNA sequences for all OTUs
│   │   ├── select_OTU_moss_509.fa                                     <-- DNA sequences of the 509 active OTUs
│   │   ├── select_OTU_moss_1.log                                      <-- Log file of phylogenetic tree construction
│   │   └── select_OTU_moss_1.treefile                                 <-- Phylogenetic tree of the 509 active OTUs
│
└── script
    ├── 1.moss_qsip_eaf_process.R                                      <-- R script for qSIP EAF calculation pipeline
    ├── 2.Extract_moss_active_OTU_IDs_for_phylogenetic_tree.R          <-- R script to extract active OTU IDs for phylogenetic analysis
    └── 3.Phylogenetic_tree_construction.sh                            <-- Shell script for phylogenetic tree construction
```


## Contents

### 1. `moss_qsip_eaf_process`
This data analysis process includes the following main steps:

1.1. **Check and load required packages**  
   - Automatically check for missing R packages, install them if necessary, and load (`qsip`, `phyloseq`, `ggplot2`, `dplyr`, etc.).

1.2. **Load metadata and qSIP dataset**  
   - Import metadata (`ide_id`, `sampleIDs`, `fraction`, `timepoint`, `isotope`, `iso_trt`, `ecosystem`, `treatment`, `replicates`).  
   - Load qSIP dataset including fraction density, qPCR, OTU abundance, and taxa data (`otu_id`, `Density.g.ml`, `avg_16S_g_soil`, `seq_abund`,  `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species` ).  
   - Exclude non-target OTUs at kingdom (Unassigned, Eukaryota, Archaea, mitochondria, chloroplast).

1.3. **Data filtering**  
   - Retain OTUs present in at least **three fractions** of **two replicates** per treatment (control or warming, light and label).  
   - Generate a filtered dataset of reproducibly detected OTUs and normalise abundances for downstream analyses.

1.4. **Calculate excess atom fraction (EAF)**  
   - Compute OTU-specific ¹⁵N EAF values by treatment and replicate.  
   - Perform calculations **with and without 1,000 bootstrap resampling iterations**.  

1.5. **Identify active OTUs (EAF > 0)**  
   - Select OTUs with significantly positive EAF values in both bootstrap and non-bootstrap datasets.  
   - Export the final active OTU table in **long** and **wide** formats for downstream tree construction and ranking analyses.

### 2. Extract moss active OTU IDs for phylogenetic tree pipeline
This step prepares the list of moss-associated active OTUs (EAF > 0) for phylogenetic tree construction.

2.1. **Load filtered dataset**  
   - Import the final filtered OTU table (`filtered_data_taxa_pro_signi_over_0.xlsx`).  
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
   - This OTU list will be used to extract sequences from `all_taxa.fa` for tree construction.

### 3. Phylogenetic tree construction
This pipeline builds phylogenetic trees of **active OTUs** by extracting the `select_OTU_moss_509` sequences identified after filtering, aligning them, and constructing a maximum likelihood tree. The resulting tree can be visualized and annotated in iTOL.

3.1. **Extract sequences**  
   - Select OTU sequences (`select_OTU_moss_509.fa`) from the full fasta file using `seqkit grep`.

3.2. **Align sequences**  
   - Perform multiple sequence alignment with MUSCLE.

3.3. **Construct a phylogenetic tree**  
   - Determine the best-fit substitution model using ModelFinder.
   - Construct a maximum likelihood phylogenetic tree in IQ-TREE using the best-fit substitution model (TIM2+F+R6) selected by ModelFinder based on the Bayesian Information Criterion (BIC). Branch support was rigorously evaluated with 1,000 ultrafast bootstrap replicates and 1,000 approximate likelihood ratio test replicates.

3.4. **Visualize and annotate**  
   - Upload the tree to [iTOL](http://itol.embl.de/) and annotate with taxonomy, EAF and ΔEAF.

## Usage
This repository contains all input data, output files, and scripts used for the qSIP analysis of active moss-associated microbial taxa.

## Citation 
Wei Zhou<sup>#</sup>, Yuxuan Bai<sup>#</sup>, Yuhong Xie, Bin Wei, Wolfgang Wanek, Kathrin Rousk, Genevieve Noyce, Dianye Zhang, Yunfeng Peng, and Yuanhe Yang<sup>*</sup>. Key role of moss in satisfying the elevated plant nitrogen demand under warming in a permafrost ecosystem. 

