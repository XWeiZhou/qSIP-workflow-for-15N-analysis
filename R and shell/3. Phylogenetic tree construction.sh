
#!/bin/bash
    # Set working directory (wd) and software/database (db) directory
    # Add environment variables and enter wd
    wd=/c/A20250331MossSIPTree_all_OTUS
    db=/c/public
    PATH=$PATH:${db}/win

# 1、Evolutionary Tree worddirectory
    cd ${wd}
    
## Match the ID column names selected after SIP-EAF filtering with all OTU.fa sequences,
## extract the corresponding gene sequences for tree construction
  seqkit grep -f select_OTU_moss_509.txt all_taxa.fa > select_OTU_moss.fa 
  
## 2. Phylogenetic tree construction
    # Input files are located in tree directory: 
    #   - otus.fa (sequences) 
    #   - annotation.txt (taxonomy and relative abundance)
    # Step 1: Perform sequence alignment with MUSCLE
    muscle -in select_OTU_moss.fa -out otus_aligned_moss.fas

    # Construct a maximum likelihood phylogenetic tree using IQ-TREE
    # -s   specify aligned sequence input file
    # -bb  set bootstrap resampling to 1000 replicates
    # -alrt set ALRT (Approximate Likelihood-Ratio Test) resampling to 1000 replicates
    # -nt  AUTO: automatically choose number of threads based on available CPUs
    # -pre specify output file prefix
    # -redo rerun analysis even if output files already exist
    
    mkdir -p iqtree
    iqtree -s otus_aligned_moss.fas \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/select_OTU_moss

## 3. Tree visualization and annotation

    # Visit http://itol.embl.de/, upload otus.nwk, 
    # then drag-and-drop the generated annotation schemes onto the tree for visualization.

    ## Scheme 1. Outer-ring colors, shapes, taxonomy, and abundance
    # annotation.txt contains OTU taxonomy annotation and abundance.
    # -a stop if input column not found (default: disabled)
    # -c convert integer columns to factors or decimals to numeric
    # -t convert ID column if label mismatch occurs
    # -w control width of color bands/areas
    # -D output directory
    # -i OTU column name
    # -l OTU display label (e.g., species/genus/family name)	
    Rscript ${db}/script/table2itol.R -a -c double -D plan1 -i OTUID -l Genus -t %s -w 0.5 annotation_moss.txt
    # Each column in the annotation file will be exported as a separate file

    ## Scheme 2. Generate abundance bar plot annotation files
    Rscript ${db}/script/table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation_moss.txt
    Rscript ${db}/script/table2itol.R -a -d -c none -D plan2_2 -b Phylum -i OTUID -l Phylum -t %s -w 0.5 annotation_moss.txt

    ## Scheme 3. Generate heatmap annotation files
    Rscript ${db}/script/table2itol.R -c keep -D plan3 -i OTUID -t %s otutab_moss.txt

    ## Scheme 4. Convert integer values to factors and generate annotation files
    Rscript ${db}/script/table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation_moss.txt

    # Visualize iqtree/otus.contree on http://itol.embl.de/ 
    # and drag files from different Plans onto the tree for annotation

    # Return to working directory
    cd ${wd}
    
###################################################################################

