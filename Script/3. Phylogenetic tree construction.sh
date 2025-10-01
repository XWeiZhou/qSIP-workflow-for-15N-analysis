
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
        
    # Return to working directory
    cd ${wd}

## 3. Tree visualisation and annotation

    Visit http://itol.embl.de/, upload **select_OTU_moss.contree**, 
    Then drag-and-drop the generated annotation schemes onto the tree for visualisation.
 
###################################################################################


