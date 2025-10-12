
#!/bin/bash
    # Set working directory (wd) and software/database (db) directory
    # Add environment variables and enter wd
    # **The following four lines must be run every time RStudio is opened**
    wd=/c/A20250331MossSIPTree_all_OTUS
    db=/c/public
    PATH=$PATH:${db}/win

# 1、Evolutionary Tree worddirectory
    cd ${wd}
    
## Match the ID column names selected after SIP-EAF filtering with all OTU.fa sequences,
## extract the corresponding gene sequences for tree construction
  seqkit grep -f select_OTU_moss_509.txt all_taxa.fa > select_OTU_moss.fa 
  
## 2. Phylogenetic tree construction
    # Input files are located in result/tree directory: 
    #   - otus.fa (sequences) 
    #   - annotation.txt (taxonomy and relative abundance)
    # Step 1: Perform sequence alignment with MUSCLE
    muscle -in select_OTU_moss.fa -out otus_aligned_moss.fas
    
    mkdir -p iqtree
	
	# Step 1: Model selection using ModelFinder
		# -s otus_aligned_moss.fas  input aligned sequences
		# -m MF     run ModelFinder to select the best substitution model
		# -nt AUTO  automatically detect the number of CPU threads
		# -pre iqtree/select_OTU_moss_modeltest : prefix for output files
	iqtree -s otus_aligned_moss.fas \
		-m MF -nt AUTO \
		-pre iqtree/select_OTU_moss_modeltest
	 	# Output model result (from IQ-TREE):
		# Akaike Information Criterion (AIC): TIM2+F+R6
		# Bayesian Information Criterion (BIC): TIM2+F+R6
		# Best-fit model according to BIC: TIM2+F+R6
		
	# Step 2: Construct the maximum likelihood tree using the best-fit model
		# -s   specify aligned sequence input file
		# -m TIM2+F+R6  use the best-fit substitution model found in Step 1
		# -bb  set bootstrap resampling to 1000 replicates
		# -alrt set ALRT (Approximate Likelihood-Ratio Test) resampling to 1000 replicates
		# -nt  AUTO automatically choose number of threads based on available CPUs
		# -pre prefix for output files of the final tree
		# -redo overwrite existing results with the same prefix
		# -seed random seed for reproducibility
    iqtree -s otus_aligned_moss.fas \
	    -m TIM2+F+R6 \
        -bb 1000 -redo -alrt 1000 -nt AUTO \
        -pre iqtree/select_OTU_moss_2 \
		-seed 13
	
## 3. Tree visualization and annotation

    # Visit http://itol.embl.de/, upload select_OTU_moss_2_13.treefile
    # then drag-and-drop the generated annotation schemes onto the tree for visualization

    # Return to working directory
    cd ${wd}
    
###################################################################################
