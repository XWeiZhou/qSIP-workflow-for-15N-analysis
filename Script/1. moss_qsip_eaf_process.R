######################################################################################################################################
## Our data analysis includes the following steps:
# 1. Check, install (if missing), and load required packages
# 2. Load qSIP data tables and choice the necessary data
# 3. Perform data filtering and quality control
# 4. calculate the EAF with and without the 1,000 bootstrap resampling iterations
# 5. Find  both with and without 1,000 bootstrap resampling EAF value all over 0

######################################################################################################################################
# 1. Check and load required packages
# Install any missing packages, then load them
check_and_install <- 
  function(packages) { for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "phyloseq") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("phyloseq")
      } else if (pkg == "qsip") {
        if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
        devtools::install_github("bramstone/qsip")
      } else {
        install.packages(pkg) }}
    library(pkg, character.only = TRUE) }}

# check，install and load the library
check_and_install(c( "data.table","qsip","ggplot2","phyloseq","dplyr","clipr","ggpointdensity","knitr","kableExtra","readxl","openxlsx", "tidyr"))

# Record the start time of the script
start_time_all <- Sys.time()

##################################################################################################
# 2.  Load qSIP data tables and choice the necessary data
# OTUs that lacked taxonomic assignment at the kingdom level, or were classified as 
# eukaryotic, archaeal, mitochondrial, chloroplast, or unassigned sequence reads were excluded from the dataset
setwd("C:/Users/Lenovo/Desktop/code/")

## Load sample metadata including sampleID, ID, Treatment, replicates, and isotope label treatment 
group_sampleID_with_control_and_rep <- read_excel("group_sampleID_with_control and rep.xlsx", sheet=1)

## Load the qSIP data and remove unnecessary columns
inital.qsip <- read.csv("qSIP_data_sheet_OTU_abundance_density.csv")[,c(-1,-3)] 

## Keep only rows corresponding to the moss ecosystem
inital.qsip <- inital.qsip[inital.qsip$ecosystem == "moss", ]
dim(inital.qsip)

## Remove the last 'Species' column
qsip <- inital.qsip[,1:dim(inital.qsip)[2]-1]

## Convert the data.frame to data.table format and inspect the structure
moss_qsip <- as.data.table(qsip)
str(moss_qsip)

## Convert sequence abundance to integer
moss_qsip$seq_abund <- as.integer(moss_qsip$seq_abund)
ssc <- c('seq_abund', 'otu_id')

# initial sequence and OTU count
cat('Initial\n')
seq_summary(moss_qsip[, ..ssc])
# We will get the summary of the initial sequences and feature
# 12,574,266 sequences
# 46,763 tax features

# Identify and remove non-target OTUs (Unassigned, Eukaryota, Archaea, mitochondrial/chloroplast) 
# to retain only bacterial OTUs for downstream analyses.
tx <- unique(moss_qsip[, c('otu_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
unassign <- tx[Kingdom == 'Unassigned', otu_id]
euk <- tx[Kingdom == 'Eukaryota', otu_id]
arch <- tx[Kingdom == 'Archaea', otu_id]
mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')][grepl('mitochond|chloroplast', tax_string, ignore.case = T), otu_id]

# Summarize non-target OTUs (Unassigned, Eukarya, Archaea, mitochondrial, chloroplast) 
# to decide on removal before downstream analyses.
cat('Unassigned\n')
seq_summary(moss_qsip[otu_id %in% unassign, ..ssc])
cat('Eukarya\n')
seq_summary(moss_qsip[otu_id %in% euk, ..ssc])
cat('Archaea\n')
seq_summary(moss_qsip[otu_id %in% arch, ..ssc])
cat('Mitochondria, Chloroplasts\n')
seq_summary(moss_qsip[otu_id %in% mito_chloro, ..ssc])

# Filter out lineages
moss_qsip <- moss_qsip[!otu_id %in% c(arch, mito_chloro, euk, unassign)]
cat('Final after quality filtering\n')
seq_summary(moss_qsip[,..ssc])
# Summary of filtered data: 12,574,266 sequences and 46,763 taxonomic features

# relativize sequence abundances (should be done after taxonomic filtering)
# normalize 16S abundances of bacteria
moss_qsip[, rel_abund := seq_abund/sum(seq_abund), by = sampleID]

###############################################################################################################
# 3. Perform data filtering
# A taxon had to be present in at least three qSIP fractions of two replicates 
# for both natural abundance and 15N-enriched treatments per control or warming treatment 

# initial sequence and OTU count
cat('Before 3 fractions filtering\n')
seq_summary(moss_qsip[, ..ssc])

# Remove taxa that occur in fewer than 3 fractions in any given replicate
moss_qsip_after3fraction <- freq_filter(moss_qsip,  min_freq = 3, filter_target = 'sampleID', tax_id = 'otu_id')
cat('After 3 fractions filtering\n')
seq_summary(moss_qsip_after3fraction[, ..ssc])
# After fraction filtering, 
# 11636347 sequences,
# 7172 tax features

# Since the data is hierarchical, extract unique OTU IDs.
# Split the data frame by sampleID into multiple subsets, then find unique OTU IDs in each subset.
# Finally, combine these unique OTU IDs into a single list.
unique_otu_ids <- lapply(split(moss_qsip_after3fraction, moss_qsip_after3fraction$sampleID), function(x) unique(x$otu_id))

df_unique_otu_ids <- 
  do.call(rbind, lapply(seq_along(unique_otu_ids), function(i) {
  data.frame(sampleID = names(unique_otu_ids)[i], otu_id = unique_otu_ids[[i]])
}))

# Count the number of OTUs per sampleID 
unique_sum <-df_unique_otu_ids %>% group_by(sampleID) %>% summarise(n =length(otu_id))
unique_sum

# Using the unique OTU IDs per treatment, identify OTUs that occur at least twice across replicates.
# Once these OTUs are identified, compare them to the full dataset to extract the final filtered information.
# A taxon had to be present in at least three qSIP fractions of two replicates 
# for both natural abundance and 15N-enriched treatments per control or warming treatment 

# Place the experimental design file (idf) in a new folder beforehand.
# Read the experimental design information so that the confirmed unique OTU ID list 
# can be matched with additional metadata for downstream analysis (using merge function).
df_unique_otu_ids_rep <- merge(df_unique_otu_ids, group_sampleID_with_control_and_rep, by="sampleID", all.x = T, sort = F)
# View(df_unique_otu_ids_rep)

# Count occurrences of each OTU ID within moss  ecosystems under control and warming treatments
otu_id_counts <- with(df_unique_otu_ids_rep, table(Treatment, category, iso_trt,  otu_id))

# Convert the result to a data frame and rename columns
df_otu_id_counts <- as.data.frame(otu_id_counts)
names(df_otu_id_counts) <- c("Treatment", "ecosystem", "iso_trt","otu_id", "count")
df_otu_id_counts$count <- as.numeric(df_otu_id_counts$count)

# Filter OTUs that appear at least twice in replicates
# The 'count' threshold can be adjusted for exploration or different filtering criteria
num_rep = 2
moss_qsip_filtered_final <- df_otu_id_counts %>% filter(count >= num_rep)
moss_qsip_filtered_final

# Count the number of OTUs per Treatment and ecosystem/category 
# to compare with manual Excel records and verify the number of OTUs after filtering for two replicates
unique_sum_2_rep <-moss_qsip_filtered_final %>% group_by(Treatment, ecosystem ,iso_trt) %>% summarise(n =length(otu_id))
unique_sum_2_rep

# After identifying OTUs present in two replicates, further identify OTUs 
# that appear in both light and label treatments(n>1 represent light and label)
moss_control <- moss_qsip_filtered_final %>% filter(Treatment=="Control" & ecosystem=="moss") %>% group_by(otu_id) %>% filter(n() > 1) %>% distinct(otu_id)
moss_warming <- moss_qsip_filtered_final %>% filter(Treatment=="Warming" & ecosystem=="moss") %>% group_by(otu_id) %>% filter(n() > 1) %>% distinct(otu_id)

# Extract the corresponding data from the three-fraction dataset 
# based on the confirmed OTU names (otu_id)
moss_control_sip <- moss_qsip_after3fraction %>% filter(treatment =="Control"& ecosystem=="moss" & otu_id %in% moss_control$otu_id) 
moss_warming_sip <- moss_qsip_after3fraction %>% filter(treatment =="Warming"& ecosystem=="moss" & otu_id %in% moss_warming$otu_id) 

# Check the output dimensions and compare with the original Excel records 
# to verify that filtering was performed correctly
dim(moss_control_sip)
dim(moss_warming_sip)

# Combine the two separate datasets into a single table
# after removing OTUs that did not appear in at least two replicates
moss_qsip_filter_data <- rbind(moss_control_sip, moss_warming_sip)
dim(moss_qsip_filter_data)
cat('After fraction and replication filtering\n')
seq_summary(moss_qsip_filter_data[, ..ssc])
# After fraction and replication filtering
# 10,765,859 sequences,
# 1,115 tax features

######################################################################################################
# 4. Calculate the EAF with and without the 1,000 bootstrap resampling iterations
moss_qsip = moss_qsip_filter_data[,1:19] 
# write.xlsx(moss_qsip,"moss_qsip.xlsx")

#Convert the data to data.table format
moss_qsip <- data.table(moss_qsip)

#convert the format
moss_qsip$treatment <- as.factor(moss_qsip$treatment)
moss_qsip$avg_16S_g_soil <- as.numeric(moss_qsip$avg_16S_g_soil)
moss_qsip$iso_trt <- factor(moss_qsip$iso_trt, levels = c("light", "label"))
moss_qsip$Density.g.ml <- as.numeric(moss_qsip$Density.g.ml )

# keep track of which columns should be kept for downstream analyses
keep_cols <- setdiff(names(moss_qsip), c('otu_id', 'sampleID', 'fraction', 'Density.g.ml','avg_16S_g_soil', 'rel_abund', 'seq_abund'))

# calculate weighted average densities
wads <- calc_wad(moss_qsip,
                 tax_id = 'otu_id', 
                 sample_id = 'sampleID', 
                 frac_id = 'fraction',
                 frac_dens = 'Density.g.ml', 
                 frac_abund = 'avg_16S_g_soil',
                 rel_abund = 'rel_abund',
                 grouping_cols = keep_cols)
wads
# calculate weighted average densities, transform to wide format
# It is particularly important to ensure that the light and label samples are ordered correctly; 
# otherwise, the EAF calculation will not work properly.
wads$treatment <- as.factor(wads$treatment)
wads$iso_trt <- factor(wads$iso_trt, levels = c("light", "label"))
str(wads)

ww <- wad_wide(wads, 
               tax_id = 'otu_id', 
               sample_id = 'sampleID',
               iso_trt = 'iso_trt', 
               isotope = 'isotope')

ww

# Make a copy of the dataset for EAF calculation
# This allows comparison between manual and automated calculations
eaf_taxa_cali <- ww

# The calculation formula is as follows:
eaf_taxa_cali[, gc_prop := (1 / 0.083506) * (light - 1.646057)][, mw_light := (0.496 * gc_prop) + 307.691][, mw_label := (((label - light) / light) + 1) * mw_light]
nat_abund_15N = 0.003663004
eaf_taxa_cali[isotope == '15N', `:=` (mw_max = mw_light + 3.517396 + (0.5024851 * gc_prop), nat_abund = nat_abund_15N)][, eaf := ((mw_label - mw_light) / (mw_max - mw_light)) * (1 - nat_abund)]
eaf_taxa_cali

# convert the format
wads$iso_trt <- factor(wads$iso_trt, levels = c("light", "label"))

# Set random seed; This ensures that any random operations produce the same result every time the code is run.
set.seed(2661)

# Taxa-specific 15N EAF was calculated at the treatment level by using 1,000 bootstrap resampling iterations
eaf_pro_1000 <- calc_excess(
  data = wads,
  wads = 'wad',
  tax_id = 'otu_id', 
  sample_id = 'sampleID', 
  iso_trt = 'iso_trt', 
  isotope = 'isotope',
  rm_outliers = F,
  correction = F,
  bootstrap = T,
  iters = 1000,
  nat_abund_15N = 0.003663004,
  grouping_cols = c('treatment'))

eaf_pro_1000

# output the  1000_bootstrap results
write.csv(eaf_pro_1000, file = "eaf_taxa_with_1000_bootstrap.csv")

# Taxa-specific 15N EAF was calculated at the replicate level without bootstrap resampling.
eaf_taxa <- calc_excess(
                    data = wads,
                    wads = 'wad',
                    tax_id = 'otu_id',
                    sample_id = 'sampleID',
                    iso_trt = 'iso_trt',
                    isotope = 'isotope',
                    correction = F,
                    rm_outliers = F,
                    nat_abund_15N = 0.003663004)

eaf_taxa
write.csv(eaf_taxa, file = "eaf_taxa_without_bootstrap.csv")

######################################################################################################
# 5. Find  both with and without 1,000 bootstrap resampling EAF value all over 0
# Create an index tag with otu_id + treatment + ecosystem,
# match records where eaf_taxa > 0 and eaf_pro > 0,
# and filter OTUs that are shared and significant in both tables
# load the data without bootstrap resampling
data_eaf_taxa <- eaf_taxa

# load the data with 1,000 bootstrap resampling iterations
data_eaf_pro <- eaf_pro_1000

# creat an index tag with otu_id + treatment + ecosystem
data_eaf_taxa[, index_column := paste(otu_id,  treatment, sep = "_")]
data_eaf_pro[, index_column := paste(otu_id, treatment, sep = "_")]

#rename the data_eaf_pro
names(data_eaf_pro) <- c("otu_id","treatment","q2.5","q50","q97.5","p_val","index_column")

# Filter data in eaf_taxa with eaf > 0, 
# and in pro_data with values > 0 and p < 0.05
filtered_eaf_taxa <- data_eaf_taxa[eaf >0]
filtered_pro_data <- data_eaf_pro[q50>0 & p_val < 0.05]

# Match by index and retain all contents from filtered_eaf_taxa 
# that have matching "index_column" entries in filtered_pro_data
filtered_data <- merge(filtered_eaf_taxa, filtered_pro_data, by = "index_column", all = FALSE)

dim(filtered_data)
filtered_data <- filtered_data[,-1]

names(filtered_data) <-  c("otu_id", "sampleID", "timepoint", "isotope", "ecosystem", 
                  "treatment", "rep", "Kingdom", "Phylum", "Class", "Order", 
                  "Family", "Genus", "wvd", "abund", "light", "label", 
                  "gc_prop", "mw_light", "mw_label", "mw_max", "nat_abund", "eaf",
                  "otu_id.y",	"treatment.y",	"q2.5",	"q50",	"q97.5",	"p_val")
# View(filtered_data)
write.xlsx(filtered_data, file = "filtered_data_taxa_signi_over_0.xlsx", rowNames = FALSE)

## Extract sampleID and eaf for creating wide_OTU_moss_result_after_eaf_over_0.csv
select_OTU_venn <- filtered_data[,c("otu_id","sampleID","eaf")]

# The purpose of these two lines is to organize the EAF data into a matrix 
# with each row representing an OTU and each column representing a sample, 
# and to fill missing values with 0 for easier downstream calculations or CSV export.
df_wide <- select_OTU_venn %>%
  pivot_wider(names_from = sampleID, values_from = eaf) %>%
  select(otu_id, W1_Moss_21, W1_Moss_22, W1_Moss_23,
         W1_Moss_24, W1_Moss_25, W1_Moss_26, 
         W1_Moss_27, W1_Moss_28, W1_Moss_29,
         W1_Moss_30, W1_Moss_31, W1_Moss_32)
df_wide[is.na(df_wide)] <- 0

#output the data file
write.csv(df_wide, "wide_OTU_moss_result_after_eaf_over_0.csv")

# Record the end time of the script
end_time_all <- Sys.time()

#total time
total_time <-  end_time_all - start_time_all
total_time

