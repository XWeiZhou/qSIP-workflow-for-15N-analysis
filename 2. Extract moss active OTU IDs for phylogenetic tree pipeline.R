# load the packages
library(dplyr)
library(clipr)
library(tidyr)
library(readxl)

#set work directory
setwd("C:/Users/Lenovo/Desktop/code/")

# Read the filtered Excel file instead of CSV
data <- read_excel("filtered_data_taxa_pro_signi_minus_0_20250330_1000.xlsx", sheet = 1)
data$eaf <- as.numeric(data$eaf)

# Convert the Treatment column to a factor variable
data$treatment <- factor(data$treatment, levels = c("Control","Warming"))

# Calculate the mean of EAF
result <- data %>%
  filter(ecosystem == "moss" & eaf > 0) %>%
  group_by(asv_id, treatment) %>%
  summarise(n =length(eaf), mean_value = mean(eaf))

# Use the pivot_wider function to convert long-format data (df) to wide-format
result <- result[,-3]
wide_result <- pivot_wider(data = result, names_from = treatment, values_from = mean_value)
wide_result
str(wide_result)

# Replace all NA values with 0
wide_result[is.na(wide_result)] <- 0

# Add a column for the EAF difference
wide_result$delta.eaf <- wide_result$Warming - wide_result$Control 
  
# load all taxa dataset
taxa <- read.csv("taxa.csv")
names(taxa) <- c("asv_id","Kingdom","Phylum","Class","Order","Family","Genus", "Species")

# Compare the two datasets and retain the table structure for OTUs that have EAF values in both
merge_result <- merge(wide_result, taxa, by = "asv_id", all.x = TRUE)
annotation <- data.frame(merge_result)
annotation

# Since Warming appears before Control, the second column is Warming instead of Control; 
# check the data carefully when naming columns
names(annotation) <- c("OTUID","Warming","Control","delta.eaf","Kingdom","Phylum","Class","Order","Family","Genus", "Species")
write.csv(annotation, "annotation_moss.csv")

# Output selected OTU names for later comparison with otus_taxa.fa, 
# to extract OTU sequences for tree construction
select_OTU_names <- as.data.frame(annotation[,1])
names(select_OTU_names) <- c("asv_id")
write.table(select_OTU_names, "select_OTU_moss_509.txt", sep = "/t", row.names = F, col.names = T, quote = F)

# Check the dimensions (rows and columns) of the data
data.moss <- read.csv("annotation_moss.csv")
dim(data.moss)

