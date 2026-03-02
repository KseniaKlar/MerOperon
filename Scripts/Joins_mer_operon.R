# Ksenia Klar
# February 17th 2026
# In this script I will remove merA only aorganisms from supplementaey data in Christakis et. al.
# Expanded Diversity and Phylogeny of mer Genes Broadens Mercury Resistance Paradigms and Reveals an Origin for MerA Among Thermophilic Archaea

# clearing the environment
rm(list=ls())

# loading the libraries
library(tidyverse)
library(dplyr)
library(dslabs)
library(Biostrings)
library(ggplot2)
library(readr)

# The dataset
mer_R <- read.csv("Data/sup_data_Christakis.csv")
View(mer_R)

# Separating by merA and merB ---------------------------------------------

mer_R$has_merA <- mer_R$merA_copies > 0
mer_R$has_merB <- (mer_R$`merB.copies` > 0) | (mer_R$`merB.like.copies` > 0)

# Keeping only the ones that have merB
mer_R <- mer_R %>%
  dplyr::filter(has_merB == TRUE)

# Removing all the genomes that only have merB copies
mer_R <- mer_R %>%
  filter(`merB.copies` > 0)

# Removing all merB like accession numbers
clean_merB_IMG <- function(df) {
  
  merB_cols <- paste0("merB.IMG.gene.ID.", 1:4)
  
  df %>%
    mutate(across(
      all_of(merB_cols),
      ~ ifelse(grepl("^[0-9]+$", .x), .x, NA)
    ))
}

mer_R <- clean_merB_IMG(mer_R)



# Extracting genome ids ---------------------------------------------------

# creating a new dataframe
merIMGid <- mer_R["IMG.Genome"]
colnames(merIMGid) <- "taxon_oid"  # new column name


# Export as tab-delimited file
write.table(merIMGid, 
            "Output/Genome ID_operon.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)    # optional: prevents quotes around text



# Extracting the IMG IDs in two batches because IMG only allows bl --------

# Extract IMG Genome column
merIMGid <- mer_R["IMG.Genome"]
colnames(merIMGid) <- "taxon_oid"

# Make sure Output folder exists
if (!dir.exists("Output")) {
  dir.create("Output")
}

# First 1000 genomes
merIMGid_1 <- merIMGid[1:1000, ]

# Remaining genomes
merIMGid_2 <- merIMGid[1001:nrow(merIMGid), ]


# Extracting in batches of 100 --------------------------------------------


# Write first batch
# Extract IMG Genome column
# merIMGid <- mer_R["IMG.Genome"]
# colnames(merIMGid) <- "taxon_oid"
# 
# # Set batch size
# batch_size <- 100
# n <- nrow(merIMGid)
# num_batches <- ceiling(n / batch_size)
# 
# # Loop over batches
# for (i in 1:num_batches) {
#   
#   # Calculate row indices for this batch
#   start_row <- (i - 1) * batch_size + 1
#   end_row <- min(i * batch_size, n)
#   
#   # Subset the batch
#   batch <- merIMGid[start_row:end_row, ]
#   
#   # File name
#   file_name <- paste0("Output/Genome_ID_batch_", i, ".tsv")
#   
#   # Write tab-delimited file
#   write.table(batch,
#               file_name,
#               sep = "\t",
#               row.names = FALSE,
#               quote = FALSE)
# }



# Combining the metabolism database with the original mer operon databas --------

mer_metabolism <- read.csv("Data/mer_metabolisms.csv")

# Joining the two dataframes
mer_R <- mer_R%>%
  left_join(mer_metabolism, 
            by = c("IMG.Genome" = "IMG.Genome.ID"))


# Adding the merR BLAST files ---------------------------------------------

# Importing the data
batch1 <- read.csv("Data/MerR/batch1.csv")
batch2 <- read.csv("Data/MerR/batch2.csv")

# combining the two dtaframes
merR <- rbind(batch1, batch2)


merR_numbered <- merR %>%
  group_by(`Genome.ID`) %>%
  mutate(hit_number = row_number()) %>%
  ungroup()

merR_wide <- merR_numbered %>%
  pivot_wider(
    names_from = hit_number,
    values_from = `Gene.ID`,
    names_prefix = "merR.ID "
  )

mer_R <- mer_R %>%
  left_join(merR_wide, by = c("IMG.Genome" = "Genome.ID"))


# Importing and cleaning up the merP and merT files -----------------------

meroperon <- list.files("Data/merO", full.names = TRUE)

merO_combined <- lapply(1:length(meroperon), function(x){
  x
  mer <- read.csv(meroperon[x], stringsAsFactors=FALSE)
  
  
  #Filepath<-paste()
  write.csv(finaldataframe, "/Gene and genome IDs for IMG screening/Outputfile.csv", row.names=FALSE)
})


# Converting the txt files to csv -----------------------------------------

# 1. Get all .txt files in your folder
txt_files <- list.files("Data/Raw data", pattern = "\\.txt$", full.names = TRUE)

# 2. Use lapply to convert each one
lapply(txt_files, function(file_path) {
  
  # Read the text file (adjust 'sep' if it's not tabs)
  # sep = "\t" for tabs, sep = "" for any whitespace
  temp_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Create the new filename by replacing .txt with .csv
  new_name <- gsub("\\.txt$", ".csv", file_path)
  
  # Save as CSV
  write.csv(temp_data, file = new_name, row.names = FALSE)
})


# Trying some things ------------------------------------------------------

# MerP

batch1_p501 <- read.csv("Data/merP/batch1_merP_Tn501.csv")
batch2_p501 <- read.csv("Data/merP/batch2_merP_Tn501.csv")

batch1_p21 <- read.csv("Data/merP/batch1_merP_Tn21.csv")
batch2_p21 <- read.csv("Data/merP/batch2_merP_Tn21.csv")

merP501 <- rbind(batch1_p501, batch2_p501)
merP21 <- rbind(batch1_p21, batch2_p21)

# MerT

batch1_t501 <- read.csv("Data/merT/batch1_merT_Tn501.csv")
batch2_t501 <- read.csv("Data/merT/batch2_merT_Tn501.csv")

batch1_t21 <- read.csv("Data/merT/batch1_merT_Tn21.csv")
batch2_t21 <- read.csv("Data/merT/batch2_merT_Tn21.csv")

merT501 <- rbind(batch1_t501, batch2_t501)
merT21 <- rbind(batch1_t21, batch2_t21)



# Cheching for non repeated values

# merP
bind_rows(merP21["Genome.ID"], merP501["Genome.ID"]) %>%
  count(Genome.ID) %>%
  filter(n == 1) %>%
  pull(Genome.ID)

# merT
bind_rows(merT21["Genome.ID"], merT501["Genome.ID"]) %>%
  count(Genome.ID) %>%
  filter(n == 1) %>%
  pull(Genome.ID)









merP <- merP %>%
  group_by(`Genome.ID`) %>%
  mutate(hit_number = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = hit_number,
    values_from = `Gene.ID`,
    names_prefix = "merP.ID "
  )

mer_R <- mer_R %>%
  left_join(merP, by = c("IMG.Genome" = "Genome.ID"))

