# Ksenia Klar
# March 2nd 2026
# This script is for extracting csv and writing csv files to be used in the IMG interface
# MUST BE USED WITH Joins_mer_operon.R code

### DO NOT CLEAR THE ENVIRONMENT OR RESTART R###

# Extracting genome IDs ---------------------------------------------------

# creating a new dataframe
merIMGid <- mer_R["IMG.Genome"]
colnames(merIMGid) <- "taxon_oid"  # new column name


# Export as tab-delimited file
write.table(merIMGid, 
            "Output/Genome ID_operon.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)



# Extracting the IMG IDs in two batches because IMG only allows 1000 genomes per blast --------

# Extract IMG Genome column
merIMGid <- mer_R["IMG.Genome"]
colnames(merIMGid) <- "taxon_oid"

# First 1000 genomes
merIMGid_1 <- merIMGid[1:1000, ]

# Remaining genomes
merIMGid_2 <- merIMGid[1001:nrow(merIMGid), ]

write.table(merIMGid_1, 
            "Output/Genome ID_operon_batch1.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)

write.table(merIMGid_2, 
            "Output/Genome ID_operon_batch2.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)



# Converting the txt files to csv -----------------------------------------

# This is to extract TSV files from each column of gene ids into spearate files

# 1. Get all .txt files in the folder
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


# Extracting merP and merT IDs to tsv for IMG -----------------------------

mer_IDs_B_AB <- function(df, gene = "merP", cols = 1:6) {
  
  # Construct gene column names
  mer_cols <- paste0(gene, ".ID.", cols)
  
  # Subset rows
  df_B_only <- df[df$has_merB == TRUE & df$has_merA == FALSE, mer_cols, drop = FALSE]
  df_BA     <- df[df$has_merB == TRUE & df$has_merA == TRUE,  mer_cols, drop = FALSE]
  
  # Flatten and remove NA
  ids_B_only <- unique(na.omit(unlist(df_B_only)))
  ids_BA     <- unique(na.omit(unlist(df_BA)))
  
  # Create output data frames with EXACT column name
  df_out_B_only <- data.frame("Gene ID" = ids_B_only, check.names = FALSE)
  df_out_BA     <- data.frame("Gene ID" = ids_BA, check.names = FALSE)
  
  # Write TSV files
  write.table(df_out_B_only,
              paste0(gene, "_orphaned.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(df_out_BA,
              paste0(gene, "_coupled.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

mer_IDs_B_AB(mer_R, gene = "merP")
mer_IDs_B_AB(mer_R, gene = "merT")


# Renaming the FASTA sequences --------------------------------------------

reheader_merP_T_fasta <- function(folder = "FASTA") {
  
  fasta_files <- list.files(folder, pattern = "\\.faa$", full.names = TRUE)
  
  lapply(fasta_files, function(file) {
    
    seqs <- Biostrings::readAAStringSet(file)
    
    df <- data.frame(
      id = sub("^([^ ]+).*", "\\1", names(seqs)),
      description = sub("^[^ ]+ ?", "", names(seqs)),
      sequence = as.character(seqs),
      stringsAsFactors = FALSE
    )
    
    # Determine prefix based on file name
    if (grepl("orphaned|MerB", file, ignore.case = TRUE)) {
      prefix <- "or"      # orphaned
    } else if (grepl("coupled|MerAB", file, ignore.case = TRUE)) {
      prefix <- "co"      # coupled
    } else {
      prefix <- ""        # fallback
    }
    
    headers <- paste0(prefix, df$id)  # only ID, with prefix
    
    out_seqs <- Biostrings::AAStringSet(df$sequence)
    names(out_seqs) <- headers
    
    # overwrite original file
    Biostrings::writeXStringSet(out_seqs, file)
    
  })
}

reheader_merP_T_fasta("Data/FASTA")


# Extracting merA IMG IDs -------------------------------------------------

# Orphaned, anaerobic merAs

orph_anae_merA <- mer_R %>% 
  filter(Oxygen.Requirement == "Anaerobe" & has_merA == T & has_merB == F) %>% 
  select(merA.IMG.gene.ID.1)

colnames(orph_anae_merA) <- "gene_oid"  # new column name

orph_anae_merA <- orph_anae_merA %>%
  filter(grepl("^[0-9]+$", gene_oid))

# Export as tab-delimited file
write.table(orph_anae_merA, 
            "Output/Gene ID_orphaned_anaerobe_merA.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)





# Orphaned other merA IMG IDs

orph_merA <- mer_R %>% 
  filter(Oxygen.Requirement != "Anaerobe" & has_merA == T & has_merB == F) %>% 
  select(merA.IMG.gene.ID.1)

colnames(orph_merA) <- "gene_oid"  # new column name

orph_merA <- orph_merA %>%
  filter(grepl("^[0-9]+$", gene_oid))

# Export as tab-delimited file
write.table(orph_merA, 
            "Output/Gene ID_orphaned_merA.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)




# Coupled merA IMG IDs

coup_merA <- mer_R %>% 
  filter(has_merA == T & has_merB == T) %>% 
  select(merA.IMG.gene.ID.1)

colnames(coup_merA) <- "gene_oid"  # new column name

coup_merA <- coup_merA %>%
  filter(grepl("^[0-9]+$", gene_oid))

# Export as tab-delimited file
write.table(coup_merA, 
            "Output/Gene ID_coupled_merA.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)


mer_R %>% 
  filter(merA.IMG.gene.ID.1 == 650971024) %>% 
  select(S)

