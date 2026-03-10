# Ksenia Klar
# March 2nd 2026
# This script is for extracting csv and writing csv files to be used in the IMG interface
# MUST BE USED WITH Joins_mer_operon.R code

# Extracting genome IDs ---------------------------------------------------

# creating a new dataframe
merIMGid <- mer_R["IMG.Genome"]
colnames(merIMGid) <- "taxon_oid"  # new column name


# Export as tab-delimited file
write.table(merIMGid, 
            "Output/Genome ID_operon.tsv", 
            sep = "\t",       # tab delimiter
            row.names = FALSE, 
            quote = FALSE)    # optional: prevents quotes around text



# Extracting the IMG IDs in two batches because IMG only allows 1000 genomes per blast --------

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


# Converting the txt files to csv -----------------------------------------

# This is to extract TSV files from each column of gene ids into spearate files

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


# Extracting merP and merT IDs to tsv for IMG -----------------------------

# creating a function that will extract each of the img id columns in a separate file

export_mer_IDs_tsv <- function(df, gene = "merP", cols = 1:6, prefix = "Gene_ID") {
  
  # Construct column names (ex: merP.ID.1, merP.ID.2, ...)
  mer_cols <- paste0(gene, ".ID.", cols)
  
  # Loop over columns
  lapply(seq_along(mer_cols), function(i) {
    
    col_name <- mer_cols[i]
    
    # Extract column
    tmp_df <- df[col_name]
    colnames(tmp_df) <- "Gene ID"
    
    # Remove NA rows
    tmp_df <- tmp_df[!is.na(tmp_df$`Gene ID`), , drop = FALSE]
    
    # File name
    tsv_file <- paste0(prefix, "_", gene, "_", cols[i], ".tsv")
    
    # Write file
    write.table(tmp_df, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  })
}

export_mer_IDs_tsv(mer_R, gene = "merP", cols = 1:6, prefix = "Gene_ID")
export_mer_IDs_tsv(mer_R, gene = "merT", cols = 1:6, prefix = "Gene_ID")


# Processing the FASTA files ----------------------------------------------
# The one created above
reheader_mer_fasta <- function(folder) {
  
  # Find all fasta files
  fasta_files <- list.files(folder, pattern = "\\.faa$", full.names = TRUE)
  
  lapply(fasta_files, function(file) {
    
    # Extract column number from filename (e.g., col_3)
    col_num <- sub(".*col_([0-9]+).*", "\\1", basename(file))
    
    seqs <- Biostrings::readAAStringSet(file)
    
    df <- data.frame(
      id = sub("^([^ ]+).*", "\\1", names(seqs)),
      description = sub("^[^ ]+ ?", "", names(seqs)),
      sequence = as.character(seqs),
      stringsAsFactors = FALSE
    )
    
    headers <- paste0("(", col_num, ")", df$id, " ", df$description)
    
    out_seqs <- Biostrings::AAStringSet(df$sequence)
    names(out_seqs) <- headers
    
    # Overwrite the original file
    Biostrings::writeXStringSet(out_seqs, file)
  })
}

reheader_mer_fasta("Data/merP FASTA")
reheader_mer_fasta("Data/merT FASTA")


# Extracting TSV based on coupled/orphaned --------------------------------

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
