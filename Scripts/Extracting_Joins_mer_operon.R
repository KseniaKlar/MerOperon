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

