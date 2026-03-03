# Ksenia Klar
# February 17th 2026
# In this script I will add information on the mer operon to the merB genomes from supplementary data in Christakis et. al.
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


# Combining the metabolism database with the original mer operon database --------

mer_metabolism <- read.csv("Data/mer_metabolisms.csv")

# Joining the two dataframes
mer_R <- mer_R%>%
  left_join(mer_metabolism, 
            by = c("IMG.Genome" = "IMG.Genome.ID"))


# Importing merP and merT files -------------------------------------------

# Structuring all file paths into one nested list
all_files_list <- list(
  merP = list(
    Tn501 = c("Data/merP/batch1_merP_Tn501.csv", "Data/merP/batch2_merP_Tn501.csv"),
    Tn21  = c("Data/merP/batch1_merP_Tn21.csv", "Data/merP/batch2_merP_Tn21.csv")
  ),
  merT = list(
    Tn501 = c("Data/merT/batch1_merT_Tn501.csv", "Data/merT/batch2_merT_Tn501.csv"),
    Tn21  = c("Data/merT/batch1_merT_Tn21.csv", "Data/merT/batch2_merT_Tn21.csv")
  )
)



# Lapply to clean merP and merT -------------------------------------------

processed_genes <- lapply(names(all_files_list), function(gene_name) {
  
# 1. Reading and binding batches from BLAST
  raw_data <- lapply(all_files_list[[gene_name]], function(files) {
    bind_rows(lapply(files, read.csv))
  })
  
# 2. Cleaning based on gene type: keeping the correctly matched ones
  if(gene_name == "merP") {
    cleaned <- lapply(raw_data, function(df) {
      filter(df, Gene.Name == "periplasmic mercuric ion binding protein") %>%
        select(Gene.Name, Genome.ID, Genome.Name, Gene.ID)
    })
  } else {
    cleaned <- lapply(raw_data, function(df) {
      filter(df, Gene.Name != "MerC mercury resistance protein") %>%
        select(Gene.Name, Genome.ID, Genome.Name, Gene.ID)
    })
  }
  
# 3. Combining the two queries from BLAST (Tn21/Tn501) keeping only unique gene IDs
  combined <- bind_rows(cleaned) %>% 
    distinct(Gene.ID, .keep_all = TRUE)
  
# 4. Pivoting wider to prepare to combine with main dataframe
  pivoted <- combined %>%
    group_by(Genome.ID) %>%
    mutate(hit_number = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      names_from = hit_number,
      values_from = Gene.ID,
      names_prefix = paste0(gene_name, ".ID.")
    ) %>%
    select(-Genome.Name)
  
  return(pivoted)
})


# Joining, adding descriptive columns and renaming ------------------------

# Assigning names so the dataframes are easily joined
names(processed_genes) <- names(all_files_list)

# Joining merP and merT with mer_R
for(gene in names(processed_genes)) {
  mer_R <- left_join(mer_R, processed_genes[[gene]], by = c("IMG.Genome" = "Genome.ID"))
}

# Adding columns with the number of copies of each gene and a TRUE/FALSE column for visualization later
mer_R <- mer_R %>%
  mutate(
    merP.copies = rowSums(!is.na(select(., starts_with("merP.ID")))),
    merT.copies = rowSums(!is.na(select(., starts_with("merT.ID")))),
    has_merP = merP.copies > 0,
    has_merT = merT.copies > 0
  ) %>% 
  mutate(across(where(is.character), ~na_if(., ""))) # Replacing all blank cells with NA for tidyness

# Renaming the gene name columns for more clarity
names(mer_R)[names(mer_R) == "Gene.Name.x"] <- "merP.comment" 
names(mer_R)[names(mer_R) == "Gene.Name.y"] <- "merT.comment"

# IDK why this doesnt work???
# %>%
#   rename(merP.comment = Gene.Name.x, 
#          merT.comment = Gene.Name.y)



# Visualizing -------------------------------------------------------------

mer_operoncontent <- mer_R %>% select(c("Domain", "IMG.Genome", "merB.copies", "merB.comment", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genome.Name...Sample.Name.x", "has_merA", "has_merB", "Metabolism", "Oxygen.Requirement", "Genome.Size....assembled", "merP.comment", "merT.comment", "merP.copies", "merT.copies", "has_merP", "has_merT"))

mer_operoncontent <- mer_operoncontent %>%
  mutate(
    quadrant = case_when(
      has_merP & has_merT  ~ "Both",
      has_merP & !has_merT ~ "merP only",
      !has_merP & has_merT ~ "merT only",
      TRUE                 ~ "Neither"
    ),
    category = case_when(
      has_merA & has_merB ~ "merAB",
      has_merB & !has_merA ~ "merB only",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category))   # keep only relevant genomes

mer_operoncontent %>% 
  ggplot(aes(x = has_merP, 
             y = has_merT, 
             color = category)) +
  geom_jitter(width = 0.5, 
              height = 0.4, 
              size = 1.5, 
              alpha = 0.6) +
  scale_x_discrete(labels = c("FALSE" = "No merP", "TRUE" = "merP")) +
  scale_y_discrete(labels = c("FALSE" = "No merT", "TRUE" = "merT")) +
  scale_color_manual(values = c("merAB" = "chocolate",
                                "merB only" = "forestgreen")) +
  theme_classic() +
  labs(x = "merP presence",
       y = "merT presence",
       color = "Genome Type",
       title = "merT and merP genes across orphaned and coupled merB genomes")


# Exploring the graph -----------------------------------------------------

# Finding out what genomes are merB, merT and merP
mer_R %>%
  filter(
    has_merB,
    has_merT,
    has_merP,
    !has_merA
  ) %>% select(IMG.Genome)

mer_operoncontent %>%
  filter(
    has_merB,
    has_merT,
    has_merP,
    !has_merA
  )

mer_R %>%
  filter(IMG.Genome %in% c(
    2791354983,
    2737471843,
    2547132255,
    2630968586,
    2728369713,
    2551306727,
    2585427827,
    2728369721
  )) %>%
  select(
    merB.IMG.gene.ID.1,
    merB.IMG.gene.ID.2,
    merB.IMG.gene.ID.3,
    merB.IMG.gene.ID.4,
    merB.IMG.gene.ID.5,
    merB.IMG.gene.ID.6,
    merB.IMG.gene.ID.7,
    merB.IMG.gene.ID.8,
    merB.IMG.gene.ID.9.12
  )
