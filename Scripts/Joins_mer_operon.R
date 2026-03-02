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


# Combining the metabolism database with the original mer operon database --------

mer_metabolism <- read.csv("Data/mer_metabolisms.csv")

# Joining the two dataframes
mer_R <- mer_R%>%
  left_join(mer_metabolism, 
            by = c("IMG.Genome" = "IMG.Genome.ID"))


# Importing merP and merT files -------------------------------------------

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



# Cleaning up the dtaframes -----------------------------------------------


# Removing all non merPs
unique(merP21$Gene.Name)

merP21 <- merP21 %>% 
  filter(Gene.Name == "periplasmic mercuric ion binding protein")
merP501 <- merP501 %>% 
  filter(Gene.Name == "periplasmic mercuric ion binding protein")

# Removing all non merTs
unique(merT501$Gene.Name)
unique(merT21$Gene.Name)

merT21 <- merT21 %>% 
  filter(Gene.Name != "MerC mercury resistance protein")
merT501 <- merT501 %>% 
  filter(Gene.Name != "MerC mercury resistance protein")

# Removing unncessary columns
# MerP21
names(merP21)

merP21 <- merP21 %>% 
  select(c("Gene.Name", "Genome.ID", "Genome.Name", "Gene.ID"))

# MerP501
names(merP501)

merP501 <- merP501 %>% 
  select(c("Gene.Name", "Genome.ID", "Genome.Name", "Gene.ID"))

# MerT21
names(merT21)

merT21 <- merT21 %>% 
  select(c("Gene.Name", "Genome.ID", "Genome.Name", "Gene.ID"))

# MerT501
names(merT501)

merT501 <- merT501 %>% 
  select(c("Gene.Name", "Genome.ID", "Genome.Name", "Gene.ID"))




# Exploring the dataframes ------------------------------------------------


# Checking for non duplicated genomes
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

# Checking for non duplicated genes

# merP
bind_rows(merP21["Gene.ID"], merP501["Gene.ID"]) %>%
  count(Gene.ID) %>%
  filter(n == 1) %>%
  pull(Gene.ID)

# merT
bind_rows(merT21["Gene.ID"], merT501["Gene.ID"]) %>%
  count(Gene.ID) %>%
  filter(n == 1) %>%
  pull(Gene.ID)



# Combining the dataframes ------------------------------------------------


merP <- bind_rows(merP21, merP501) %>%
  distinct(Gene.ID, .keep_all = TRUE)
merT <- bind_rows(merT21, merT501) %>%
  distinct(Gene.ID, .keep_all = TRUE)


# Pivoting wider
# merP
merP <- merP %>%
  group_by(`Genome.ID`) %>%
  mutate(hit_number = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = hit_number,
    values_from = `Gene.ID`,
    names_prefix = "merP.ID."
  )

# merT
merT <- merT %>%
  group_by(`Genome.ID`) %>%
  mutate(hit_number = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = hit_number,
    values_from = `Gene.ID`,
    names_prefix = "merT.ID."
  )


# Combining with the master mer dataframe ---------------------------------



# Importing and cleaning up the merP and merT files -----------------------

meroperon <- list.files("Data/merO", full.names = TRUE)

merO_combined <- lapply(1:length(meroperon), function(x){
  x
  mer <- read.csv(meroperon[x], stringsAsFactors=FALSE)
  
  
  #Filepath<-paste()
  write.csv(finaldataframe, "/Gene and genome IDs for IMG screening/Outputfile.csv", row.names=FALSE)
})

# Adding the merR BLAST files ---------------------------------------------

mer_R <- mer_R %>%
  left_join(merR_wide, by = c("IMG.Genome" = "Genome.ID"))
