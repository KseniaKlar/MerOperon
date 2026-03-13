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
library(gridExtra)
library(tibble)
library(msa)
library(pwalign)
library(magrittr)


# The dataset
mer_R <- read.csv("Data/sup_data_Christakis.csv")
View(mer_R)

# Separating by merA and merB ---------------------------------------------

mer_R$has_merA <- mer_R$merA_copies > 0
mer_R$has_merB <- (mer_R$`merB.copies` > 0) | (mer_R$`merB.like.copies` > 0)


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

mer_metabolism <- read.csv("Data/mer_metabolisms_all.csv")

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



# Joining -----------------------------------------------------------------

# Assigning names so the dataframes are easily joined
names(processed_genes) <- names(all_files_list)

# Joining merP and merT with mer_R
for(gene in names(processed_genes)) {
  mer_R <- left_join(mer_R, processed_genes[[gene]], by = c("IMG.Genome" = "Genome.ID"))
}



# Tidying -----------------------------------------------------------------
names(mer_R)

# Adding columns with the number of copies of each gene and a TRUE/FALSE column for visualization later
mer_R <- mer_R %>%
  mutate(
    merP.copies = rowSums(!is.na(select(., starts_with("merP.ID")))),
    merT.copies = rowSums(!is.na(select(., starts_with("merT.ID")))),
    has_merP = merP.copies > 0,
    has_merT = merT.copies > 0
  )
# Renaming the gene name columns for more clarity
names(mer_R)[names(mer_R) == "Gene.Name.x"] <- "merP.comment" 
names(mer_R)[names(mer_R) == "Gene.Name.y"] <- "merT.comment"

mer_R <- mer_R %>% 
  mutate(across(where(is.character), ~na_if(., ""))) %>% # Replacing all blank cells with NA for tidyness
  select(-c("merA.comment", "merA.IMG.gene.ID.1", "merA.IMG.gene.ID.2", "merA.IMG.gene.ID.3", "merA.IMG.gene.ID.4", "merA.IMG.gene.ID.5", "merA.IMG.gene.ID.6", "merA.IMG.gene.ID.7", "merA.IMG.gene.ID.8", "merB.like.copies", "Content.in.mer.genes", "merB.comment", "Genome.Name...Sample.Name.x", "merB.IMG.gene.ID.1", "merB.IMG.gene.ID.2", "merB.IMG.gene.ID.3", "merB.IMG.gene.ID.4", "merB.IMG.gene.ID.5", "merB.IMG.gene.ID.6", "merB.IMG.gene.ID.7", "merB.IMG.gene.ID.8", "merB.IMG.gene.ID.9.12")) %>% # Removing unncessary columns
  dplyr::rename(Sample.Name = `Genome.Name...Sample.Name.y`) %>% 
  mutate(Oxygen.Requirement = case_when(
    Oxygen.Requirement %in% c("Aerobe", "Obligate aerobe") ~ "Aerobe",
    Oxygen.Requirement == "Microaerophilic" ~ "Microaerophilic",
    Oxygen.Requirement %in% c("Facultative", "Facultative anaerobe") ~ "Facultative",
    Oxygen.Requirement %in% c("Anaerobe", "Obligate anaerobe") ~ "Anaerobe",
    TRUE ~ NA_character_ #Tidying the oxygen requirement column to have consistent data
  ))

# Running an alignment ----------------------------------------------------

# Read protein sequences
merB_seqs_orphaned <- readAAStringSet("Data/orphanedMerB.faa")
merB_seqs_coupled <- readAAStringSet("Data/coupledMerB.faa")
# Run MUSCLE alignment
merB_alignment_orphaned <- msa(merB_seqs_orphaned, method = "Muscle")
merB_alignment_coupled <- msa(merB_seqs_coupled, method = "Muscle")

# Getting the % identity

as(merB_alignment_orphaned, "AAStringSet") %>% 
  pwalign::stringDist(method = "hamming") %>% 
  as.matrix() %>% 
  # Use ncol() instead of width() for alignment objects
  {1 - (. / ncol(merB_alignment_orphaned)) } %>% 
  .[lower.tri(.)] %>% 
  mean() %>% 
  { . * 100 } %>% 
  round(2) %>% 
  { print(paste0("Total % Pairwise Identity: ", ., "%")) }


as(merB_alignment_coupled, "AAStringSet") %>% 
  pwalign::stringDist(method = "hamming") %>% 
  as.matrix() %>% 
  # Use ncol() instead of width() for alignment objects
  {1 - (. / ncol(merB_alignment_coupled)) } %>% 
  .[lower.tri(.)] %>% 
  mean() %>% 
  { . * 100 } %>% 
  round(2) %>% 
  { print(paste0("Total % Pairwise Identity: ", ., "%")) }



# Visualizing -------------------------------------------------------------

p1 <- mer_R %>%
  filter(!is.na(Oxygen.Requirement)) %>%
  mutate(
    GenomeType = case_when(
      has_merA & has_merB ~ "merA and merB",
      !has_merA & has_merB ~ "merB only",
      has_merA & !has_merB ~ "merA only",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(GenomeType)) %>%
  count(GenomeType, Oxygen.Requirement) %>%
  group_by(GenomeType) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "", y = prop, fill = Oxygen.Requirement)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~GenomeType, strip.position = "bottom") +
  labs(
    title = "Oxygen Requirements\nAmong merA and merB Genomes",
    fill = "Oxygen Requirement"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    strip.placement = "outside",
    strip.text = element_text(size = 12)
  )

p1


p2 <- mer_R %>% mutate(
  quadrant = case_when(
    has_merP & has_merT  ~ "Both",
    has_merP & !has_merT ~ "merP only",
    !has_merP & has_merT ~ "merT only",
    TRUE                 ~ "Neither"
  ),
  category = case_when(
    has_merA & has_merB ~ "coupled",
    has_merB & !has_merA ~ "orphaned",
    TRUE ~ NA_character_
  )
) %>%
  filter(!is.na(category)) %>% 
  ggplot(aes(x = has_merP, 
             y = has_merT, 
             color = category)) +
  geom_jitter(width = 0.5, 
              height = 0.4, 
              size = 1.5, 
              alpha = 0.6) +
  scale_x_discrete(labels = c("FALSE" = "No merP", "TRUE" = "merP")) +
  scale_y_discrete(labels = c("FALSE" = "No merT", "TRUE" = "merT")) +
  scale_color_manual(values = c("coupled" = "chocolate",
                                "orphaned" = "forestgreen")) +
  theme_classic() +
  labs(x = "merP presence",
       y = "merT presence",
       color = "merB",
       title = "merT and merP genes across\norphaned and coupled merB genomes") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5),
    legend.title = element_text(
      face = "bold",
      size = 14),
    legend.text = element_text(
      size = 12),
    axis.ticks = element_blank()
  ) 
p2

# Lines separating the quandrents
# +
#   geom_vline(xintercept = 1.5, color = "gray40") +
#   geom_hline(yintercept = 1.5, color = "gray40")

ggsave("Output/MerT and merP presence.png", plot = p2, bg = "white")


p3 <- mer_R %>%
  mutate(
    GenomeType = case_when(
      has_merA & has_merB ~ "merA and merB",
      !has_merA & has_merB ~ "merB only",
      has_merA & !has_merB ~ "merA only",
      TRUE ~ NA_character_
    )
  ) %>%
  ggplot(aes(x = Phylum, fill = GenomeType)) +
  geom_bar() +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5),
    legend.title = element_text(
      face = "bold",
      size = 14),
    legend.text = element_text(
      size = 12)
  ) +
  labs(
    x = "Phylum",
    y = "Count",
    fill = "Genome Type",
    title = "Distribution of Orphaned and Coupled merB by Phylum"
  ) +
  scale_fill_manual(values = c("merA and merB" = "deeppink3", "merB only" = "#56B4E9", "merA only" = "forestgreen"))

p3



p4 <- cor(
  mer_R[,c("merA_copies","merB.copies","merT.copies","merP.copies")],
  method = "spearman",
  use = "complete.obs"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene1") %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "rho") %>%
  ggplot(aes(gene1, gene2, fill = rho)) +
  geom_tile() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5),
    legend.title = element_text(
      face = "bold",
      size = 14),
    legend.text = element_text(
      size = 12),
    axis.ticks = element_blank()) +
  labs(x = NULL,
       y = NULL,
       title = "Correlation of mer Gene Copy Presence",
       fill = "Spearman ρ")

p4


p5 <- mer_R %>%
  filter(!is.na(Oxygen.Requirement)) %>%
  mutate(operon_status = ifelse(has_merA & has_merB & has_merT & has_merP, "Complete", "Incomplete")) %>% 
  ggplot(aes(x = Oxygen.Requirement, fill = operon_status)) +
  geom_bar() +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5),
    legend.title = element_text(
      face = "bold",
      size = 14),
    legend.text = element_text(
      size = 12)
  ) +
  labs(
    x = "Oxygen requirement",
    y = "Count",
    fill = "Operon",
    title = "Distribution of oxygen requirments\nby operon completeness"
  )

p5


p6 <- mer_R %>%
  mutate(
    operon_category = case_when(
      has_merB & !has_merA & !has_merT & !has_merP ~ "B",
      has_merB & !has_merA & (has_merT | has_merP) & !(has_merT & has_merP) ~ "BT/BP",
      has_merB & !has_merA & (has_merT & has_merP) ~ "BPT",
      has_merB & has_merA & !has_merT & !has_merP ~ "BA",
      has_merB & has_merA & (has_merT | has_merP) & !(has_merT & has_merP) ~ "BAP/BAT",
      has_merB & has_merA & has_merT & has_merP ~ "BAPT"
    ),
    operon_category = factor(
      operon_category,
      levels = c("B", "BT/BP", "BPT", "BA", "BAP/BAT", "BAPT")
    )
  ) %>% 
  ggplot(aes(x = operon_category, y = Genome.Size....assembled)) +
  geom_boxplot(fill = "olivedrab4") +
  stat_summary(fun = mean, 
               geom = "point", 
               shape = 20, 
               size = 3, 
               color = "deeppink3") +  # add mean points
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 16, 
                              hjust = 0.5),
  ) +
  labs(
    x = "Operon completeness",
    y = "Genome size",
    title = "Genome size by operon completeness"
  )

p6

p1
p2
p3
p4
p5
p6
a1
a2

a1 <- merB_alignment_orphaned %>%
  consensusMatrix(as.prob = TRUE) %>%
  {
    data.frame(
      Position = 1:ncol(.),
      # Calculate Variation (100% minus the most common AA %)
      Variation = (1 - apply(., 2, max)) * 100,
      AA = apply(., 2, function(x) names(which.max(x)))
    )
  } %>%
  ggplot(aes(x = Position, 
             y = Variation)) +
  geom_col(aes(fill = Variation), 
           show.legend = FALSE) +
  scale_fill_gradient(low = "lightgrey", 
                      high = "forestgreen") +
  # Make room for the labels at the bottom
  scale_y_continuous(limits = c(NA, 100), 
                     expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5)
  ) +
  labs(
    title = "Orphaned merB Positional Variation",
    y = "% Consensus",
    x = "Alignment Position"
  )

a1

a2 <- merB_alignment_coupled %>%
  consensusMatrix(as.prob = TRUE) %>%
  {
    data.frame(
      Position = 1:ncol(.),
      # Calculate Variation (100% minus the most common AA %)
      Variation = (1 - apply(., 2, max)) * 100,
      AA = apply(., 2, function(x) names(which.max(x)))
    )
  } %>%
  ggplot(aes(x = Position, 
             y = Variation)) +
  geom_col(aes(fill = Variation), 
           show.legend = FALSE) +
  scale_fill_gradient(low = "lightgrey", 
                      high = "steelblue3") +
  # Make room for the labels at the bottom
  scale_y_continuous(limits = c(NA, 100), 
                     expand = c(0, 0)) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5)
  ) +
  labs(
    title = "Coupled merB Positional Variation",
    y = "% Consensus",
    x = "Alignment Position"
  )

a2



# Exploring the graphs -----------------------------------------------------

# Finding out what genomes are merB, merT and merP
mer_R %>%
  filter(
    has_merB,
    has_merT,
    has_merP,
    !has_merA
  ) %>% select(IMG.Genome)


# Statistical analysis ----------------------------------------------------

# Checking in there is a statistical correlation for p6
mer_R %>%
  mutate(
    operon_category = case_when(
      has_merB & !has_merA & !has_merT & !has_merP ~ "b",
      has_merB & !has_merA & (has_merT | has_merP) & !(has_merT & has_merP) ~ "bt/bp",
      has_merB & !has_merA & (has_merT & has_merP) ~ "btp",
      has_merB & has_merA & !has_merT & !has_merP ~ "ba",
      has_merB & has_merA & (has_merT | has_merP) & !(has_merT & has_merP) ~ "bap/bat",
      has_merB & has_merA & has_merT & has_merP ~ "bapt",
      TRUE ~ NA_character_
    ),
    operon_category = factor(
      operon_category,
      levels = c("b", "bt/bp", "btp", "ba", "bap/bat", "bapt")
    )
  ) %>%
  filter(!is.na(Genome.Size....assembled) & !is.na(operon_category)) %>%
  {pairwise.wilcox.test(.$Genome.Size....assembled, .$operon_category, p.adjust.method = "BH")}

# Checing if there are any correlations amongst mer copies represented in p4
cor(
  mer_R[,c("merA_copies","merB.copies","merT.copies","merP.copies")],
  method = "spearman",
  use = "complete.obs"
)



ids <- c(
  2510630952,
  2510646086,
  2870236048,
  2637850024,
  2699273517
)

mer_R %>%
  filter(if_any(starts_with("merT.ID"), ~ . %in% ids))


mer_R %>% 
  filter(Oxygen.Requirement == "Anaerobe" & has_merA == TRUE) %>% 
  select(Sample.Name)

