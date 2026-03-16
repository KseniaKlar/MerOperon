# Ksenia Klar
# March 13th 2026
# In this script I will learn to visulaize my tree from GToTree

# clearing the environment
rm(list=ls())

# loading the libraries
library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)


# Load your tree
tree <- read.newick("GTT/tree.tre")
ko_hits <- read_tsv("GTT/KO_counts.tsv")

ko_hits <- ko_hits %>% 
  select(assembly_id, K00520, K00221) %>%
  column_to_rownames("assembly_id")

ko_hits <- ko_hits %>%
  rownames_to_column(var = "accession")

species_map <- data.frame(
  accession = c(
    "GCA_000015365.1",
    "GCA_000213595.1",
    "GCA_000702605.1",
    "GCA_000709065.1",
    "GCA_000710085.1",
    "GCA_000787015.1",
    "GCA_000960735.1",
    "GCA_001312365.1",
    "GCA_001316705.2",
    "GCA_001530215.1",
    "GCA_001591465.1",
    "GCA_001619405.1",
    "GCA_001715785.1",
    "GCA_002159175.1",
    "GCA_900106065.1"
  ),
  species = c(
    "Marinobacter nauticus VT8",
    "Desmospora sp. 8437",
    "Exiguobacterium acetylicum DSM 20416",
    "Cupriavidus metallidurans NE12",
    "Pseudomonas syringae pv. atrofaciens LMG 5095",
    "Acinetobacter baumannii 269",
    "Listeria innocua 12KSM",
    "Variovorax sp. RA8 JCM 16519",
    "Enterobacter kobei ST54:941713674",
    "Burkholderia vietnamiensis HI3908",
    "Cytobacillus firmus NBRC 15306",
    "Bacillus cereus B4155",
    "Brochothrix thermosphacta Bth-7816",
    "Flavonifractor sp. An92",
    "Pseudomonas mandelii LMG 21607"
  ),
  stringsAsFactors = FALSE
)

ko_hits <- ko_hits %>%
  left_join(species_map, by = "accession")



# Prepare gene info
gene_info <- ko_hits %>%
  mutate(
    # Factor for tip color with levels
    tip_color = factor(
      case_when(
        K00520 > 0 & K00221 > 0 ~ "merA & merB",
        K00520 > 0 & K00221 == 0 ~ "Only merA",
        K00520 == 0 & K00221 > 0 ~ "Only merB",
        TRUE ~ "Neither"
      ),
      levels = c("merA & merB", "Only merA", "Only merB", "Neither")
    ),
    # Repeat letters for copy numbers
    label_AB = paste0(
      ifelse(K00520 > 0, strrep("A", K00520), ""),
      ifelse(K00221 > 0, strrep("B", K00221), "")
    )
  )

gene_info <- gene_info %>%
  mutate(label = str_replace(species, "^((?:\\S+\\s){2})", "\\1\n")) %>% 
  select(-species) %>% 
  dplyr::rename(species = label)


# Color mapping for legend
color_mapping <- c(
  "merA & merB" = "forestgreen",
  "Only merA" = "deeppink3",
  "Only merB" = "lightblue3",
  "Neither" = "black"
)

# Circular tree
ggtree(tree, layout = "circular") %<+% gene_info +
  geom_tiplab(aes(label = species), 
              size = 4, 
              align = TRUE,
              show.legend = FALSE) +  # tip labels colored
  geom_tippoint(aes(color = tip_color), size = 2.5) +
  geom_text2(aes(label = label_AB), 
             nudge_x = 0.03,
             nudge_y = -0.21, 
             size = 2,
             color = "red") +  # repeated letters
  scale_color_manual(values = color_mapping, 
                     name = "Gene Presence") +   # legend
  theme_tree2() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = c(0.9, 0.56)
  ) +
  labs(title = NULL) +
  xlim_tree(2)


ggsave(
  "mer_tree.png",
  p,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave("tree.png",
       p + coord_fixed() + xlim_tree(0, 1.2),
       width = 12, height = 12, dpi = 300)


