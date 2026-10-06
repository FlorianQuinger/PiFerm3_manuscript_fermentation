library(here)

source(here("50_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

meta <- readRDS("clean/meta1.RDS")

# function to transform to relative abundances

to_rel_abd <- function(input) {
  output <- input %>%
    group_by(sampleid) %>%
    mutate(rel_abd = counts/sum(counts)*100) %>%
    ungroup() %>%
    dplyr::select(-counts) %>%
    relocate(rel_abd, .after = sampleid)
  return(output)
}


#############################
### For Ferment only
#############################

# load data 

taxonomy <- read_tsv("data/16S_ferment/taxonomy.tsv")
rarefied <- read_tsv("data/16S_ferment/rarefied_table.tsv", skip = 1) %>%
  dplyr::rename("asv" = "#OTU ID") %>%
  dplyr::rename_all(., ~gsub("IS-Pig-", "", .))
feature <- read_tsv("data/16S_ferment/feature_table.tsv", skip = 1) %>%
  dplyr::rename("asv" = "#OTU ID")  %>%
  dplyr::rename_all(., ~gsub("IS-Pig-", "", .))

# format taxonomy

taxonomy2 <- taxonomy %>%
  dplyr::select(asv = "Feature ID", taxon = Taxon) %>%
  separate(taxon, into = c("domain", "phylum", "class", "order", "family", "genus", "species") ,sep = ";") %>%
  mutate(domain = str_remove(domain, "d__"),
         phylum = str_remove(phylum, "p__"),
         class = str_remove(class, "c__"),
         order = str_remove(order, "o__"),
         family = str_remove(family, "f__"),
         genus = str_remove(genus, "g__"),
         species = str_remove(species, "s__")) %>% 
  mutate(phylum = ifelse(is.na(phylum), paste("unclassified", domain), phylum),
         class = ifelse(is.na(class), 
                        ifelse(str_detect(phylum, "unclassified"), phylum, paste("unclassified", phylum)), class),
         order = ifelse(is.na(order), 
                        ifelse(str_detect(class, "unclassified"), class, paste("unclassified", class)), order),
         family = ifelse(is.na(family), 
                         ifelse(str_detect(order, "unclassified"), order, paste("unclassified", order)), family),
         genus = ifelse(is.na(genus), 
                        ifelse(str_detect(family, "unclassified"), family, paste("unclassified", family)), genus),
         species = ifelse(is.na(species), 
                          ifelse(str_detect(genus, "unclassified"), genus, paste("unclassified", genus)), species)) %>%
  filter(!class == "Cyanobacteriia") %>%
  dplyr::select(asv, R1 = domain, P = phylum, C = class, O = order, "F" = family, G = genus, S = species)

saveRDS(taxonomy2, "clean/5_16s_ferment_asv_taxonomy.RDS")

# to long format

rarefied_long_counts <- rarefied %>%
  pivot_longer(-asv, names_to = "sampleid", values_to = "counts") %>%
  inner_join(taxonomy2, by = "asv")

saveRDS(rarefied_long_counts, "clean/5_16s_ferment_taxa_rarefied_raw_counts.RDS")

rarefied_long_counts_filtered <- rarefied_long_counts %>%
  filter_frequency_and_abundance(all_matrices = TRUE, frequency_cutoff = 0)

saveRDS(rarefied_long_counts_filtered, "clean/5_16s_ferment_taxa_rarefied_filtered_counts.RDS")

feature_long_counts <- feature %>%
  pivot_longer(-asv, names_to = "sampleid", values_to = "counts")  %>%
  inner_join(taxonomy2, by = "asv") 

saveRDS(feature_long_counts, "clean/5_16s_ferment_taxa_feature_raw_counts.RDS")

feature_long_counts_filtered <- feature_long_counts %>%
  filter_frequency_and_abundance(all_matrices = TRUE, frequency_cutoff = 0)

saveRDS(feature_long_counts_filtered, "clean/5_16s_ferment_taxa_feature_filtered_counts.RDS")

# to relative abundance

rarefied_long_rel_abd <- rarefied_long_counts %>%
  to_rel_abd()

saveRDS(rarefied_long_rel_abd, "clean/5_16s_ferment_taxa_rarefied_raw_rel_abd.RDS")

rarefied_long_rel_abd_filtered <- rarefied_long_counts_filtered %>%
  to_rel_abd()

saveRDS(rarefied_long_rel_abd_filtered, "clean/5_16s_ferment_taxa_rarefied_filtered_rel_abd.RDS")

feature_long_rel_abd <- feature_long_counts %>%
  to_rel_abd()

saveRDS(feature_long_rel_abd, "clean/5_16s_ferment_taxa_feature_raw_rel_abd.RDS")

feature_long_rel_abd_filtered <- feature_long_counts_filtered %>%
  to_rel_abd()

saveRDS(feature_long_rel_abd_filtered, "clean/5_16s_ferment_taxa_feature_filtered_rel_abd.RDS")

# to ktable style format

rarefied_ktable_counts <- calculate_rank_abundance(rarefied_long_counts, sum_column = "counts")
saveRDS(rarefied_ktable_counts, "clean/5_16s_ferment_taxonomy_rarefied_raw_counts.RDS")

rarefied_ktable_counts_filtered <- calculate_rank_abundance(rarefied_long_counts_filtered, sum_column = "counts")
saveRDS(rarefied_ktable_counts_filtered, "clean/5_16s_ferment_taxonomy_rarefied_filtered_counts.RDS")

feature_ktable_counts <- calculate_rank_abundance(feature_long_counts, sum_column = "counts")
saveRDS(feature_ktable_counts, "clean/5_16s_ferment_taxonomy_feature_raw_counts.RDS")

feature_ktable_counts_filtered <- calculate_rank_abundance(feature_long_counts_filtered, sum_column = "counts")
saveRDS(feature_ktable_counts_filtered, "clean/5_16s_ferment_taxonomy_feature_filtered_counts.RDS")

rarefied_ktable_rel_abd <- calculate_rank_abundance(rarefied_long_rel_abd, sum_column = "rel_abd")
saveRDS(rarefied_ktable_rel_abd, "clean/5_16s_ferment_taxonomy_rarefied_raw_rel_abd.RDS")

rarefied_ktable_rel_abd_filtered <- calculate_rank_abundance(rarefied_long_rel_abd_filtered, sum_column = "rel_abd")
saveRDS(rarefied_ktable_rel_abd_filtered, "clean/5_16s_ferment_taxonomy_rarefied_filtered_rel_abd.RDS")

feature_ktable_rel_abd <- calculate_rank_abundance(feature_long_rel_abd, sum_column = "rel_abd")
saveRDS(feature_ktable_rel_abd, "clean/5_16s_ferment_taxonomy_feature_raw_rel_abd.RDS")

feature_ktable_rel_abd_filtered <- calculate_rank_abundance(feature_long_rel_abd_filtered, sum_column = "rel_abd")
saveRDS(feature_ktable_rel_abd_filtered, "clean/5_16s_ferment_taxonomy_feature_filtered_rel_abd.RDS")

