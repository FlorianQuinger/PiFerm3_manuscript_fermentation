
library(here)
library(NormalyzerDE)
library(DIMAR)

source(here("50_omics_functions.R"))

# function to transform to relative abundances

to_rel_abd <- function(input) {
  output <- input %>%
    group_by(sampleid) %>%
    mutate(rel_abd = intensity/sum(intensity)*100) %>%
    ungroup() %>%
    dplyr::select(-intensity) %>%
    relocate(rel_abd, .after = sampleid)
  return(output)
}

# function to transform to log 2 (+1)

to_log2 <- function(input) {
  output <- input %>%
    mutate(log2 = log2(intensity+1)) %>%
    dplyr::select(-intensity) %>%
    relocate(log2, .after = sampleid)
  return(output)
}

# load metadata

meta <- readRDS("clean/meta1.RDS")

########################
### Fermented Pea
#######################

#folder <- here("data/MetaLab_ferment_dual/") # previous Output used, but apparently with pfind 3.2.3 and without sample 405
folder <- here("data/MetaLab1.1_ferment_dual_3.2.0/")

# summary file
summary <- read_tsv(here(folder, "final_summary.tsv")) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  mutate(sampleid = str_extract(str_remove(raw_file, "Florian_"), "\\d{3}")) %>%
  filter(!is.na(sampleid)) %>%
  dplyr::select(sampleid, 'ms/ms', 'ms/ms_identified', 'ms/ms_identified_[%]', 'peptide_sequences_identified')

saveRDS(summary, "clean/6_ferment_summary_wide.RDS")

summary_long <- summary %>%
  pivot_longer(-sampleid, names_to = "parameter", values_to = "value") %>%
  saveRDS("clean/6_ferment_summary_long.RDS")

# filter Final proteins file

nrow(read_tsv(here(folder, "final_proteins.tsv"))) # ca. 3400 proteins

proteins <- read_tsv(here(folder, "final_proteins.tsv")) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  filter(`peptide_counts_(all)` > 1) %>% # remove one hit wonders
  dplyr::select(-intensity) %>% # remove total intensity column
  pivot_longer(starts_with("intensity"), names_to = "sampleid", values_to = "intensity") %>%
  mutate(sampleid = str_extract(str_remove(sampleid, "intensity_intensity_florian_"), "\\d{3}")) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(proteingroup = protein_ids) %>% # keep all ids as proteingroups
  separate(protein_ids, into = "proteinid", sep = ";") %>% # extract first protein 
  dplyr::select(proteingroup, proteinid, sampleid, intensity) # right order for processing

add_origin <- function(input_df) {
  output <- input_df %>%
    mutate(origin = case_when(str_detect(proteinid, "MGYG") ~ "micro",
                              str_detect(proteinid, "PIG") ~ "pig",
                              str_detect(proteinid, "PEA") ~ "pea",
                              str_detect(proteinid, "WHEAT") ~ "wheat",
                              str_detect(proteinid, "HORVV") ~ "barley",
                              str_detect(proteinid, "BRANA") ~ "rape",
                              .default = "error")) %>% # add origin by leading protein
    dplyr::select(proteinid, proteingroup, sampleid, intensity, origin) # select relevant columns
  return(output)
}

proteins_unfiltered <- filter_frequency_and_abundance(proteins, feature_columns = c("proteingroup", "proteinid"),
                                                    abundance_cutoff = 0, frequency_cutoff = 0, filter_df = FALSE) %>% 
  add_origin() # only filtered by frequency
saveRDS(proteins_unfiltered, "clean/6_ferment_proteins_long_unfiltered_intensity.RDS")


proteins_filtered <- filter_frequency_and_abundance(proteins, feature_columns = c("proteingroup", "proteinid"),
                                                    abundance_cutoff = 0, filter_df = FALSE) %>% 
  add_origin() # only filtered by frequency
saveRDS(proteins_filtered, "clean/6_ferment_proteins_long_raw_intensity.RDS")

proteins_filtered2 <- filter_frequency_and_abundance(proteins, feature_columns = c("proteingroup", "proteinid"),
                                                     abundance_cutoff = 5e-4, filter_df = FALSE) %>% 
  add_origin() # also filtered by abundance
saveRDS(proteins_filtered2, "clean/6_ferment_proteins_long_raw_intensity_filtered.RDS")

# transform to relative abundances

proteins_rel_abd <- to_rel_abd(proteins_filtered)
saveRDS(proteins_rel_abd, "clean/6_ferment_proteins_long_raw_rel_abd.RDS")

proteins_rel_abd_filtered <- to_rel_abd(proteins_filtered2)
saveRDS(proteins_rel_abd_filtered, "clean/6_ferment_proteins_long_raw_rel_abd_filtered.RDS")

# transform to log2

proteins_log2 <- to_log2(proteins_filtered)
saveRDS(proteins_log2, "clean/6_ferment_proteins_long_raw_log2.RDS")

proteins_log2_filtered <- to_log2(proteins_filtered2)
saveRDS(proteins_log2_filtered, "clean/6_ferment_proteins_long_raw_log2_filtered.RDS")

# Normalize Protein data using NormalyzerDE
# from Ileum # using filtered2 proteins

proteins_to_normalize <- proteins_filtered %>%
  dplyr::select(proteinid, sampleid, intensity) %>%
  mutate(sampleid = str_c("s", sampleid)) %>% # workaround to get rid of unmatched colnames due to sci notation
  pivot_wider(names_from = "sampleid", values_from = "intensity") %>%
  write_tsv("temp/ferment_proteins_to_normalize.tsv")
design_to_normalize <- proteins_filtered %>%
  left_join(meta, by = "sampleid") %>%
  mutate(diet = ifelse(diet == "P3D7", "3", diet)) %>% # workaround: set P3D7 to 3, so it is replicated
  mutate(sampleid = str_c("s", sampleid)) %>% # workaround to get rid of unmatched colnames due to sci notation
  dplyr::select(sample = sampleid, group = diet) %>%
  distinct() %>%
  write_tsv("temp/ferment_design_to_normalize.tsv")

# normalyzer(jobName = "proteins_ferment_normalized",
#            designPath = "temp/ferment_design_to_normalize.tsv",
#            dataPath = "temp/ferment_proteins_to_normalize.tsv",
#            outputDir = "temp/") # CycLoess best, median also good

proteins_normalized <- as.data.frame(read_tsv("temp/proteins_ferment_normalized/CycLoess-normalized.txt"))

# Impute with dimar

proteins_to_impute <- as.matrix(proteins_normalized[,-1])
rownames(proteins_to_impute) <- proteins_normalized$proteinid

dimarPlotHeatmap(proteins_to_impute)

proteins_to_impute_test <- proteins_to_impute[sample(1:nrow(proteins_to_impute), 1000),1:7]
# dimar(mtx = proteins_to_impute_test,
#             group = read_tsv("temp/ferment_design_to_normalize.tsv")$group) # too less observations

proteins_imputed <- DIMAR::dimarDoOptimalImputation(proteins_to_impute, method = "impSeqRob")[[1]]

dimarPlotHeatmap(proteins_imputed)

# clean file

proteins_norm_imp <- proteins_imputed %>%
  as_tibble(rownames = "proteinid") %>%
  pivot_longer(-proteinid, names_to = "sampleid", values_to = "log2") %>%
  mutate(sampleid = str_extract(sampleid, "\\d{3}")) %>% # remove "s" from normalization
  left_join(distinct(dplyr::select(proteins_filtered, proteinid, proteingroup, origin)), by = "proteinid") %>%
  dplyr::select(proteinid, proteingroup, sampleid, log2, origin)

saveRDS(proteins_norm_imp, "clean/6_ferment_proteins_long_norm_imp_log2.RDS")


# as intensities

proteins_norm_imp_intensity <- proteins_norm_imp %>%
  mutate(intensity = (2^log2)-1) %>%
  dplyr::select(proteinid, proteingroup, sampleid, intensity, origin)

saveRDS(proteins_norm_imp_intensity, "clean/6_ferment_proteins_long_norm_imp_intensity.RDS")

proteins_norm_imp_intensity_filtered <- filter_frequency_and_abundance(proteins_norm_imp_intensity, 
                                                                       feature_columns = c("proteinid", "proteingroup"), 
                                                                       abundance_cutoff = 5e-4, filter_df = FALSE)

saveRDS(proteins_norm_imp_intensity_filtered, "clean/6_ferment_proteins_long_norm_imp_intensity_filtered.RDS")

# as relative abundances

proteins_norm_imp_rel_abd <- to_rel_abd(proteins_norm_imp_intensity)

saveRDS(proteins_norm_imp_rel_abd, "clean/6_ferment_proteins_long_norm_imp_rel_abd.RDS")

proteins_norm_imp_rel_abd_filtered <- filter_frequency_and_abundance(proteins_norm_imp_rel_abd, 
                                                                     feature_columns = c("proteinid", "proteingroup"),
                                                                     abundance_cutoff = 5e-4, filter_df = FALSE)

saveRDS(proteins_norm_imp_rel_abd_filtered, "clean/6_ferment_proteins_long_norm_imp_rel_abd_filtered.RDS")

# to log2

proteins_norm_imp_log2_filtered <- to_log2(proteins_norm_imp_intensity_filtered)

saveRDS(proteins_norm_imp_log2_filtered, "clean/6_ferment_proteins_long_norm_imp_log2_filtered.RDS")

# to taxa

#prepare genomes file
genomes <- read_tsv(here(folder, "Genome.tsv")) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  dplyr::select(bin = genome, starts_with("gtdb"), -c("gtdb_name", "gtdb_rank", "gtdb_kingdom")) %>%
  dplyr::rename("R1" = gtdb_superkingdom, P = gtdb_phylum, C = gtdb_class, O = gtdb_order,
                "F" = gtdb_family, "G" = gtdb_genus, "S" = gtdb_species) # 19 taxa


# for unnormalised unimputed file 
taxa_raw_intensity <- calculate_taxa_intensity(proteins_filtered, genomes = genomes)

saveRDS(taxa_raw_intensity, "clean/6_ferment_taxa_long_raw_intensity.RDS")

taxa_raw_intensity_filtered <- calculate_taxa_intensity(proteins_filtered2, genomes = genomes)

saveRDS(taxa_raw_intensity_filtered, "clean/6_ferment_taxa_long_raw_intensity_filtered.RDS")

# for imputed files
taxa_norm_imp_intensity <- calculate_taxa_intensity(proteins_norm_imp_intensity, genomes = genomes)

saveRDS(taxa_norm_imp_intensity, "clean/6_ferment_taxa_long_norm_imp_intensity.RDS")

taxa_norm_imp_intensity_filtered <- calculate_taxa_intensity(proteins_norm_imp_intensity_filtered, genomes = genomes)

saveRDS(taxa_norm_imp_intensity_filtered, "clean/6_ferment_taxa_long_norm_imp_intensity_filtered.RDS")

# as relative abundances

taxa_raw_rel_abd <- to_rel_abd(taxa_raw_intensity)
saveRDS(taxa_raw_rel_abd, "clean/6_ferment_taxa_long_raw_rel_abd.RDS")

taxa_raw_rel_abd_filtered <- to_rel_abd(taxa_raw_intensity_filtered)
saveRDS(taxa_raw_rel_abd_filtered, "clean/6_ferment_taxa_long_raw_rel_abd_filtered.RDS")

taxa_norm_imp_rel_abd <- to_rel_abd(taxa_norm_imp_intensity)
saveRDS(taxa_norm_imp_rel_abd, "clean/6_ferment_taxa_long_norm_imp_rel_abd.RDS")

taxa_norm_imp_rel_abd_filtered <- to_rel_abd(taxa_norm_imp_intensity_filtered)
saveRDS(taxa_norm_imp_rel_abd_filtered, "clean/6_ferment_taxa_long_norm_imp_rel_abd_filtered.RDS")

# as log2

taxa_raw_log2 <- to_log2(taxa_raw_intensity)
saveRDS(taxa_raw_log2, "clean/6_ferment_taxa_long_raw_log2.RDS")

taxa_raw_log2_filtered <- to_log2(taxa_raw_intensity_filtered)
saveRDS(taxa_raw_log2_filtered, "clean/6_ferment_taxa_long_raw_log2_filtered.RDS")

taxa_norm_imp_log2 <- to_log2(taxa_norm_imp_intensity)
saveRDS(taxa_norm_imp_log2, "clean/6_ferment_taxa_long_norm_imp_log2.RDS")

taxa_norm_imp_log2_filtered <- to_log2(taxa_norm_imp_intensity_filtered)
saveRDS(taxa_norm_imp_log2_filtered, "clean/6_ferment_taxa_long_norm_imp_log2_filtered.RDS")

# agglomerate at different levels # removed filtering since filtered proteins are used

taxonomy_raw_intensity <- calculate_rank_abundance(taxa_raw_intensity[,-1], sum_column = "intensity")
#taxonomy_raw_intensity_filtered <- filter_frequency_and_abundance_ktable(taxonomy_raw_intensity)
saveRDS(taxonomy_raw_intensity, "clean/6_ferment_taxonomy_long_raw_intensity.RDS")

taxonomy_raw_intensity_filtered <- calculate_rank_abundance(taxa_raw_intensity_filtered[,-1], sum_column = "intensity")
saveRDS(taxonomy_raw_intensity_filtered, "clean/6_ferment_taxonomy_long_raw_intensity_filtered.RDS")

taxonomy_norm_imp_intensity <- calculate_rank_abundance(taxa_norm_imp_intensity[,-1], sum_column = "intensity")
#taxonomy_norm_imp_intensity_filtered <- filter_frequency_and_abundance_ktable(taxonomy_norm_imp_intensity)
saveRDS(taxonomy_norm_imp_intensity, "clean/6_ferment_taxonomy_long_norm_imp_intensity.RDS")

taxonomy_norm_imp_intensity_filtered <- calculate_rank_abundance(taxa_norm_imp_intensity_filtered[,-1], sum_column = "intensity")
saveRDS(taxonomy_norm_imp_intensity_filtered, "clean/6_ferment_taxonomy_long_norm_imp_intensity_filtered.RDS")

# for rel_abd

taxonomy_raw_rel_abd <- calculate_rank_abundance(taxa_raw_rel_abd[,-1], sum_column = "rel_abd")
#taxonomy_raw_rel_abd_filtered <- filter_frequency_and_abundance_ktable(taxonomy_raw_rel_abd)
saveRDS(taxonomy_raw_rel_abd, "clean/6_ferment_taxonomy_long_raw_rel_abd.RDS")

taxonomy_raw_rel_abd_filtered <- calculate_rank_abundance(taxa_raw_rel_abd_filtered[,-1], sum_column = "rel_abd")
saveRDS(taxonomy_raw_rel_abd_filtered, "clean/6_ferment_taxonomy_long_raw_rel_abd_filtered.RDS")

taxonomy_norm_imp_rel_abd <- calculate_rank_abundance(taxa_norm_imp_rel_abd[,-1], sum_column = "rel_abd")
#taxonomy_norm_imp_rel_abd_filtered <- filter_frequency_and_abundance_ktable(taxonomy_norm_imp_rel_abd)
saveRDS(taxonomy_norm_imp_rel_abd, "clean/6_ferment_taxonomy_long_norm_imp_rel_abd.RDS")

taxonomy_norm_imp_rel_abd_filtered <- calculate_rank_abundance(taxa_norm_imp_rel_abd_filtered[,-1], sum_column = "rel_abd")
saveRDS(taxonomy_norm_imp_rel_abd_filtered, "clean/6_ferment_taxonomy_long_norm_imp_rel_abd_filtered.RDS")

# to functions

# protein functions
func <- read_tsv(here(folder, "functions.tsv")) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  dplyr::select(-starts_with("Intensity"), -group_id, -peptide_count, -psm_count, -pep) %>%
  dplyr::rename(proteinid = name, gos = gene_ontology_id)

saveRDS(func, "clean/6_ferment_functions_long_raw.RDS")

func_meta <- func %>%
  dplyr::select(-proteinid, -protein_name) %>% # protein_name differs between seed_orthologs
  distinct() 

# function to calculate function intensity 

calculate_function_intensity <- function(input) { # without dividing intensity of proteingroups
  output <- input %>% 
    dplyr::select(-proteinid) %>% #remove to mitigate duplicates
    separate_rows(proteingroup, sep = ";") %>%
    filter(str_detect(proteingroup, "^MGYG")) %>% # filter out host proteins
    dplyr::rename(proteinid = proteingroup) %>%
    group_by(sampleid, proteinid) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>% # sum up intensities per sample and protein
    pivot_wider(names_from = sampleid, values_from = intensity) %>%
    pivot_longer(-proteinid, names_to = "sampleid", values_to = "intensity") %>%
    mutate(intensity = ifelse(is.na(intensity), 0, intensity)) %>% # replace NAs with 0
    inner_join(func, by = "proteinid") %>% # lost some proteins
    dplyr::select(seed_ortholog, sampleid, intensity) %>% 
    group_by(sampleid, seed_ortholog) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>%
    pivot_wider(names_from = sampleid, values_from = intensity) %>%
    pivot_longer(-seed_ortholog, names_to = "sampleid", values_to = "intensity") %>%
    mutate(intensity = ifelse(is.na(intensity), 0, intensity)) %>%
    inner_join(func_meta, by = "seed_ortholog")
  return(output)
}

# for untransformed data

functions_intensity <- calculate_function_intensity(proteins_filtered)
saveRDS(functions_intensity, "clean/6_ferment_functions_long_raw_intensity.RDS")

functions_intensity_filtered <- calculate_function_intensity(proteins_filtered2)
saveRDS(functions_intensity_filtered, "clean/6_ferment_function_long_raw_intensity_filtered.RDS")

# for transformed data

functions_norm_imp_intensity <- calculate_function_intensity(proteins_norm_imp_intensity)
saveRDS(functions_norm_imp_intensity, "clean/6_ferment_functions_long_norm_imp_intensity.RDS")

functions_norm_imp_intensity_filtered <- calculate_function_intensity(proteins_norm_imp_intensity_filtered)
saveRDS(functions_norm_imp_intensity_filtered, "clean/6_ferment_functions_long_norm_imp_intensity_filtered.RDS")

# to rel_abd

functions_rel_abd <- to_rel_abd(functions_intensity)
saveRDS(functions_rel_abd, "clean/6_ferment_functions_long_raw_rel_abd.RDS")

functions_rel_abd_filtered <- to_rel_abd(functions_intensity_filtered)
saveRDS(functions_rel_abd_filtered, "clean/6_ferment_functions_long_raw_rel_abd_filtered.RDS")

functions_norm_imp_rel_abd <- to_rel_abd(functions_norm_imp_intensity)
saveRDS(functions_norm_imp_rel_abd, "clean/6_ferment_functions_long_norm_imp_rel_abd.RDS")

functions_norm_imp_rel_abd_filtered <- to_rel_abd(functions_norm_imp_intensity_filtered)
saveRDS(functions_norm_imp_rel_abd_filtered, "clean/6_ferment_functions_long_norm_imp_rel_abd_filtered.RDS")

# to log2

functions_log2 <- to_log2(functions_intensity)
saveRDS(functions_log2, "clean/6_ferment_functions_long_raw_log2.RDS")

functions_log2_filtered <- to_log2(functions_intensity_filtered)
saveRDS(functions_log2_filtered, "clean/6_ferment_functions_long_raw_log2_filtered.RDS")

functions_norm_imp_log2 <- to_log2(functions_norm_imp_intensity)
saveRDS(functions_norm_imp_log2, "clean/6_ferment_functions_long_norm_imp_log2.RDS")

functions_norm_imp_log2_filtered <- to_log2(functions_norm_imp_intensity_filtered)
saveRDS(functions_norm_imp_log2_filtered, "clean/6_ferment_functions_long_norm_imp_log2_filtered.RDS")

# to Kegg functions

kegg_intensity <- calculate_kegg_abundance(functions_intensity, abundance_column = "intensity")
saveRDS(kegg_intensity, "clean/6_ferment_kegg_long_raw_intensity.RDS")

kegg_intensity_filtered <- calculate_kegg_abundance(functions_intensity_filtered, abundance_column = "intensity")
saveRDS(kegg_intensity_filtered, "clean/6_ferment_kegg_long_raw_intensity_filtered.RDS")

kegg_norm_imp_intensity <- calculate_kegg_abundance(functions_norm_imp_intensity, abundance_column = "intensity")
saveRDS(kegg_norm_imp_intensity, "clean/6_ferment_kegg_long_norm_imp_intensity.RDS")

kegg_norm_imp_intensity_filtered <- calculate_kegg_abundance(functions_norm_imp_intensity_filtered, abundance_column = "intensity")
saveRDS(kegg_norm_imp_intensity_filtered, "clean/6_ferment_kegg_long_norm_imp_intensity_filtered.RDS")

# as rel_abd

kegg_rel_abd <- to_rel_abd(kegg_intensity)
saveRDS(kegg_rel_abd, "clean/6_ferment_kegg_long_raw_rel_abd.RDS")

kegg_rel_abd_filtered <- to_rel_abd(kegg_intensity_filtered)
saveRDS(kegg_rel_abd_filtered, "clean/6_ferment_kegg_long_raw_rel_abd_filtered.RDS")

kegg_norm_imp_rel_abd <- to_rel_abd(kegg_norm_imp_intensity)
saveRDS(kegg_norm_imp_rel_abd, "clean/6_ferment_kegg_long_norm_imp_rel_abd.RDS")

kegg_norm_imp_rel_abd_filtered <- to_rel_abd(kegg_norm_imp_intensity_filtered)
saveRDS(kegg_norm_imp_rel_abd_filtered, "clean/6_ferment_kegg_long_norm_imp_rel_abd_filtered.RDS")

# as log2

kegg_log2 <- to_log2(kegg_intensity)
saveRDS(kegg_log2, "clean/6_ferment_kegg_long_raw_log2.RDS")

kegg_log2_filtered <- to_log2(kegg_intensity_filtered)
saveRDS(kegg_log2_filtered, "clean/6_ferment_kegg_long_raw_log2_filtered.RDS")

kegg_norm_imp_log2 <- to_log2(kegg_norm_imp_intensity)
saveRDS(kegg_norm_imp_log2, "clean/6_ferment_kegg_long_norm_imp_log2.RDS")

kegg_norm_imp_log2_filtered <- to_log2(kegg_norm_imp_intensity_filtered)
saveRDS(kegg_norm_imp_log2_filtered, "clean/6_ferment_kegg_long_norm_imp_log2_filtered.RDS")

# to GO

go_intensity <- calculate_go_abundance(functions_intensity, abundance_column = "intensity")
saveRDS(go_intensity, "clean/6_ferment_go_long_raw_intensity.RDS")

go_intensity_filtered <- calculate_go_abundance(functions_intensity_filtered, abundance_column = "intensity")
saveRDS(go_intensity_filtered, "clean/6_ferment_go_long_raw_intensity_filtered.RDS")

go_norm_imp_intensity <- calculate_go_abundance(functions_norm_imp_intensity, abundance_column = "intensity")
saveRDS(go_norm_imp_intensity, "clean/6_ferment_go_long_norm_imp_intensity.RDS")

go_norm_imp_intensity_filtered <- calculate_go_abundance(functions_norm_imp_intensity_filtered, abundance_column = "intensity")
saveRDS(go_norm_imp_intensity_filtered, "clean/6_ferment_go_long_norm_imp_intensity_filtered.RDS")

# to rel_abd

go_rel_abd <- to_rel_abd(go_intensity)
saveRDS(go_rel_abd, "clean/6_ferment_go_long_raw_rel_abd.RDS")

go_rel_abd_filtered <- to_rel_abd(go_intensity_filtered)
saveRDS(go_rel_abd_filtered, "clean/6_ferment_go_long_raw_rel_abd_filtered.RDS")

go_norm_imp_rel_abd <- to_rel_abd(go_norm_imp_intensity)
saveRDS(go_norm_imp_rel_abd, "clean/6_ferment_go_long_norm_imp_rel_abd.RDS")

go_norm_imp_rel_abd_filtered <- to_rel_abd(go_norm_imp_intensity_filtered)
saveRDS(go_norm_imp_rel_abd_filtered, "clean/6_ferment_go_long_norm_imp_rel_abd_filtered.RDS")

# to log2

go_log2 <- to_log2(go_intensity)
saveRDS(go_log2, "clean/6_ferment_go_long_raw_log2.RDS")

go_log2_filtered <- to_log2(go_intensity_filtered)
saveRDS(go_log2_filtered, "clean/6_ferment_go_long_raw_log2_filtered.RDS")

go_norm_imp_log2 <- to_log2(go_norm_imp_intensity)
saveRDS(go_norm_imp_log2, "clean/6_ferment_go_long_norm_imp_log2.RDS")

go_norm_imp_log2_filtered <- to_log2(go_norm_imp_intensity_filtered)
saveRDS(go_norm_imp_log2_filtered, "clean/6_ferment_go_long_norm_imp_log2_filtered.RDS")

# to COGs

cog_intensity <- calculate_cog_abundance(functions_intensity, abundance_column = "intensity")
saveRDS(cog_intensity, "clean/6_ferment_cog_long_raw_intensity.RDS")

cog_intensity_filtered <- calculate_cog_abundance(functions_intensity_filtered, abundance_column = "intensity")
saveRDS(cog_intensity_filtered, "clean/6_ferment_cog_long_raw_intensity_filtered.RDS")

cog_norm_imp_intensity <- calculate_cog_abundance(functions_norm_imp_intensity, abundance_column = "intensity")
saveRDS(cog_norm_imp_intensity, "clean/6_ferment_cog_long_norm_imp_intensity.RDS")

cog_norm_imp_intensity_filtered <- calculate_cog_abundance(functions_norm_imp_intensity_filtered, abundance_column = "intensity")
saveRDS(cog_norm_imp_intensity_filtered, "clean/6_ferment_cog_long_norm_imp_intensity_filtered.RDS")

# to rel_abd

cog_rel_abd <- to_rel_abd(cog_intensity)
saveRDS(cog_rel_abd, "clean/6_ferment_cog_long_raw_rel_abd.RDS")

cog_rel_abd_filtered <- to_rel_abd(cog_intensity_filtered)
saveRDS(cog_rel_abd_filtered, "clean/6_ferment_cog_long_raw_rel_abd_filtered.RDS")

cog_norm_imp_rel_abd <- to_rel_abd(cog_norm_imp_intensity)
saveRDS(cog_norm_imp_rel_abd, "clean/6_ferment_cog_long_norm_imp_rel_abd.RDS")

cog_norm_imp_rel_abd_filtered <- to_rel_abd(cog_norm_imp_intensity_filtered)
saveRDS(cog_norm_imp_rel_abd_filtered, "clean/6_ferment_cog_long_norm_imp_rel_abd_filtered.RDS")

# to log2

cog_log2 <- to_log2(cog_intensity)
saveRDS(cog_log2, "clean/6_ferment_cog_long_raw_log2.RDS")

cog_log2_filtered <- to_log2(cog_intensity_filtered)
saveRDS(cog_log2_filtered, "clean/6_ferment_cog_long_raw_log2_filtered.RDS")

cog_norm_imp_log2 <- to_log2(cog_norm_imp_intensity) 
saveRDS(cog_norm_imp_log2, "clean/6_ferment_cog_long_norm_imp_log2.RDS")

cog_norm_imp_log2_filtered <- to_log2(cog_norm_imp_intensity_filtered)
saveRDS(cog_norm_imp_log2_filtered, "clean/6_ferment_cog_long_norm_imp_log2_filtered.RDS")

# to Cazy # changed for cazy reannotation

cazy_intensity <- calculate_cazy_abundance(functions_intensity, abundance_column = "intensity", cazy_column = "cazy2")
saveRDS(cazy_intensity, "clean/6_ferment_cazy_long_raw_intensity.RDS")

cazy_intensity_filtered <- calculate_cazy_abundance(functions_intensity_filtered, abundance_column = "intensity", cazy_column = "cazy2")
saveRDS(cazy_intensity_filtered, "clean/6_ferment_cazy_long_raw_intensity_filtered.RDS")

cazy_norm_imp_intensity <- calculate_cazy_abundance(functions_norm_imp_intensity, abundance_column = "intensity", cazy_column = "cazy2")
saveRDS(cazy_norm_imp_intensity, "clean/6_ferment_cazy_long_norm_imp_intensity.RDS")

cazy_norm_imp_intensity_filtered <- calculate_cazy_abundance(functions_norm_imp_intensity_filtered, abundance_column = "intensity", cazy_column = "cazy2")
saveRDS(cazy_norm_imp_intensity_filtered, "clean/6_ferment_cazy_long_norm_imp_intensity_filtered.RDS")

# to rel_abd

cazy_rel_abd <- to_rel_abd(cazy_intensity)
saveRDS(cazy_rel_abd, "clean/6_ferment_cazy_long_raw_rel_abd.RDS")

cazy_rel_abd_filtered <- to_rel_abd(cazy_intensity_filtered)
saveRDS(cazy_rel_abd_filtered, "clean/6_ferment_cazy_long_raw_rel_abd_filtered.RDS")

cazy_norm_imp_rel_abd <- to_rel_abd(cazy_norm_imp_intensity)
saveRDS(cazy_norm_imp_rel_abd, "clean/6_ferment_cazy_long_norm_imp_rel_abd.RDS")

cazy_norm_imp_rel_abd_filtered <- to_rel_abd(cazy_norm_imp_intensity_filtered)
saveRDS(cazy_norm_imp_rel_abd_filtered, "clean/6_ferment_cazy_long_norm_imp_rel_abd_filtered.RDS")

# to log2

cazy_log2 <- to_log2(cazy_intensity)
saveRDS(cazy_log2, "clean/6_ferment_cazy_long_raw_log2.RDS")

cazy_log2_filtered <- to_log2(cazy_intensity_filtered)
saveRDS(cazy_log2_filtered, "clean/6_ferment_cazy_long_raw_log2_filtered.RDS")

cazy_norm_imp_log2 <- to_log2(cazy_norm_imp_intensity) 
saveRDS(cazy_norm_imp_log2, "clean/6_ferment_cazy_long_norm_imp_log2.RDS")

cazy_norm_imp_log2_filtered <- to_log2(cazy_norm_imp_intensity_filtered)
saveRDS(cazy_norm_imp_log2_filtered, "clean/6_ferment_cazy_long_norm_imp_log2_filtered.RDS")

# to EC

ec_intensity <- calculate_ec_abundance(functions_intensity, abundance_column = "intensity", ec_column = "ec_id2")
saveRDS(ec_intensity, "clean/6_ferment_ec_long_raw_intensity.RDS")

ec_intensity_filtered <- calculate_ec_abundance(functions_intensity_filtered, abundance_column = "intensity", ec_column = "ec_id2")
saveRDS(ec_intensity_filtered, "clean/6_ferment_ec_long_raw_intensity_filtered.RDS")

ec_norm_imp_intensity <- calculate_ec_abundance(functions_norm_imp_intensity, abundance_column = "intensity", ec_column = "ec_id2")
saveRDS(ec_norm_imp_intensity, "clean/6_ferment_ec_long_norm_imp_intensity.RDS")

ec_norm_imp_intensity_filtered <- calculate_ec_abundance(functions_norm_imp_intensity_filtered, abundance_column = "intensity", ec_column = "ec_id2")
saveRDS(ec_norm_imp_intensity_filtered, "clean/6_ferment_ec_long_norm_imp_intensity_filtered.RDS")

# to rel_abd

ec_rel_abd <- to_rel_abd(ec_intensity)
saveRDS(ec_rel_abd, "clean/6_ferment_ec_long_raw_rel_abd.RDS")

ec_rel_abd_filtered <- to_rel_abd(ec_intensity_filtered)
saveRDS(ec_rel_abd_filtered, "clean/6_ferment_ec_long_raw_rel_abd_filtered.RDS")

ec_norm_imp_rel_abd <- to_rel_abd(ec_norm_imp_intensity)
saveRDS(ec_norm_imp_rel_abd, "clean/6_ferment_ec_long_norm_imp_rel_abd.RDS")

ec_norm_imp_rel_abd_filtered <- to_rel_abd(ec_norm_imp_intensity_filtered)
saveRDS(ec_norm_imp_rel_abd_filtered, "clean/6_ferment_ec_long_norm_imp_rel_abd_filtered.RDS")

# to log2

ec_log2 <- to_log2(ec_intensity)
saveRDS(ec_log2, "clean/6_ferment_ec_long_raw_log2.RDS")

ec_log2_filtered <- to_log2(ec_intensity_filtered)
saveRDS(ec_log2_filtered, "clean/6_ferment_ec_long_raw_log2_filtered.RDS")

ec_norm_imp_log2 <- to_log2(ec_norm_imp_intensity) 
saveRDS(ec_norm_imp_log2, "clean/6_ferment_ec_long_norm_imp_log2.RDS")

ec_norm_imp_log2_filtered <- to_log2(ec_norm_imp_intensity_filtered)
saveRDS(ec_norm_imp_log2_filtered, "clean/6_ferment_ec_long_norm_imp_log2_filtered.RDS")

# Pea protein function

# create tsv for mapping
pea_ids <- proteins_filtered %>%
  filter(origin %in% c("pea")) %>%
  dplyr::select(proteinid) %>%
  distinct() %>%
  mutate(proteinid = str_extract(proteinid, "[A-Za-z0-9]+_[A-Z]+$"))

write_tsv(pea_ids, file = here(folder, "pea_ids.tsv"))

# prepare feed functions

pea_func <- read_tsv(here(folder, "idmapping_2026_07_31.tsv")) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  dplyr::select(-from) %>%
  dplyr::rename(protein = entry, kegg_ko = kegg, gos = gene_ontology_ids) %>%
  mutate(kegg_ko = str_remove_all(kegg_ko, "ssc:|;$|\\s"),
         gos = str_remove_all(gos, "$;|\\s"),
         kegg_ko = str_replace_all(kegg_ko, ";", ","),
         gos = str_replace_all(gos, ";", ","),
         kegg_ko = ifelse(is.na(kegg_ko), "-", kegg_ko),
         gos = ifelse(is.na(gos), "-", gos)) # prepare for go and kegg functions

saveRDS(pea_func, "clean/6_ferment_pea_functions_long_raw.RDS")
