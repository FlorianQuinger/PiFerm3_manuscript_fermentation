library(here)

source(here("50_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#save switch
#save = TRUE
save = FALSE

meta <- readRDS("clean/meta1.RDS")

# load data

rarefied_counts <- readRDS("clean/5_16s_ferment_taxa_rarefied_raw_counts.RDS")

rarefied_rel_abd <- readRDS("clean/5_16s_ferment_taxa_rarefied_raw_rel_abd.RDS")
rarefied_rel_abd_filtered <- readRDS("clean/5_16s_ferment_taxa_rarefied_raw_rel_abd.RDS") 

rarefied_ktable_rel_abd <- readRDS("clean/5_16s_ferment_taxonomy_rarefied_raw_rel_abd.RDS") 
rarefied_ktable_rel_abd_filtered <- readRDS("clean/5_16s_ferment_taxonomy_rarefied_filtered_rel_abd.RDS") 

feature_counts <- readRDS("clean/5_16s_ferment_taxa_feature_raw_counts.RDS") 
feature_counts_filtered <- readRDS("clean/5_16s_ferment_taxa_feature_filtered_counts.RDS") 

feature_rel_abd <- readRDS("clean/5_16s_ferment_taxa_feature_raw_rel_abd.RDS")
feature_rel_abd_filtered <- readRDS("clean/5_16s_ferment_taxa_feature_filtered_rel_abd.RDS")

feature_ktable_counts <- readRDS("clean/5_16s_ferment_taxonomy_feature_raw_counts.RDS")
feature_ktable_counts_filtered <- readRDS("clean/5_16s_ferment_taxonomy_feature_filtered_counts.RDS")

feature_ktable_rel_abd <- readRDS("clean/5_16s_ferment_taxonomy_feature_raw_rel_abd.RDS")
feature_ktable_rel_abd_filtered <- readRDS("clean/5_16s_ferment_taxonomy_feature_filtered_rel_abd.RDS")

taxonomy <- readRDS("clean/5_16s_ferment_asv_taxonomy.RDS")


#######################
# Alpha diversity
#######################

# qiime2 files

alpha_observed <- read_tsv("data/16S_ferment/alpha-diversity_observed.tsv")
alpha_evenness <- read_tsv("data/16S_ferment/alpha-diversity_evenness.tsv")
alpha_shannon <- read_tsv("data/16S_ferment/alpha-diversity_shannon.tsv")
alpha_faith <- read_tsv("data/16S_ferment/alpha-diversity_faith_pd.tsv")

alpha_q <- alpha_observed %>% 
  inner_join(alpha_evenness, by = "...1") %>%
  inner_join(alpha_shannon, by = "...1") %>%
  inner_join(alpha_faith, by = "...1") %>%
  dplyr::rename(observed = observed_features, evenness = pielou_evenness, 
                shannon = shannon_entropy, faith = faith_pd, sampleid = `...1`) %>%
  mutate(sampleid = str_remove(sampleid, "IS-Pig-")) %>%
  inner_join(meta, by = "sampleid")

ggplot(pivot_longer(alpha_q, c(observed, evenness, shannon, faith), names_to = "index", values_to = "value"), 
       aes(x = diet, y = value)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  facet_wrap(vars(matrix, index), scales = "free_y")

# comparison 

comp_obs_fe <- combined_comparison_ferment(select_response(alpha_q, "observed"), transformation = "test")
comp_evenness_fe <- combined_comparison_ferment(select_response(alpha_q, "evenness"), transformation = "test")
comp_shannon_fe <- combined_comparison_ferment(select_response(alpha_q, "shannon"), transformation = "test")
comp_faith_fe <- combined_comparison_ferment(select_response(alpha_q, "faith"), transformation = "test")


table_obs_fe <- create_results_table(select_response(alpha_q, "observed"), comp_obs_fe, response = "observed_fe")
table_evenness_fe <- create_results_table(select_response(alpha_q, "evenness"), comp_evenness_fe, response = "evenness_fe")
table_shannon_fe <- create_results_table(select_response(alpha_q, "shannon"), comp_shannon_fe, response = "shannon_fe")
table_faith_fe <- create_results_table(select_response(alpha_q, "faith"), comp_faith_fe, response = "faith_fe")


table_alpha_q <- table_obs_fe %>%
  left_join(table_evenness_fe, by = "diet") %>%
  left_join(table_shannon_fe, by = "diet") %>%
  left_join(table_faith_fe, by = "diet") 

write_tsv(table_alpha_q, "tables/57_16S_ferment_alpha_div.txt")


########################
# beta diversity
########################

# Bray curtis from qiime2

bray <- read_tsv("data/16S_ferment/distance-matrix.tsv") %>%
  mutate(`...1` = str_remove(`...1`, "IS-Pig-")) %>%
  rename_all(., ~gsub("IS-Pig-", "",.)) %>%
  column_to_rownames("...1") 

bray_fe <- bray %>%
  as.dist()

do_ordination(bray_fe, region = "ferment", title = "Taxonomy 16S ferment", save_name = "57_pcoa_ferment_16s", save = save)

# Beta diversity with rarefied rel abundances

do_ordination(rarefied_rel_abd, region = "ferment", title = "Taxonomy 16S ferment", save_name = "57_pcoa_ferment_16s", save = save)


###############################
## Taxa barplots
##############################

# ferment

taxa_barplot_from_ktable(rarefied_ktable_rel_abd_filtered, meta = meta, selected_rank = "P", selected_matrix = "ferment", 
                         title = "Phyla ferment", save_name = "57_taxa_barplot_ferment_16s_phylum", save = save)

taxa_barplot_from_ktable(rarefied_ktable_rel_abd_filtered, meta = meta, selected_rank = "C", selected_matrix = "ferment", 
                         title = "Class ferment", save_name = "57_taxa_barplot_ferment_16s_class", save = save)

taxa_barplot_from_ktable(rarefied_ktable_rel_abd_filtered, meta = meta, selected_rank = "O", selected_matrix = "ferment", 
                         title = "Orders ferment", save_name = "57_taxa_barplot_ferment_16s_order", save = save)

taxa_barplot_from_ktable(rarefied_ktable_rel_abd_filtered, meta = meta, selected_rank = "F", selected_matrix = "ferment", 
                         title = "Families ferment", save_name = "57_taxa_barplot_ferment_16s_family", save = save)

taxa_barplot_from_ktable(rarefied_ktable_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "ferment", 
                         title = "Genera ferment", save_name = "57_taxa_barplot_ferment_16s_genus", save = save)

# list all genera in ferment

genera_ferment <- rarefied_ktable_rel_abd_filtered %>%
  filter_ferment() %>%
  filter(rank == "G") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd))

asv_ferment <- rarefied_rel_abd %>%
  filter_ferment() %>%
  group_by(asv) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  filter(rel_abd > 0) %>%
  arrange(rel_abd) %>%
  rename(ferment = rel_abd)


###########################
# differential abundance
###########################

# ferment

fe_ancom_P <- perform_ancombc_and_plot(ktable = feature_ktable_counts_filtered, meta = meta, selected_rank = "P", 
                                       selected_matrix = "ferment")
fe_ancom_C <- perform_ancombc_and_plot(ktable = feature_ktable_counts_filtered, meta = meta, selected_rank = "C", 
                                       selected_matrix = "ferment")
fe_ancom_O <- perform_ancombc_and_plot(ktable = feature_ktable_counts_filtered, meta = meta, selected_rank = "O", 
                                       selected_matrix = "ferment")
fe_ancom_F <- perform_ancombc_and_plot(ktable = feature_ktable_counts_filtered, meta = meta, selected_rank = "F", 
                                       selected_matrix = "ferment")
fe_ancom_G <- perform_ancombc_and_plot(ktable = feature_ktable_counts_filtered, meta = meta, selected_rank = "G", 
                                       selected_matrix = "ferment")

fe_ancom_low <- perform_ancombc_and_plot(ktable = feature_counts_filtered, meta = meta, selected_rank = "low", 
                                         selected_matrix = "ferment")

