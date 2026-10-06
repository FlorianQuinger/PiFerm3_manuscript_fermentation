library(here)
library(venn)

source(here("50_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#save switch
#save = TRUE
save = FALSE

# load meta

meta <- readRDS("clean/meta1.RDS")

ec_table <- read_tsv("data/ec_meta_new.txt") %>%
  mutate(Substrate = ifelse(str_detect(Substrate, "galactosides"), "alpha-galactosides", Substrate)) %>%
  dplyr::rename(EC = 'Enzyme Commission number', substrate = Substrate)

# load summary file

summary_long <- readRDS("clean/6_ferment_summary_long.RDS") %>%
  filter(!parameter == "ms/ms") %>%
  left_join(dplyr::select(meta, sampleid, matrix, animal, period, diet), by = "sampleid") %>%
  mutate(parameter = case_when(parameter == "ms/ms_identified" ~ "MS identified",
                               parameter == "ms/ms_identified_[%]" ~ "MS identified %",
                               parameter == "peptide_sequences_identified" ~ "Peptides"))

# load raw protein file 

proteins <- readRDS("clean/6_ferment_proteins_long_raw_intensity.RDS")

proteins_rel_abd <- readRDS("clean/6_ferment_proteins_long_raw_rel_abd.RDS")

proteins_norm_imp_rel_abd <- readRDS("clean/6_ferment_proteins_long_norm_imp_rel_abd_filtered.RDS")

proteins_norm_imp_intensity <- readRDS("clean/6_ferment_proteins_long_norm_imp_intensity_filtered.RDS")

taxa_raw <- readRDS("clean/6_ferment_taxa_long_raw_intensity.RDS")

taxa_norm_imp_rel_abd <- readRDS("clean/6_ferment_taxa_long_norm_imp_rel_abd_filtered.RDS") 

taxa_norm_imp_intensity <- readRDS("clean/6_ferment_taxa_long_norm_imp_intensity_filtered.RDS")

taxonomy_norm_imp_rel_abd <- readRDS("clean/6_ferment_taxonomy_long_norm_imp_rel_abd_filtered.RDS")

taxonomy_norm_imp_intensity <- readRDS("clean/6_ferment_taxonomy_long_norm_imp_intensity_filtered.RDS")

kegg_norm_imp_rel_abd <- readRDS("clean/6_ferment_kegg_long_norm_imp_rel_abd_filtered.RDS")

kegg_norm_imp_intensity <- readRDS("clean/6_ferment_kegg_long_norm_imp_intensity_filtered.RDS")

cog_norm_imp_intensity <- readRDS("clean/6_ferment_cog_long_norm_imp_intensity_filtered.RDS")

cazy_norm_imp_intensity <- readRDS("clean/6_ferment_cazy_long_norm_imp_intensity_filtered.RDS")

ec_norm_imp_intensity <- readRDS("clean/6_ferment_ec_long_norm_imp_intensity_filtered.RDS")

ec_norm_imp_rel_abd <- readRDS("clean/6_ferment_ec_long_norm_imp_rel_abd_filtered.RDS")

# objects without sample 409

proteins_norm_imp_rel_abd_filtered <- filter(proteins_norm_imp_rel_abd, sampleid != "409")

proteins_norm_imp_intensity_filtered <- filter(proteins_norm_imp_intensity, sampleid != "409")

taxa_raw_filtered <- filter(taxa_raw, sampleid != "409")

taxa_norm_imp_rel_abd_filtered <- filter(taxa_norm_imp_rel_abd, sampleid != "409")

taxa_norm_imp_intensity_filtered <- filter(taxa_norm_imp_intensity, sampleid != "409")

taxonomy_norm_imp_rel_abd_filtered <- filter(taxonomy_norm_imp_rel_abd, sampleid != "409")

taxonomy_norm_imp_intensity_filtered <- filter(taxonomy_norm_imp_intensity, sampleid != "409")

kegg_norm_imp_rel_abd_filtered <- filter(kegg_norm_imp_rel_abd, sampleid != "409")

kegg_norm_imp_intensity_filtered <- filter(kegg_norm_imp_intensity, sampleid != "409")

cog_norm_imp_intensity_filtered <- filter(cog_norm_imp_intensity, sampleid != "409")

cazy_norm_imp_intensity_filtered <- filter(cazy_norm_imp_intensity, sampleid != "409")

ec_norm_imp_intensity_filtered <- filter(ec_norm_imp_intensity, sampleid != "409")

ec_norm_imp_rel_abd_filtered <- filter(ec_norm_imp_rel_abd, sampleid != "409")

# functions file for joining

functions <- readRDS("clean/6_ferment_functions_long_raw.RDS")

pea_functions <- readRDS("clean/6_ferment_pea_functions_long_raw.RDS")

# filter microbial proteins

proteins_micro_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

proteins_micro_rel_abd_filtered <- proteins_norm_imp_rel_abd_filtered %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

proteins_pea_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin == "pea") %>%
  recalculate_rel_abd()

proteins_pea_rel_abd_filtered <- proteins_norm_imp_rel_abd_filtered %>%
  filter(origin == "pea") %>%
  recalculate_rel_abd()

proteins_micro_intensity <- proteins_norm_imp_intensity %>%
  filter(origin == "micro")

proteins_micro_intensity_filtered <- proteins_norm_imp_intensity_filtered %>%
  filter(origin == "micro") 

proteins_pea_intensity <- proteins_norm_imp_intensity %>%
  filter(origin == "pea") 

proteins_pea_intensity_filtered <- proteins_norm_imp_intensity_filtered %>%
  filter(origin == "pea") 

functions_filtered <- functions %>%
  filter(proteinid %in% proteins_micro_intensity_filtered$proteinid)

# plot

ggplot(summary_long, aes(x = sampleid, y = value, fill = matrix)) +
  geom_bar(stat = "identity", show.legend = F) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(rows = vars(parameter), scales = "free_y") +
  labs(y = "") +
  scale_fill_manual(values = c("grey50", "grey70")) + 
  theme(axis.text.x = element_text(size = 8))

save_big("62_ms_summary_bar")

# plot as boxplots

ggplot(summary_long, aes(x = diet, y = value)) +
  geom_violin(draw_quantiles = c(0.5), fill = "snow2") +
  geom_quasirandom(aes(shape = diet, color = period), width = 0.3, size = 2) +
  #geom_jitter(aes(shape = period, color = diet), width = 0.3, size = 2) +
  facet_wrap(~parameter, scales = "free_y") +
  scale_color_manual(values = colors) +
  labs(x = "", y = "") + 
  scale_x_discrete(guide = guide_axis(angle = 30))

save_big("62_ms_summary_violin")

# plot number of proteingroups per sample

proteingroups_per_sample <- proteins %>%
  mutate(abu = ifelse(intensity > 0, 1, 0)) %>%
  group_by(sampleid) %>%
  summarize(abu = sum(abu)) %>%
  left_join(meta, by = "sampleid")

ggplot(proteingroups_per_sample, aes(x = diet, y = abu)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.3, size = 2) +
  labs(x = "matrix", y = "number of protein groups", title = "Number of protein groups quantified")

save_big("62_proteingroups_boxplot")

# plot number of proteingroups of different origin per sample

proteingroups_per_sample_origin <- proteins %>%
  mutate(abu = ifelse(intensity > 0, 1, 0)) %>%
  group_by(sampleid, origin) %>%
  summarize(abu = sum(abu)) %>%
  left_join(meta, by = "sampleid")

ggplot(proteingroups_per_sample_origin, aes(x = sampleid, y = abu, fill = origin)) +
  geom_bar(stat = "identity") +
  labs(y = "number of protein groups", title = "Number of protein groups quantified") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 2700), expand = c(0,0))

save_big("62_proteingroups_origin")

proteingroups_origin <- proteins_norm_imp_intensity %>%
  group_by(proteinid, origin) %>%
  summarize(intensity = sum(intensity), .groups = "drop") %>%
  mutate(abu = ifelse(intensity > 0, 1, 0)) %>%
  group_by(origin) %>%
  summarise(abu = sum(abu), .groups = "drop")

# plot pea - micro ratio

proteins_origin <- proteins %>%
  filter(origin %in% c("pea", "micro")) %>% # only host and microbial proteins
  group_by(sampleid, origin) %>%
  summarize(intensity = sum(intensity), .groups = "drop") %>%
  group_by(sampleid) %>%
  mutate(rel_abd = intensity/sum(intensity)) %>%
  ungroup() %>%
  left_join(meta, by = "sampleid")


proteins_origin %>%
  group_by(diet, origin) %>%
  summarize(rel_abd = mean(rel_abd), .groups = "drop") %>%
  ggplot(aes(x = diet, y = rel_abd, fill = origin)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_boxplot(data = filter(proteins_origin, origin == "pea"), 
               aes(x = diet, y = rel_abd), 
               alpha = 0, width = 0.5, lwd = 1,
               inherit.aes = F) + 
  geom_quasirandom(data = filter(proteins_origin, origin == "pea"), 
                   aes(x = diet, y = rel_abd),
                   width = 0.2, alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "relative abundance", title = "Proportion of microbial and host proteins")

save_big("62_micro_pea_ratio_bar")

proteins_origin %>%
  filter(sampleid != "409") %>%
  filter(origin == "micro") %>%
  pull(rel_abd) %>%
  mean()

# taxonomic annotation of proteingroups

proteingroups_taxonomy <- proteins_norm_imp_intensity_filtered %>%
  filter(str_detect(proteinid, "MGYG")) %>%
  mutate(bin = str_extract(proteinid, "[A-Za-z0-9]+")) %>%
  left_join(taxa_raw %>% dplyr::select(bin, G) %>% distinct(), by = "bin") %>%
  dplyr::select(-sampleid, -intensity) %>%
  distinct() %>%
  group_by(G) %>%
  summarise(n = length(proteinid), .groups = "drop") # 128 Companilactobacillus proteins / 2271 genes of Companilactobacillus kimchii = 5.6%

ggplot(proteingroups_taxonomy, aes(x = 1, y = n, fill = G)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors)



####### Ordination

do_ordination(proteins_norm_imp_rel_abd, title = "All proteins ferment", save_name = "proteingroups_all", save = save,
              second_indicator = "period")

do_ordination(proteins_norm_imp_rel_abd_filtered, title = "All proteins ferment", save_name = "proteingroups_all", save = save,
              second_indicator = "period", region = "ferment")

do_ordination(proteins_micro_rel_abd, title = "Microbial proteins ferment", save_name = "proteingroups_micro", save = save,
              second_indicator = "period")

do_ordination(proteins_micro_rel_abd_filtered, title = "Microbial proteins ferment", save_name = "proteingroups_micro", save = save,
              second_indicator = "period", region = "ferment")

do_ordination(taxa_norm_imp_rel_abd, title = "Microbial taxa ferment", save_name = "taxa", save = save,
              second_indicator = "period")

do_ordination(taxa_norm_imp_rel_abd_filtered, title = "Microbial taxa ferment", save_name = "taxa", save = save,
              second_indicator = "period", region = "ferment")

do_ordination(kegg_norm_imp_rel_abd, title = "Microbial functions ferment", save_name = "function_kegg", save = save,
              second_indicator = "period")

do_ordination(kegg_norm_imp_rel_abd_filtered, title = "Microbial functions ferment", save_name = "function_kegg", save = save,
              second_indicator = "period", region = "ferment")

do_ordination(proteins_pea_rel_abd, title = "Pea proteins ferment", save_name = "proteingroups_pea", save = save,
              second_indicator = "period")

do_ordination(proteins_pea_rel_abd_filtered, title = "Pea proteins ferment", save_name = "proteingroups_pea", save = save,
              second_indicator = "period", region = "ferment")


############ Taxa barplot

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "P", selected_matrix = "ferment", 
                         title = "Phyla ferment", save_name = "62_taxa_barplot_ferment_phylum", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "C", selected_matrix = "ferment", 
                         title = "Class ferment", save_name = "62_taxa_barplot_ferment_class", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "O", selected_matrix = "ferment", 
                         title = "Order ferment", save_name = "62_taxa_barplot_ferment_order", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "F", selected_matrix = "ferment", 
                         title = "Family ferment", save_name = "62_taxa_barplot_ferment_family", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "ferment", 
                         title = "Genus ferment", save_name = "62_taxa_barplot_ferment_genus", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "ferment", 
                         title = "Genus ferment", save_name = "62_taxa_barplot_ferment_genus", save = save, threshold = 0)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "ferment", 
                         title = "Species ferment", save_name = "62_taxa_barplot_ferment_species", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S", selected_matrix = "ferment", 
                         title = "Species ferment",save = F)

taxa_list_400 <- taxa_norm_imp_rel_abd %>%
  filter(sampleid != "409") %>%
  group_by(bin, R1, P, C, O, `F`, G, S) %>%
  summarize(rel_abd = mean(rel_abd))

taxa_list_400_G <- taxa_norm_imp_rel_abd_filtered %>%
  group_by(R1, P, C, O, `F`, G) %>%
  summarize(rel_abd = mean(rel_abd))

taxa_list_409 <- taxa_norm_imp_rel_abd %>%
  filter(sampleid == "409") %>%
  group_by(bin, R1, P, C, O, `F`, G, S) %>%
  summarize(rel_abd = mean(rel_abd))

######### Barplots for proteins

taxa_barplot_from_ktable(dplyr::select(proteins_micro_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ferment", title = "Microbial proteins ferment",
                         save_name = "62_barplot_micro_proteins_ferment", save = save)

taxa_barplot_from_ktable(dplyr::select(proteins_pea_rel_abd_filtered, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ferment", title = "Pea proteins ferment",
                         save_name = "62_barplot_pea_proteins_ferment", save = save)

taxa_barplot_from_ktable(dplyr::select(kegg_norm_imp_rel_abd, name = kegg_ko, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ferment", title = "Functions KEGG ferment",
                         save_name = "62_barplot_functions_kegg_ferment", save = save)

taxa_barplot_from_ktable(dplyr::select(ec_norm_imp_rel_abd_filtered, name = ec_id2, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ferment", title = "Functions EC ferment",
                         save_name = "62_barplot_functions_ec_ferment", save = save)

# Alpha diversity

## observed species

observed_species <- taxa_raw_filtered %>%
  mutate(rel_abd = ifelse(intensity > 0, 1, 0)) %>% # set to absence presence
  group_by(sampleid) %>%
  summarise(observed = sum(rel_abd))

## shannon

#matrix <- taxa_raw_rel_abd %>%
matrix <- taxa_norm_imp_rel_abd_filtered %>% # differences by using imputed file
  dplyr::select(bin, sampleid, rel_abd) %>%
  pivot_wider(names_from = "bin", values_from = "rel_abd") %>%
  column_to_rownames("sampleid") %>%
  mutate_all(~ifelse(is.na(.), 0, .))

shannon <- tibble(sampleid = rownames(matrix),
                  shannon = vegan::diversity(matrix, index = "shannon"),
                  simpson = vegan::diversity(matrix, index = "invsimpson"))

# plot

alpha <- observed_species %>%
  inner_join(shannon, by = "sampleid") %>%
  inner_join(meta, by = "sampleid")

ggplot(pivot_longer(alpha, c(shannon, observed, simpson), names_to = "index", values_to = "value"), aes(x = diet, y = value)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  facet_wrap(vars(index), scales = "free_y")

#comparison

comp_shannon <- combined_comparison(select_response(alpha, "shannon"), transformation = "test", 
                                    model_style = "ferment") 
# comp_obs <- combined_comparison(select_response(alpha, "observed"), transformation = "test", 
#                                    model_style = "ferment")
comp_simpson <- combined_comparison(select_response(alpha, "simpson"), transformation = "test", 
                                       model_style = "ferment")

#table_obs <- create_results_table(select_response(alpha, "observed"), comp_obs, response = "observed_il")
table_shannon <- create_results_table(select_response(alpha, "shannon"), comp_shannon, response = "shannon_il")
table_simpson <- create_results_table(select_response(alpha, "simpson"), comp_simpson, response = "simpson_il")

table_alpha <- table_simpson %>%
  inner_join(table_shannon, by = "diet")

write_tsv(table_alpha, "tables/62_alpha_diversity.txt")

# differential abundance analysis

ancom_P <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "P",
                                       selected_matrix = "ferment")

ancom_C <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "C",
                                       selected_matrix = "ferment")

ancom_O <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "O",
                                       selected_matrix = "ferment")

ancom_F <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "F",
                                       selected_matrix = "ferment")

ancom_G <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "G",
                                       selected_matrix = "ferment")

ancom_S <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "S",
                                       selected_matrix = "ferment")

ancom_low <- perform_ancombc_and_plot(ktable = taxa_norm_imp_intensity_filtered, meta = meta, selected_rank = "low",
                                         selected_matrix = "ferment")

# DEA for proteins

proteins_results <- perform_edger_and_plot(proteins_norm_imp_intensity_filtered, selected_matrix = "ferment", save = save,
                                                 abundance_column = "intensity", save_name = "proteingroups_all")

proteins_micro_results <- perform_edger_and_plot(proteins_micro_intensity_filtered, selected_matrix = "ferment", save = save,
                                                    abundance_column = "intensity", save_name = "proteingroups_micro")

proteins_pea_results <- perform_edger_and_plot(proteins_pea_intensity_filtered, selected_matrix = "ferment", save = save,
                                                 abundance_column = "intensity", save_name = "proteingroups_pea")

kegg_micro_results <- perform_edger_and_plot(kegg_norm_imp_intensity_filtered, selected_matrix = "ferment", save = save,
                                               abundance_column = "intensity", save_name = "kegg_micro")

cog_micro_results <- perform_edger_and_plot(cog_norm_imp_intensity_filtered, selected_matrix = "ferment", save = save,
                                             abundance_column = "intensity", save_name = "cog_micro")

cazy_micro_results <- perform_edger_and_plot(cazy_norm_imp_intensity_filtered, selected_matrix = "ferment", save = save,
                                             abundance_column = "intensity", save_name = "cazy_micro")

proteins_micro_rel_abd_filtered %>%
  left_join(meta, by = "sampleid") %>%
  filter(proteinid == "MGYGMiFoDB36_01138") %>%
  ggplot(aes(x = diet, y = rel_abd)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  labs(title = "Limosilactobacillus reuteri_I pyruvate dehydrogenase")

proteins_micro_rel_abd_filtered %>%
  left_join(meta, by = "sampleid") %>%
  filter(proteinid == "MGYGMiFoDB56_01726") %>%
  ggplot(aes(x = diet, y = rel_abd)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  labs(title = "Limosilactobacillus vaginalis \n xylulose-5-phosphate phosphoketolase")

# Pea proteins

functions_pea_rel_abd_filtered <- proteins_pea_rel_abd_filtered %>%
  group_by(proteinid) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  mutate(entry_name = str_remove(proteinid, "^[a-z]+\\|[A-Z0-9]+\\|")) %>%
  inner_join(pea_functions, by = "entry_name") %>%
  group_by(protein_names) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  add_column(sampleid = "401") %>%
  dplyr::rename(name = protein_names) # alpha galactosidase present

taxa_barplot_from_ktable(functions_pea_rel_abd_filtered, meta = meta, selected_rank = "low", selected_matrix = "ferment",
                         save = F)

# Protein function lists

proteins_micro_rel_abd_400 <- proteins_micro_rel_abd %>%
  filter(sampleid != "409") %>%
  group_by(proteinid) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  left_join(functions, by = "proteinid")

library(pathview)

nmr <- readRDS("clean/7_nmr_ferment_fm.RDS") %>% 
  pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
  filter(sampleid != "409") %>%
  group_by(metabolite) %>%
  summarise(concentration = mean(concentration), .groups = "drop") %>%
  left_join(read_tsv("data/nmr_metabolite_meta.txt"), by = "metabolite") %>%
  filter(!is.na(kegg_compound)) %>%
  dplyr::select(metabolite, kegg_compound)

nutrition_metabolites <- read_tsv("data/nutrition_metabolite_meta.txt")

compounds <- rbind(nmr, nutrition_metabolites)

# C. kimchii

Ckimchi_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGPF11641")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

taxa_barplot_from_ktable(dplyr::select(Ckimchi_proteins_400, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "")

save_big("62_barplot_ckimchi_proteins")

Ckimchi_cog_400 <- Ckimchi_proteins_400 %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Ckimchi_cog_400, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_ckimchi_cog_cat")

Ckimchi_kegg_400 <- Ckimchi_proteins_400 %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- annotate_keggs(Ckimchi_kegg_400$kegg_ko)

Ckimchi_kegg_400 <- Ckimchi_kegg_400 %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Ckimchi_kegg_400, name = kegg_ko, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

filter(Ckimchi_kegg_400, str_detect(pathways, "Glycolysis / Gluconeogenesis"))
cat(filter(Ckimchi_kegg_400, str_detect(pathways, "Fructose and mannose metabolism"))$kegg_ko)
filter(Ckimchi_kegg_400, str_detect(pathways, "Pyruvate metabolism"))
filter(Ckimchi_kegg_400, str_detect(pathways, "Two-component system"))
filter(Ckimchi_kegg_400, str_detect(pathways, "Starch and sucrose metabolism"))
cat(filter(Ckimchi_kegg_400, str_detect(pathways, "Pentose phosphate pathway"))$kegg_ko)
filter(Ckimchi_kegg_400, str_detect(pathways, "Galactose metabolism"))
filter(Ckimchi_kegg_400, str_detect(pathways, "ABC transporters"))

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00010", species = "ko", 
         kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
         plot.col.key = F, out.suffix = "Ckimchi_glycolysis")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00051", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Ckimchi_mannose")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00052", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Ckimchi_galactose")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "02010", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Ckimchi_transporters")

# Weissella

Wconfusa_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGPF1022")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

taxa_barplot_from_ktable(dplyr::select(Wconfusa_proteins_400, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = 1)  +
  labs(fill = "protein", x = "")

save_big("62_barplot_wconfusa_proteins", width = 25)

Wconfusa_cog_400 <- Wconfusa_proteins_400 %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Wconfusa_cog_400, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_wconfusa_cog_cat")

Wconfusa_kegg_400 <- Wconfusa_proteins_400 %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- annotate_keggs(Wconfusa_kegg_400$kegg_ko)

Wconfusa_kegg_400 <- Wconfusa_kegg_400 %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Wconfusa_kegg_400, name = kegg_ko, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

filter(Wconfusa_kegg_400, str_detect(pathways, "Arginine biosynthesis"))$kegg_ko
filter(Wconfusa_kegg_400, str_detect(pathways, "Glycolysis / Gluconeogenesis"))$kegg_ko %>% cat()
filter(Wconfusa_kegg_400, str_detect(pathways, "ABC transporter"))$kegg_ko
filter(Wconfusa_kegg_400, str_detect(pathways, "Pentose phosphate pathway"))$kegg_ko %>% cat()
filter(Wconfusa_kegg_400, str_detect(pathways, "Butanoate"))$kegg_ko
filter(Wconfusa_kegg_400, str_detect(pathways, "Starch and sucrose metabolism"))$kegg_ko
filter(Wconfusa_kegg_400, str_detect(pathways, "Pyruvate metabolism"))$kegg_ko %>% cat

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00010", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Wconfusa_glycolysis")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00220", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Wconfusa_arginine")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00030", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Wconfusa_ppp")

pv <- pathview(gene.data = keggs$kegg_ko, cpd.data = compounds$kegg_compound, pathway.id = "00620", species = "ko", 
               kegg.dir = "temp/", gene.idtype = "KEGG", cpd.idtype = "kegg", kegg.native = T, res = 1000, cex = 10e-10,
               plot.col.key = F, out.suffix = "Wconfusa_pyruvate")


# Limosilactobacillus vaginalis

Lvaginalis_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGMiFoDB56")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

taxa_barplot_from_ktable(dplyr::select(Lvaginalis_proteins_400, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = 2)

save_big("62_barplot_lvaginalis_proteins")

Lvaginalis_kegg_400 <- Lvaginalis_proteins_400 %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- annotate_keggs(Lvaginalis_kegg_400$kegg_ko)

Lvaginalis_kegg_400 <- Lvaginalis_kegg_400 %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Lvaginalis_kegg_400, name = kegg_ko, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

# Kosakonia cowanii

Kcowanii_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGMiFoDB1")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

taxa_barplot_from_ktable(dplyr::select(Kcowanii_proteins_400, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = 2)

Kcowanii_kegg_400 <- Kcowanii_proteins_400 %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- annotate_keggs(Kcowanii_kegg_400$kegg_ko)

Kcowanii_kegg_400 <- Kcowanii_kegg_400 %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Kcowanii_kegg_400, name = kegg_ko, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

# for 409

proteins_micro_rel_abd_409 <- proteins_micro_rel_abd %>%
  filter(sampleid == "409") %>%
  group_by(proteinid) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  left_join(functions, by = "proteinid")

# Companilactobacillus kimchii
# does the function of C. kimchii differ in 409?

Ckimchi_proteins <- proteins_micro_rel_abd %>%
  filter(str_detect(proteinid, "MGYGPF11641")) %>%
  recalculate_rel_abd() %>%
  left_join(functions, by = "proteinid")

taxa_barplot_from_ktable(dplyr::select(Ckimchi_proteins, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "") # no real difference

save_big("62_barplot_ckimchi_proteins_comparison")

Ckimchi_cog <- Ckimchi_proteins %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Ckimchi_cog, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_ckimchi_cog_cat_comparison")

Ckimchi_kegg <- Ckimchi_proteins %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- Ckimchi_kegg %>% 
  select(kegg_ko) %>% 
  distinct() %>%
  pull(kegg_ko) %>%
  annotate_keggs()

Ckimchi_kegg <- Ckimchi_kegg %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Ckimchi_kegg, name = name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

# Enterobacter mori

Emori_proteins <- proteins_micro_rel_abd %>%
  filter(str_detect(proteinid, "MGYGcFMD1006")) %>%
  recalculate_rel_abd() %>%
  left_join(functions, by = "proteinid")

taxa_barplot_from_ktable(dplyr::select(Emori_proteins, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "") # different from 400

save_big("62_barplot_Emori_proteins_comparison")

Emori_cog <- Emori_proteins %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Emori_cog, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_Emori_cog_cat_comparison")

Emori_kegg <- Emori_proteins %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- Emori_kegg %>% 
  select(kegg_ko) %>% 
  distinct() %>%
  pull(kegg_ko) %>%
  annotate_keggs()

Emori_kegg <- Emori_kegg %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Emori_kegg, name = name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

save_big("62_barplot_Emori_kegg_comparison")

Emori_kegg_409 <- filter(Emori_kegg, sampleid == "409")

filter(Emori_kegg_409, str_detect(pathways, "Glycolysis / Gluconeogenesis"))$kegg_ko
filter(Emori_kegg_409, str_detect(pathways, "Pyruvate metabolism"))$kegg_ko

# Kosakonia kowanii

Kkowanii_proteins <- proteins_micro_rel_abd %>%
  filter(str_detect(proteinid, "MGYGMiFoDB1")) %>%
  recalculate_rel_abd() %>%
  left_join(functions, by = "proteinid")

taxa_barplot_from_ktable(dplyr::select(Kkowanii_proteins, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "") # different from 400

save_big("62_barplot_Kkowanii_proteins_comparison")

Kkowanii_cog <- Kkowanii_proteins %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Kkowanii_cog, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_Kkowanii_cog_cat_comparison")

Kkowanii_kegg <- Kkowanii_proteins %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- Kkowanii_kegg %>% 
  select(kegg_ko) %>% 
  distinct() %>%
  pull(kegg_ko) %>%
  annotate_keggs()

Kkowanii_kegg <- Kkowanii_kegg %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Kkowanii_kegg, name = name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

save_big("62_barplot_Kkowanii_kegg_comparison")

Kkowanii_kegg_409 <- filter(Kkowanii_kegg, sampleid == "409")

filter(Kkowanii_kegg_409, str_detect(pathways, "Glycolysis / Gluconeogenesis"))$kegg_ko %>% cat()
filter(Kkowanii_kegg_409, str_detect(pathways, "Pyruvate metabolism"))$kegg_ko


# Clostridium beijerinckii

Cbeiji_proteins <- proteins_micro_rel_abd %>%
  filter(str_detect(proteinid, "MGYGMiFoDB3")) %>%
  recalculate_rel_abd() %>%
  left_join(functions, by = "proteinid")

taxa_barplot_from_ktable(dplyr::select(Cbeiji_proteins, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "") # different from 400

save_big("62_barplot_Cbeiji_proteins_comparison")

Cbeiji_cog <- Cbeiji_proteins %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Cbeiji_cog, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_Cbeiji_cog_cat_comparison")

Cbeiji_kegg <- Cbeiji_proteins %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- Cbeiji_kegg %>% 
  select(kegg_ko) %>% 
  distinct() %>%
  pull(kegg_ko) %>%
  annotate_keggs()

Cbeiji_kegg <- Cbeiji_kegg %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Cbeiji_kegg, name = name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

save_big("62_barplot_Cbeiji_kegg_comparison")

Cbeiji_kegg_409 <- filter(Cbeiji_kegg, sampleid == "409")

filter(Cbeiji_kegg_409, str_detect(pathways, "Butanoate metabolism"))$kegg_ko %>% cat()

#Lactococcus petauri

Lpetauri_proteins <- proteins_micro_rel_abd %>%
  filter(str_detect(proteinid, "MGYGcFMD438")) %>%
  recalculate_rel_abd() %>%
  left_join(functions, by = "proteinid")

taxa_barplot_from_ktable(dplyr::select(Lpetauri_proteins, name = protein_name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "protein", x = "") # different from 400

save_big("62_barplot_Lpetauri_proteins_comparison")

Lpetauri_cog <- Lpetauri_proteins %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category")

taxa_barplot_from_ktable(dplyr::select(Lpetauri_cog, name = cat_description, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5) +
  labs(fill = "COG", x = "")

save_big("62_barplot_Lpetauri_cog_cat_comparison")

Lpetauri_kegg <- Lpetauri_proteins %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")%>%
  recalculate_rel_abd()

keggs <- Lpetauri_kegg %>% 
  select(kegg_ko) %>% 
  distinct() %>%
  pull(kegg_ko) %>%
  annotate_keggs()

Lpetauri_kegg <- Lpetauri_kegg %>%
  left_join(keggs, by = "kegg_ko")

taxa_barplot_from_ktable(dplyr::select(Lpetauri_kegg, name = name, rel_abd, sampleid), 
                         meta = meta, selected_rank = "low", selected_matrix = "ferment", threshold = .5)

save_big("62_barplot_Lpetauri_kegg_comparison")

Lpetauri_kegg_409 <- filter(Lpetauri_kegg, sampleid == "409")

