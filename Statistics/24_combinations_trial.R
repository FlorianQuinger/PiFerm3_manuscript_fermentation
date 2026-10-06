library(here)

source("0_general_functions.R")
source("20_pretrial_functions.R")

save = FALSE
#save = TRUE

# load data
combinations <- read_tsv("data/24_combinations.txt") %>%
  mutate(combination = as.character(combination), 
         Bezeichnung = factor(Bezeichnung, levels= c("blank", "A10", "P4", "P2", "A10,P4", "A10,P2", "P4,P2", "A10,P4,P2")))

meta <- read_tsv("data/24_meta.txt") %>%
  mutate(block = as.character(block),
         box = as.character(box),
         sampleid = as.character(sampleid),
         combination = as.character(combination)) %>%
  inner_join(dplyr::select(combinations, combination, Bezeichnung), by = "combination") %>%
  dplyr::rename(bacteria = "Bezeichnung", species = combination) %>%
  mutate(enzyme = ifelse(enzyme == "none", "none", "enzyme"))

ph <- read_tsv("data/24_ph.txt") %>%
  mutate(sampleid = as.character(sampleid))

ammonia <- read_tsv("data/24_ammonia.txt") %>%
  mutate(sampleid = as.character(sampleid), 
         mg_nh4n_100g = as.numeric(mg_nh4n_100g)) %>%
  filter(!is.na(mg_nh4n_100g))

acids <- read_tsv("data/24_acids.txt") %>%
  mutate(sampleid = as.character(sampleid)) %>%
  filter(lactate > 0) %>%
  mutate(lactate = ifelse(sampleid == "137", NA, lactate))

# combine with meta

ph_combined <- ph %>%
  inner_join(meta, by = "sampleid")

ammonia_combined <- ammonia %>%
  inner_join(meta, by = "sampleid")

acids_combined <- acids %>%
  inner_join(meta, by = "sampleid") 


# comparison of ph at different timepoints one factorial

ph0 <- combined_comparison_and_results(df = ph_combined, selected_response = "0", transformation = "none",
                                       type = "h", response_name = "pH 0", y_axis = "pH", save_name = "ph_0_pea",
                                       save = save, digits = 2)

ph2 <- combined_comparison_and_results(df = ph_combined, selected_response = "2", transformation = "none",
                                       type = "h", response_name = "pH 2", y_axis = "pH", save_name = "ph_2_pea",
                                       save = save, digits = 2)

ph4 <- combined_comparison_and_results(df = ph_combined, selected_response = "4", transformation = "none",
                                       type = "h", response_name = "pH 4", y_axis = "pH", save_name = "ph_4_pea",
                                       save = save, digits = 2)

ph6 <- combined_comparison_and_results(df = ph_combined, selected_response = "6", transformation = "none",
                                       type = "h", response_name = "pH 6", y_axis = "pH", save_name = "ph_6_pea",
                                       save = save, digits = 2)

ph8 <- combined_comparison_and_results(df = ph_combined, selected_response = "8", transformation = "none",
                                       type = "h", response_name = "pH 8", y_axis = "pH", save_name = "ph_8_pea",
                                       save = save, digits = 2)

ph12 <- combined_comparison_and_results(df = ph_combined, selected_response = "12", transformation = "test",
                                       type = "h", response_name = "pH 12", y_axis = "pH", save_name = "ph_12_pea",
                                       save = save, digits = 2)

ph24 <- combined_comparison_and_results(df = ph_combined, selected_response = "24", transformation = "test",
                                        type = "h", response_name = "pH 24", y_axis = "pH", save_name = "ph_24_pea",
                                        save = save, digits = 2)

# combine results

ph_table_pea <- cbind(ph0, "pH 2" = ph2[,2], "pH 4" = ph4[,2], "pH 6" = ph6[,2], "pH 8" = ph8[,2],
                      "pH 12" = ph12[,2], "pH 24" = ph24[,2])

write_tsv(ph_table_pea, "tables/24_combinations_ph.tsv")
save_results_table(ph_table_pea, title = "pH combinations", "24_ph_table_pea")

# plot ph

ph_mean <- ph_combined %>%
  pivot_longer(-c(sampleid, species, enzyme, block, box, bacteria), names_to = "time", values_to = "ph") %>%
  group_by(bacteria, time) %>%
  summarize(mean = mean(ph, na.rm=T), n = length(ph),
            sd = sd(ph, na.rm= T), .groups = "drop") %>%
  mutate(time = as.numeric(time))

ggplot(ph_mean, aes(x = time, y = mean)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1) +
  geom_line(aes(group = bacteria, color = bacteria), size = 1.5) +
  geom_point(aes(shape = bacteria), size = 3) +
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_color_manual(values = colors) +
  labs(title = "pH curve combinations", x = "time (h)", y = "pH", shape = "combinations", color = "combinations")

save_big("24_ph_curve_combinations")

# compare ammonia 

ammo <- combined_comparison_and_results(df = ammonia_combined, selected_response = "mg_nh4n_100g", 
                                                      transformation = "test",
                                                      type = "s", response_name = "ammonia", y_axis = "mg NH4-N/100g", save_name = "ammonia",
                                                      save = save)

# compare acetate and lactate 

lac <- combined_comparison_and_results(df = acids_combined, selected_response = "lactate", 
                                                     transformation = "test",
                                                     type = "h", response_name = "Lactate", y_axis = "mmol/kg", save_name = "lactate",
                                                     save = save)

ace <- combined_comparison_and_results(df = acids_combined, selected_response = "acetate", 
                                                     transformation = "test",
                                                     type = "h", response_name = "Acetate", y_axis = "mmol/kg", save_name = "acetate",
                                                     save = save)


# save table

ammo_acids_table <- cbind(ammo, "lactate" = lac[,2], "acetate" = ace[,2])

write_tsv(ammo_acids_table, "tables/24_combinations_acids_ammonia.tsv")
save_results_table(ammo_acids_table, title = "combinations", "24_ammo_ace_lac_table")

# Multivariate analysis

pca_df <- ph %>%
  dplyr::select(sampleid, ph_8 = "8", ph_24 = "24") %>%
  inner_join(ammonia, by = "sampleid") %>%
  inner_join(acids, by = "sampleid") %>%
  filter(!is.na(lactate)) %>%
  dplyr::rename(ammonia = mg_nh4n_100g) %>%
  inner_join(meta, by = "sampleid") %>%
  filter(bacteria != "blank")

pca <- prcomp(~ph_24+ph_8+ammonia+acetate+lactate, data = pca_df, scale. = TRUE, center = TRUE)
rownames(pca$x) <- pca_df$bacteria
biplot(pca)
rownames(pca$x) <- pca_df$sampleid

pca_out <- pca$x %>%
  as_tibble(rownames = "sampleid") %>%
  inner_join(meta, by = "sampleid")

ggplot(pca_out, aes(x = PC1, y = PC2, color = bacteria)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  scale_color_manual(values = colors)

save_big("24_pca")

pca_mean_df <- ph %>%
  dplyr::select(sampleid, ph_8 = "8", ph_24 = "24") %>%
  inner_join(ammonia, by = "sampleid") %>%
  inner_join(acids, by = "sampleid") %>%
  dplyr::rename(ammonia = mg_nh4n_100g) %>%
  inner_join(meta, by = "sampleid") %>%
  filter(bacteria != "blank") %>%
  group_by(bacteria) %>%
  summarize(ph_8 = mean(ph_8),
            ph_24 = mean(ph_24), 
            ammonia = mean(ammonia),
            lactate = mean(lactate, na.rm=T),
            acetate = mean(acetate))

pca_mean <- prcomp(~ph_24+ph_8+ammonia+acetate+lactate, data = pca_mean_df, scale. = TRUE, center = TRUE)
rownames(pca_mean$x) <- pca_mean_df$bacteria

jpeg(filename="plots/24_pca_mean.jpg", width = 15, height = 15, units = "cm", res = 100)
biplot(pca_mean)
dev.off()

screeplot(pca_mean)
