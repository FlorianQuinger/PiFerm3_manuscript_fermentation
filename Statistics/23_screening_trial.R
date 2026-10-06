library(here)

source("0_general_functions.R")
source("20_pretrial_functions.R")

save = FALSE
#save = TRUE

# load data

meta <- read_tsv("data/23_meta.txt") %>%
  mutate(block = as.character(block))

bacteria <- read_tsv("data/23_bacteria.txt")

ph <- read_tsv("data/23_ph.txt") %>%
  mutate(`48` = as.numeric(`48`))

ammonia <- read_tsv("data/23_ammonia.txt") %>%
  mutate(sampleid = str_remove(sampleid, " \\(frisch\\)"))

acids <- read_tsv("data/23_acids.txt")%>%
  mutate(total = acetate + lactate) # total amount of acids produced

ammonia_combined <- ammonia %>%
  left_join(meta, by = "sampleid") %>%
  left_join(bacteria, by = "bacteria") %>%
  mutate(substrate = ifelse(is.na(substrate), sampleid, substrate),
         bacteria = ifelse(is.na(bacteria), sampleid, bacteria),
         species = ifelse(is.na(species), sampleid, species),
         species_short = ifelse(is.na(species_short), sampleid, species_short)) %>%
  dplyr::select(-block)

acids_combined <- acids %>%
  left_join(meta, by = "sampleid") %>%
  left_join(bacteria, by = "bacteria") %>%
  mutate(substrate = ifelse(is.na(substrate), sampleid, substrate),
         bacteria = ifelse(is.na(bacteria), sampleid, bacteria),
         species = ifelse(is.na(species), sampleid, species),
         species_short = ifelse(is.na(species_short), sampleid, species_short)) %>%
  dplyr::select(-block)

ph_combined <- ph%>%
  inner_join(meta, by = "sampleid") %>%
  inner_join(bacteria, by = "bacteria")

ph_pea <- ph_combined %>%
  filter(substrate == "Erbse") %>%
  mutate(bacteria = ifelse(species == "Lactiplantibacillus plantarum", paste(species_short, bacteria), species_short),
         bacteria = factor(bacteria, levels = c("Control", sort(unique(bacteria)[which(unique(bacteria)!="Control")]))))

# comparison of pH at different timepoints

ph0 <- combined_comparison_and_results(df = ph_pea, selected_response = "0", transformation = "none",
                                       type = "h", response_name = "pH 0", y_axis = "pH", save_name = "ph_0_pea",
                                       save = save, digits = 2)

ph2 <- combined_comparison_and_results(df = ph_pea, selected_response = "2", transformation = "none",
                                       type = "h", response_name = "pH 2", y_axis = "pH", save_name = "ph_2_pea",
                                       save = save, digits = 2)

ph4 <- combined_comparison_and_results(df = ph_pea, selected_response = "4", transformation = "none",
                                       type = "h", response_name = "pH 4", y_axis = "pH", save_name = "ph_4_pea",
                                       save = save, digits = 2)

ph6 <- combined_comparison_and_results(df = ph_pea, selected_response = "6", transformation = "test",
                                       type = "h", response_name = "pH 6", y_axis = "pH", save_name = "ph_6_pea",
                                       save = save, digits = 2)

ph8 <- combined_comparison_and_results(df = ph_pea, selected_response = "8", transformation = "test",
                                       type = "h", response_name = "pH 8", y_axis = "pH", save_name = "ph_8_pea",
                                       save = save, digits = 2)

ph12 <- combined_comparison_and_results(df = ph_pea, selected_response = "12", transformation = "test",
                                       type = "h", response_name = "pH 12", y_axis = "pH", save_name = "ph_12_pea",
                                       save = save, digits = 2)

ph24 <- combined_comparison_and_results(df = ph_pea, selected_response = "24", transformation = "test",
                                        type = "h", response_name = "pH 24", y_axis = "pH", save_name = "ph_24_pea",
                                        save = save, digits = 2)

# combine results

ph_table_pea <- cbind(ph0, "pH 2" = ph2[,2], "pH 4" = ph4[,2], "pH 6" = ph6[,2], "pH 8" = ph8[,2],
                      "pH 12" = ph12[,2], "pH 24" = ph24[,2]) %>%
  left_join(dplyr::select(bacteria, bacteria, species_short), by = "bacteria") %>%
  mutate(bacteria = ifelse(is.na(species_short), bacteria, species_short)) %>%
  dplyr::select(-species_short)

write_tsv(ph_table_pea, "tables/23_screening_ph.tsv")
save_results_table(ph_table_pea, title = "pH pea", "23_ph_table_pea")

# plot curves for pea

ph_mean_pea <- ph %>%
  dplyr::select(-`48`) %>%
  pivot_longer(-sampleid, names_to = "time", values_to = "ph") %>%
  inner_join(meta, by = "sampleid") %>%
  filter(substrate == "Erbse") %>%
  group_by(bacteria, time) %>%
  summarize(mean = mean(ph, na.rm=T), n = length(ph),
            sd = sd(ph, na.rm= T), .groups = "drop") %>%
  inner_join(bacteria, by = "bacteria") %>%
  mutate(bacteria = ifelse(species == "Lactiplantibacillus plantarum", paste(species_short, bacteria), species_short),
         bacteria = factor(bacteria, levels = c("Control", "L. plantarum A10", "C. kimchii", "L. spicheri", "W. confusa", "B. licheniformis", "L. plantarum A2", "P. pentosaceous", "L. suionicum", "B. subtilis"))) %>%
  mutate(time = as.numeric(time))

ggplot(ph_mean_pea, aes(x = time, y = mean)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1) +
  geom_line(aes(group = bacteria, color = fermentation), size = 1.5) +
  geom_point(aes(shape = bacteria), size = 3) +
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10)) +
  scale_color_manual(values = colors) +
  labs(title = "pH curve pea", x = "time (h)", y = "pH", shape = "species")

save_big("23_ph_curve_pea_24h")


#### Ammonia plotting

ammonia_pea <- ammonia_combined %>%
  filter(substrate == "Erbse") %>%
  mutate(species_short = ifelse(species_short == "Erbse", "unfermented", species_short)) %>%
  mutate(species_short = ifelse(species == "Lactiplantibacillus plantarum", 
                                paste(species_short, bacteria), species_short),
  species_short = factor(species_short, levels=c("unfermented", "Control", "L. plantarum A10", "C. kimchii", "L. spicheri", "W. confusa", "B. licheniformis", "L. plantarum A2", "P. pentosaceous", "L. suionicum", "B. subtilis")))

write_tsv(ammonia_pea, "tables/23_screening_ammonia.tsv")

# plots for pea

ggplot(ammonia_pea, aes(x = species_short, y = mg_NH4N_100g)) +
  geom_bar(stat = "identity", width = 0.4)+
  scale_y_continuous(limits = c(min(0, ammonia_pea$mg_NH4N_100g),1.1*max(ammonia_pea$mg_NH4N_100g)), expand = expansion(mult = c(0, .1))) +
  labs(x = "bacteria", y = "mg NH4-N in 100g ferment", title = "NH4-N in fermented pea") +
  guides(x=guide_axis(angle = 30))

save_big("23_ammonia_pea_mg_100g")

ammonia_pea %>%
  filter(mg_NH4N_formation != 0) %>%
  ggplot(aes(x = species_short, y = mg_NH4N_formation, fill = fermentation)) +
  geom_bar(stat = "identity", width = 0.4)+
  scale_y_continuous(limits = c(min(0, ammonia_pea$mg_NH4N_formation),1.1*max(ammonia_pea$mg_NH4N_formation)), expand = expansion(mult = c(0, .1))) +
  labs(x = "bacteria", y = "mg NH4-N formation in 100g ferment", title = "NH4-N formation in fermented pea") +
  guides(x=guide_axis(angle = 30)) +
  scale_fill_manual(values = colors)

save_big("23_ammonia_pea_mg_formation_100g")

ammonia_pea %>%
  filter(N_loss_percent != 0) %>%
  ggplot(aes(x = species_short, y = N_loss_percent, fill = fermentation)) +
  geom_bar(stat = "identity", width = 0.4)+
  scale_y_continuous(limits = c(min(0, ammonia_pea$N_loss_percent),1.1*max(ammonia_pea$N_loss_percent)), expand = expansion(mult = c(0, .1))) +
  labs(x = "bacteria", y = "NH4-N formation of total N (%)", title = "relative NH4-N formatian of total N") +
  guides(x=guide_axis(angle = 30)) +
  scale_fill_manual(values = colors)

save_big("23_ammonia_pea_relative")

# Acids plots

acids_pea <- acids_combined %>%
  filter(substrate == "Erbse") %>%
  mutate(species_short = ifelse(species == "Lactiplantibacillus plantarum", 
                                paste(species_short, bacteria), species_short),
         species_short = factor(species_short, levels=c("Control", "L. plantarum A10", "C. kimchii", "L. spicheri", "W. confusa", "B. licheniformis", "L. plantarum A2", "P. pentosaceous", "L. suionicum", "B. subtilis"))) %>%
  pivot_longer(c(acetate, lactate), names_to = "acid", values_to = "concentration")

write_tsv(acids_pea, "tables/23_screening_acids.tsv")


ggplot(acids_pea, aes(x = species_short, y = concentration)) +
  geom_bar(stat = "identity", width = 0.4, position = "dodge") +
  labs(x = "bacteria", y = "concentration (mmol/kg)", title = "Acid concentrations pea") +
  guides(x=guide_axis(angle = 30)) +
  scale_fill_manual(values = colors) +
  facet_grid(cols = vars(acid), scales = "free_y") +
  theme(legend.position = "bottom")

save_big("23_acids_pea")