source("0_general_functions.R")

library(here)

source("0_general_functions.R")
source("20_pretrial_functions.R")
source("50_omics_functions.R")
source("E:/R/source/ggplot2_theme_Fermentation.R")




# Figure 1 fermentation strains annotation

# import files

ec_table <- read_tsv("data/ec_meta.txt")
ec_table <- read_tsv("data/ec_meta_new.txt") %>%
  mutate(Substrate = ifelse(str_detect(Substrate, "galactosides"), "alpha-galactosides", Substrate)) %>%
  dplyr::rename(EC = 'Enzyme Commission number', substrate = Substrate)
strain_meta = read_tsv("data/strain_meta.txt") %>%
  dplyr::select(species, strain) %>%
  filter(!is.na(strain)) %>%
  mutate(strain = tolower(strain)) %>%
  mutate(species = paste(species, strain)) %>%
  mutate(species = str_replace(species, "dsm", "DSM ")) %>%
  mutate(species = str_replace(species, "a3_flongle", "A3")) %>%
  mutate(species = str_replace(species, "a10_flongle", "A10")) %>%
  mutate(species = str_replace(species, "a2_flongle", "A2")) %>%
  mutate(species = str_replace(species, "i32_dragonflye", "I32"))

files_prokka = list.files("data/fermentation_strains3/")
files_deep = list.files("data/fermentation_strains4/")

# bind strain files

for (i in files_prokka) {
  strain_name = str_remove(i, ".tsv")
  file_name = paste0("prokka_", str_remove(i, ".tsv"))
  file = read_tsv(paste0("data/fermentation_strains3/", i)) %>% 
    dplyr::select("locus" = locus_tag, "ec_number" = EC_number) %>%
    add_column(strain = strain_name)
  file_name = paste0("prokka_", str_remove(i, ".tsv"))
  assign(file_name, file)
}

for (i in files_deep) {
  strain_name = str_remove(i, ".tsv")
  file = read_tsv(paste0("data/fermentation_strains4/", i)) %>% 
    dplyr::select("locus" = sequence_ID, "ec_number" = prediction) %>%
    add_column(strain = strain_name)
  file_name = paste0("deep_", str_remove(i, ".tsv"))
  assign(file_name, file)
}

# bind in one df

prokka_strains <- do.call(rbind, mget(ls(pattern = "prokka_dsm|prokka_A|prokka_I"))) %>%
  mutate(strain = tolower(strain)) %>%
  inner_join(strain_meta, by = "strain")

deep_strains <- do.call(rbind, mget(ls(pattern = "deep_dsm|deep_A|deep_I"))) %>%
  mutate(strain = tolower(strain)) %>%
  inner_join(strain_meta, by = "strain") %>%
  separate_rows(ec_number, sep = ";") %>%
  mutate(ec_number = str_remove(ec_number, "EC:"))

combined_strains <- rbind(prokka_strains, deep_strains) %>%
  distinct() %>%
  filter(!is.na(ec_number)) %>%
  filter(ec_number != "None")

# find ecs

ec <- unique(ec_table$EC)

create_reaction_matrix <- function(strain_df, ec_list) {
  strain_ec_matrix <- matrix(nrow = length(unique(strain_df$species)), ncol = length(ec_list))
  colnames(strain_ec_matrix) <- ec_list
  rownames(strain_ec_matrix) <- unique(strain_df$species)
  
  for (i in 1:nrow(strain_ec_matrix)) {
    for (j in 1:ncol(strain_ec_matrix)) {
      strain_name = rownames(strain_ec_matrix)[i]
      symbol_name = colnames(strain_ec_matrix)[j]
      check <- strain_df %>%
        filter(species == strain_name) %>%
        filter(ec_number == symbol_name)
      value = nrow(check) 
      strain_ec_matrix[i,j] <- value
    }
  }
  print(heatmap(strain_ec_matrix))
  return(strain_ec_matrix)
}

# create ec matix

ec_matrix <- create_reaction_matrix(combined_strains, ec)


ec_plot_selection <- ec_matrix %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% c("Lactiplantibacillus plantarum A10",
                        "Lactiplantibacillus plantarum A2", 
                        "Pediococcus pentosaceus A3",
                        "Bacillus licheniformis I32",
                        "Weissella confusa DSM 20196",
                        "Levilactobacillus spicheri DSM 15429",
                        "Leuconostoc suionicum DSM 20241",
                        "Companilactobacillus kimchii DSM 13961",
                        "Bacillus subtilis DSM 3257")) %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  group_by(ec_number) %>%
  filter(sum(value)>0) %>%
  ungroup() %>%
  mutate(ec_number = factor(ec_number, levels = unique(ec_table$EC))) %>%
  ggplot(aes(x = ec_number, y = species, fill = value)) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = T) +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_markdown(family = "arial", size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_gradient(low = "white", high = "black", breaks = c(0,3,6,9,12)) +
  labs(fill = "Number of genes") +
  scale_y_discrete(labels = c("*Bacillus licheniformis* I32",
                              "*Bacillus subtilis* DSM 3257",
                              "*Companilactobacillus kimchii* DSM 13961",
                              "*Lactiplantibacillus plantarum* A10",
                              "*Lactiplantibacillus plantarum* A2",
                              "*Leuconostoc suionicum* DSM 20241",
                              "*Levilactobacillus spicheri* DSM 15429",
                              "*Pediococcus pentosaceus* A3",
                              "*Weissella confusa* DSM 20196"))

legend <- ec_matrix %>% 
  as_tibble(rownames = "species") %>%
  filter(species %in% c("Lactiplantibacillus plantarum A10",
                        "Lactiplantibacillus plantarum A2", 
                        "Pediococcus pentosaceus A3",
                        "Bacillus licheniformis I32",
                        "Weissella confusa DSM 20196",
                        "Levilactobacillus spicheri DSM 15429",
                        "Leuconostoc suionicum DSM 20241",
                        "Companilactobacillus kimchii DSM 13961",
                        "Bacillus subtilis DSM 3257")) %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  group_by(ec_number) %>%
  filter(sum(value)>0) %>%
  ungroup() %>%
  mutate(ec_number = factor(ec_number, levels = ec_table$EC)) %>%
  separate(substrate, into = c("substrate"), sep = ",") %>%
  mutate(substrate = factor(substrate, levels = unique(substrate))) %>%
  dplyr::select(ec_number, substrate) %>%
  distinct()%>%
  ggplot(aes(x = ec_number, y = 1, fill = substrate)) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = T) +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.key.size = unit(.4, "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(position = "bottom", ncol = 2, title = "Substrate", theme(legend.title.position = "top")))

fig1 <- legend / ec_plot_selection +
  plot_layout(heights = c(1,12), guides = "collect") &
  theme(legend.position = "bottom")
print(fig1)

ggsave(filename= "95_fig1.jpeg",
       plot = fig1,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 12,
       scale=1,
       dpi=600)

# Figure 2 of screening trial


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

ph_mean_pea <- ph %>%
  pivot_longer(-sampleid, names_to = "time", values_to = "ph") %>%
  inner_join(meta, by = "sampleid") %>%
  filter(substrate == "Erbse") %>%
  group_by(bacteria, time) %>%
  summarize(mean = mean(ph, na.rm=T), n = length(ph),
            sd = sd(ph, na.rm= T), .groups = "drop") %>%
  inner_join(bacteria, by = "bacteria") %>%
  mutate(bacteria = ifelse(species == "Lactiplantibacillus plantarum", paste(species_short, bacteria), species_short),
         bacteria = factor(bacteria, levels = c(c("Control", "B. licheniformis", "B. subtilis", "C. kimchii", 
                                                  "L. plantarum A10",  "L. plantarum A2",  "L. suionicum",
                                                  "L. spicheri", "P. pentosaceous", "W. confusa")))) %>%
  mutate(fermentation = ifelse(fermentation %in% c("Bacillus", "none"), "other", fermentation)) %>%
  mutate(time = as.numeric(time))

fig2a <- ggplot(ph_mean_pea, aes(x = time, y = mean)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 0.5) +
  geom_line(aes(group = bacteria, color = fermentation), size = 1) +
  geom_point(aes(shape = bacteria), size = 2) +
  scale_x_continuous(breaks = c(0,2,4,6,8,12,24,48)) +
  scale_shape_manual(values = c(4,2,3,1,5,6,7,8,9,10),
                     labels = c("Control", 
                                expression(italic("B. licheniformis")~"I32"),
                                expression(italic("B. subtilis")~"DSM 3257"),
                                expression(italic("C. kimchii")~"DSM 13961"), 
                                expression(italic("L. plantarum")~"A10"),
                                expression(italic("L. plantarum")~"A2"),
                                expression(italic("L. suionicum")~"DSM 20241"),
                                expression(italic("L. spicheri")~"DSM 15429"),
                                expression(italic("P. pentosaceus")~"A3"), 
                                expression(italic("W. confusa")~"DSM 20196"))) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "grey")) +
  labs(x = "Time (h)", y = "pH", shape = "Bacterial strain", color = "Type") +
  guides(shape = guide_legend(position = "right", ncol = 1, title = "Bacterial strain",
                              theme(legend.title.position = "top")),
         color = guide_legend(position = "bottom", ncol = 3, title = "Type",
                              theme(legend.title.position = "left"))) +
  theme(legend.justification.bottom = "left",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


# acids and ammonia

acids_pea <- acids_combined %>%
  filter(substrate == "Erbse") %>%
  mutate(species_short = ifelse(species == "Lactiplantibacillus plantarum", 
                                paste(species_short, bacteria), species_short),
         species_short = factor(species_short, levels = c(c("Control", "B. licheniformis", "B. subtilis", "C. kimchii", 
                                                  "L. plantarum A10",  "L. plantarum A2",  "L. suionicum",
                                                  "L. spicheri", "P. pentosaceous", "W. confusa"))))

ammonia_pea <- ammonia_combined %>%
  filter(substrate == "Erbse") %>%
  filter(!bacteria %in% c("Erbse"))  %>%
  mutate(mmol_NH4N_kg = mg_NH4N_100g / 14.007 * 10) %>%
  mutate(species_short = ifelse(species == "Lactiplantibacillus plantarum", 
                                paste(species_short, bacteria), species_short),
         species_short = factor(species_short, levels = c(c("Control", "B. licheniformis", "B. subtilis", "C. kimchii", 
                                                            "L. plantarum A10",  "L. plantarum A2",  "L. suionicum",
                                                            "L. spicheri", "P. pentosaceous", "W. confusa"))))

metabolites <- acids_pea %>%
  left_join(ammonia_pea, by = "species_short") %>%
  pivot_longer(c(acetate, lactate, mmol_NH4N_kg), names_to = "metabolite", values_to = "value")

fig2b <- ggplot(metabolites, aes(x = species_short, y = value)) +
  geom_bar(stat = "identity", width = 0.4, position = "dodge") +
  labs(y = "Concentration (mmol/kg)") +
  guides(x=guide_axis(angle = 45)) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = c("Control",
                                expression(italic("B. licheniformis")~"I32"),
                                expression(italic("B. subtilis")~"DSM 3257"),
                                expression(italic("C. kimchii")~"DSM 13961"),
                                expression(italic("L. plantarum")~"A10"),
                                expression(italic("L. plantarum")~"A2"),
                                expression(italic("L. suionicum")~"DSM 20241"),
                                expression(italic("L. spicheri")~"DSM 15429"),
                                expression(italic("P. pentosaceus")~"A3"),
                                expression(italic("W. confusa")~"DSM 20196"))) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        strip.text = element_markdown(family = "arial", size = 14)) +
  facet_wrap(~metabolite, scales = "free_y", labeller = labeller(metabolite = c(acetate = "Acetic acid", 
                                                                                lactate = "Lactic acid",
                                                                                mmol_NH4N_kg = "Ammonium")))

# fig2b <- ggplot(acids_pea, aes(x = species_short, y = acetate)) +
#   geom_bar(stat = "identity", width = 0.4, position = "dodge") +
#   labs(y = "Concentration (mmol/kg)") +
#   guides(x=guide_axis(angle = 30)) +
#   scale_fill_manual(values = colors) +
#   scale_x_discrete(labels = c("Control", 
#                                 expression(italic("B. licheniformis")~"I32"),
#                                 expression(italic("B. subtilis")~"DSM 3257"),
#                                 expression(italic("C. kimchii")~"DSM 13961"), 
#                                 expression(italic("L. plantarum")~"A10"),
#                                 expression(italic("L. plantarum")~"A2"),
#                                 expression(italic("L. suionicum")~"DSM 20241"),
#                                 expression(italic("L. spicheri")~"DSM 15429"),
#                                 expression(italic("P. pentosaceus")~"A3"), 
#                                 expression(italic("W. confusa")~"DSM 20196"))) +
#   scale_y_continuous(limits = c(0,45), expand = c(0,0)) +
#   theme(axis.title.x = element_blank())
# 
# fig2c <- ggplot(acids_pea, aes(x = species_short, y = lactate)) +
#   geom_bar(stat = "identity", width = 0.4, position = "dodge") +
#   labs(y = "Concentration (mmol/kg)") +
#   guides(x=guide_axis(angle = 30)) +
#   scale_fill_manual(values = colors) +
#   scale_x_discrete(labels = c("Control", 
#                               expression(italic("B. licheniformis")~"I32"),
#                               expression(italic("B. subtilis")~"DSM 3257"),
#                               expression(italic("C. kimchii")~"DSM 13961"), 
#                               expression(italic("L. plantarum")~"A10"),
#                               expression(italic("L. plantarum")~"A2"),
#                               expression(italic("L. suionicum")~"DSM 20241"),
#                               expression(italic("L. spicheri")~"DSM 15429"),
#                               expression(italic("P. pentosaceus")~"A3"), 
#                               expression(italic("W. confusa")~"DSM 20196"))) +
#   scale_y_continuous(limits = c(0,150), expand = c(0,0)) +
#   theme(axis.title.x = element_blank())
# 
# fig2d <- ggplot(ammonia_pea, aes(x = species_short, y = mmol_NH4N_kg)) +
#   geom_bar(stat = "identity", width = 0.4, position = "dodge") +
#   labs(y = "Concentration (mmol/kg)") +
#   guides(x=guide_axis(angle = 30)) +
#   scale_fill_manual(values = colors) +
#   scale_x_discrete(labels = c("Control", 
#                               expression(italic("B. licheniformis")~"I32"),
#                               expression(italic("B. subtilis")~"DSM 3257"),
#                               expression(italic("C. kimchii")~"DSM 13961"), 
#                               expression(italic("L. plantarum")~"A10"),
#                               expression(italic("L. plantarum")~"A2"),
#                               expression(italic("L. suionicum")~"DSM 20241"),
#                               expression(italic("L. spicheri")~"DSM 15429"),
#                               expression(italic("P. pentosaceus")~"A3"), 
#                               expression(italic("W. confusa")~"DSM 20196"))) +
#   scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
#   theme(axis.title.x = element_blank())

fig2 <- fig2a /
  fig2b +
  plot_layout(design = "
              AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
              AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
              AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
              BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
              BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
              ") &
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(filename= "95_fig2.jpeg",
       plot = fig2,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 18,
       scale=1,
       dpi=600)

# Figure 3 of combinations trial

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

# plot ph

ph_mean <- ph_combined %>%
  pivot_longer(-c(sampleid, species, enzyme, block, box, bacteria), names_to = "time", values_to = "ph") %>%
  group_by(bacteria, time) %>%
  summarize(mean = mean(ph, na.rm=T), n = length(ph),
            sd = sd(ph, na.rm= T), .groups = "drop") %>%
  mutate(time = as.numeric(time))

fig3a <- ggplot(ph_mean, aes(x = time, y = mean)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 0.5) +
  geom_line(aes(group = bacteria, color = bacteria), size = 1, show.legend = F) +
  geom_point(aes(shape = bacteria, color = bacteria), size = 3) +
  scale_x_continuous(breaks = c(0,2,4,6,8,12,24)) +
  scale_color_manual(values = c("grey", colors[-c(1,2,3,4,5,11)]),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  scale_shape_manual(values = c(4, 16,16,16, 15,15,15, 17),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  labs(x = "Time (h)", y = "pH", shape = "Bacterial strains", color = "Bacterial strains") +
  # guides(shape = guide_legend(position = "right", ncol = 1, title = "Bacterial strains",
  #                             theme(legend.title.position = "top")),
  #        color = guide_legend(position = "bottom", ncol = 3, title = "Type",
  #                             theme(legend.title.position = "left"))) +
  theme(legend.justification.bottom = "left",
        legend.text = element_markdown(family = "arial", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(.4, "cm"))

# acids and ammonia

ammo <- combined_comparison_and_results(df = ammonia_combined, selected_response = "mg_nh4n_100g", 
                                        transformation = "test",
                                        type = "s", response_name = "ammonia", y_axis = "mg NH4-N/100g", save_name = "ammonia",
                                        save = save) %>%
  pivot_wider(names_from = bacteria, values_from = ammonia) %>% 
  pivot_longer(-c("Pooled SEM", "P-value"), values_to = "value", names_to = "bacteria") %>%
  add_column(metabolite = "ammonia")

lac <- combined_comparison_and_results(df = acids_combined, selected_response = "lactate", 
                                       transformation = "test",
                                       type = "h", response_name = "Lactate", y_axis = "mmol/kg", save_name = "lactate",
                                       save = save) %>%
  pivot_wider(names_from = bacteria, values_from = Lactate) %>% 
  pivot_longer(-c("Pooled SEM", "P-value"), values_to = "value", names_to = "bacteria") %>%
  add_column(metabolite = "lactate")

ace <- combined_comparison_and_results(df = acids_combined, selected_response = "acetate", 
                                       transformation = "test",
                                       type = "h", response_name = "Acetate", y_axis = "mmol/kg", save_name = "acetate",
                                       save = save) %>%
  pivot_wider(names_from = bacteria, values_from = Acetate) %>% 
  pivot_longer(-c("Pooled SEM", "P-value"), values_to = "value", names_to = "bacteria") %>%
  add_column(metabolite = "acetate")

metabolites <- rbind(ammo, lac, ace) %>%
  mutate(cld = str_extract(value, "[a-z]+"),
         value = str_remove(value, "[a-z]+"),
         value = as.numeric(value),
         `Pooled SEM` = as.numeric(`Pooled SEM`),
         bacteria = factor(bacteria, levels = levels(meta$bacteria)),
         metabolite = factor(metabolite, levels = c("acetate", "lactate", "ammonia"))) %>%
  mutate(value = ifelse(metabolite == "ammonia", value / 14.007 * 10, value)) %>% # transform to mmol/kg
  group_by(metabolite) %>%
  mutate(max = (max(value) + `Pooled SEM`)*1.1)

fig3b <- ggplot(metabolites, aes(x = bacteria, y = value)) +
  geom_errorbar(aes(ymax = value + `Pooled SEM`, ymin = value - `Pooled SEM`), width = 0.25, size = 0.5) +
  geom_bar(aes(fill = bacteria), color = "grey20", stat = "identity", width = 0.4, position = "dodge", show.legend = F) +
  geom_point(aes(x = bacteria, y = max, shape = bacteria, color = bacteria), show.legend = F, size = 3) +
  geom_blank(aes(x = bacteria, y = max*1.03)) +
  geom_text(aes(x = bacteria, y = value + `Pooled SEM` + max * 0.04 , label = cld), family = "arial", size = 4) +
  labs(y = "Concentration (mmol/kg)") +
  guides(x=guide_axis(angle = 45)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("grey", colors[-c(1,2,3,4,5,11)]),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  scale_color_manual(values = c("grey", colors[-c(1,2,3,4,5,11)]),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  scale_shape_manual(values = c(4, 16,16,16, 15,15,15, 17),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        strip.text = element_markdown(family = "arial", size = 14)) +
  facet_wrap(~metabolite, scales = "free_y", labeller = labeller(metabolite = c(acetate = "Acetic acid", 
                                                                                lactate = "Lactic acid",
                                                                                ammonia = "Ammonium")))

legend <- ggplot(metabolites, aes(x = bacteria, y = 1)) +
  geom_point(aes(x = bacteria, y = 0, shape = bacteria, color = bacteria), size = 2, show.legend = F) +
  scale_color_manual(values = c("grey", colors[-c(1,2,3,4,5,11)]),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  scale_shape_manual(values = c(4, 16,16,16, 15,15,15, 17),
                     labels = c("Control",
                                "*L. plantarum* A10",
                                "*C. kimchii* DSM 13961",
                                "*L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961",
                                "*L. plantarum* A10 +<br> *L. spicheri* DSM 15429",
                                "*C. kimchii* DSM 13961 +<br> *L. spicheri* DSM 15429",
                                "*L. plantarum* A10 +<br> *C. kimchii* DSM 13961 +<br> *L.spicheri* DSM 15429"
                     )) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~metabolite)

fig3b_ <- fig3b /
  legend +
  plot_layout(heights = c(7, 0.2), guides = "collect")

fig3 <- fig3a /
  fig3b /
  legend +
  plot_layout(design = "
              AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
              AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA#
              BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
              BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
              ") &
  plot_annotation(tag_levels = list(c("a", "b", " "))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))


ggsave(filename= "95_fig3.jpeg",
       plot = fig3,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 18,
       scale=1,
       dpi=600)

# Figure 4 panels

# load nutrition data

pea_analysis <- read_tsv("data/nutrition_pea_pre_post.txt") %>%
  dplyr::select(-CA) %>% # exclude due to duplication and no data at this point
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  mutate(Fermentation = factor(c("unfermented", "unfermented", "fermented", "fermented"), 
                            levels = c("unfermented", "fermented"))) %>%
  mutate_all(., ~as.character(.)) %>%
  pivot_longer(-c(Fermentation, diet, description), names_to = "nutrient", values_to = "concentration")

nutrients_selected <- pea_analysis %>%
  filter(nutrient %in% c("tdf", "arg", "insp6...43", "sucrose", "raffinose", "galactose*/maltitol")) %>%
  mutate(nutrient = case_when(nutrient == "tdf" ~ "Total dietary fibre",
                              nutrient == "arg" ~ "Arginine", 
                              nutrient == "insp6...43" ~ "InsP~6~",
                              nutrient == "sucrose" ~ "Sucrose",
                              nutrient == "raffinose" ~ "Raffinose",
                              nutrient == "galactose*/maltitol" ~ "Galactose/Maltitol")) %>%
  mutate(Fermentation = factor(Fermentation, levels = c("fermented", "unfermented"))) %>%
  mutate(concentration = as.numeric(concentration)) %>%
  mutate(concentration = ifelse(nutrient %in% c("Sucrose", "Raffinose", "Galactose/Maltitol"), 
                                concentration/1000, concentration)) %>%
  group_by(nutrient, Fermentation) %>%
  mutate(min = min(concentration),
         max = max(concentration),
         mean = mean(concentration)) %>%
  ungroup() %>%
  dplyr::select(-c(diet, concentration, description)) %>%
  distinct() %>%
  mutate(nutrient = factor(nutrient, levels = c("Total dietary fibre","Galactose/Maltitol", "Sucrose",
                                                "Raffinose","InsP~6~","Arginine")))

fig4a1 <- ggplot(filter(nutrients_selected, nutrient == "Total dietary fibre"),
              aes(y = nutrient, x = mean, fill = Fermentation)) +
  geom_bar(stat = "identity", position = "dodge", width = .6, color = "black", show.legend = F) +
  geom_errorbar(aes(xmin = min, xmax = max), position = position_dodge(width = .6), width = .2) +
  scale_fill_manual(values = colors[c(4,1)]) +
  scale_x_continuous(name = "g/kg DM",
                     expand = c(0,0), limits = c(0,200)) +
  labs(fill = "") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))

fig4a2 <- ggplot(filter(nutrients_selected, nutrient != "Total dietary fibre"),
       aes(y = nutrient, x = mean, fill = Fermentation)) +
  geom_bar(stat = "identity", position = "dodge", width = .6, color = "black") +
  geom_errorbar(aes(xmin = min, xmax = max), position = position_dodge(width = .6), width = .2) +
  scale_fill_manual(values = colors[c(4,1)]) +
  scale_x_continuous(name = "g/kg DM", 
                     expand = c(0,0), limits = c(0,30)) +
  scale_y_discrete(labels = c("Galactose/Maltitol", "Sucrose", "Raffinose",
                              expression("InsP"[6]),
                              "Arginine")) +
  labs(fill = "") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

# metabolites

nmr_ferment <- readRDS("clean/7_nmr_ferment_fm.RDS") %>%
  filter(sampleid != "409") 

nmr_ferment_mean <- nmr_ferment %>%
  pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
  group_by(metabolite) %>%
  mutate(min = min(concentration),
         max = max(concentration),
         mean = mean(concentration)) %>%
  ungroup() %>%
  dplyr::select(-concentration, -sampleid) %>%
  distinct() %>%
  filter(metabolite %in% c("lactate", "ethanol", "acetate", "propionate", "butyrate")) %>%
  mutate(metabolite = case_when(metabolite == "lactate" ~ "Lactate",
                                metabolite == "acetate" ~ "Acetate",
                                metabolite == "ethanol" ~ "Ethanol",
                                metabolite == "propionate" ~ "Propionate",
                                metabolite == "butyrate" ~ "Butyrate")) %>%
  arrange(mean) %>%
  mutate(metabolite = factor(metabolite, levels = metabolite))

fig4b <- ggplot(nmr_ferment_mean,
                aes(y = metabolite, x = mean)) +
  geom_bar(stat = "identity", position = "dodge", width = .6, color = "black", show.legend = F) +
  geom_errorbar(aes(xmin = min, xmax = max), position = position_dodge(width = .6), width = .2) +
  scale_x_continuous(expand = c(0,0), limits = c(0,225)) +
  labs(x = "mmol/kg") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))



# taxa barplots for ferment

# load metaproteomics and 16S data

meta <- readRDS("clean/meta1.RDS")

feature_ktable_rel_abd_filtered <- readRDS("clean/5_16s_ferment_taxonomy_feature_filtered_rel_abd.RDS") 

taxonomy_norm_imp_rel_abd_filtered <- readRDS("clean/6_ferment_taxonomy_long_norm_imp_rel_abd_filtered.RDS") %>% 
  filter(sampleid != "409")

filtered_table_fe_g <- filter_ktable(feature_ktable_rel_abd_filtered, meta = meta, selected_rank = "G",
                                     selected_matrix = "ferment") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  add_column(technique = "16S rRNA gene sequencing")

filtered_table_fe_p <- filter_ktable(taxonomy_norm_imp_rel_abd_filtered, meta = meta, selected_rank = "G",
                                     selected_matrix = "ferment") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  add_column(technique = "Metaproteomics")

plotting_table <- rbind(filtered_table_fe_g, filtered_table_fe_p) %>%
  group_by(name) %>%
  mutate(max = max(rel_abd)) %>%
  ungroup() %>%
  mutate(name = ifelse(max < .5, "other", name)) %>%
  group_by(technique, name, max) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  arrange(desc(max)) %>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

fig4c <- ggplot(plotting_table, aes(x = technique, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "", y = "Relative abundance (%)", fill = "Genera") +
  scale_fill_manual(values = c(colors[1:(length(unique(plotting_table$name))-1)], "grey")) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0)) + 
  scale_x_discrete(labels = c("16S rRNA",
                              "Proteins")) +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "plain"),
        legend.key.size = unit(.4, "cm")) 

# Microbial proteins

proteins_norm_imp_rel_abd <- readRDS("clean/6_ferment_proteins_long_norm_imp_rel_abd_filtered.RDS")

functions <- readRDS("clean/6_ferment_functions_long_raw.RDS")

proteins_norm_imp_rel_abd_filtered <- filter(proteins_norm_imp_rel_abd, sampleid != "409")

proteins_micro_rel_abd_filtered <- proteins_norm_imp_rel_abd_filtered %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

proteins_micro_rel_abd_400 <- proteins_micro_rel_abd_filtered %>%
  group_by(proteinid) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  left_join(functions, by = "proteinid")

Ckimchi_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGPF11641")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

# Wconfusa

Wconfusa_proteins_400 <- proteins_micro_rel_abd_400 %>%
  filter(str_detect(proteinid, "MGYGPF1022")) %>%
  mutate(sampleid = "401") %>% # for compatibility reasons
  recalculate_rel_abd() 

# COG

Ckimchi_cog_400 <- Ckimchi_proteins_400 %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category") %>%
  add_column(genus = "Companilactobacillus")

Wconfusa_cog_400 <- Wconfusa_proteins_400 %>%
  calculate_cog_cat_abundance(abundance_column = "rel_abd") %>%
  recalculate_rel_abd() %>%
  left_join(cog_cat_table, by = "cog_category") %>%
  add_column(genus = "Weissella")

plotting_table <- rbind(Ckimchi_cog_400, Wconfusa_cog_400) %>%
  group_by(cat_description) %>%
  mutate(max = max(rel_abd)) %>%
  ungroup() %>%
  mutate(cat_description = ifelse(max < 1, "other", cat_description)) %>%
  group_by(genus, cat_description) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(cat_description = factor(cat_description, levels = c(unique(cat_description)[-which(unique(cat_description)=="other")], "other")))

fig4d <- ggplot(plotting_table, aes(x = genus, y = rel_abd, fill = cat_description)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "", y = "Relative abundance (%)", fill = "COG category") +
  scale_fill_manual(values = c(colors[1:(length(unique(plotting_table$cat_description))-1)], "grey")) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0)) + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), 
        legend.key.size = unit(.4, "cm")) 


fig4 <- fig4a1 /
  fig4a2 /
  ((free(fig4b) | fig4c) + plot_layout(widths = c(7,10))) /
  free(fig4d, side = "lr") +
  plot_layout(heights = c(1,4,7.5,7.5)) +
  plot_annotation(tag_levels = list(c("a", "", "b", "c", "d"))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))


fig4 <- fig4b /
  ((free((fig4a1 / fig4a2) + plot_layout(heights = c(1,5))) | fig4c) + plot_layout(widths = c(10,8))) /
  free(fig4d, side = "lr") +
  plot_layout(heights = c(3,7.5,7.5)) +
  plot_annotation(tag_levels = list(c("a", "b", "", "c", "d"))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(filename= "95_fig4.jpeg",
       plot = fig4,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=1,
       dpi=600)

# Figure S2 16S

# Bray curtis from qiime2

meta <- readRDS("clean/meta1.RDS")

bray <- read_tsv("data/16S_ferment/distance-matrix.tsv") %>%
  mutate(`...1` = str_remove(`...1`, "IS-Pig-")) %>%
  rename_all(., ~gsub("IS-Pig-", "",.)) %>%
  column_to_rownames("...1") 

bray_fe <- bray %>%
  as.dist()

figS2a <- do_ordination(bray_fe, region = "ferment", title = "", save_name = "57_pcoa_ferment_16s", save = save) +
  theme(plot.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = colors, labels = c("no", "yes")) +
  scale_shape_manual(values = c(15,16), labels = c("no", "yes")) +
  labs(color = "Carbohydrase addition", shape = "Carbohydrase addition")

# ancombc2 at asv level

feature_counts_filtered <- readRDS("clean/5_16s_ferment_taxa_feature_filtered_counts.RDS") 

filtered_table <- filter_ktable(feature_counts_filtered, meta = meta, selected_rank = "low", 
                                selected_matrix = "ferment")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "ferment") %>%
  column_to_rownames("sampleid")

# perform ancombc

output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + period", 
                     p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                     taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) # lib_cut = 1000 

res_pair_list <- extract_pairwise_taxa_to_list(output = output, selected_matrix = "ferment")

figS2b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  theme(plot.title = element_blank(), 
        strip.text = element_blank()) +
  scale_size_manual(values = c(2))


figS2 <- figS2a /
  figS2b +
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(filename= "95_figS2.jpeg",
       plot = figS2,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 20,
       scale=1,
       dpi=600)

# Figure S3 metaproteomics

proteins_norm_imp_rel_abd <- readRDS("clean/6_ferment_proteins_long_norm_imp_rel_abd_filtered.RDS")

proteins_norm_imp_rel_abd_filtered <- filter(proteins_norm_imp_rel_abd, sampleid != "409")

proteins_micro_rel_abd_filtered <- proteins_norm_imp_rel_abd_filtered %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

figS3a <- do_ordination(proteins_micro_rel_abd_filtered, title = "Microbial proteins ferment", save_name = "proteingroups_micro", save = save, region = "ferment") +
  theme(plot.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = colors, labels = c("no", "yes")) +
  scale_shape_manual(values = c(15,16), labels = c("no", "yes")) +
  labs(color = "Carbohydrase addition", shape = "Carbohydrase addition")

# ancombc at genus level

taxonomy_norm_imp_intensity <- readRDS("clean/6_ferment_taxonomy_long_norm_imp_intensity_filtered.RDS")

taxonomy_norm_imp_intensity_filtered <- filter(taxonomy_norm_imp_intensity, sampleid != "409")



filtered_table <- filter_ktable(taxonomy_norm_imp_intensity_filtered, meta = meta, selected_rank = "G", 
                                selected_matrix = "ferment")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "ferment") %>%
  column_to_rownames("sampleid")

# perform ancombc

output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) # lib_cut = 1000 

res_pair_list <- extract_pairwise_taxa_to_list(output = output, selected_matrix = "ferment")

figS3b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  theme(plot.title = element_blank(), 
        strip.text = element_blank()) +
  scale_size_manual(values = c(2))

figS3 <- figS3a /
  figS3b +
  plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(filename= "95_figS3.jpeg",
       plot = figS3,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 20,
       scale=1,
       dpi=600)



# ckimchii proteins

plotting_table <- Ckimchi_proteins_400 %>%
  group_by(protein_name) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(protein_name = ifelse(rel_abd < 1, "other", protein_name)) %>%
  group_by(protein_name) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  arrange(desc(rel_abd)) %>%
  mutate(protein_name = factor(protein_name, levels = c(unique(protein_name)[-which(unique(protein_name)=="other")], "other")))

figS4a <- ggplot(plotting_table, aes(x = 1, y = rel_abd, fill = protein_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "", y = "Relative abundance (%)", fill = "Protein", title = "Companilactobacillus") +
  scale_fill_manual(values = c(colors[1:(length(unique(plotting_table$protein_name))-1)], "grey")) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0)) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "italic"),
        legend.key.size = unit(.4, "cm")) 

# wconfusa proteins

plotting_table <- Wconfusa_proteins_400 %>%
  group_by(protein_name) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(protein_name = ifelse(rel_abd < 2, "other", protein_name)) %>%
  group_by(protein_name) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  arrange(desc(rel_abd)) %>%
  mutate(protein_name = factor(protein_name, levels = c(unique(protein_name)[-which(unique(protein_name)=="other")], "other")))

figS4b <- ggplot(plotting_table, aes(x = 1, y = rel_abd, fill = protein_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "", y = "Relative abundance (%)", fill = "Protein", title = "Weissella") +
  scale_fill_manual(values = c(colors[1:(length(unique(plotting_table$protein_name))-1)], "grey")) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0)) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "italic"),
        legend.key.size = unit(.4, "cm")) 

# pea proteins

proteins_pea_rel_abd_filtered <- proteins_norm_imp_rel_abd_filtered %>%
  filter(origin == "pea") %>%
  recalculate_rel_abd()

pea_functions <- readRDS("clean/6_ferment_pea_functions_long_raw.RDS")


functions_pea_rel_abd_filtered <- proteins_pea_rel_abd_filtered %>%
  group_by(proteinid) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  mutate(entry_name = str_remove(proteinid, "^[a-z]+\\|[A-Z0-9]+\\|")) %>%
  inner_join(pea_functions, by = "entry_name") %>%
  group_by(protein_names) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  add_column(sampleid = "401") %>%
  dplyr::rename(name = protein_names) # alpha galactosidase present

plotting_table <- functions_pea_rel_abd_filtered %>%
  group_by(name) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  mutate(name = ifelse(rel_abd < 1, "other", name)) %>%
  group_by(name) %>%
  summarize(rel_abd = sum(rel_abd), .groups = "drop") %>%
  arrange(desc(rel_abd)) %>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

figS4c <- ggplot(plotting_table, aes(x = 1, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "", y = "Relative abundance (%)", fill = "Protein", title = "Pea proteins") +
  scale_fill_manual(values = c(colors[1:(length(unique(plotting_table$name))-1)], "grey")) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0)) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "plain"),
        legend.key.size = unit(.4, "cm")) 


figS4 <- figS4a /
  figS4b /
  figS4c +
  plot_layout(heights = c(1,1,1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(filename= "95_figS4.jpeg",
       plot = figS4,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 18,
       height = 20,
       scale=1,
       dpi=600)

