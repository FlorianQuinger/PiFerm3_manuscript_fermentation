source("0_general_functions.R")

library(patchwork)
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

write_tsv(as_tibble(ec_matrix, rownames = "species"), "tables/14_ec_matrix.txt")

ec_matrix %>% as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  group_by(ec_number) %>%
  summarise(sum = sum(value), .groups = "drop") %>%
  mutate(abu = ifelse(sum > 0, 1, 0)) %>%
  pull(abu) %>%
  sum()

ec_plot <- ec_matrix %>% 
  as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  group_by(ec_number) %>%
  filter(sum(value)>0) %>%
  ungroup() %>%
  mutate(ec_number = factor(ec_number, levels = unique(ec_table$EC))) %>%
  ggplot(aes(x = ec_number, y = species, fill = value)) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = T) +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical") +
  scale_fill_gradient(low = "white", high = "black", breaks = c(0,3,6,9,12)) +
  labs(fill = "Number of genes") 

save_big("14_ec_matrix_print", width = 18, height = 22)

legend <- ec_matrix %>% as_tibble(rownames = "species") %>%
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
        legend.key.size = unit(.4, "cm")) +
  guides(fill = guide_legend(position = "bottom", ncol = 2, title = "Substrate", theme(legend.title.position = "top")))

p <- legend / ec_plot +
  plot_layout(heights = c(.5,22), guides = "collect") &
  theme(legend.position = "bottom")
print(p)

save_big("14_ec_matrix_print", width = 18, height = 22)

##################

ec_matrix %>% as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  mutate(ec_number = paste(substrate, ec_number)) %>%
  group_by(ec_number) %>%
  filter(sum(value)>0) %>%
  ungroup() %>%
  ggplot(aes(x = factor(ec_number), y = species, fill = ifelse(value >1, log2(value), value))) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 5)) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "", y = "")

save_big("13_ec_matrix_combined")

ec_matrix %>% as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value")

ranking <- ec_matrix %>% as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  #filter(!substrate %in% c("starch", "maltose", "dextrin")) %>%
  filter(str_detect(substrate,
                     "alpha-galactosides|sucrose|cellulose|cellobiose|cellodextrin|pectin|pectate|digalacturonate|xylan")) %>%
  group_by(species) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  add_column(substrate = "score")

ec_matrix %>% as_tibble(rownames = "species") %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  group_by(species, substrate) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  #add_row(ranking) %>%
  ggplot(aes(x = factor(substrate, levels = c(unique(ec_table$substrate), "score")), y = species, fill = ifelse(value >1, log2(value), value))) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 5)) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "", y = "")

save_big("13_ec_matrix_substrates_combined")

ec_matrix %>% as_tibble(rownames = "species") %>%
  filter(species %in% c("Lactiplantibacillus plantarum a10_flongle",
                        "Lactiplantibacillus plantarum a2_flongle", 
                        "Pediococcus pentosaceus a3_flongle",
                        "Bacillus licheniformis i32_dragonflye",
                        "Weissella confusa dsm20196",
                        "Levilactobacillus spicheri dsm15429",
                        "Leuconostoc suionicum dsm20241",
                        "Companilactobacillus kimchii dsm13961",
                        "Bacillus subtilis dsm3257")) %>%
  mutate(species = factor(species, levels = c("Bacillus subtilis dsm3257", "Leuconostoc suionicum dsm20241", "Pediococcus pentosaceus a3_flongle", "Lactiplantibacillus plantarum a2_flongle", "Bacillus licheniformis i32_dragonflye",  "Weissella confusa dsm20196", "Levilactobacillus spicheri dsm15429", "Companilactobacillus kimchii dsm13961", "Lactiplantibacillus plantarum a10_flongle"))) %>%
  pivot_longer(-species, names_to = "ec_number", values_to = "value") %>%
  left_join(distinct(dplyr::select(ec_table, EC, substrate)), by = c("ec_number" = "EC")) %>%
  group_by(species, substrate) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  #add_row(ranking) %>%
  ggplot(aes(x = factor(substrate, levels = c(unique(ec_table$substrate), "score")), y = species, fill = ifelse(value >1, log2(value), value))) +
  geom_tile(color = "black", lwd = 0.25, linetype = 1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 10)) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "", y = "")

save_big("13_ec_matrix_substrates_combined_strains")
