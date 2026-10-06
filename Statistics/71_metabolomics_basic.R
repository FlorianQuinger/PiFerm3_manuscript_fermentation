library(here)
library(venn)

source(here("0_general_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

meta <- readRDS("clean/meta1.RDS") %>%
  add_row(sampleid = "409", diet = "X")

# load nmr files

nmr_ferment <- readRDS("clean/7_nmr_ferment_dm.RDS")

# save flag
save = F
#save = T

# create venn diagramms

ferment_abu <- nmr_ferment %>%
  pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
  group_by(metabolite) %>%
  summarise(abu = 1, .groups = "drop") %>%
  mutate(matrix = "ferment")

# plots

# barplot of metabolite concentrations

nmr_ferment %>%
  pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
  group_by(metabolite) %>%
  mutate(mean = mean(concentration),
         metabolite = ifelse(mean <= 10, "other", metabolite)) %>%
  ungroup() %>%
  ggplot(aes(x = sampleid, y = concentration, fill = metabolite)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0,1990), expand = c(0,0)) +
  labs(title = "Metabolite concentrations in ferment samples", y = "mmol/kg DM")

save_big("71_nmr_ferment_stacked_barplot")

# barplot

nmr_ferment %>%
  pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
  left_join(dplyr::select(meta, sampleid, diet), by = "sampleid") %>%
  group_by(metabolite, diet) %>%
  mutate(mean_con = mean(concentration), 
          sd_con = sd(concentration),
         upper = mean_con + sd_con, 
         lower = mean_con - sd_con) %>%
  ungroup() %>%
  dplyr::select(-sampleid, -concentration) %>%
  distinct() %>%
  ggplot(aes(x = metabolite, y = mean_con, group = diet)) +
  geom_bar(aes(fill = diet), stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9), width = 0.1) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_log10() +
  scale_fill_manual(values = colors) +
  labs(title = "Metabolite concentrations in ferments", y = "log10(mmol/kg DM)", x = "") +
  theme(axis.text.x = element_text(size = 10))

save_big("71_nmr_ferment_dodged_barplot")
