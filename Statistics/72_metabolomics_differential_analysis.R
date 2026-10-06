library(here)
library(venn)

source(here("50_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

dry_matter <- read_tsv("data/dry_matter.txt") %>%
  mutate(sampleid = as.character(sampleid))

meta <- readRDS("clean/meta1.RDS") %>%
  left_join(dry_matter, by = "sampleid")

# load nmr files

nmr_ferment_dm <- readRDS("clean/7_nmr_ferment_dm.RDS")
nmr_ferment_fm <- readRDS("clean/7_nmr_ferment_fm.RDS")


# saveswitch

save = FALSE
# save = TRUE

# prepare dfs for analysis function

fer_dm <- prepare_metabolomics(nmr_ferment_dm)
fer_fm <- prepare_metabolomics(nmr_ferment_fm)


# comparison of two ferments


fer_dm_table <- loop_comparison_metabolomics_ferment(fer_dm, y_axis = "mmol/kg DM", table_title = "NMR ferment DM", 
                                             save_name = "ferment_dm", save = save)
write_tsv(fer_dm_table, "tables/72_metabolomics_ferment_dm.txt")

fer_fm_table <- loop_comparison_metabolomics_ferment(filter(fer_fm, diet != "P3D7"),
                                                     y_axis = "mmol/kg FM", table_title = "NMR ferment FM", 
                                                     save_name = "ferment_fm", save = save, digits = 2)
write_tsv(fer_fm_table, "tables/72_metabolomics_ferment_fm.txt")

