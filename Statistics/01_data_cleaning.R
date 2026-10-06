
library(here)

source(here("0_general_functions.R"))

# load data

meta <- read_tsv("data/meta.txt")
nmr_meta <- read_tsv("data/nmr_meta.txt") %>%
  mutate(sampleid = as.character(sampleid))

# clean meta


meta1 <- meta %>%
  mutate(across(c(sampleid, square, animal, period, diet), as.character)) %>%
  left_join(dplyr::select(nmr_meta, sampleid, ph), by = "sampleid") %>% #add ph to meta
  saveRDS("clean/meta1.RDS")

################
#NMR cleaning
################

# load nmr data

nmr_id_to_sampleid <- read_tsv("data/nmr_id_to_sampleid.txt") %>%
  dplyr::select(sampleid = "Kuh-Id", nmrid = "...6") %>%
  mutate(sampleid = str_remove(sampleid, "AP3_"),
         nmrid = str_replace(nmrid, "maeu", "maeu0"))
dry_matter <- read_tsv("data/dry_matter.txt") %>%
  mutate(sampleid = as.character(sampleid))

nmr_ferment <- read_tsv("data/nmr_ferment.txt", skip = 2) %>%
  dplyr::select(-c())


# clean metabolomics data
clean_nmr <- function(input) {
  output <- input %>%
    rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
    dplyr::rename(nmrid = ...1, "dss-d6" = "dss-d6_(chemical_shape_indicator)") %>%
    mutate(nmrid = str_remove(nmrid, "_ex1.cnx")) %>%
    inner_join(nmr_id_to_sampleid, by = "nmrid") %>%
    relocate(sampleid) %>%
    dplyr::select(-c(nmrid, "dss-d6")) %>% # remove TSP standard
    arrange(sampleid)
  return(output)
}

nmr_ferment_clean <- clean_nmr(nmr_ferment)

# calculate per FM
# original unit = mM = mmol/l = mmol/kg
# account for 0.5 mM concentration of TSP, but scaling to 5 mM -> /10
# account for addition of 60 µl buffer to 54o µl sample -> /0.9
# to calculate per kg FM -> dilution = (weighin+1200µl)/weighin
# mmol/kg * dilution
# -> mmol/kg FM

nmr_to_fm <- function(input) {
  output <- input %>%
    pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
    left_join(dplyr::select(nmr_meta, sampleid, weigh_in), by = "sampleid") %>%
    mutate(concentration = ifelse(is.na(concentration), 0, concentration), # impute missing values with 0
           concentration = concentration/10, # correct for wrong TSP concentration
           concentration = concentration/0.9, # account for buffer addition
           dilution = (1200+weigh_in)/weigh_in,
           concentration = concentration * dilution) %>%
    dplyr::select(-c(dilution, weigh_in)) %>%
    pivot_wider(names_from = "metabolite", values_from = "concentration")
  return(output)
}

nmr_ferment_fm <- nmr_to_fm(nmr_ferment_clean) %>%
  dplyr::select(-citrate)
saveRDS(nmr_ferment_fm, "clean/7_nmr_ferment_fm.RDS")
write_tsv(nmr_ferment_fm, "tables/7_nmr_ferment_fm.txt")

# calculate per dm
# /(DM%/100)
# -> mmol/g DM

nmr_to_dm <- function(input) {
  output <- input %>%
    pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
    left_join(dry_matter, by = "sampleid") %>%
    mutate(concentration = concentration / (dry_matter/100)) %>%
    dplyr::select(-dry_matter) %>%
    pivot_wider(names_from = "metabolite", values_from = "concentration")
  return(output)
}

nmr_ferment_dm <- nmr_to_dm(nmr_ferment_fm)
saveRDS(nmr_ferment_dm, "clean/7_nmr_ferment_dm.RDS")
write_tsv(nmr_ferment_dm, "tables/7_nmr_ferment_dm.txt")
