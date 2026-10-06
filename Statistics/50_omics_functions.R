library(vegan)
library(ANCOMBC)
library(patchwork)
library(ROTS)
library(edgeR)
library(MicrobiomeProfiler)
library(clusterProfiler)
library(org.Ss.eg.db)
library(KEGGREST)
library(ggrepel)
library(pheatmap)
library(SmCCNet)
library(mixOmics)
library(RCy3)
library(igraph)
library(DESeq2)
library(psych)
library(corrplot)
library(ggplotify)

source(here("0_general_functions.R"))

##################################
# Data cleaning functions
##################################

# sum up to the respective ranks function for data cleaning, to create kreport style table from mmseqs2

calculate_rank_abundance <- function(df, sum_column) {
  ranks = c("R1", "P", "C", "O", "F", "G", "S")
  rank_list <- list()
  for (i in 1:length(ranks)) {
    rank <- ranks[i]
    filtered_df <- df %>%
      dplyr::select(sampleid, !!sym(rank), !!sym(sum_column)) %>%
      filter(!is.na(!!sym(rank))) %>%
      group_by(sampleid, !!sym(rank)) %>%
      summarize(!!sym(sum_column) := sum(!!sym(sum_column)), .groups = "drop") %>%
      pivot_longer(-c(sampleid, !!sym(sum_column)), names_to = "rank", values_to = "name") %>%
      pivot_wider(names_from = sampleid, values_from = !!sym(sum_column)) %>% # introducing NAs that will be replaced by zeros
      pivot_longer(-c(rank, name), names_to = "sampleid", values_to = quo_name(sum_column)) %>%
      mutate(!!sym(sum_column) := ifelse(is.na(!!sym(sum_column)), 0, !!sym(sum_column))) %>% # replace NA with 0
      dplyr::select(name, rank, sampleid, !!sym(sum_column)) %>% # reorder similar to k2report
      {if (sum_column == "rel_abd") group_by(., sampleid) %>% mutate(., rel_abd=rel_abd/sum(rel_abd)*100) %>% ungroup(.) else .} # if rel_abd is choosen, recalculate relative abundances
    
    rank_list[[i]] <- filtered_df
  }
  kstyle_df <- bind_rows(rank_list)
  
  return(kstyle_df)
}

# filter out features with low frequency, takes dataframe with sampleid, abundance, features, and optional rank
# grouping cols should be cols used for frequency calculation -> features and optional rank
# abundance column is rel_abd or reads

filter_by_frequency <- function(df, grouping_cols, abundance_column, frequency_cutoff = 1/3) {
  df_freq <- df %>%
    mutate(abu = ifelse(!!sym(abundance_column) > 0, 1, 0)) %>%
    group_by(!!!syms(grouping_cols)) %>% # syms accepts vector or single string
    mutate(frequency = sum(abu)/length(abu)) %>%
    ungroup() 
  
  p <- ggplot(df_freq, aes(x = frequency)) +
    geom_histogram() +
    #scale_x_log10() +
    geom_vline(xintercept = frequency_cutoff, col = "red", lwd = 1) +
    labs(title = "distribution of frequency values", x = "relative abundance [%/100]")
  print(p)
  
  df_freq_filtered <- df_freq%>%
    filter(frequency > frequency_cutoff) %>%
    dplyr::select(-c(abu, frequency))
  
  print(str_c("Frequency filtering kept", length(unique(df_freq_filtered[[1]])), "of", 
              length(unique(df_freq[[1]])), "features, discarding",
              sum(filter(df_freq, frequency < frequency_cutoff)[[abundance_column]]) 
              / sum(df_freq[[abundance_column]]) *100, 
              "% of", abundance_column,
              sep = " "))
  
  return(df_freq_filtered)
}

# filter out low abundant features, takes dataframe with sampleid, abundance, features, and optional rank
# grouping cols should be cols used for abundance average calculation -> features and optional rank
# abundance column is rel_abd or reads


filter_by_abundance <- function(df, grouping_cols, abundance_column, rank_column = NULL, abundance_cutoff = 0.01) {
  abundance_cutoff <- abundance_cutoff / 100 # convert % into number
  
  df_abu <- df %>%
    {if (is.null(rank_column)) group_by(., sampleid) else group_by(., sampleid, !!sym(rank_column)) } %>%
    mutate(sum = sum(!!sym(abundance_column))) %>%
    ungroup() %>%
    mutate(rel_abd2 = !!sym(abundance_column) / sum) %>%
    group_by(!!!syms(grouping_cols)) %>%
    mutate(avg = mean(rel_abd2)) %>%
    ungroup()
  
  p <- ggplot(df_abu, aes(x = avg)) +
    geom_histogram() +
    scale_x_log10() +
    geom_vline(xintercept = abundance_cutoff, col = "red", lwd = 1) +
    labs(title = "distribution of relative abundance values", x = "relative abundance [%/100]")
  print(p)
  
  df_abu_filtered <- df_abu %>%
    filter(avg > abundance_cutoff) %>% # filtering out features that are in avg lower than the cutoff relative to the sum of abundance at this level
    dplyr::select(-c(sum, rel_abd2, avg))
  
  print(str_c("Abundance filtering kept", length(unique(df_abu_filtered[[1]])), "of", 
              length(unique(df_abu[[1]])), "features, discarding",
              sum(filter(df_abu, avg < abundance_cutoff)$avg) / sum(df_abu$avg) *100, 
              "% of", abundance_column,
              sep = " "))
  
  return(df_abu_filtered)
}

# recalculates rel_abd after some features were filtered out, always calculates to 100 % at each level
# rank column if ktable style is included for grouping

recalculate_rel_abd <- function(df, rank_column = NULL) {
  df_recalc <- df %>%
    {if (is.null(rank_column)) group_by(., sampleid) else group_by(., sampleid, !!sym(rank_column)) } %>%
    mutate(rel_abd = rel_abd / sum(rel_abd) * 100) %>%
    ungroup()
  return(df_recalc)
}

# function combining abundance and frequency filtering for ktable style table and optionally recalculates rel_abd

filter_frequency_and_abundance_ktable <- function(ktable, frequency_cutoff = 1/3, abundance_cutoff = 1e-3) {
  # detect columns: should be features, rank, sampleid, abundance
  feature_column <- colnames(ktable)[1]
  rank_column <- colnames(ktable)[2]
  abundance_column <- colnames(ktable)[4]
  
  print(paste("Feature column is", feature_column))
  print(paste("Rank column is called", rank_column)) 
  print(paste("Abundance column is called", abundance_column))
  
  # process separately for ileum and faeces
  ktable_il <- filter_ileum(ktable)
  ktable_fa <- filter_faeces(ktable)
  
  # filter by frequency
  ktable_il_freq <- filter_by_frequency(df = ktable_il, grouping_cols = c(feature_column, rank_column),
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  ktable_fa_freq <- filter_by_frequency(df = ktable_fa, grouping_cols = c(feature_column, rank_column),
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  
  #filter out low abundant features
  ktable_il_abu <- filter_by_abundance(df = ktable_il_freq, grouping_cols = c(feature_column, rank_column),
                                       abundance_column = abundance_column, rank_column = rank_column,
                                       abundance_cutoff = abundance_cutoff)
  ktable_fa_abu <- filter_by_abundance(df = ktable_fa_freq, grouping_cols = c(feature_column, rank_column),
                                       abundance_column = abundance_column, rank_column = rank_column,
                                       abundance_cutoff = abundance_cutoff)
  
  # recalculate rel_abd if not reads
  if (abundance_column == "rel_abd") {
    ktable_il_out <- recalculate_rel_abd(df = ktable_il_abu, rank_column = rank_column)
    ktable_fa_out <- recalculate_rel_abd(df = ktable_fa_abu, rank_column = rank_column)
  } else {
    ktable_il_out <- ktable_il_abu
    ktable_fa_out <- ktable_fa_abu
  }
  
  ktable_out <- rbind(ktable_il_out, ktable_fa_out)
  
  return(ktable_out)
}

# combines frequency and abundance filtering for table without ranks
# filter if ileum and faeces should be separated and processed separately, False for ferment

filter_frequency_and_abundance <- function(df, feature_columns = NULL,
                                           frequency_cutoff = 1/3, abundance_cutoff = 1e-3,
                                           filter_df = TRUE, all_matrices = FALSE) {
  # detect columns: should be features, sampleid, abundance
  if (is.null(feature_columns)) { # if not manually defined, try to find 
    feature_column <- colnames(df)[1]
    abundance_column <- colnames(df)[3]
  } else {
    feature_column <- feature_columns
    abundance_column <- colnames(df)[length(feature_columns)+2] # number of feature columns determines position of abundance
  }
  
  print(paste("Feature column is", feature_column))
  print(paste("Abundance column is called", abundance_column))
  
  if (isTRUE(all_matrices)) {
    # process separately for ileum and faeces and ferment
    df_il <- filter_ileum(df)
    df_fa <- filter_faeces(df)
    df_fe <- filter_ferment(df)
    
    # filter by frequency
    df_il_freq <- filter_by_frequency(df = df_il, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    df_fa_freq <- filter_by_frequency(df = df_fa, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    df_fe_freq <- filter_by_frequency(df = df_fe, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    
    #filter out low abundant features
    df_il_abu <- filter_by_abundance(df = df_il_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    df_fa_abu <- filter_by_abundance(df = df_fa_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    df_fe_abu <- filter_by_abundance(df = df_fe_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    
    # recalculate rel_abd if not reads or intensity
    if (abundance_column == "rel_abd") {
      df_il_out <- recalculate_rel_abd(df_il_abu)
      df_fa_out <- recalculate_rel_abd(df_fa_abu)
      df_fe_out <- recalculate_rel_abd(df_fe_abu)
    } else {
      df_il_out <- df_il_abu
      df_fa_out <- df_fa_abu
      df_fe_out <- df_fe_abu
    }
    
    df_out <- rbind(df_il_out, df_fa_out, df_fe_out)
  } else if (isTRUE(filter_df)) {
    # process separately for ileum and faeces
    df_il <- filter_ileum(df)
    df_fa <- filter_faeces(df)
    
    # filter by frequency
    df_il_freq <- filter_by_frequency(df = df_il, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    df_fa_freq <- filter_by_frequency(df = df_fa, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    
    #filter out low abundant features
    df_il_abu <- filter_by_abundance(df = df_il_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    df_fa_abu <- filter_by_abundance(df = df_fa_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    
    # recalculate rel_abd if not reads or intensity
    if (abundance_column == "rel_abd") {
      df_il_out <- recalculate_rel_abd(df_il_abu)
      df_fa_out <- recalculate_rel_abd(df_fa_abu)
    } else {
      df_il_out <- df_il_abu
      df_fa_out <- df_fa_abu
    }
    
    df_out <- rbind(df_il_out, df_fa_out)
    
  } else {
    # filter by frequency
    df_freq <- filter_by_frequency(df = df, grouping_cols = feature_column,
                                      abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
    
    #filter out low abundant features
    df_abu <- filter_by_abundance(df = df_freq, grouping_cols = feature_column,
                                     abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
    
    # recalculate rel_abd if not reads or intensity
    if (abundance_column == "rel_abd") {
      df_out <- recalculate_rel_abd(df_abu)
    } else {
      df_out <- df_abu
    }
  }
  
  return(df_out)
}

# calculate abundances of single ko's from an eggnog style table

calculate_kegg_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), kegg_ko) %>%
    filter(kegg_ko != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(kegg_ko, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(kegg_ko)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, kegg_ko) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    mutate(kegg_ko = str_remove(kegg_ko, "ko:")) %>%
    dplyr::select(kegg_ko, sampleid, !!sym(abundance_column))
  return(output)
}

# calculate abundances of single go's from an eggnog style table, go column should be named "gos"

calculate_go_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), gos) %>%
    filter(gos != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(gos, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(gos)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, gos) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(go = gos, sampleid, !!sym(abundance_column))
  return(output)
}

# caclulate abundance of single cogs from an eggnog style table, cog column should be called cog_accession

calculate_cog_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), cog_accession) %>%
    filter(!is.na(cog_accession)) %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(cog_accession, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(cog_accession)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, cog_accession) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(cog_accession, sampleid, !!sym(abundance_column))
  return(output)
}

# caclulate abundance of cog categories from an eggnog style table, cog column should be called cog_accession

calculate_cog_cat_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), cog_category) %>%
    filter(!is.na(cog_category)) %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(cog_category, sep = "") %>% # split after each character
    filter(str_detect(cog_category, "[A-Z]")) %>% #remove empty rows
    group_by(id) %>%
    mutate(factor = 1 / length(cog_category)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, cog_category) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(cog_category, sampleid, !!sym(abundance_column))
  return(output)
}

cog_cat_table <- tribble(
  ~cog_category, ~description,
  "J", "Translation, ribosomal structure and biogenesis",
  "A", "RNA processing and modification",
  "K", "Transcription",
  "L", "Replication, recombination and repair",
  "B", "Chromatin structure and dynamics",
  "D", "Cell cycle control, cell division, chromosome partitioning",
  "Y", "Nuclear structure",
  "V", "Defense mechanisms",
  "T", "Signal transduction mechanisms",
  "M", "Cell wall/membrane/envelope biogenesis",
  "N", "Cell motility",
  "Z", "Cytoskeleton",
  "W", "Extracellular structures",
  "U", "Intracellular trafficking, secretion, and vesicular transport",
  "O", "Posttranslational modification, protein turnover, chaperones",
  "C", "Energy production and conversion",
  "G", "Carbohydrate transport and metabolism",
  "E", "Amino acid transport and metabolism",
  "F", "Nucleotide transport and metabolism",
  "H", "Coenzyme transport and metabolism",
  "I", "Lipid transport and metabolism",
  "P", "Inorganic ion transport and metabolism",
  "Q", "Secondary metabolites biosynthesis, transport and catabolism",
  "R", "General function prediction only",
  "S", "Function unknown"
) %>%
  mutate(cat_description = paste(cog_category, "-", description))


# calculate abundance of cazymes from eggnog style table

calculate_cazy_abundance<- function(input, abundance_column, cazy_column = "cazy") {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), !!sym(cazy_column)) %>%
    filter(!!sym(cazy_column) != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(!!sym(cazy_column), sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(!!sym(cazy_column))) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, !!sym(cazy_column)) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(!!sym(cazy_column), sampleid, !!sym(abundance_column))
  return(output)
}

# calculate abundance of EC's from eggnog style table

calculate_ec_abundance<- function(input, abundance_column, ec_column = "ec_id") {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), !!sym(ec_column)) %>%
    filter(!!sym(ec_column) != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(!!sym(ec_column), sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(!!sym(ec_column))) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, !!sym(ec_column)) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(!!sym(ec_column), sampleid, !!sym(abundance_column))
  return(output)
}

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

####################################
# data analysis functions
##################################

# function to do pcoa

create_distance_matrix <- function(input_df, region = NULL, dissimilarity_index = "bray", 
                                   abundance_column = "rel_abd") {
  first_column <- colnames(input_df)[1]
  
  if (is.null(region)) {
    print("Warning uncommon columns will be excluded")
  } else if (region == "ileum") {
    input_df <- filter(input_df, str_detect(sampleid, "^1"))
  } else if (region == "faeces") {
    input_df <- filter(input_df, str_detect(sampleid, "^3"))
  } else if (region == "ferment") {
    input_df <- filter(input_df, str_detect(sampleid, "^4"))
  }
  
  mat <- input_df %>%
    dplyr::select(!!sym(first_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = first_column, values_from = abundance_column) %>%
    dplyr::select(-sampleid) %>%
    as.matrix()
  
  labels <- input_df %>%
    dplyr::select(!!sym(first_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = first_column, values_from = abundance_column) %>%
    pull(sampleid)
  
  rownames(mat) <- labels
  
  bray <- vegdist(mat, method = dissimilarity_index, na.rm = T)
  
  return(bray)
}

create_pcoa <- function(bray, meta = meta, title = "", second_indicator = "diet") {
  
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  
  pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
  pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")
  
  x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
    inner_join(meta, by = "sampleid")
  
  # calculate centroids
  x_centers <- x %>%
    group_by(matrix, diet) %>%
    summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")
  
  x_star <- x %>%
    group_by(diet) %>%
    mutate(c1 = mean(V1), c2 = mean(V2)) %>%
    ungroup()
  
  p <- ggplot(x) +
    geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                               #color = diet, 
                                               shape = diet, 
                                               fill = matrix),
               show.legend = FALSE, size = 8, alpha = 0.5) +
    #geom_segment(data = x_star, mapping = aes(x = V1, y = V2,
    #                                         xend = c1, yend = c2),
    #           show.legend = F, lwd = 1.5, alpha = 0.2) +
    geom_point(aes(x = V1, y = V2, shape = diet, color = !!sym(second_indicator)), size = 4) +
    #geom_label(aes(label = sampleid)) +
    #geom_label(aes(label = sampling_time)) +
    #geom_label(aes(label = animal)) +
    stat_ellipse(aes(x = V1, y = V2, shape = diet, color = diet, fill = matrix),
                 level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
    labs(x = pco1, y = pco2, title = title) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(15,16,17,18,21,22,23,24))
  
  return(p)
}

# function to create permutation matrix, needed when observations are unbalanced, expects meta with period column which is used for permuting

create_permutation_matrix <- function(bray, meta = meta, region, treatment_col, design = c("double_latin_square", "row_column")) {
  period_no <- length(unique(meta$period))
  # meta with all sampleids of one region
  if (region == "ileum") {
    filtered_meta <- filter_ileum(meta)
  } else if (region == "faeces") {
    filtered_meta <- filter_faeces(meta)
  } else if (region == "ferment") {
    filtered_meta <- filter_ferment(meta) %>%
      filter(sampleid != "409")
  }
  # meta with sampleids in bray
  reduced_meta <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  
  perms <- allPerms(period_no, control = how(observed = T)) # possible permutations if only periods are permuted
  double_perms <- expand.grid(square1 = seq_len(nrow(perms)),square2 = seq_len(nrow(perms)))[-1,] # excluding observed perm
  perms_double <- matrix(ncol = period_no*2, nrow = nrow(double_perms)) # all possible combinations of perms for double latin square
  for (i in 1:nrow(double_perms)) {
    double_perm <- double_perms[i,]
    perm <- c(perms[double_perm[[1]],], perms[double_perm[[2]],]) # concatenate single perms for double perms
    perms_double[i,] <- perm
  }
  
  if (design == "row_column" || region == "ferment") {
    perms <- perms[-1,] # exclude observed column from perms
  } else if (design == "double_latin_square") {
    perms <- perms_double
  }
  
  
  permutations <- matrix(ncol = nrow(reduced_meta), nrow = nrow(perms))
  for (i in 1:nrow(perms)) {
    # permute period according to perms
    if (design == "row_column" || region == "ferment") {
      # permute identical for both squares
      permutation <- c(which(filtered_meta$period == perms[i, 1]),
                       which(filtered_meta$period == perms[i, 2]),
                       which(filtered_meta$period == perms[i, 3]),
                       which(filtered_meta$period == perms[i, 4])) 
    } else if (design == "double_latin_square") {
      # normal order is s1p1, s2p1, s1p2, s2p2, s1p3, s2p3, s1p4, s2p4
      permutation <- c(which(filtered_meta$period == perms[i, 1] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 5] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 2] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 6] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 3] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 7] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 4] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 8] & filtered_meta$square == "2"))
    } else {print("no design selected")}
    if (nrow(filtered_meta) != nrow(reduced_meta)) {
      missing_rows <- setdiff(filtered_meta, reduced_meta)
      for (j in nrow(missing_rows)) {
        missing_row <- missing_rows[j,]
        missing_position <- which(filtered_meta$sampleid == missing_row$sampleid) # find position where row was removed
        permutation <- permutation[-missing_position] # remove removed position from permutation
        permuted_position <- which(permutation == missing_position) # detect where missing position was permuted
        missing_label <- filtered_meta[[treatment_col]][missing_position] # detect which label missing position had
        permutation[permuted_position] <- setdiff(which(filtered_meta[[treatment_col]] == missing_label), missing_position)[1] # replace permuted position with identical label from other position
        permutation <- ifelse(permutation > missing_position, permutation - 1, permutation)# substract 1 from each position > than missing position 
      }
    }
    permutations[i,] <- permutation
  }
  
  # for row column design and balanced design, matrix should be the same as
  #ctrl <- with( #main block, never permute between blocks, only within
  #test_meta, how(plots = Plots(strata = animal),  #plot, permutation defined by within argument
  #                           within = Within(type = "free", constant = T))) #constant leads to same within plot permutation for all plots
  
  return(permutations)
}

# do a betadisper test

do_betadisper <- function(bray, meta = meta, permutation_matrix) {
  #filter meta
  x <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  
  # perform betadisper
  betadis <- betadisper(bray, group= x$diet, type = "centroid")
  permbetadis <- permutest(betadis, permutations = permutation_matrix, pairwise = T)
  
  return(list(betadis,permbetadis))
}

# plot betadispersion as boxplots

plot_betadisper <- function(betadis, meta = meta, title = NULL) {
  # extract permuation p values
  treatment_comparisons <- betadis[[2]]$pairwise$permuted %>%
    as_tibble(rownames = "comparison") %>%
    separate(comparison, into = c("comp1", "comp2"))
  
  # create p matrix
  p_mat <- matrix(nrow = length(unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))),
                  ncol = length(unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))))
  colnames(p_mat) <- unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))
  rownames(p_mat) <- unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))
  diag(p_mat) <- 1
  
  for (i in 1:nrow(treatment_comparisons)) {
    p_mat[treatment_comparisons$comp1[i], treatment_comparisons$comp2[i]] <- treatment_comparisons$value[i]
    p_mat[treatment_comparisons$comp2[i], treatment_comparisons$comp1[i]] <- treatment_comparisons$value[i]
  }
  
  # get cld
  cld <- generate_cld(p_mat) %>%
    dplyr::rename(diet = treatments) %>%
    mutate(cld = str_remove(treatments_cld, "\\d"))
  
  # get distances
  dispersion <- betadis[[1]]$distances %>%
    as_tibble(rownames = "sampleid") %>%
    left_join(meta, by = "sampleid") %>%
    left_join(cld, by = "diet") %>%
    group_by(diet) %>%
    mutate(mean_dis = mean(value)) %>%
    ungroup()
  
  P <- betadis[[2]]$tab$`Pr(>F)`[1]
  Print <- round_p_value(P)
  
  p <- ggplot(dispersion, aes(x = diet, y = value, fill = diet)) +
    geom_boxplot(width = .5, outliers = F, show.legend = F, alpha = .4) +
    geom_quasirandom(width = .2, size=2, show.legend = F, pch = 21, alpha = .5) +
    {if (P < 0.05) geom_text(aes(y = mean_dis, label = cld))} +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(y = "Distance to centroid", title = title, 
         subtitle = paste("Permuation test: P", Print))
  
  return(p)
}

do_permanova <- function(bray, meta = meta, permutation_matrix) {
  x <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  # perform PERMANOVA
  
  permanova <- adonis2(bray~diet, data = x, permutations = permutation_matrix)
  
  
  ############ permanova with full formula and without permutation design
  #permanova <- adonis2(bray ~ diet + animal + period, data = x, permutations = 9999, by = "margin")
  
  return(permanova)
}

#generate cld out of matrix and output treatments with cld letters
generate_cld <- function(p_matrix) {
  # insert and absorb algorithm #
  ###############################
  # Get treatment names
  treatments <- rownames(p_matrix)
  
  # find significant differences
  sig <- list()
  counter = 1
  for (i in 1:length(treatments)) { # the treatment to compare with
    for (j in 1:length(treatments)) { # the treatment of interest
      if (p_matrix[i,j] < 0.05) { # if significant different, add to group
        sig[[counter]] <- sort(c(treatments[[i]], treatments[j]))
        counter = counter + 1
      }
    }
  }
  
  # remove identical comparisons
  sig <- unique(sig)
  
  # initialize dataframe
  df <- data.frame(treatments = treatments, column1 = 1)
  
  if (length(sig) > 0) {
    for (i in 1:length(sig)) {
      t1 <- sig[[i]][1]
      t2 <- sig[[i]][2]
      pos_t1 <- which(treatments == t1)
      pos_t2 <- which(treatments == t2)
      # insert
      df_inserted <- df[1]
      for (j in 2:ncol(df)) {
        col = df[,j]
        if (col[pos_t1] == 1 & col[pos_t2] == 1) {
          col1 = col
          col1[pos_t1] <- 0 # replace with 0
          col2 = col
          col2[pos_t2] <- 0 # replace with 0
          df_inserted <- cbind(df_inserted, col1, col2)
        } else {
          df_inserted <- cbind(df_inserted, col) # if no second column needs to be created
        }
      }
      # absorb
      df_absorbed <- df[1]
      for (j in 2:ncol(df_inserted)) {
        contained = 0
        col1 = df_inserted[,j]
        for (k in 2:ncol(df_inserted)) {
          if (j == k) {
            next
          } else {
            col2 = df_inserted[,k]
            # check if col1 is contained in col2
            contained_mini = 1
            for (l in 1:length(col1)) { # iterate over column
              if (col1[l] == 1) { # check only 1 positions
                if (col2[l] == 0) { # if col2 is not 1 at the same position
                  contained_mini = 0 # contained is set FALSE
                }
              }
            }
            # if column is contained, set global contained to TRUE
            if (contained_mini == 1) {
              contained = 1
            }
          }
        }
        if (contained == 0) {
          df_absorbed = cbind(df_absorbed, col1)
        }
      }
      df <- df_absorbed # to start loop again
    }
  }
  
  names(df)[2:ncol(df)] <- sort(letters[1:(ncol(df))-1], decreasing = T)
  
  cld <- df[, c(1,sort(c(2:ncol(df)), decreasing = T))]
  # replace 0 and 1 by letters
  for (i in 2:ncol(cld)) {
    colname = names(cld)[i]
    for (j in 1:nrow(cld)) {
      cld[j,i] <- ifelse(cld[j,i] == 1, colname, "")
    }
  }
  
  # combine treatments with letters
  treatments_cld <- c()
  for (i in 1:nrow(cld)) {
    concat <- str_c(cld[i,],collapse = "")
    treatments_cld[i] <- concat
  }
  
  return(data.frame(treatments, treatments_cld))
}

pairwise_permanova <- function(bray, meta = meta, treatment = "diet", region, design) {
  meta_filtered <- filter(meta, sampleid %in% rownames(bray)) # remove rows with other matrices and other treatments
  treatments <- sort(unique(meta_filtered[[which(colnames(meta_filtered) == treatment)]]))
  
  #create frame to save values
  p_mat <- matrix(data = 1, nrow = length(treatments), ncol = length(treatments), 
                  dimnames = list(treatments,treatments))
  
  #generate different comparison possibilities
  comparisons <- tibble("comb1" = as.character(), "comb2" = as.character())
  for (i in 1:length(treatments)) {
    comb1 <- treatments[i]
    for (j in 1:length(treatments)) {
      comb2 <- treatments[j]
      # "only upper side of the matrix"...more or less
      if (j > i) {
        comp <- tibble("comb1" = comb1, "comb2" = comb2)
        comparisons <- comparisons %>%
          add_row(comp)
      }
    }
  }
  
  bray_matrix <- as.matrix(bray)
  
  comparisons <- add_column(comparisons, p_value = NA, p_adj = NA)
  
  for (i in 1:nrow(comparisons)) {
    current_comp1 <- comparisons$comb1[i]
    current_comp2 <- comparisons$comb2[i]
    
    # filter relevant sampleids
    if (region == "ileum") {
      meta_filtered <- filter_ileum(meta) %>%
        filter(!!sym(treatment) %in% c(current_comp1, current_comp2))
    } else if (region == "faeces") {
      meta_filtered <- filter_faeces(meta) %>%
        filter(!!sym(treatment) %in% c(current_comp1, current_comp2))
    }
    
    bray_filtered <- as.dist(bray_matrix[rownames(bray_matrix) %in% meta_filtered$sampleid,
                                         colnames(bray_matrix) %in% meta_filtered$sampleid])
    
    permutation_matrix <- create_permutation_matrix(bray = bray_filtered, meta = meta_filtered, region = region,
                                                    treatment_col = treatment, design = design) # not in use right now
    
    meta_reduced <- as_tibble(as.matrix(bray_filtered), rownames = "sampleid") %>%
      dplyr::select(sampleid) %>%
      inner_join(meta_filtered, by = "sampleid")
    
    permanova <- do_permanova(bray_filtered, meta = meta_reduced, permutation_matrix = 
                                with(meta_reduced, how(plots = Plots(strata = animal), 
                                                       within = Within(type = "free", constant = F))))
    
    comparisons$p_value[i] <- permanova$`Pr(>F)`[1]
  }
  
  # adjust p value using BH
  #comparisons$p_adj <- p.adjust(comparisons$p_value, method = "BH")
  comparisons$p_adj <- comparisons$p_value # unadjusted due to errors
  
  # fill matrix
  for (i in 1:nrow(comparisons)) {
    current_comp1 <- comparisons$comb1[i]
    current_comp2 <- comparisons$comb2[i]
    p_mat[which(rownames(p_mat) == current_comp1),
          which(colnames(p_mat) == current_comp2)] <- comparisons$p_adj[i]
  }
  
  # mirror matrix
  p_mat[lower.tri(p_mat)] <- t(p_mat)[lower.tri(p_mat)] # transpose to ensure right ordering
  print(p_mat)
  
  # add cld
  cld <- generate_cld(p_mat)
  
  return(cld)
}

do_ordination <- function(input_df, region = NULL, title = "", save = FALSE, save_name = NULL,
                          second_indicator = "diet", dissimilarity_index = "bray", abundance_column = "rel_abd") {
  
  # check if input_df is already a distance matrix
  if ("dist" %in% class(input_df)) {
    bray <- input_df
  } else {
    bray <- create_distance_matrix(input_df, region = region, dissimilarity_index = dissimilarity_index,
                                   abundance_column = abundance_column)
  }
  
  if (!is.null(region)){
    permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = region, 
                                                    treatment_col = "diet", design = "double_latin_square")
    
    # perform betadisper
    betadis <- do_betadisper(bray, meta = meta, permutation_matrix = permutation_matrix)
    print(betadis)
    
    p <- plot_betadisper(betadis = betadis, meta = meta, title = title)
    print(p)
    
    if (isTRUE(save)) {
      script_no <- get_script_number()
      save_name2 <- paste0(script_no, "_betadis_", save_name)
      save_big(save_name2)
    }
    
    # perform permanova
    permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)
    print(permanova)
    
    # perform pairwise permanova if significant
    if (permanova$`Pr(>F)`[1] < 0.05) {
      cld <- pairwise_permanova(bray, meta = meta, region = region, design = "double_latin_square")
      # replace diet groups with cld's for plotting
      meta_cld <- meta %>%
        inner_join(cld, by = c("diet"="treatments")) %>%
        dplyr::select(-diet) %>%
        dplyr::rename(diet = treatments_cld)
    } else {meta_cld = meta}
  } else {meta_cld = meta}
  
  p <- create_pcoa(bray, meta = meta_cld, title = title, second_indicator = second_indicator) 
  if (!is.null(region)) {
    Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), paste("=", round(permanova$`Pr(>F)`[1],3)))
    R2 <- format(round(permanova$R2[1], 3),nsmall = 3)
    p <- p + labs(subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print)))
  }
  return(p)
  
  if (isTRUE(save)) {
    script_no <- get_script_number()
    save_name2 <- paste0(script_no, "_pcoa_", save_name)
    save_big(save_name2)
  }
}

# create taxa barplots from kraken structured table

filter_ktable <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                          selected_matrix = c("ileal digesta", "faeces", "ferment")) {
  filtered_table <- ktable %>%
    inner_join(dplyr::select(meta, sampleid, matrix), by = "sampleid") %>%
    {if (selected_rank != "low") filter(., rank == selected_rank) else .} %>%
    {if (selected_rank != "low") dplyr::select(., -rank) else .} %>%
    filter(matrix == selected_matrix) %>%
    dplyr::select(-matrix)
  return(filtered_table)
}

aggregate_low_abundant_taxa <- function(filtered_table, threshold = 1) {
  aggregated_table <- filtered_table %>%
    group_by(name) %>%
    mutate(mean = mean(rel_abd)) %>%
    ungroup() %>%
    mutate(name = ifelse(mean < threshold, "other", name)) %>%
    group_by(sampleid, name) %>%
    summarize(rel_abd = sum(rel_abd), .groups = "drop")
  return(aggregated_table)
}

create_taxa_barplot <- function(table, meta, title = "") {
  plotting_table <- table %>%
    inner_join(meta, by = "sampleid") %>%
    group_by(diet, name) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop")
  if ("other" %in% plotting_table$name) { 
    plotting_table <- mutate(plotting_table, name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))
    } # fits other at the end of the legend
  
  ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    labs(x = "Diet", y = "Relative abundance [%]", fill = "Taxa", title = title) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0,100.01), expand = c(0,0))
}

# plotting function to create individual barplots for each samples, grouped by a grouping factor
# grouping factor should be a vector c() with a grouping factor for the facetting, and one for the x axis label

create_taxa_barplot_individual <- function(table, meta, title = "", grouping_factors) {
  grouping_factor <- grouping_factors[1]
  x_axis <- grouping_factors[2]
  plotting_table <- table %>%
    inner_join(meta, by = "sampleid") %>%
    mutate(sampleid = str_c(!!sym(grouping_factor), sampleid, sep = "_")) %>%
    mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other"))) # fits other at the end of the legend
  
  ggplot(plotting_table, aes(x = !!sym(x_axis), y = rel_abd, fill = name)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    labs(y = "Relative abundance [%]", fill = "Taxa", title = title, subtitle = grouping_factor) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0,100.1), expand = c(0,0)) + # 100.1 to not exclude rows due to out of range
    facet_grid(cols= vars(!!sym(grouping_factor))) +
    theme(plot.subtitle = element_text(hjust = 0.5))
}

# if grouping_factors is not defined, a barplot will be created on means per diet, otherwise a bar for each sample will be plotted and ordered according to grouping variables

taxa_barplot_from_ktable <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S"), 
                                     selected_matrix = c("ileal digesta", "faeces"), title = "", threshold = 1,
                                     save_name, save = F, grouping_factors = NULL) {
  filtered_table <- filter_ktable(ktable = ktable, meta = meta, selected_rank = selected_rank, 
                                  selected_matrix = selected_matrix)
  
  aggregated_table <- aggregate_low_abundant_taxa(filtered_table = filtered_table, threshold = threshold)
  
  if (is.null(grouping_factors)){
    p <- create_taxa_barplot(aggregated_table, meta = meta, title = title)
  } else {
    p <- create_taxa_barplot_individual(aggregated_table, meta = meta, title = title, 
                                        grouping_factors = grouping_factors)
  }

  if (isTRUE(save)) {
    save_big(name = save_name)
  }
    
  return(p)
}

# Differential abundance analysis with ancombc

# columns name, rank, sampleid, abundance expected
prepare_table_for_ancombc <- function(input_df) {
  
  feature_column <- colnames(input_df)[1]
  abundance_column <- colnames(input_df)[3]
  
  matrix <- input_df %>%
    dplyr::select(feature_column, sampleid, abundance_column) %>% # added due to duplicate rownames error for "low"
    pivot_wider(names_from = feature_column, values_from = abundance_column) %>%
    arrange(sampleid) %>%
    column_to_rownames("sampleid")
  return(matrix)
}

filter_meta <- function(meta, selected_matrix = c("ileal digesta", "faeces")) {
  filtered_meta <- meta %>%
    filter(matrix == selected_matrix)
  return(filtered_meta)
}

# extract taxa from output and change to long format

extract_pairwise_taxa_to_list <- function(output, selected_matrix) {
  res_pair_list <- list()
  if (selected_matrix == "ferment") {
    res <- output$res # normal output for ileum and faeces
    comparisons <- c("diet4")
    comparisons_renamed <- c("diet4-diet3")
  } else {
    res <- output$res_pair # pairwise for ileum and faeces
    comparisons <- c("diet2", "diet3", "diet4", "diet3_diet2", "diet4_diet2", "diet4_diet3")
    comparisons_renamed <- c("diet2-diet1", "diet3-diet1", "diet4-diet1", "diet3-diet2", "diet4-diet2", "diet4-diet3")
  }
  
  for (i in 1:length(comparisons)) {
    current_comparison <- comparisons[i]
    current_lfc <- str_c("lfc_", current_comparison)
    current_se <- str_c("se_", current_comparison)
    current_p <- str_c("p_", current_comparison)
    current_q <- str_c("q_", current_comparison)
    current_ss <- str_c("passed_ss_", current_comparison)
    
    res_pair_filtered <- res %>%
      dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                    q = !!sym(current_q), ss = !!sym(current_ss)) %>%
      mutate(comparison = comparisons_renamed[i])
    res_pair_list[[i]] <- res_pair_filtered
    names(res_pair_list)[i] <- comparisons_renamed[i]
  }
  return(res_pair_list)
}

# waterfall plots of significant taxa for each comparison

create_waterfall_plot_from_list <- function(res_pair_list, q_threshold = .05, lfc_threshold = 0) {
  res_pair_plot_list <- list()
  for (i in 1:length(res_pair_list)) {
    current_df <- res_pair_list[[i]] %>%
      filter(q < q_threshold) %>%
      mutate(indicator_color = case_when(lfc > lfc_threshold & ss == 1 ~ "#33a02c", # color according to lfc and if passed sensitivity test
                                         lfc > lfc_threshold & ss == 0 ~ "#b2df8a",
                                         lfc < -1*lfc_threshold & ss == 1 ~ "#e31a1c",
                                         lfc < -1*lfc_threshold & ss == 0 ~ "#fb9a99",
                                         .default = "#000000")) %>%
      arrange(lfc) %>%
      mutate(taxon = factor(taxon, levels = taxon))
    
    p <- ggplot(current_df, aes(x = lfc, y = taxon, fill = indicator_color)) +
      geom_bar(stat = "identity", show.legend = T, color = "grey20", width = 0.7) +
      geom_errorbar(aes(xmin = lfc - se, xmax = lfc + se), color = "grey20", width = 0.3) +
      labs(x = "log fold change", y = "", title = names(res_pair_list)[i]) + 
      theme(axis.text.y = element_text(size = 8)) +
      scale_fill_identity() # use colors in column
    
    res_pair_plot_list[[i]] <- p
  }
  
  if (length(res_pair_plot_list) > 1) { # if there are more than one plot (ileal digesta and faeces)
    p_new <- (res_pair_plot_list[[1]] | res_pair_plot_list[[2]] | res_pair_plot_list[[3]]) /
      (res_pair_plot_list[[4]] | res_pair_plot_list[[5]] | res_pair_plot_list[[6]])
    
    print(p_new)
  } else {
    print(p)
  }
  
  return(res_pair_plot_list)
}

# volcano plots for each comparison
create_volcano_plot_from_list <- function(res_pair_list, title, q_threshold = .05, lfc_threshold = 0) {
  res_pair_long <- bind_rows(res_pair_list) %>%
    mutate(indicator_color = case_when(q < q_threshold & lfc > lfc_threshold & ss == 1 ~ "#33a02c", # color according to lfc and if passed sensitivity test
                                       q < q_threshold & lfc > lfc_threshold & ss == 0 ~ "#b2df8a",
                                       q < q_threshold & lfc < -1*lfc_threshold & ss == 1 ~ "#e31a1c",
                                       q < q_threshold & lfc < -1*lfc_threshold & ss == 0 ~ "#fb9a99",
                                       .default = "#000000"))
  
  p <- ggplot(res_pair_long, aes(x = lfc, y = -log10(p), color = indicator_color, size = indicator_color)) +
    geom_point(show.legend = F) +
    geom_vline(xintercept = 0 + lfc_threshold, linetype = 2, color = "grey30") +
    geom_vline(xintercept = 0 - lfc_threshold, linetype = 2, color = "grey30") +
    scale_color_identity() +
    scale_size_manual(values = c(.2,2,2,2,2)) +
    facet_wrap(~comparison, ncol = 3) +
    labs(x = expression(Log[2]~fold~change), 
         y = expression("Significance ("*-log[10]~italic(P)*-value*")"), title = title)
  
  print(p)
  
  return(p)
}

# combined function to perform ancombc


perform_ancombc_and_plot <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                     selected_matrix = c("ileal digesta", "faeces", "ferment")) {
  
  filtered_table <- filter_ktable(ktable, meta = meta, selected_rank = selected_rank, 
                                  selected_matrix = selected_matrix)
  
  matrix <- prepare_table_for_ancombc(filtered_table)
  
  filtered_meta <- filter_meta(meta = meta, selected_matrix = selected_matrix) %>%
    column_to_rownames("sampleid")
  
  # perform ancombc
  if (selected_matrix == "ferment") {
    output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + period", 
                       p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                       taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) # lib_cut = 1000 
  } else {
    output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                       #rand_formula = "(1|animal) + (1|period)", # not converging with random effects structure
                       p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                       taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) # lib_cut = 1000
  }
 
  
  res_pair_list <- extract_pairwise_taxa_to_list(output = output, selected_matrix = selected_matrix)
  
  # build save names
  script_no <- get_script_number()
  save_name <- str_c(script_no, "")
  
  # plot
  waterfall_list <- create_waterfall_plot_from_list(res_pair_list)
  save_big(name = str_c(script_no, "waterfall", str_replace(selected_matrix, " ", "_"), selected_rank, sep = "_"))
  
  create_volcano_plot_from_list(res_pair_list, title = str_c(selected_matrix, selected_rank, sep = " "))
  save_big(name = str_c(script_no, "volcano", str_replace(selected_matrix, " ", "_"), selected_rank, sep = "_"))
  
  # create output object
  out <- list()
  out[["input"]] <- filtered_table
  out[["meta"]] <- filtered_meta
  out[["ancombc_output"]] <- output
  out[["res_pair_list"]] <- res_pair_list
  out[["waterfall_list"]] <- waterfall_list
  
  return(out)
}


# helper function to create pval matrix from dataframe with p vals

create_p_matrix_from_df <- function(taxon_pvals) {
  p_matrix <- matrix(nrow = 4, ncol = 4, dimnames = list(c("diet1", "diet2", "diet3", "diet4"),
                                                         c("diet1", "diet2", "diet3", "diet4")))
  diag(p_matrix) <- 1
  
  for (i in 1:nrow(taxon_pvals)) {
    row <- taxon_pvals[i,] %>%
      separate(comparison, into = c("treat1", "treat2"), sep = "-")
    p_matrix[which(rownames(p_matrix) == row$treat1), which(colnames(p_matrix) == row$treat2)] <- row$q
    p_matrix[which(rownames(p_matrix) == row$treat2), which(colnames(p_matrix) == row$treat1)] <- row$q
  }
  
  return(p_matrix)
}

# create filtered plot from ancombc output

create_plot_for_taxon <- function(output_object, df_rel_abd, meta, 
                                  selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                  selected_matrix = c("ileal digesta", "faeces"), selected_taxon,
                                  save = F) {
  
  filtered_df <- filter_ktable(df_rel_abd, meta = meta, selected_rank = selected_rank, 
                               selected_matrix = selected_matrix)
  
  p_adj <- filter(output_object$ancombc_output$res_global, taxon == selected_taxon)$q_val
  Print <- round_p_value(p_adj)
  
  taxon_pvals <- bind_rows(output_object$res_pair_list) %>%
    filter(taxon == selected_taxon)
  
  p_matrix <- create_p_matrix_from_df(taxon_pvals)
  
  cld <- generate_cld(p_matrix) %>%
    mutate(cld = treatments_cld,
           cld = str_remove(cld, treatments),
           diet = str_remove(treatments, "diet")) 
  
  taxon_df <- filtered_df %>%
    filter(name == selected_taxon) %>%
    left_join(dplyr::select(meta, sampleid, diet), by = "sampleid") %>%
    left_join(cld, by = "diet") %>%
    group_by(diet) %>%
    mutate(mean = mean(rel_abd)) %>%
    ungroup()
  
  p <- ggplot(taxon_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(outliers = F, fill = "snow2", width = 0.3) +
    geom_beeswarm() +
    geom_text(aes(y = mean, label = cld), size = 5, position = position_nudge(x = 0.2)) +
    annotate("text", x = 2.5, y = max(taxon_df$rel_abd), label = paste("q", Print), size = 6) +
    labs(x = "Diet", y = "Relative abundance (%)", title = selected_taxon)
  print(p)
  
  if (isTRUE(save)) {
    # build save names
    script_no <- get_script_number()
    
    save_big(name = str_c(script_no, "taxon_boxplot", str_replace(selected_matrix, " ", "_"), selected_rank, 
                          selected_taxon, sep = "_"))
  }
}

# loop taxa plot function ofer a vector of significant taxa

create_plots_for_taxa <- function(output_object, df_rel_abd, meta, 
                                  selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                  selected_matrix = c("ileal digesta", "faeces"), 
                                  taxa_vector, save) {
  for (i in 1:length(taxa_vector)) {
    create_plot_for_taxon(output_object = output_object, df_rel_abd = df_rel_abd, meta = meta,
                          selected_rank = selected_rank, selected_matrix = selected_matrix,
                          selected_taxon = taxa_vector[i], save = save)
  }
}


# edgeR

# creates edgeR object from table with sampleid column, feature column in first place, and custom intensity or reads column

create_edger_object <- function(input_df, abundance_column) {
  feature_column <- colnames(input_df)[1]
  
  output <- input_df %>%
    dplyr::select(!!sym(feature_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = sampleid, values_from = !!sym(abundance_column))
  mat <- as.matrix(output[,-1])
  rownames(mat) <- output[[1]]
  
  groups <- tibble("sampleid" = colnames(mat)) %>%
    inner_join(meta, by = "sampleid") %>%
    dplyr::select(sampleid, animal, period, diet)
  # combine to edger object
  edger_object <- DGEList(counts = mat, samples=groups)
  
  return(edger_object)
}

edger_analysis <- function(edger_object, title = "", feature_column, selected_matrix) {
  #create design matrix
  if (selected_matrix == "ferment") {
    design <- model.matrix(~diet+period, data = edger_object$samples) # reduced formula for ferment
  } else { # for ileum and faeces
    design <- model.matrix(~diet+animal+period, data = edger_object$samples)
  }
  rownames(design) <- edger_object$samples$sampleid
  #estimate dispersion
  disp <- estimateDisp(edger_object, design)
  plotBCV(disp)
  #fit glm
  fit <- glmFit(disp, design = design)
  
  # do anova
  if (selected_matrix == "ferment") {
    anova <- glmLRT(fit, coef = 2) # for ferment only second coef should be tested (diet4)
  } else {
    anova <- glmLRT(fit, coef = 2:4)
  }

  anova_adj <- decideTests(anova)
  anova_adj_sum <- summary(anova_adj)
  
  # continue if significant
  if (selected_matrix == "ferment") {
    if ((anova_adj_sum[1] + anova_adj_sum[3]) == 0) {
      print("No significant differences.")
      plotMD(anova)
      continue = F #dummy
    } else {
      contrasts <- design[0,] 
      contrasts <- rbind(contrasts,
                         "diet4-diet3" = c(0,1,0,0,0))
      results <- vector("list", length = nrow(contrasts))
      names(results) <- rownames(contrasts)
      
      continue = T #dummy
    }
  } else {
    if (anova_adj_sum[2] == 0) {
      print("No significant differences.")
      plotMD(anova)
      
      continue = F #dummy
    } else {
      contrasts <- design[0,] 
      contrasts <- rbind(contrasts,
                         "diet2-diet1" = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                         "diet3-diet1" = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                         "diet4-diet1" = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                         "diet3-diet2" = c(0,-1,1,0,0,0,0,0,0,0,0,0,0,0),
                         "diet4-diet2" = c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0),
                         "diet4-diet3" = c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0))
      results <- vector("list", length = nrow(contrasts))
      names(results) <- rownames(contrasts)
      
      continue = T #dummy
    }
  }
  
  if (isTRUE(continue)) {
    for (i in 1:nrow(contrasts)) {
      contrast <- as.vector(contrasts[i,])
      result <- glmTreat(fit, contrast = contrast)#, lfc = 1) # test + filters lfc > 1
      sig <- decideTests(result, adjust.method = "BH") # adjust P and filter according to 0.05
      table <- topTags(result, n=nrow(fit$counts), adjust.method = "BH", sort.by = "none", p.value = 1)$table # output all unsorted
      print(summary(sig))
      #save in output
      results[[i]][["result"]] <- result
      results[[i]][["significant"]] <- sig
      results[[i]][["table"]] <- table
    }
  }
  
  return(results)
}

# takes edger output list and filters every table for significant results
# stores output table in edger output object

create_table_sig_from_edger <- function(edger_results, feature_column) {
  for (i in 1:length(edger_results)) {
    table <- edger_results[[i]][["table"]] %>%
      as_tibble(rownames = feature_column) %>%
      mutate(contrast = names(edger_results)[i],
             sig = ifelse(FDR < 0.05, "sig", "not sig"),
             sig = ifelse(is.na(sig), "not sig", sig))
    edger_results[[i]][["table_sig"]] <- table # safe table to list for later
  }
  
  return(edger_results)
}

# concatenates edger outputs and creates a volcano plot

create_volcano_plot_from_edger <- function(edger_results, title = "", no_label = 0) {
  
  # concatenate all results and create plot
  for (i in 1:length(edger_results)) {
    table <- edger_results[[i]][["table_sig"]]
    if (i == 1) {
      all_results <- table 
    } else { # from second iteration onwards combine previous tables with latest table
      all_results <- rbind(all_results, table)
    }
  }
  
  first_column <- colnames(all_results)[1]
  # add labels for ggrepel
  all_results <- all_results %>%
    group_by(contrast, sig) %>%
    mutate(rank_asc = row_number(logFC),
           rank_desc = row_number(-logFC),
           label = ifelse(sig == "sig" & ((rank_desc <= no_label & logFC > 0 )| (rank_asc <= no_label & logFC < 0)), 
                          !!sym(first_column), "")) %>%
    ungroup()
  
  p <- ggplot(all_results, aes(x = logFC, y = -log10(PValue), color = sig, size = sig)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey30") +
    geom_point(show.legend = T) +
    facet_wrap(~contrast, ncol = 3) +
    scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
    scale_size_manual(values = c(0.1,1), labels = c("q ≥ 0.05", "q < 0.05")) +
    labs(x = expression(Log[2]~fold~change), 
         y = expression("Significance ("*-log[10]~italic(P)*-value*")"), title = title, 
         color = "Significance", size = "Significance") +
    theme(legend.position = "bottom")
  
  if (no_label != 0) {
    p <- p +
      geom_text_repel(aes(x = logFC, y = -log10(PValue), label = label), color = "black", show.legend = F, 
                      max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, max.time = 3, size = 3, 
                      inherit.aes = F, segment.color = "grey30") +
      scale_x_continuous(expand = expansion(mult = 0.2)) +
      scale_y_continuous(expand = expansion(mult = c(0,0.2)))
  }
  return(p)
}

perform_edger_and_plot <- function(input_df, selected_matrix = c("ileum", "faeces", "ferment"), abundance_column, save_name, save = F) {
  
  # filter
  if (selected_matrix == "ileum") {
    filtered_df <- filter_ileum(input_df)
  } else if (selected_matrix == "faeces") {
    filtered_df <- filter_faeces(input_df)
  } else if (selected_matrix == "ferment") {
    filtered_df <- filter_ferment(input_df)
  }
  
  edger_object <- create_edger_object(filtered_df, abundance_column = abundance_column)
  
  edger_results <- edger_analysis(edger_object = edger_object, title = str_c(save_name, selected_matrix, sep = " "),
                                  feature_column = colnames(input_df)[1], selected_matrix = selected_matrix)
  
  edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = colnames(input_df)[1])
  
  p <- create_volcano_plot_from_edger(edger_results = edger_results, title = str_c(save_name, selected_matrix, sep = " "),
                                      no_label = 10)
  print(p)
  
  # build save names
  if (isTRUE(save)) {
    script_no <- get_script_number()
    save_big(name = str_c(script_no, save_name, str_replace(selected_matrix, " ", "_"), sep = "_"))
  }
  
  return(edger_results)
}

# add functions to edger result object based on meta file with all functions
# join_by defines column to join, must be same in results object and functions_object

add_functions <- function(results_object, functions_object, join_by) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      #filter(sig == "sig") %>%
      left_join(functions_object, by = join_by)
    results_object[[i]][["table_sig"]] <- table_func # save as table sig
  }
  return(results_object)
}

# add taxonomy to edger results object based on meta file with all taxa assignments
# join_by defines join column, but should be bins

add_taxonomy <- function(results_object, taxonomy_object, join_by) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      mutate(bin = proteinid) %>%
      separate(bin, into = "bin", sep = "_") %>%
      #filter(sig == "sig") %>%
      left_join(taxonomy_object, by = join_by)
    results_object[[i]][["table_sig"]] <- table_func # save as table sig
  }
  return(results_object)
}

# counts the number of proteins from each taxa and plots the 5 most abundant significant taxa
# needs a annotated edger object

create_taxa_count_plot_from_edger <- function(edger_results) {
  plot_list <- list()
  
  for (i in 1:length(edger_results)) {
    current_df <- edger_results[[i]]$table_sig %>%
      filter(sig == "sig") %>%
      mutate(direction = ifelse(logFC > 0, 1, -1)) %>%
      mutate(indicator_color = ifelse(logFC > 0, "#33a02c", "#e31a1c")) %>%
      group_by(contrast, direction, indicator_color, G) %>%
      summarise(count = length(G), .groups = "drop") %>%
      mutate(count = direction * count) %>%
      mutate(rank_desc = row_number(desc(count)),
             rank_asc = row_number(count)) %>%
      filter(direction == -1 & rank_asc <= 5 | 
               direction == 1 & rank_desc <= 5) %>%
      arrange(count) %>%
      mutate(rank_asc = factor(rank_asc, levels = rank_asc)) %>%
      mutate(text_pos = ifelse(direction == 1, G, ""),
             text_neg = ifelse(direction == -1, G, "")) %>%
      mutate(count_pos = ifelse(direction == 1, count, ""),
             count_neg = ifelse(direction == -1, abs(count), ""))
    
    p <- ggplot(current_df, aes(x = count, y = rank_asc, fill = indicator_color))  +
      geom_bar(stat = "identity", show.legend = F, color = "grey20", width = 0.7) +
      labs(x = "# proteins", y = "genus", title = names(edger_results)[i]) + 
      geom_text(aes(label = text_pos, y = rank_asc, x = 0), hjust = 1, size = 3) +
      geom_text(aes(label = text_neg, y = rank_asc, x = 0), hjust = 0, size = 3) +
      geom_text(aes(label = count_pos, y = rank_asc, x = count), hjust = 0, fontface = "bold") +
      geom_text(aes(label = count_neg, y = rank_asc, x = count), hjust = 1, fontface = "bold") +
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank()) +
      scale_x_continuous(limits = c(-max(abs(current_df$count)), max(abs(current_df$count))), expand = c(.2)) +
      scale_fill_identity() # use colors in column
    
    plot_list[[i]] <- p
  }
  
  if (length(plot_list) > 1) { # if there are more than one plot (ileal digesta and faeces)
    p_new <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
      (plot_list[[4]] | plot_list[[5]] | plot_list[[6]])
    
    print(p_new)
  } else {
    print(p)
  }
  
  return(plot_list)
}

# from protein or gene differential abundance analysis with edger
# relocates the function column to be used for enrichment analysis

transform_to_function <- function(results_object, function_name) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      filter(!is.na(!!sym(function_name))) %>%
      filter(!!sym(function_name) != "-") %>%
      mutate(!!sym(function_name) := str_remove(!!sym(function_name), ";$")) %>%
      separate_rows(!!sym(function_name), sep = ",") %>%
      separate_rows(!!sym(function_name), sep = ";") %>%    
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "^ssc:")) %>%
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "^ko:")) %>%
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "\\s")) %>%
      dplyr::relocate(!!sym(function_name))
    
    results_object[[i]][["table_sig"]] <- table_func
  }
  return(results_object)
}

# Enrichment analysis from edger results object
# input should be an edger results object with kegg ko's in first column of table_sig 

enrich_kegg <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    feature_name <- names(table_sig)[1] # kegg ko's must be the first column
    universe <- table_sig[[feature_name]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(feature_name))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(feature_name))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichKO(overexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(name = paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(paste("Overexpressed KEGGs for", contrast))
        cat("\n", overexpressed, "\n")
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichKO(underexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(paste("Underexpressed KEGGs for", contrast))
        cat("\n", underexpressed, "\n")
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# enrichment function for kegg modules
# still a bit buggy
enrich_mkegg <- function(results_object, save_name, save) { 
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    feature_name <- names(table_sig)[1] # kegg ko's must be the first column
    universe <- table_sig[[feature_name]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(feature_name))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(feature_name))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichModule(overexpressed, universe = universe, minGSSize = 6) # minGsize not recognized correct
      if (!is.null(enrichment)) {
        if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
          p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
          print(p)
          if (isTRUE(save)) {
            save_big(name = paste0(save_name, "_", contrast, "_overexpressed"))
          }
          print(paste("Overexpressed KEGGs for", contrast))
          cat("\n", overexpressed, "\n")
        } else {
          print(paste("No enriched KEGGs found for", contrast))
        }
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichModule(underexpressed, universe = universe)
      if (!is.null(enrichment)) {
        if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
          p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
          print(p)
          if (isTRUE(save)) {
            save_big(paste0(save_name, "_", contrast, "_underexpressed"))
          }
          print(paste("Underexpressed KEGGs for", contrast))
          cat("\n", underexpressed, "\n")
        } else {
          print(paste("No enriched KEGGs found for", contrast))
        }
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# creates a heatmap from an edger object and a object containing only significant proteins incl. rel_abd per sample
# heatmap contains cld

create_heatmap_from_edger <- function(results_object, sig_df, function_df = NULL, heatmap_y = "protein_names",
                                      top = NULL, save_name, save) {
  
  first_column <- colnames(sig_df)[1]
  
  sig_proteins <- unique(sig_df[[sym(first_column)]])
  
  #aggregate results in one table
  for (i in 1:length(results_object)) {
    protein_pval <- results_object[[i]]$table_sig %>%
      dplyr::select(first_column, comparison = contrast, q = FDR, logFC) # adapt to previous functions
    if (i == 1) {
      protein_pvals <- protein_pval
    } else {
      protein_pvals <- rbind(protein_pvals, protein_pval)
    }
  }
  
  protein_pvals <- protein_pvals %>%
    filter(!!sym(first_column) %in% sig_proteins)
  
  if (!is.null(top)) { # filter top proteins by logfc
    protein_logfc <- protein_pvals %>%
      group_by(!!sym(first_column)) %>%
      summarise(logFC = max(logFC), .groups = "drop") %>%
      mutate(rank = row_number(desc(abs(logFC)))) %>%
      filter(rank <= top) %>%
      pull(!!sym(first_column))
    protein_pvals <- protein_pvals %>%
      filter(!!sym(first_column) %in% protein_logfc)
    sig_df <- sig_df %>%
      filter(!!sym(first_column) %in% protein_logfc)
    sig_proteins <- protein_logfc # replace for next loop
  }
  
  for (i in 1:length(sig_proteins)) {
    protein_pvals_filtered <- protein_pvals %>%
      filter(!!sym(first_column) == sig_proteins[i])
    
    p_matrix <- create_p_matrix_from_df(protein_pvals_filtered)
    
    cld <- generate_cld(p_matrix) %>%
      mutate(cld = treatments_cld,
             cld = str_remove(cld, treatments),
             diet = str_remove(treatments, "diet")) %>%
      add_column(!!sym(first_column) := sig_proteins[i]) %>%
      dplyr::select(first_column, diet, cld)
    
    if (i == 1) {
      protein_cld <- cld
    } else {
      protein_cld <- rbind(protein_cld, cld)
    }
  }
  
  if (!is.null(function_df)) {
    if ("character" %in% class(function_df) && function_df == "annotate_keggs") { # annotate keggs and save as function_df
      function_df <- annotate_keggs(unique(sig_proteins)) %>%
        rename(protein_names = name, !!sym(first_column) := kegg_ko)
    }
  }

  heatmap <- sig_df %>%
    left_join(meta, by = "sampleid") %>%
    group_by(!!sym(first_column), diet) %>%  # identical protein-names will be averaged
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    {if (!is.null(function_df)) left_join(., function_df, by = first_column) else .} %>% 
    left_join(protein_cld, by = c(first_column, "diet")) %>%
    group_by(!!sym(first_column)) %>%
    mutate(value = scale(rel_abd)) %>%
    ungroup()
  
  p <- ggplot(heatmap, aes(x = diet, y = !!sym(heatmap_y), fill = value)) +
    geom_tile(show.legend = F) +
    geom_text(aes(label = cld), size = 2) +
    coord_fixed() +
    scale_y_discrete(position = "right") +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 8)) +
    labs(y = "")
  print(p)
  
  if (isTRUE(save)) {
    save_big(name = save_name)
  }
}


##############################
# Functions from Metaproteomics
##############################

# calculate taxa from protein intensities
# takes protein df and genomes df

calculate_taxa_intensity <- function(input, genomes) {
  output <- input %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(proteingroup, sep = ";") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(proteingroup)) %>% # calculating factor
    ungroup() %>%
    mutate(intensity = intensity * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    filter(str_detect(proteingroup, "^MGYG")) %>% # filter out host proteins
    separate(proteingroup, into = "bin", sep = "_") %>%
    group_by(sampleid, bin) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>% # sum up intensities per sample and genome
    left_join(genomes, by = "bin") %>%
    dplyr::select(bin, sampleid, intensity, "R1", P, C, O, "F", G, S)
  return(output)
}


# enrichment analysis of COGs for microbial proteins
# takes edger results object
# first column of sig df should be COGs -> either from cog dea analysis or translated to cog

enrich_cog <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[["first_column"]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichCOG(overexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched COGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed COGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichCOG(underexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched COGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed COGs for", contrast))
    }
  }
}

# Enrichtment analysis for host keggs
# takes edger results object
# first column of sig df should be keggs -> either from cog dea analysis or translated to kegg
# default organism is sus scrofa, "psat" for pisum sativum

enrich_kegg_host <- function(results_object, save_name, save, organism = "ssc") {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[[first_column]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichKEGG(overexpressed, universe = universe, organism = organism)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichKEGG(underexpressed, universe = universe, organism = organism)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# Enrichtment analysis for host gos
# takes edger results object
# first column of sig df should be go -> either from cog dea analysis or translated to go

enrich_go_host <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[[first_column]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichGO(overexpressed, "org.Ss.eg.db",universe = universe, keyType = "GO", ont = "BP")
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched GOs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichGO(underexpressed, "org.Ss.eg.db",universe = universe, keyType = "GO", ont = "BP")
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched GOs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# output significant proteins etc for correlation purposes

get_significant <- function(venn_input, input_df) {
  first_column <- colnames(input_df)[1]
  output <- venn_input %>%
    as_tibble(rownames = first_column) %>%
    dplyr::select(!!sym(first_column)) %>%
    inner_join(input_df, by = first_column) 
  return(output)
}

# Venn diagramms
# takes results object of ancombc or edgeR and input_df with rel_abd
# selected_rank and selected_matrix need to be specified in case of a ktable
# selected_matrix must be specified to filter for matrix

create_venn_for_sig <- function(results_object, input_df, selected_matrix = c("ileal digesta", "faeces"),
                                selected_rank = NULL, save_name, save) {
  # pre test for ancombc objects
  if (names(results_object)[1] == "input") { # indicator for ancombc object
    results_object <- results_object[["res_pair_list"]] # change results object to res_pair_list
  }
  
  for (i in 1:length(results_object)) {
    first_column <- colnames(input_df)[1]
    contrast <- names(results_object)[i]
    table <- results_object[[i]]
    #ancombc vs edgeR block
    if (is.data.frame(table)) { # detect if table is already dataframe (ancombc) or list (edgeR)
      table <- table %>%
        mutate(sig = ifelse(q < 0.05, "sig", NA)) %>% # if ancombc add sig column
        dplyr::rename(!!sym(first_column) := taxon) # unify name of first column
    } else {
      table <- table[["table_sig"]] # choose table_sig if edgeR
    }
     table_sig <- table %>%
      dplyr::select(!!sym(first_column), sig) %>%
      mutate(sig = ifelse(sig == "sig", 1, 0),
             sig = ifelse(is.na(sig), 0, sig)) %>% # replace NAs separately
      dplyr::rename(!!sym(contrast):=sig)
    if (i == 1) {
      all_sig <- table_sig
    } else { # from second iteration onwards combine previous tables with latest table
      all_sig <- all_sig %>%
        left_join(table_sig, by = first_column)
    }
  }
  sig_df <- as.data.frame(all_sig[,-1])
  rownames(sig_df) <- all_sig[[first_column]]
  sig_df <- sig_df[rowSums(sig_df) > 0,]
  
  #create venn
  save_venn(sig_df)
  
  # create sig output table
  if (is.null(selected_rank)) { # if no rank is specified -> no ktable
    output_df <- filter_ktable(input_df, meta = meta, selected_rank = "low", 
                              selected_matrix = selected_matrix)
  } else {
    output_df <- filter_ktable(input_df, meta = meta, selected_rank = selected_rank, 
                              selected_matrix = selected_matrix)
  }
  
  output <- get_significant(sig_df, output_df) # return tibble with rel_abd for correlating
  
  # create save name and save
  if (isTRUE(save)) {
    no <- get_script_number()
    matrix <- str_extract(selected_matrix, ".{2}")
    saveRDS(output, file = str_c("clean/", str_c(no, "sig", matrix, selected_rank, save_name, sep = "_"), ".RDS"))
  }
  
  return(output)
}

# Function to plot a barplot for one specific function, showing the taxa origin for each diet
# needs a df with function and taxa combined
# type of function needs to be specified e.g. cazy and which functions should be filtered

plot_func_tax_origin <- function(func_taxa_df = func_tax_combined, function_filter, function_type = "kegg_ko",
                                 region = c("ileal digesta", "faeces", "ferment"), taxa_level = "G") {
  first_column <- colnames(func_taxa_df)[1]
  
  if (region == "ileal digesta") {
    prefiltered_df <- filter_ileum(func_taxa_df)
  } else if (region == "faeces") {
    prefiltered_df <- filter_faeces(func_taxa_df)
  } else if (region == "ferment") {
    prefiltered_df <- filter_ferment(func_taxa_df)
  }
  
  filtered_df <- prefiltered_df %>%
    left_join(meta %>% dplyr::select(-description), by = "sampleid") %>%
    filter(str_detect(!!sym(function_type), function_filter)) %>%
    group_by(!!sym(first_column), !!sym(taxa_level), diet) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    group_by(!!sym(taxa_level), diet) %>%
    summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
    group_by(!!sym(taxa_level)) %>% 
    mutate(mean_rel_abd = mean(rel_abd)) %>%
    ungroup() %>%
    mutate(!!sym(taxa_level) := ifelse(mean_rel_abd < 0.00025*sum(rel_abd), "other", !!sym(taxa_level)))
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    labs(title = function_filter)
}

# function that plots the protein abundance as boxplots
# works with everything that is based on first column
# has to have a rel abd column
# dont use | in protein

plot_protein_abundance <- function(input_df, region, protein) {
  first_column <- colnames(input_df)[1]
  
  if (region == "ileal digesta") {
    prefiltered_df <- filter_ileum(input_df)
  } else if (region == "faeces") {
    prefiltered_df <- filter_faeces(input_df)
  } else if (region == "ferment") {
    prefiltered_df <- filter_ferment(input_df)
  }
  
  filtered_df <- prefiltered_df %>%
    filter(str_detect(!!sym(first_column), protein)) %>%
    left_join(meta, by = "sampleid")
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = protein)
}

# function that plots the taxa abundance as boxplots
# works with a ktable
# has to have a rel abd column
# name has to be in first column
# dont use | in protein

plot_taxa_abundance <- function(ktable, region, taxa_level, taxa) {
  first_column <- colnames(ktable)[1]
  
  if (region == "ileal digesta") {
    prefiltered_df <- filter_ileum(ktable)
  } else if (region == "faeces") {
    prefiltered_df <- filter_faeces(ktable)
  } else if (region == "ferment") {
    prefiltered_df <- filter_ferment(ktable)
  }
  
  filtered_df <- prefiltered_df %>%
    filter(rank == taxa_level) %>%
    dplyr::select(-rank) %>%
    filter(str_detect(!!sym(first_column), taxa)) %>%
    left_join(meta, by = "sampleid")
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = taxa)
}

##############################
# functions from metabolomics
#############################

# function for preparation

prepare_metabolomics <- function(input) {
  output <- input %>%
    inner_join(dplyr::select(meta, sampleid, square, animal, period, diet)) %>%
    relocate(animal, period, square, diet) %>%
    dplyr::select(-sampleid)
  return(output)
}

# analysis function

loop_comparison_metabolomics <- function(input_df, y_axis, table_title, save_name, save = save) {
  shapiro_table <- data.frame("metabolite" = as.character(), "p_value"=as.numeric(), transformation = as.character()) # to check for bad transformations
  pb <- txtProgressBar(min = 1, max = ncol(input_df), style = 3) # initialize progessbar
  for (i in 5:ncol(input_df)) { #assuming 4 meta columns in the df
    setTxtProgressBar(pb, i) # print progessbar
    
    number <- get_script_number()
    
    current_response <- colnames(input_df)[i]
    # filter df for i'th response
    filtered_df <- select_response(input_df, current_response)
    # perform reml
    out <- combined_comparison(df = filtered_df, transformation = "test")
    # save shapiro results
    shapiro_table <- rbind(shapiro_table, c(current_response, round(out$shapiro$p.value, 3), out$transformation[1]))
    # print plot
    plot_pairwise(filtered_df = filtered_df, output = out, selected_response = current_response)
    if (isTRUE(save)) {
      save_big(paste0(number, "_", save_name, "_", current_response))
    }
    # print and save plot
    create_results_plot(input_df = filtered_df, comparison_object = out, response_name = current_response, y_axis = y_axis)
    if (isTRUE(save)) {
      save_big(paste0(number, "_result_", save_name, "_", current_response))
    }
    table <- create_results_table(input_df = filtered_df, comparison_object = out, response_name = current_response, digits = 1)
    table <- rbind(table, c("exact_p", out$ano[which(rownames(out$ano) == "diet"), 
                                                             which(colnames(out$ano) == "Pr(>F)")])) # for exact p value adjustment
    
    if (i == 5) { # create table if first iteration, join afterwards
      final_table <- table
    } else {
      final_table <- final_table %>%
        left_join(table, by = "diet")
    }
    
  }
  final_table <- rbind(final_table[-7,], 
                       c("P-adj", 
                         round(p.adjust(final_table[7, 2:ncol(final_table)], method = "BH"), 3)))
  if (isTRUE(save)) {
    final_table_gt <- save_results_table(final_table, title = table_title, name =  paste0(number, "_table_", save_name))
    gtsave(final_table_gt, filename = paste0(number, "_table_", save_name, ".png"), "plots", vwidth = 3000)
  }

  close(pb) # close progessbar
  print(shapiro_table)
  return(final_table)
}


#####################
# Metabolomics functions for ferment
#######################

# test model for ferment

test_models_ferment <- function(df) {
  model <- lmer(response ~ diet + (1|period), data = df)
  if (VarCorr(model)["period"] == 0) {
    best_model <- aov(response ~ diet + period, data = df)
  } else {
    best_model <- model
  }
  print(best_model)
  return(best_model)
}

# own tukey function with ferment model

transformTukey2_ferment <- function(df, start = -10, end = 10, int = 0.25) { 
  result_frame <- data.frame(transformation = seq(start, end, int),
                             W = 0)
  for (i in 1:nrow(result_frame)) {
    transformation_factor <- result_frame$transformation[i]
    transformed_df <- df
    # test transformations
    if(transformation_factor > 0) {
      transformed <- df$response^(transformation_factor)
      transformed_df$response <- transformed
    } else if (transformation_factor == 0) {
      transformed <- log(df$response)
      transformed_df$response <- transformed
    } else if (transformation_factor < 0) {
      transformed <- -1 * df$response^transformation_factor
      transformed_df$response <- transformed
    }
    if (any(is.infinite(transformed_df$response)) == FALSE & any(is.nan(transformed_df$response)) == FALSE) { # only validate if no NA or infinities in transformed data
      x <- suppressMessages(suppressWarnings(
        { capture.output(model <- test_models_ferment(transformed_df)) }
      ))
      test_statistic <- shapiro.test(resid(model))
      if (deviance(model) < sqrt(.Machine$double.eps)) { # if precision gets to low (similar to test in car::Anova)
        result_frame$W[i] <- 0 # do not use this transformation
      } else {
        result_frame$W[i] <- test_statistic$statistic
      }
    }
  }
  plot(x = result_frame$transformation, y = result_frame$W)
  W <- result_frame$transformation[which.max(result_frame$W)] 
  print(W)
  # generate output
  if(W > 0) {
    transformed <- df$response^(W)
  } else if (W == 0) {
    transformed <- log(df$response)
  } else if (W < 0) {
    transformed <- -1 * df$response^W
  }
  transformed_list <- list(transformed, W)
  return(transformed_list)
}

# function to perform transformations on concentration column with tukey2 for ferment

transform_data_ferment <- function(df, transformation = c("gauss", "tukey", "tukey2", "log", "logit"), type = c("h", "hh", "s")) {
  if (transformation == "gauss") {
    transformed <- Gaussianize(df$response, type = type, return.tau.mat = T)
    
  } else if (transformation == "tukey") {
    
    transformed1 <- transformTukey(df$response)
    transformed2 <- transformTukey(df$response, returnLambda = T)
    transformed <- list(transformed1, transformed2)
  } else if (transformation == "tukey2") { # leads to choosing the worst model
    transformed <- transformTukey2_ferment(df)
  } else if (transformation == "log") {
    transformed1 <- log(df$response)
    transformed2 <- 1
    transformed <- list(transformed1, transformed2)
  } else if (transformation == "logit") {
    transformed1 <- log((df$response/100)/(1-(df$response/100)))
    transformed2 <- 1
    transformed <- list(transformed1, transformed2)
  }
  df$response <- transformed[[1]][1:length(transformed[[1]])]
  transformation_factor <- transformed[[2]]
  transformation <- list(df, transformation_factor)
  return(transformation)
}

# function to test different transformations and find the best one according to shapiro wilk test with transformation function for ferment

test_transformations_ferment <- function(df) {
  result_frame <- data.frame(transformation = c("none", "tukey", "tukey2", "log", "logit", "gauss", "gauss", "gauss"),
                             type = c(NA, NA, NA, NA, NA, "h", "hh", "s"),
                             W = 0)
  for (i in 1:nrow(result_frame)) {
    transformation <- result_frame$transformation[i]
    type <- result_frame$type[i]
    if (transformation == "none") {
      transformed_df <- df
    } else {
      transformed_df <- transform_data_ferment(df, transformation = transformation, type = type)[[1]]
    }
    
    if (any(is.infinite(transformed_df$response)) == FALSE & any(is.nan(transformed_df$response)) == FALSE) {
      x <- suppressMessages(suppressWarnings(
        { capture.output(model <- test_models_ferment(transformed_df)) }
      ))
      test_statistic <- shapiro.test(resid(model))
      if (deviance(model) < sqrt(.Machine$double.eps)) { # if precision gets to low (similar to test in car::Anova)
        result_frame$W[i] <- 0 # do not use this transformation
      } else {
        result_frame$W[i] <- test_statistic$statistic
      }
    }
    
  }
  best_transformation <- result_frame$transformation[which.max(result_frame$W)] # returning only first 
  best_type <- result_frame$type[which.max(result_frame$W)]
  return(list(best_transformation = best_transformation, best_type = best_type))
}

# combined comparison function for ferments

combined_comparison_ferment <- function(df,transformation = "none", 
                                type = c("h", "hh", "s"), model_style = "klein") { #, correction_filter
  if (transformation == "test") {
    best_transformation = test_transformations_ferment(df)
    transformation <- best_transformation$best_transformation
    type <- best_transformation$best_type
    print(paste("Best transformation:", transformation))
  }
  
  if (transformation == "none") {
    transformation_factor = NA # compatible for output
    
    model <- test_models_ferment(df)
    
    n <- create_n_table(df)
    
    shapiro <- evaluate_model(model)
    print(shapiro)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    
    cld2 <- cld # for compatibility with transformation (original means)
    
  } else {
    
    transformed <- transform_data_ferment(df, transformation = transformation, type = type)
    transformed_df <- transformed[[1]]
    transformation_factor <- transformed[[2]]
    
    model <- test_models_ferment(transformed_df)
    
    n <- create_n_table(df)
    
    shapiro <- evaluate_model(model)
    print(shapiro)
    
    ano <- do_anova(model)
    
    cld <- pairwise_comparisons(model)
    
    cld <- backtransform_data(cld, transformation = transformation, transformation_factor = transformation_factor)
    
    model_untransformed <- test_models_ferment(df)
    cld2 <- pairwise_comparisons(model_untransformed)
    
  }
  
  if (ano[which(rownames(ano) == "diet"), which(colnames(ano) == "Pr(>F)")] > .05) {
    cld <- cld %>%
      mutate(.group = "")
  } # if not significant, remove contrasts (if calculated)
  
  return(list(ano = ano, cld = cld, cld2 = cld2, n = n, 
              transformation = c(transformation, transformation_factor, type),
              shapiro = shapiro))
}

# plotting function for ferment

plot_pairwise_ferment <- function(filtered_df, output, selected_response, save = F) { #, correction_filter
  cld <- output$cld
  P <- output$ano[which(rownames(output$ano) == "diet"), which(colnames(output$ano) == "Pr(>F)")]
  Print <- round_p_value(P)
  p <- ggplot(cld, aes(x = diet)) +
    geom_boxplot(data = filtered_df, aes(y = response), width = 0.2, position = position_nudge(x = -0.1)) +
    geom_point(data = filtered_df, aes(y = response, color = period), 
               size = 2, position = position_nudge(x = -0.1), stroke = 2) +
    {if(P < 0.05)geom_text(aes(y = emmean, label = str_trim(.group)), size = 5, position = position_nudge(x = 0.2))}+
    labs(x = "diet", y = "value", title = paste(selected_response)) + #, correction_filter
    annotate("label", label = paste("diet\n P", Print), x = 2.25, y = max(filtered_df$response), size = 5) +
    #scale_color_manual(values = c("black", "grey20", "grey40", "grey60")) +
    scale_color_manual(values = c("red", "blue", "green", "yellow")) +
    scale_shape_manual(values = c(16,17,15,18,21,24,22,23))
  print(p)
}

# analysis function

loop_comparison_metabolomics_ferment <- function(input_df, y_axis, table_title, save_name, save = save, digits = 1) {
  shapiro_table <- data.frame("metabolite" = as.character(), "p_value"=as.numeric(), transformation = as.character()) # to check for bad transformations
  pb <- txtProgressBar(min = 1, max = ncol(input_df), style = 3) # initialize progessbar
  for (i in 5:ncol(input_df)) { #assuming 4 meta columns in the df
    setTxtProgressBar(pb, i) # print progessbar
    number <- get_script_number()
    
    current_response <- colnames(input_df)[i]
    # filter df for i'th response
    filtered_df <- select_response(input_df, current_response)
    # remove rows with NA (for functional activity)
    filtered_df <- filter(filtered_df, !is.na(response))
    # perform reml
    out <- combined_comparison_ferment(df = filtered_df, transformation = "test")
    # save shapiro results
    shapiro_table <- rbind(shapiro_table, c(current_response, round(out$shapiro$p.value, 3), out$transformation[1]))
    # print plot
    plot_pairwise_ferment(filtered_df = filtered_df, output = out, selected_response = current_response)
    if (isTRUE(save)) {
      save_big(paste0(number, "_", save_name, "_", current_response))
    }
    # print and save plot
    create_results_plot(input_df = filtered_df, comparison_object = out, response_name = current_response, y_axis = y_axis)
    if (isTRUE(save)) {
      save_big(paste0(number, "_result_", save_name, "_", current_response))
    }
    table <- create_results_table(input_df = filtered_df, comparison_object = out, response_name = current_response, digits = digits)
    table <- rbind(table, c("exact_p", out$ano[which(rownames(out$ano) == "diet"), 
                                               which(colnames(out$ano) == "Pr(>F)")])) # for exact p value adjustment
    
    if (i == 5) { # create table if first iteration, join afterwards
      final_table <- table
    } else {
      final_table <- final_table %>%
        left_join(table, by = "diet")
    }
    
  }
  final_table <- rbind(final_table[-7,], 
                       c("P-adj", 
                         round(p.adjust(final_table[5, 2:ncol(final_table)], method = "BH"), 3)))
  if (isTRUE(save)) {
    final_table_gt <- save_results_table(final_table, title = table_title, name =  paste0(number, "_table_", save_name))
    gtsave(final_table_gt, filename = paste0(number, "_table_", save_name, ".png"), "plots", vwidth = 3000)
  }
  
  close(pb) # close progessbar
  print(shapiro_table)
  return(final_table)
}

# add information to kegg_kos 
# takes a vector of kegg_kos

annotate_keggs <- function(kegg_kos) {
  pb <- txtProgressBar(min = 1, max = length(kegg_kos), style = 3) # initialize progessbar
  out_list <- list()
  for (i in 1:length(kegg_kos)) {
    ko <- kegg_kos[i]
    setTxtProgressBar(pb, i) # print progressbar
    
    out <- tryCatch(
      {
        kegg_info <- keggGet(ko)[[1]]
        tibble(kegg_ko = ko,
               name = kegg_info$NAME,
               pathways = str_c(kegg_info$PATHWAY, collapse = ","),
               pathway_ids = str_c(names(kegg_info$PATHWAY), collapse = ","),
               modules = str_c(kegg_info$MODULE, collapse = ","),
               module_ids = str_c(names(kegg_info$MODULE),collapse= ","),
               reactions = str_c(kegg_info$REACTION, collapse = ","),
               reaction_ids = str_c(names(kegg_info$REACTION), collapse = ","),
               brite = str_c(kegg_info$BRITE, collapse = ","))
      },
      error = function(e) {
        tibble(kegg_ko = ko,
               name = NA,
               pathways = NA,
               pathway_ids = NA,
               modules = NA,
               module_ids = NA,
               reactions = NA,
               reaction_ids = NA,
               brite = NA)
      }
    )
    
    out_list[[i]] <- out
  }
  out <- do.call("rbind", out_list)
  close(pb) # close progessbar
  
  return(out)
}


##############################
# functions from metabolomics
#############################

# function for preparation

prepare_metabolomics <- function(input) {
  output <- input %>%
    inner_join(dplyr::select(meta, sampleid, square, animal, period, diet)) %>%
    relocate(animal, period, square, diet) %>%
    dplyr::select(-sampleid)
  return(output)
}

# analysis function

loop_comparison_metabolimics <- function(input_df, y_axis, table_title, save_name, save = save) {
  shapiro_table <- data.frame("metabolite" = as.character(), "p_value"=as.numeric(), transformation = as.character()) # to check for bad transformations
  pb <- txtProgressBar(min = 1, max = ncol(input_df), style = 3) # initialize progessbar
  for (i in 5:ncol(input_df)) { #assuming 4 meta columns in the df
    setTxtProgressBar(pb, i) # print progressbar
    
    number <- get_script_number()
    
    current_response <- colnames(input_df)[i]
    # filter df for i'th response
    filtered_df <- select_response(input_df, current_response)
    # remove rows with NA (for functional activity)
    filtered_df <- filter(filtered_df, !is.na(response))
    # perform reml
    out <- combined_comparison(df = filtered_df, transformation = "test")
    # save shapiro results
    shapiro_table <- rbind(shapiro_table, c(current_response, round(out$shapiro$p.value, 3), out$transformation[1]))
    # print plot
    plot_pairwise(filtered_df = filtered_df, output = out, selected_response = current_response)
    if (isTRUE(save)) {
      save_big(paste0(number, "_", save_name, "_", current_response))
    }
    # print and save plot
    create_results_plot(input_df = filtered_df, comparison_object = out, response_name = current_response, y_axis = y_axis)
    if (isTRUE(save)) {
      save_big(paste0(number, "_result_", save_name, "_", current_response))
    }
    table <- create_results_table(input_df = filtered_df, comparison_object = out, response_name = current_response, digits = 1)
    table <- rbind(table, c("exact_p", out$ano[which(rownames(out$ano) == "diet"), 
                                               which(colnames(out$ano) == "Pr(>F)")])) # for exact p value adjustment
    
    if (i == 5) { # create table if first iteration, join afterwards
      final_table <- table
    } else {
      final_table <- final_table %>%
        left_join(table, by = "diet")
    }
    
  }
  final_table <- rbind(final_table[-7,], 
                       c("P-adj", 
                         round(p.adjust(final_table[7, 2:ncol(final_table)], method = "BH"), 3)))
  if (isTRUE(save)) {
    final_table_gt <- save_results_table(final_table, title = table_title, name =  paste0(number, "_table_", save_name))
    gtsave(final_table_gt, filename = paste0(number, "_table_", save_name, ".png"), "plots", vwidth = 3000)
  }
  
  close(pb) # close progessbar
  print(shapiro_table)
  return(final_table)
}




###############################
# Functions for functional activity
#######################################

# function to join a gene and a protein frame by sampleid and join_column2
# df's should be filtered for taxa level beforehand
# gets filtered per matrix for missingness < than missingness and imputes 0s by a small value

create_integrated_df <- function(df1, df2, join_column2, missingness = 0) {
  integrated <- inner_join(df1, df2, by = c("sampleid", join_column2))
  
  # filtering for ileum
  integrated_il <- filter_ileum(integrated) %>%
    group_by(!!sym(join_column2)) %>%
    filter((sum(rel_abd_prot == 0)/length(rel_abd_prot)) <= missingness & 
             (sum(rel_abd_gene == 0)/length(rel_abd_gene)) <= missingness) %>%
    ungroup()
  
  if (nrow(integrated_il) > 0) {
    # calculate discarded abundance
    integrated_il_discarded <- integrated_il %>%
      group_by(sampleid) %>%
      summarise(rel_abd_prot_disc = 100-sum(rel_abd_prot),
                rel_abd_gene_disc = 100-sum(rel_abd_gene), .groups = "drop") %>%
      pivot_longer(-sampleid, names_to = "parameter", values_to = "value") %>%
      left_join(meta, by = "sampleid")
    
    p <- ggplot(integrated_il_discarded, aes(x = diet, y = value)) +
      geom_boxplot() +
      geom_quasirandom() +
      facet_wrap(~parameter) +
      labs(title = "Discarded abundance ileal digesta")
    print(p)
  }

  
  # filtering for faeces
  integrated_fa <- filter_faeces(integrated) %>%
    group_by(!!sym(join_column2)) %>%
    filter((sum(rel_abd_prot == 0)/length(rel_abd_prot)) <= missingness & 
             (sum(rel_abd_gene == 0)/length(rel_abd_gene)) <= missingness) %>%
    ungroup()
  
  if (nrow(integrated_fa) > 0) {
    # calculate discarded abundance
    integrated_fa_discarded <- integrated_fa %>%
      group_by(sampleid) %>%
      summarise(rel_abd_prot_disc = 100-sum(rel_abd_prot),
                rel_abd_gene_disc = 100-sum(rel_abd_gene), .groups = "drop") %>%
      pivot_longer(-sampleid, names_to = "parameter", values_to = "value") %>%
      left_join(meta, by = "sampleid")
    
    p <- ggplot(integrated_fa_discarded, aes(x = diet, y = value)) +
      geom_boxplot() +
      geom_quasirandom() +
      facet_wrap(~parameter) +
      labs(title = "Discarded abundance faeces")
    print(p)
  }
  
  # filtering for ferment
  integrated_fe <- filter_ferment(integrated) %>%
    group_by(!!sym(join_column2)) %>%
    filter((sum(rel_abd_prot == 0)/length(rel_abd_prot)) <= missingness & 
             (sum(rel_abd_gene == 0)/length(rel_abd_gene)) <= missingness) %>%
    ungroup()
  
  if (nrow(integrated_fe) > 0) {
    # calculate discarded abundance
    integrated_fe_discarded <- integrated_fe %>%
      group_by(sampleid) %>%
      summarise(rel_abd_prot_disc = 100-sum(rel_abd_prot),
                rel_abd_gene_disc = 100-sum(rel_abd_gene), .groups = "drop") %>%
      pivot_longer(-sampleid, names_to = "parameter", values_to = "value") %>%
      left_join(meta, by = "sampleid")
    
    p <- ggplot(integrated_fe_discarded, aes(x = diet, y = value)) +
      geom_boxplot() +
      geom_quasirandom() +
      facet_wrap(~parameter) +
      labs(title = "Discarded abundance ileal digesta")
    print(p)
  }
  
  # no imputation but removal of Inf
  integrated <- rbind(integrated_il, integrated_fa, integrated_fe) %>%
    # imputation with a small value
    # mutate(rel_abd_prot = ifelse(rel_abd_prot == 0, 1e-8, rel_abd_prot), # impute with a low value to mitigate Inf
    #        rel_abd_gene = ifelse(rel_abd_gene == 0, 1e-8, rel_abd_gene)) %>%
    # imputation with minimum value
    # mutate(rel_abd_prot = ifelse(rel_abd_prot == 0, NA, rel_abd_prot), # replace 0 with NA
    #        rel_abd_gene = ifelse(rel_abd_gene == 0, NA, rel_abd_gene)) %>%
    # group_by(!!sym(join_column2)) %>%
    # mutate(rel_abd_prot = ifelse(is.na(rel_abd_prot), min(rel_abd_prot, na.rm = T), rel_abd_prot), # impute minimum without NA
    #        rel_abd_gene = ifelse(is.na(rel_abd_gene), min(rel_abd_gene, na.rm = T), rel_abd_gene)) %>%
    # ungroup() %>%
    mutate(func_act = log2(rel_abd_prot/rel_abd_gene)) %>%
    filter(func_act != Inf) %>%
    filter(func_act != -Inf)
  
  return(integrated)
}

# function creating a scatterplot plotting the protein vs the DNA abundance
# takes an integrated df with rel_abd prot and gene columns
# averages over all samples included

create_func_act_scatter_plot <- function(integrated_df) {
  feature_column <- colnames(integrated_df)[1]
  
  aggregated_df <- integrated_df %>%
    group_by(!!sym(feature_column)) %>%
    summarise(mean_rel_abd_prot = mean(rel_abd_prot),
              sd_rel_abd_prot = sd(rel_abd_prot),
              mean_rel_abd_gene = mean(rel_abd_gene),
              sd_rel_abd_gene = sd(rel_abd_gene), .groups = "drop") %>%
    mutate(rank_prot = row_number(desc(mean_rel_abd_prot)),
           rank_gene = row_number(desc(mean_rel_abd_gene)),
           label = ifelse(rank_prot <= 10 | rank_gene <= 10, !!sym(feature_column), ""))
  
  p <- ggplot(aggregated_df, aes(x = mean_rel_abd_gene, y = mean_rel_abd_prot)) +
    # geom_errorbar(aes(ymin = mean_rel_abd_prot - sd_rel_abd_prot, ymax = mean_rel_abd_prot + sd_rel_abd_prot),
    #               width = .05, lwd = 1) +
    # geom_errorbar(aes(xmin = mean_rel_abd_gene - sd_rel_abd_gene, xmax = mean_rel_abd_gene + sd_rel_abd_gene),
    #               width = .05, lwd = 1) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = label), color = "black", show.legend = F, 
                    max.overlaps = Inf, box.padding = 0.5,
                    min.segment.length = 0, max.time = 3, size = 3, nudge_y = -.1) +
    labs(x = "relative abundance DNA (%)", y = "relative abundance proteins (%)") +
    scale_x_log10(limits = c(min(aggregated_df$mean_rel_abd_prot, aggregated_df$mean_rel_abd_gene),
                             max(aggregated_df$mean_rel_abd_prot, aggregated_df$mean_rel_abd_gene)), 
                  expand = c(.1,.1)) +
    scale_y_log10(limits = c(min(aggregated_df$mean_rel_abd_prot, aggregated_df$mean_rel_abd_gene),
                             max(aggregated_df$mean_rel_abd_prot, aggregated_df$mean_rel_abd_gene)), 
                  expand = c(.1,.1)) +
    # scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
    # scale_x_continuous(limits = c(0,100), expand = c(0,0)) +
    geom_abline(intercept = 0, slope = 1) 
  return(p)
}

# function for filtering and plotting
# if for ferment (region) only diets 3 and 4 are used for pivotting

filter_results_table <- function(results_table, p_threshold, feature_name, region = c("ileum", "faeces", "ferment")) {
  if (region == "ferment") {
    diets = c("3", "4")
  } else {
    diets = c("1", "2", "3", "4")
  }
  
  filtered_table <- results_table %>%
    column_to_rownames("diet") %>%
    t() %>%
    as_tibble(rownames = feature_name) %>%
    mutate(`P-adj` = as.numeric(`P-adj`),
           `P-value` = as.numeric(`P-value`)) %>%
    filter(`P-adj` < p_threshold) %>%
    pivot_longer(diets, names_to = "diet", values_to = "letter") %>%
    mutate(mean = as.numeric(str_remove_all(letter, "[a-z]+$")),
           letter = str_extract_all(letter, "[a-z]+$"))
  return(filtered_table)
}

plot_functional_activity <- function(integrated_df, sig_results, region = c("ileum", "faeces"), name_column = NULL) {
  feature_column <- colnames(sig_results)[1]
  
  # if name_column not defined, define it as feature column
  if (is.null(name_column)) {
    name_column <- feature_column
  }
  
  if (region == "ferment") {
    filtered_df <- filter_ferment(integrated_df)
  } else if (region == "ileum") {
    filtered_df <- filter_ileum(integrated_df)
  } else {
    filtered_df <- filter_faeces(integrated_df)
  }
  
  plotting_df <- filtered_df %>%
    inner_join(dplyr::select(meta, sampleid, square, animal, period, diet), by = "sampleid") %>%
    inner_join(sig_results, by = c(feature_column, "diet"))
  
  p <- ggplot(plotting_df, aes(x = diet)) +
    facet_wrap(facets = vars(!!sym(name_column)), scales = "free_y") +
    geom_boxplot(aes(y = func_act), fill = "#A4BED5FF", outliers = F, alpha = .5) +
    geom_quasirandom(aes(y = func_act, size = rel_abd_prot), color = "#023743FF", width = 0.2, alpha = .7) +
    geom_quasirandom(aes(y = func_act, size = rel_abd_gene), color = "#FED789FF", show.legend = F, width = 0.2, 
                     alpha = .5) +
    labs(size = "rel. abd. %", y = "log2 protein/gene ratio") +
    geom_text(aes(x = ((length(unique(diet))+1) /2), y = max(func_act), label = paste("q = ", `P-adj`))) +
    geom_text(aes(y = mean, label = letter))
  print(p)
  return(dplyr::select(plotting_df, !!sym(feature_column), sampleid, func_act))
}

# function taking an integrated frame with columns for rel_abd_gene and rel_abd_prot and precalculated func_act
# frames should be filtered for region
# testing dna agains protein rel_abd with wilcoxon test and creating a volcano plot based on mean func_act and adj p

test_protein_against_dna <- function(input) {
  # remove Inf values
  input <- input %>%
    filter(func_act != Inf)
  
  first_column <- colnames(input)[1]
  
  iterations <- unique(input[[first_column]])
  
  output <- tibble(!!sym(first_column) := iterations, func_act = NA, p_value = NA, p_adj = NA, 
                   rel_abd_prot = NA, rel_abd_gene = NA)
  
  for (i in 1:length(iterations)) {
    iteration <- iterations[i]
    
    comparison <- input %>%
      filter(!!sym(first_column) == iteration)
    
    wilcox <- wilcox.test(comparison$rel_abd_gene, comparison$rel_abd_prot, paired = T)
    
    p <- wilcox$p.value
    
    output$p_value[i] <- p
    
    output$func_act[i] <- median(comparison$func_act) # changed to median
    output$rel_abd_prot[i] <- median(comparison$rel_abd_prot)
    output$rel_abd_gene[i] <- median(comparison$rel_abd_gene)
    
  }
  
  # adjust p_value
  
  output$p_adj <- p.adjust(output$p_value, method = "holm")
  
  output <- output %>%
    mutate(sig = ifelse(p_adj < 0.05, "sig", "not sig")) %>%
    mutate(rank_asc = row_number(func_act),
           rank_desc = row_number(-func_act),
           label = ifelse(sig == "sig" & (rank_desc <= 10 | rank_asc <= 10), !!sym(first_column), ""))
  
  p <- ggplot(output, aes(x = func_act, y = -log10(p_value))) +
    geom_vline(xintercept = 0, color = "grey20") +
    geom_hline(yintercept = 0, color = "grey20") + 
    geom_point(aes(size = rel_abd_prot), show.legend = T, color = "#023743FF", alpha = .7) +
    geom_point(aes(size = rel_abd_gene), show.legend = F, color = "#FED789FF") +
    geom_text_repel(aes(label = label), color = "black", show.legend = F, max.overlaps = Inf, box.padding = 0.5,
                    min.segment.length = 0, max.time = 3, size = 3, nudge_y = -.1) +
    scale_color_manual(values = c("grey20", "red")) +
    #scale_size_manual(values = c(1,2)) +
    labs(x = "log2 fold change", y = "-log10 P-value", size = "rel. abd. %") 
  
  print(p)
  
  return(output)
}

# functional distance between metagenomic and metaproteomic frame
# func dist function

calculate_func_dist <- function(integrated_df) {
  
  func_dist <- tibble(sampleid = unique(integrated_df$sampleid),
                      func_dist = NA)
  
  for (i in 1:nrow(func_dist)) {
    current_sampleid <- func_dist$sampleid[i]
    
    filtered_matrix <- integrated_df %>%
      filter(sampleid == current_sampleid) %>%
      dplyr::select(rel_abd_gene, rel_abd_prot) %>%
      as.matrix() %>%
      t()
    
    bray <- vegdist(filtered_matrix)
    
    func_dist$func_dist[i] <- bray
    
  }
  
  func_dist_out <- func_dist %>%
    inner_join(meta, by = "sampleid")
  
  p <- ggplot(func_dist_out, aes(x = diet, y = func_dist))+
    geom_boxplot(outliers = F) +
    geom_quasirandom(aes(color = animal)) +
    facet_grid(~matrix) +
    scale_color_manual(values = colors)
  print(p)
  
  return(func_dist_out)
}

###########################################
#### Correlation function
###########################################

# takes one df and correlates each column with each other

correlate_sym <- function(df, method = "pearson", save_name, save = FALSE) {
  if (!is.matrix(df)) {
    matrix <- as.matrix(df[,-1])
  } else {
    matrix <- df
  }
  corr_matrix <- corr.test(matrix, method = method, adjust = "fdr")
  annotation <- corr_matrix$p
  colnames(annotation) <- c(1:ncol(annotation))
  annotation[which(lower.tri(annotation, diag = T))] <- 1
  annotation_list <- annotation %>%
    as_tibble(rownames = "y") %>%
    mutate(y = nrow(annotation) - row_number() + 1) %>%
    pivot_longer(-y, names_to = "x", values_to = "value") %>%
    mutate(value = case_when(value < 0.001 ~ "***",
                             value < 0.01 ~ "**",
                             value < 0.05 ~ "*",
                             value < 0.1 ~ ".",
                             .default = " "))
  if (isTRUE(save)) jpeg(filename = paste0("plots/", save_name, ".jpeg"), 
                         width =15, height = 10, unit="cm", res = 1000)
  corrplot.mixed(corr_matrix$r, upper = "shade", tl.pos = "lt", diag = "n", tl.col = "black", tl.srt = 45)
  text(x = annotation_list$x, y = annotation_list$y, label = annotation_list$value, cex = 1)
  if (isTRUE(save)) dev.off()
  corr_out <- corr_matrix$r
  corr_out[corr_matrix$p>=0.05] <- "ns"
  corr_out <- as_tibble(corr_out, rownames = "pea variety")
  return(corr_out)
}

# takes df's and checks for sampleno differences

correlate_unsym <- function(df1, df2, method = "pearson", save_name, save = FALSE, meta = NULL, exclude_diet =NULL,
                            p.adjust = "fdr", filter_significant = T) {
  
  if (!is.null(exclude_diet)) {
    samplenos <- meta %>% 
      filter(!diet %in% c(exclude_diet)) %>% # exclude diet X and pull sampleids
      pull(sampleno)
    
    df1 <- filter(df1, sampleno %in% samplenos) # filter dfs to only contain filtered sampleids
    df2 <- filter(df2, sampleno %in% samplenos)
  }
  
  sampleno_intersect <- intersect(df1$sampleno, df2$sampleno) # to account for missing samples
  
  if (!is.matrix(df1)) {
    matrix1 <- df1 %>%
      filter(sampleno %in% sampleno_intersect) %>%
      column_to_rownames("sampleno") %>%
      as.matrix()
  } else {
    matrix1 <- df1
  }
  if (!is.matrix(df2)) {
    matrix2 <- df2 %>%
      filter(sampleno %in% sampleno_intersect) %>%
      column_to_rownames("sampleno") %>%
      as.matrix()
  } else {
    matrix2 <- df2
  }
  
  corr_matrix <- corr.test(matrix1, matrix2, method = method, adjust = p.adjust)
  annotation <- corr_matrix$p.adj
  
  # filter to only significant, if true
  if (isTRUE(filter_significant)) {
    corr_matrix$r <- corr_matrix$r[c(which(apply(annotation, 1, FUN = min, na.rm = T) < 0.05)), 
                                   c(which(apply(annotation, 2, FUN = min, na.rm = T) < 0.05)), drop = FALSE]
    annotation <- annotation[which(apply(annotation, 1, FUN = min, na.rm = T) < 0.05), 
                             which(apply(annotation, 2, FUN = min, na.rm = T) < 0.05), drop = FALSE]
  }
  
  colnames(annotation) <- c(1:ncol(annotation))
  annotation_list <- annotation %>%
    as_tibble(rownames = "y") %>%
    mutate(y = nrow(annotation) - row_number() + 1) %>%
    pivot_longer(-y, names_to = "x", values_to = "value") %>%
    mutate(value = case_when(value < 0.001 ~ "***",
                             value < 0.01 ~ "**",
                             value < 0.05 ~ "*",
                             value < 0.1 ~ ".",
                             .default = " "))
  if (isTRUE(save)) jpeg(filename = paste0("plots/", save_name, ".jpeg"), 
                         width =15, height = 10, unit="cm", res = 1000)
  corrplot(corr_matrix$r, method = "shade", tl.col = "black", tl.srt = 45)
  text(x = annotation_list$x, y = annotation_list$y-0.2, label = annotation_list$value, cex = 1)
  if (isTRUE(save)) dev.off()
  corr_out <- corr_matrix$r
  corr_out[annotation >=0.05] <- "ns"
  corr_out <- as_tibble(corr_out, rownames = "pea variety")
  return(corr_out)
}

# function to perform correlations dietwise

correlate_unsym_multiple <- function(df1, df2, meta, by = "diet", method = "pearson", save_name, save = FALSE) {
  levels <- sort(unique(meta[[by]])) # define levels
  
  par(mfrow = c(1,length(levels))) # set plotting window
  
  matrix_list <- list()
  
  for (i in 1:length(levels)) {
    samplenos <- meta %>%
      filter(!!sym(by) == i) %>%
      pull(sampleno)
    
    df1sub <- filter(df1, sampleno %in% samplenos)
    df2sub <- filter(df2, sampleno %in% samplenos)
    
    matrix <- correlate_unsym(df1 = df1sub, df2 = df2sub, method = method, save = FALSE, p.adjust = "none")
    matrix_list[[i]] <- matrix %>%
      mutate_all(~ifelse(.=="ns", 0, .)) %>% 
      column_to_rownames("pea variety") %>%
      mutate_all(as.numeric) %>%
      as.matrix()
  }
  
  par(mfrow = c(1,1))
  
  # multiply matrices
  for (i in 2:length(levels)) {
    if (i == 2) {
      matrix_out = matrix_list[[1]] * matrix_list[[i]]
    } else {
      matrix_out = matrix_out * matrix_list[[i]]
    }
  }
  matrix_out <- matrix_out[which(rowSums(matrix_out)>0), which(colSums(matrix_out)>0), drop = FALSE]
  
  corrplot(matrix_out)
  return(matrix_out)
}

#create a single scatterplot with labels for diet and maybe other factors

single_correlation_plot <- function(df1, df2, variable1, variable2, meta, second_indicator = "animal") {
  sampleno_intersect <- intersect(df1$sampleno, df2$sampleno) # to account for missing samples
  
  df1 <- dplyr::select(df1, sampleno, !!sym(variable1))
  df2 <- dplyr::select(df2, sampleno, !!sym(variable2))
  
  df <- inner_join(df1, df2, by = "sampleno") %>%
    #dplyr::select(sampleno, !!sym(variable1), !!sym(variable2)) %>%
    left_join(meta, by = "sampleno")
  
  cor <- cor.test(df[[variable1]], df[[variable2]])
  print(cor)
  
  ggplot(df, aes(x = !!sym(variable1), y = !!sym(variable2))) +
    geom_point(aes(shape = diet, color = !!sym(second_indicator)), size = 5) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(16,17,15,18,21,24,22,23)) +
    geom_smooth(method = "lm", se = F)
}


###########################################
### Multiomics functions
############################################

# merge nutrition tables with meta and select relevant columns

prepare_nutrition_for_correlation <- function(tibble) {
  tibble_prepared <- meta_nutrition %>%
    inner_join(tibble, join_by(animal, period)) %>%
    dplyr::select(-c(animal, period, square, diet, description)) %>%
    dplyr::select(where(~ !is.numeric(.) || n_distinct(.) > 1)) %>% # exclude numeric columns with all the same value
    arrange(sampleno)
  return(tibble_prepared)
}

# transforms long format df into a wide format df with features as column names and sampleno as first column
# sampleid is transformed into sampleno and sorted

prepare_omics_for_correlation <- function(input_df, abundance_column) {
  first_column <- colnames(input_df)[1] # 1st column should be feature column
  output <- input_df %>%
    mutate(sampleno = str_extract(sampleid, "[0-9]{2}$")) %>%
    dplyr::select(sampleno, !!sym(first_column), !!sym(abundance_column)) %>%
    pivot_wider(names_from = !!sym(first_column), values_from = !!sym(abundance_column)) %>%
    arrange(sampleno)
  return(output)
}

# preparing data for mixomics

prepare_for_mixomics <- function(df, sampleno_intersect) {
  mat <- df %>%
    filter(sampleno %in% sampleno_intersect) %>%
    arrange(sampleno) %>%
    column_to_rownames("sampleno") %>%
    as.matrix()
  mat <- mat[,which(!colSums(mat) == 0)]
  return(mat)
}

# filter function for binary outcome

filter_input_diablo <- function(omics_list, meta, include) {
  meta_filtered <- meta %>%
    filter(diet %in% include) %>%
    mutate(diet = as.factor(diet)) #as factor for PLS
  
  list_filtered <- omics_list
  for (i in 1:length(list_filtered)) {
    list_filtered[[i]] <- list_filtered[[i]][rownames(list_filtered[[i]]) %in% meta_filtered$sampleno,]
  }
  
  output_list <- list()
  output_list[[1]] <- list_filtered
  output_list[[2]] <- meta_filtered
  return(output_list)
}

# detect pairwise correlations
diablo_pairwise_correlations <- function(omics_list) {
  cor_matrix = matrix(ncol = length(omics_list), nrow = length(omics_list), 
                  dimnames = list(names(omics_list), names(omics_list)))
  diag(cor_matrix) <- 0
  for (i in 1:length(omics_list)) {
    X <- omics_list[[i]]
    for (j in 1:length(omics_list)) {
      Y <- omics_list[[j]]
      if (j > i) {
        print(paste(names(omics_list)[i], "vs", names(omics_list)[j]))
        pls <- spls(X, Y, keepX = c(min(ncol(X), 20),min(ncol(X), 20)), 
                    keepY = c(min(ncol(Y), 20),min(ncol(Y), 20)))
        plotVar(pls, var.names = T, cutoff = 0.5)
        print(cor(pls$variates$X, pls$variates$Y))
        cor_value <- cor(pls$variates$X, pls$variates$Y)[1,1]
        cor_matrix[i,j] <- cor_value
        cor_matrix[j,i] <- cor_value
        Sys.sleep(2)
      }
    }
  }

  return(cor_matrix)
}

# set up and test model, return final model

diablo_auto_model <- function(omics_list, meta, design_value, ncomp_value = NULL) {
  # measure execution time
  start_time <- Sys.time()
  
  #create design matrix # just take the given matrix
  # design = matrix(design_value, ncol = length(omics_list), nrow = length(omics_list), 
  #                 dimnames = list(names(omics_list), names(omics_list)))
  # diag(design) = 0 # set diagonal to 0s
  design <- design_value
  
  # initialise 
  
  ncomp_init <- ifelse(length(meta) > 10, 10, 3) # for ferment a lower number of comps is necessary
  diablo_init <- block.splsda(X = omics_list, Y = meta, ncomp = ncomp_init, design = design)
  
  print(diablo_init)
  
  # test for best number of comps
  diablo_perf <- perf(diablo_init, validation = "loo",
                      progressBar = T, seed = 1112, nrepeat = 5) # nrepeat to obtain optimal ncomp numbers
  
  par(mfrow = c(1,1))
  plot(diablo_perf)
  
  print(diablo_perf$choice.ncomp)
  
  if (is.null(ncomp_value)) {
    #ncomp = min(unlist(diablo_perf$choice.ncomp))
    ncomp = as.integer(readline(prompt = "Enter number of components"))
    print(paste("Choice of components:", ncomp))
  } else {
    ncomp = ncomp_value
  }
  
  
  # test different numbers of features
  test.keepX <- list()
  for (i in 1:length(omics_list)) {
    test.keepX[[names(omics_list)[i]]] = seq(5, min(ncol(omics_list[[i]]), 45), 10)
  }
  
  diablo_tune <- tune.block.splsda(X = omics_list, Y = meta,
                                   ncomp = ncomp, test.keepX = test.keepX,
                                   design = design, validation = "loo",
                                   dist = "mahalanobis.dist",
                                   BPPARAM = BiocParallel::SnowParam(workers = 8), # 8 threads
                                   progressBar = T, seed = 1112)
  plot(diablo_tune)
  
  list.keepX <- diablo_tune$choice.keepX
  print("Choice of features:")
  print(list.keepX)
  
  # final model
  
  diablo_final <- block.splsda(X = omics_list, Y = meta, ncomp = ncomp,
                               keepX = list.keepX, design = design)
  
  # measure execution time
  end_time <- Sys.time()
  print(end_time-start_time)
  
  return(diablo_final)
}

# Create every plot with final diablo model

plot_diablo <- function(final_model, interactive_network = FALSE, cutoff_value = 0.8, save_name = NULL, 
                        return_network = FALSE, return_full_matrix = FALSE) {
  plotDiablo(final_model, ncomp = min(final_model$ncomp[1], 2))
  if (final_model$ncomp[1] > 1) {
    plotIndiv(final_model)
    plotArrow(final_model)
    plotVar(final_model)
  }
  
  par(mfrow = c(1,1))
  
  # helper to select use components for network
  comps <- names(final_model$ncomp)[1:(length(final_model$ncomp)-1)]
  list_comps <- list()
  for (i in 1:length(comps)) {
    list_comps[[i]] <- seq(1, final_model$ncomp[1])
    names(list_comps)[i] <- comps[i]
  }
  
  matrix <- circosPlot(final_model, cutoff = cutoff_value, size.variables = 1, line = T, size.labels = 1.5)
  network <- network(final_model, cutoff = cutoff_value, blocks = seq(1, length(final_model$names$blocks)-1),
                     comp = list_comps, interactive = interactive_network,
                     size.node = 0.1, cex.node.name = 0.5, lwd.edge = 1, symkey = T,
                     color.edge = c(color.jet(100)[1:25], color.jet(100)[76:100]),
                     save = if (is.null(save_name)) save_name else "jpeg",
                     name.save = if (is.null(save_name)) NULL else { 
                       paste0("plots/", get_script_number(), "_diablo_network_", save_name)})
  
  if (!isTRUE(return_full_matrix)) {
    matrix[which(abs(matrix)<cutoff_value)] <- 0
    matrix <- matrix[which(rowSums(matrix) != 0),
                     which(colSums(matrix) != 0)]
  }
  
  pheatmap(matrix, col = color.jet(100), scale = "none", fontsize = 8)
  
  plotLoadings(final_model, comp = 1)
  plotLoadings(final_model, comp = 2)
  #cimDiablo(final_model, legend.position = "left", transpose = T)
  if (isTRUE(return_network)) {
    return(network)
  } else if (isTRUE(return_full_matrix)) {
    return(circosPlot(final_model, cutoff = 0))
  } else {
    return(matrix)
  }
}

# filter correlation matrix to a subnetwork given a vector of names and create a heatmap
# if a component is given, the diablo model has to be supplied and the vector is created based on the features in the component

filter_submatrix <- function(matrix, vector, component = NULL, final_model) {
  if (!is.null(component)) {
    vector <- c()
    for (i in 1:(length(final_model$loadings)-1)) {
      loadings <- final_model$loadings[[i]][,component, drop = F]
      names <- rownames(loadings[which(abs(loadings) > 0),,drop = F])
      vector <- c(vector, names) # append to vector
    }
  }
  not_found <- vector[!vector %in% rownames(matrix)]
  if (length(not_found) > 0){
    print(paste("Not found:", not_found))
  }
  
  index_vector <- which(rownames(matrix) %in% vector)
  submatrix <- matrix[index_vector, index_vector]
  
  # remove features not correlated
  diag(submatrix) <- diag(submatrix) / 100 # helper to mitigate filtering due to self correlation
  submatrix[which(abs(submatrix)<.2)] <- 0 # .2 should be below every threshold
  submatrix <- submatrix[which(rowSums(submatrix) != 0),
                   which(colSums(submatrix) != 0)]
  #diag(submatrix) <- diag(submatrix) * 100
  
  pheatmap(submatrix)
  return(submatrix)
}

# create a igraph from a adjacency matrix

create_igraph_from_matrix <- function(matrix, save_name = NULL) {
  g <- graph_from_adjacency_matrix(matrix, mode = "undirected", weighted = T, diag = F)
  layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
  # assign colors to weight values
  color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
  E(g)$color <- color.jet(100)[color_index]
  E(g)$width <- .5
  node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
  color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                            str_detect(V(g)$name, "^g_") ~ 2,
                            str_detect(V(g)$name, "^p_") ~ 3,
                            str_detect(V(g)$name, "^ssc_") ~ 4,
                            str_detect(V(g)$name, "_PEA$") ~ 5,
                            str_detect(V(g)$name, "^m_") ~ 6,
                            .default = 0)
  V(g)$color <- node_color[color_vector]
  V(g)$label.family <- "sans"
  V(g)$label.font <- 1
  V(g)$label.color <- "black"
  V(g)$label.cex <- .15
  if (!is.null(save_name)) jpeg(filename = paste0("plots/", get_script_number(), "_igraph_", save_name, ".jpeg"),
                                width = 4000, height = 4000, unit="px", res = 1200, pointsize = 20)
  par(mar = c(0,0,0,0))
  plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = TRUE,
       vertex.shape = "rectangle", 
       vertex.size = str_width(V(g)$name)*2,
       vertex.size2 = 4,
       vertex.frame.width = .1)
  legend("topright",
         legend = c(1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
         fill = rev(color.jet(11)),
         border = NA,
         y.intersp = 0.5,
         cex = .2, text.font = 1)
  legend("bottomright",
         legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
         fill = node_color,
         cex = .2)
  if (!is.null(save_name)) dev.off() 
}



##################
# Functional reduncancy
####################

# Function that calculates functional redundancy
# takes a table with proteinid, respective taxonomic classification (G), functional classification (kegg_cog), sampleid, and relative abundance
# takes table (on genus level) with name, sampleid, and relative abundance
# meta with sampleid needed

calculate_functional_redundancy <- function(pro_gen_cog, genus_table, meta) {
  # build PCN
  
  for (i in unique(meta$sampleid)) {
    #print(i)
    pro_gen_cog_sam_sum_wide <- pro_gen_cog %>%
      dplyr::select(G, kegg_cog, sampleid, rel_abd) %>%
      filter(sampleid == i) %>%
      filter(!is.na(kegg_cog)) %>% 
      group_by(G, kegg_cog) %>%
      summarise(sum = sum(rel_abd)) %>%
      pivot_wider(names_from = G, values_from = sum)
    assign(paste0("pcn", i), pro_gen_cog_sam_sum_wide)
  }
  
  # build 01 PCNs
  # 
  # for (i in unique(meta$sampleid)) {
  #   print(i)
  #   PCN_table_single_filtered_sum_col <- get(paste0("pcn", i)) %>%
  #     replace(is.na(.), 0) %>%
  #     pivot_longer(-kegg_cog, names_to = "G", values_to = "value") %>%
  #     group_by(G) %>%
  #     filter(sum(value) > 0) %>%
  #     ungroup() %>%
  #     group_by(kegg_cog) %>%
  #     filter(sum(value) > 0) %>%
  #     ungroup() %>%
  #     mutate(value = ifelse(value == 0, 0, 1)) %>%
  #     pivot_wider(names_from = kegg_cog, values_from = value)
  #   assign(paste0("pcn01", i), PCN_table_single_filtered_sum_col)
  # }
  
  # dij calculation
  
  dij <- function(mat, x, y) {
    min <- sum(pmin(mat[which(rownames(mat) == x),], mat[which(rownames(mat) == y),]))
    max <- sum(pmax(mat[which(rownames(mat) == x),], mat[which(rownames(mat) == y),]))
    dij <- 1 - (min / max)
    return(dij)
  }
  
  
  pb <- txtProgressBar(min = 1, max = length(unique(meta$sampleid)), style = 3) # initialize progessbar
  for (idx in 1:length(unique(meta$sampleid))) {
    i <- unique(meta$sampleid)[idx]
    setTxtProgressBar(pb, idx) # print progessbar
    PCN <- get(paste0("pcn", i)) %>%
      replace(is.na(.), 0) %>%    
      pivot_longer(-kegg_cog, names_to = "G", values_to = "value") %>%
      group_by(kegg_cog) %>%
      filter(sum(value) > 0) %>%
      ungroup() %>%
      group_by(G) %>%
      filter(sum(value) > 0) %>%
      mutate(norm = value / sum(value)) %>%
      ungroup() %>%
      dplyr::select(-value) %>%
      pivot_wider(names_from = kegg_cog, values_from = norm)
    PCNmat <- as.matrix(PCN[,-1])
    rownames(PCNmat) <- PCN$G
    dijmat <- matrix(nrow = nrow(PCNmat), ncol = nrow(PCNmat))
    rownames(dijmat) <- PCN$G
    colnames(dijmat) <- PCN$G
    
    for (j in rownames(PCNmat)) {
      #print(j)
      for (k in rownames(PCNmat)) {
        #print(k)
        dijmat[j,k] <- dij(mat = PCNmat, x = j, y = k)
      }
    }
    assign(paste0("dij", i), dijmat)
  }
  close(pb) # close progessbar
  
  # FR calculation
  
  genus_table_p_matrix <- genus_table %>%
    group_by(sampleid) %>%
    mutate(rel_abd = rel_abd/sum(rel_abd)) %>%
    pivot_wider(names_from = sampleid, values_from = rel_abd)
  
  output <- matrix(ncol = 4, nrow = length(unique(meta$sampleid))) # matrix to store output
  rownames(output) <- unique(meta$sampleid)
  colnames(output) <- c("FR", "nFR", "GSI", "FD")
  
  for (i in unique(meta$sampleid)) {
    dijmat <- get(paste0("dij", i))
    genus_table <- genus_table_p_matrix %>% # discarding genera not in dij
      filter(name %in% rownames(dijmat)) %>%
      mutate(name = factor(name, levels = rownames(dijmat))) %>%
      arrange(name)
    for (z in 1:length(genus_table$name)) {
      if (genus_table$name[z] != rownames(dijmat)[z]) {
        print("Names do not match!")
      }
    }
    genus_table$pipj <- 0
    genus_table$dijpipj <- 0
    
    for (j in row_number(genus_table)) {
      genus_table$pipj[j] <- 
        sum(genus_table[[i]][j] * 
              genus_table[[i]][-j])
      
      genus_table$dijpipj[j] <- 
        sum(genus_table[[i]][j] * genus_table[[i]][-j] * 
              (1 - dijmat[which(rownames(dijmat) == genus_table$name[j]),
                          which(rownames(dijmat) != genus_table$name[j])]))
    }
    GSI <- sum(genus_table$pipj) 
    FR <- sum(genus_table$dijpipj)
    output[which(rownames(output) == i), "FR"] <- FR
    output[which(rownames(output) == i), "nFR"] <- FR / GSI
    output[which(rownames(output) == i), "GSI"] <- GSI
    output[which(rownames(output) == i), "FD"] <- GSI - FR
  }
  
  return(output)
}

