# PiFerm3_manuscript_fermentation

## Overview

This repository contains the code for the paper "Selection and Evaluation of a Starter Culture for the Fermenta-tion of Field Pea Grains in Swine Nutrition".
The code is separated in three directories. The directory "16S" contains a jupyter notebook with the workflow used for processing of raw reads using Qiime2. The directory "DB_construction" contains shell scripts that were used to assemble and annotate the protein database for metaproteomics. The directory "Statistics" contains all R scripts used for data analysis including data wrangling, statistical analysis, and visualization.

### `16S`

| File | Description |
|---|---|
| `20260826_Qiime2_ferment.ipynb` | QIIME 2 workflow for processing paired-end 16S sequencing data. |

### `DB_construction`

| File | Description |
|---|---|
| `081_ferment_database_drep.sh` | Builds dereplicated genome databases with dRep. |
| `082_ferment_prokka.sh` | Annotates dereplicated genomes with Prokka. |
| `083_ferment_protein_eggnog.sh` | Assigns protein functions using eggNOG-mapper. |
| `084_ferment_db_creation.sh` | Packages genome sequences and annotations into analysis databases. |
| `085_ferment_gtdbtk.sh` | Classifies genomes taxonomically with GTDB-Tk. |

### `Statistics`

| File | Description |
|---|---|
| `0_general_functions.R` | Shared plotting, modelling, transformation, and utility functions. |
| `01_data_cleaning.R` | Cleans metadata and NMR/metabolomics data. |
| `05_data_cleaning_16S.R` | Processes 16S taxonomy, filters features, calculates abundances, and creates aggregated tables. |
| `06_data_cleaning_metaproteomics.R` | Cleans, normalizes, imputes, and annotates metaproteomics data. |
| `14_fermentation_strains_prokka+deep.R` | Compares predicted enzyme activities across fermentation strains. |
| `20_pretrial_functions.R` | Provides reusable ANOVA, transformation, comparison, plotting, and table-generation functions. |
| `21_pretrial1.R` | Performs statistical analysis and generates plots for the first fermentation experiment. |
| `23_screening_trial.R` | Analyzes pH, ammonia, acetate, and lactate measurements from the second fermentation experiment. |
| `24_combinations_trial.R` | Analyzes bacterial/enzyme combinations and performs principal component analysis for the third fermentation experiment. |
| `50_omics_functions.R` | Provides shared functions for omics data processing, statistics, visualization, and annotation. |
| `57_16S_ferment.R` | Performs 16S diversity, ordination, taxonomy, and differential-abundance analyses for fermentation samples. |
| `62_ferment.R` | Analyzes fermentation metaproteomics, taxonomy, functional annotations, diversity, and differential abundance. |
| `71_metabolomics_basic.R` | Creates basic summaries and plots of NMR metabolites. |
| `72_metabolomics_differential_analysis.R` | Performs differential metabolomics comparisons and writes result tables. |
| `95_Fermentation_plots.R` | Generates the main and supplementary figures for the manuscript. |
