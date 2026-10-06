# PiFerm3_manuscript_fermentation

## Overview

This repository contains the code for the paper "Selection and Evaluation of a Starter Culture for the Fermenta-tion of Field Pea Grains in Swine Nutrition".
The code is separated in three directories. The directory "16S" contains a jupyter notebook with the workflow used for processing of raw reads using Qiime2. The directory "DB_construction" contains shell scripts that were used to assemble and annotate the protein database for metaproteomics. The directory "Statistics" contains all R scripts used for data analysis including data wrangling, statistical analysis, and visualization.

### 16S
  20260826_Qiime2_ferment.ipynb   QIIME 2 workflow for paired-end 16S sequencing

### DB_construction
  081_ferment_database_drep.sh     Builds dereplicated genome databases with dRep
  082_ferment_prokka.sh            Annotates dereplicated genomes with Prokka
  083_ferment_protein_eggnog.sh    Assigns protein functions with eggNOG-mapper
  084_ferment_db_creation.sh       Packages genome and annotation files into databases
  085_ferment_gtdbtk.sh            Classifies genomes taxonomically with GTDB-Tk

### Statistics
  0_general_functions.R            Shared plotting, modelling, transformation, and utility functions
  01_data_cleaning.R               Metadata and NMR/metabolomics cleaning
  05_data_cleaning_16S.R           16S taxonomy, filtering, abundance, and aggregation
  06_data_cleaning_metaproteomics.R Metaproteomics cleaning, normalization, imputation, and annotation
  14_fermentation_strains_prokka+deep.R
                                    Compares predicted enzyme activities across fermentation strains
  20_pretrial_functions.R          Reusable ANOVA, transformation, comparison, plotting, and table functions
  21_pretrial1.R                   Statistical analysis and plots for the first fermentation experiment
  23_screening_trial.R             Analysis of pH, ammonia, acetate, and lactate for the second fermentation experiment
  24_combinations_trial.R          Analysis of bacterial/enzyme combinations and PCA for the third fermentation experiment
  50_omics_functions.R             Shared omics-analysis functions
  57_16S_ferment.R                 Diversity, ordination, taxonomy, and differential-abundance analysis of 16S data
  62_ferment.R                     Metaproteomics, taxonomy, functions, diversity, and differential analysis for ferments
  71_metabolomics_basic.R          Basic NMR metabolite summaries and plots
  72_metabolomics_differential_analysis.R
                                    Differential metabolomics comparisons and output tables
                                    
  95_Fermentation_plots.R          Generates the main and supplementary manuscript figures
