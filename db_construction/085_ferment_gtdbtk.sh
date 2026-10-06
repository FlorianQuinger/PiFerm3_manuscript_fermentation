#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=500g
#SBATCH --job-name=085_gtdbtk
#SBATCH --output=logs/085_gtdbtk.out
#SBATCH --error=logs/085_gtdbtk.err
cd /pfs/10/work/ho_yogau97-PiFerm3/PiFerm3_metagenomics

module load devel/miniforge
conda activate gtdbtk_v2.5.2

#environment variable for database
conda env config vars set GTDBTK_DATA_PATH="/pfs/10/project/db/gtdb/release226"

#environment variables
DB=/pfs/10/project/db/gtdb/release226

DBdrepMiFoDB=00_data/08_db_drep_MiFoDB/dereplicated_genomes
DBdrepcFMD=00_data/08_db_drep_cFMD/dereplicated_genomes
DBdrepdual=00_data/08_db_drep_dual/dereplicated_genomes

OUTMiFoDB=00_data/08_gtdbtk_MiFoDB
OUTcFMD=00_data/08_gtdbtk_cFMD
OUTdual=00_data/08_gtdbtk_dual

OUTMiFoDBreport=02_report/08_gtdbtk_MiFoDB
OUTcFMDreport=02_report/08_gtdbtk_cFMD
OUTdualreport=02_report/08_gtdbtk_dual

#create directory
mkdir -p $OUTMiFoDB
mkdir -p $OUTcFMD
mkdir -p $OUTdual

mkdir -p $OUTMiFoDBreport
mkdir -p $OUTcFMDreport
mkdir -p $OUTdualreport

# reactivate environment
conda activate gtdbtk_v2.5.2

#run gtdb classify workflow on different dbs
gtdbtk classify_wf --genome_dir $DBdrepMiFoDB --out_dir $OUTMiFoDB --cpus 32 -x fna --skip_ani_screen #--mash_db $DB/gtdbtk_data.msh 
gtdbtk classify_wf --genome_dir $DBdrepcFMD --out_dir $OUTcFMD --cpus 32 -x fna --skip_ani_screen #--mash_db $DB/gtdbtk_data.msh 
gtdbtk classify_wf --genome_dir $DBdrepdual --out_dir $OUTdual --cpus 32 -x fna --skip_ani_screen #--mash_db $DB/gtdbtk_data.msh 

# copy output files to reports

cp $OUTMiFoDB/*.summary.tsv $OUTMiFoDBreport/
cp $OUTcFMD/*.summary.tsv $OUTcFMDreport/
cp $OUTdual/*.summary.tsv $OUTdualreport/

