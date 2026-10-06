#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=081_db_drep
#SBATCH --output=logs/081_db_drep.out
#SBATCH --error=logs/081_db_drep.err
cd /pfs/10/work/ho_yogau97-PiFerm3/PiFerm3_metagenomics

module load devel/miniforge
conda activate bbmap_v39.33

#environment variables
DBown=02_report/08_mags_own
DBMiFoDB=02_report/08_mags_MiFoDB
DBcFMD=02_report/08_mags_cFMD

DBallMiFoDB=00_data/08_db_raw_MiFoDB
DBallcFMD=00_data/08_db_raw_cFMD
DBalldual=00_data/08_db_raw_dual

DBdrepMiFoDB=00_data/08_db_drep_MiFoDB
DBdrepcFMD=00_data/08_db_drep_cFMD
DBdrepdual=00_data/08_db_drep_dual

DBoutMiFoDB=02_report/08_db_drep_MiFoDB
DBoutcFMD=02_report/08_db_drep_cFMD
DBoutdual=02_report/08_db_drep_dual

#create directory
mkdir -p $DBallMiFoDB
mkdir -p $DBallcFMD
mkdir -p $DBalldual

mkdir -p $DBdrepMiFoDB
mkdir -p $DBdrepcFMD
mkdir -p $DBdrepdual

mkdir -p $DBoutMiFoDB
mkdir -p $DBoutcFMD
mkdir -p $DBoutdual

# copy bins to respective folders # renaming due to contig name length limit of prokka

cp $DBown/* $DBallMiFoDB
cp $DBown/* $DBallcFMD
cp $DBown/* $DBalldual

i=1
for MAG in $DBMiFoDB/*; do
	newNAME=MGYGMiFoDB$i
	rename.sh in=$MAG out=$DBallMiFoDB/$newNAME.fna prefix=$newNAME
	rename.sh in=$MAG out=$DBalldual/$newNAME.fna prefix=$newNAME

	echo "renamed $MAG to $newNAME"

 	((i++))
done

i=1
for MAG in $DBcFMD/*; do
	newNAME=MGYGcFMD$i
	rename.sh in=$MAG out=$DBallcFMD/$newNAME.fna prefix=$newNAME
	rename.sh in=$MAG out=$DBalldual/$newNAME.fna prefix=$newNAME

	echo "renamed $MAG to $newNAME"

	((i++))
done

# activate drep
conda activate drep_v3.6.2

# drep own bins + MiFoDB
dRep dereplicate $DBdrepMiFoDB -g $DBallMiFoDB/* -sa 0.95 -nc 0.3 -p 32 # tresholds as mgnify dbs

# drep own bins + cFMD
dRep dereplicate $DBdrepcFMD -g $DBallcFMD/* -sa 0.95 -nc 0.3 -p 32 # tresholds as mgnify dbs

# drep own bins + MiFoDB + cFMD
dRep dereplicate $DBdrepdual -g $DBalldual/* -sa 0.95 -nc 0.3 -p 32 # tresholds as mgnify dbs

# move drep plots to reports folder

cp -r $DBdrepMiFoDB/figures $DBoutMiFoDB
cp -r $DBdrepcFMD/figures $DBoutcFMD
cp -r $DBdrepdual/figures $DBoutdual
